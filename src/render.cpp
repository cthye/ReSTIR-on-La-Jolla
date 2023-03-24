#include "render.h"
#include "intersection.h"
#include "material.h"
#include "parallel.h"
#include "path_tracing.h"
#include "vol_path_tracing.h"
#include "restir_path_tracing.h"
// #include "pcg.h"
#include "progress_reporter.h"
#include "scene.h"
#include "debug.h"


/// Render auxiliary buffers e.g., depth.
Image3 aux_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    parallel_for([&](const Vector2i &tile) {
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Ray ray = sample_primary(scene.camera, Vector2((x + Real(0.5)) / w, (y + Real(0.5)) / h));
                RayDifferential ray_diff = init_ray_differential(w, h);
                if (std::optional<PathVertex> vertex = intersect(scene, ray, ray_diff)) {
                    Real dist = distance(vertex->position, ray.org);
                    Vector3 color{0, 0, 0};
                    if (scene.options.integrator == Integrator::Depth) {
                        color = Vector3{dist, dist, dist};
                    } else if (scene.options.integrator == Integrator::ShadingNormal) {
                        // color = (vertex->shading_frame.n + Vector3{1, 1, 1}) / Real(2);
                        color = vertex->shading_frame.n;
                    } else if (scene.options.integrator == Integrator::MeanCurvature) {
                        Real kappa = vertex->mean_curvature;
                        color = Vector3{kappa, kappa, kappa};
                    } else if (scene.options.integrator == Integrator::RayDifferential) {
                        color = Vector3{ray_diff.radius, ray_diff.spread, Real(0)};
                    } else if (scene.options.integrator == Integrator::MipmapLevel) {
                        const Material &mat = scene.materials[vertex->material_id];
                        const TextureSpectrum &texture = get_texture(mat);
                        auto *t = std::get_if<ImageTexture<Spectrum>>(&texture);
                        if (t != nullptr) {
                            const Mipmap3 &mipmap = get_img3(scene.texture_pool, t->texture_id);
                            Vector2 uv{modulo(vertex->uv[0] * t->uscale, Real(1)),
                                       modulo(vertex->uv[1] * t->vscale, Real(1))};
                            // ray_diff.radius stores approximatedly dpdx,
                            // but we want dudx -- we get it through
                            // dpdx / dpdu
                            Real footprint = vertex->uv_screen_size;
                            Real scaled_footprint = max(get_width(mipmap), get_height(mipmap)) *
                                                    max(t->uscale, t->vscale) * footprint;
                            Real level = log2(max(scaled_footprint, Real(1e-8f)));
                            color = Vector3{level, level, level};
                        }
                    }
                    img(x, y) = color;
                } else {
                    img(x, y) = Vector3{0, 0, 0};
                }
            }
        }
    }, Vector2i(num_tiles_x, num_tiles_y));

    return img;
}

Image3 path_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);
    std::cout << "light number: " << size(scene.lights) << std::endl;

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        // Use a different rng stream for each thread.
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);

        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {               
                Spectrum radiance = make_zero_spectrum();
                int spp = scene.options.samples_per_pixel;
                for (int s = 0; s < spp; s++) {
                    radiance += path_tracing(scene, x, y, rng);
                }
                img(x, y) = radiance / Real(spp);
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    return img;
}


const std::vector<Vector2> &direction = {{0, -1}, {0, 1}, {-1, 0}, {1, 0}};

Image3 restir_path_render(const Scene &scene) {
    std::cout << "restir's M: " << scene.options.ris_samples << std::endl;
    std::cout << "restir's unbiased method: " << scene.options.unbiased << std::endl;
    std::cout << "light number: " << size(scene.lights) << std::endl;


    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);
    ImageReservoir imgReservoir(w, h);
    ImageReservoir prevImgReservoir(w, h);


    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    //* want to reuse the random number seed per tile
    std::vector<std::vector<pcg32_state> > rngs(num_tiles_x);
    for (int rng_i = 0; rng_i < num_tiles_x; rng_i++) {
        rngs[rng_i] = std::vector<pcg32_state>(num_tiles_y);

        for (int rng_j = 0; rng_j < num_tiles_y; rng_j++) {
            rngs[rng_i][rng_j] = init_pcg32(rng_j * num_tiles_x + rng_i);
        }
    }
    int spp = scene.options.samples_per_pixel;    
    ProgressReporter reporter(spp);
    for (int s = 0; s < spp; s++) {
        //! how to wait for all the threads finish for each round...
        parallel_for([&](const Vector2i &tile) {
            // Use a different rng stream for each thread.
            //pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
            // pcg32_state first_rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);

            int x0 = tile[0] * tile_size;
            int x1 = min(x0 + tile_size, w);
            int y0 = tile[1] * tile_size;
            int y1 = min(y0 + tile_size, h);
            for (int y = y0; y < y1; y++) {
                for (int x = x0; x < x1; x++) {
                    Reservoir rsv = init_reservoir();
                    std::vector<Reservoir> reservoirs;
                    // reservoirs.push_back(prevImgReservoir(x, y));

                    // std::cout << "x " << x << "y " << y << std::endl;
                    if(s > 0) {
                        for(int d = 0; d < 4; d++) {
                            int neighbour_x = x + 30 * (1 -  next_pcg32_real<Real>(rngs[tile[0]][tile[1]]) * 2);
                            int neighbour_y = y + 30 * (1 - next_pcg32_real<Real>(rngs[tile[0]][tile[1]]) * 2);
                            // if(neighbour_x < 0 || neighbour_x >= w || neighbour_y < 0 || neighbour_y >= h) {
                            //     continue;
                            // }
                            // if(neighbour_x < x0 || neighbour_x >= x1 || neighbour_y < y0 || neighbour_y >= y1) {
                            //     continue;
                            // }
                            if(neighbour_x < 0 || neighbour_x >= w || neighbour_y < 0 || neighbour_y >= h) {
                                continue;
                            }
                            reservoirs.push_back(prevImgReservoir(neighbour_x, neighbour_y));
                            if(debug(x,y)) {
                                std::cout << "======segment fault location log======" << std::endl;
                            }
                        }
                    }
                    if(debug(x,y) && s) {
                        std::cout << "vector size " << size(reservoirs) << std::endl;
                    }
                    bool reuse = !(s==0);
                    // bool reuse = false;

                    Spectrum radiance = restir_path_tracing(scene, x, y, rngs[tile[0]][tile[1]], rsv, reuse, reservoirs);
                    // Spectrum radiance = restir_path_tracing(scene, x, y, rng, rsv, reuse, reservoirs);
                    
                    // std::cout << "radiance " << radiance << std::endl;
                    img(x, y) += radiance / Real(spp);
                    if(debug(x,y)) {
                        std::cout << "img(x,y) " << img(x,y) << std::endl;
                    }
                    imgReservoir(x, y) = rsv;
                    
                }
            }
        }, Vector2i(num_tiles_x, num_tiles_y));

        reporter.update(1);
        prevImgReservoir = imgReservoir;
    }

    reporter.done();
    return img;
}

Image3 vol_path_render(const Scene &scene) {
    int w = scene.camera.width, h = scene.camera.height;
    Image3 img(w, h);

    constexpr int tile_size = 16;
    int num_tiles_x = (w + tile_size - 1) / tile_size;
    int num_tiles_y = (h + tile_size - 1) / tile_size;

    auto f = vol_path_tracing;
    if (scene.options.vol_path_version == 1) {
        f = vol_path_tracing_1;
    } else if (scene.options.vol_path_version == 2) {
        f = vol_path_tracing_2;
    } else if (scene.options.vol_path_version == 3) {
        f = vol_path_tracing_3;
    } else if (scene.options.vol_path_version == 4) {
        f = vol_path_tracing_4;
    } else if (scene.options.vol_path_version == 5) {
        f = vol_path_tracing_5;
    } else if (scene.options.vol_path_version == 6) {
        f = vol_path_tracing;
    }

    ProgressReporter reporter(num_tiles_x * num_tiles_y);
    parallel_for([&](const Vector2i &tile) {
        // Use a different rng stream for each thread.
        pcg32_state rng = init_pcg32(tile[1] * num_tiles_x + tile[0]);
        int x0 = tile[0] * tile_size;
        int x1 = min(x0 + tile_size, w);
        int y0 = tile[1] * tile_size;
        int y1 = min(y0 + tile_size, h);
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                Spectrum radiance = make_zero_spectrum();
                int spp = scene.options.samples_per_pixel;
                for (int s = 0; s < spp; s++) {
                    Spectrum L = f(scene, x, y, rng);
                    if (isfinite(L)) {
                        // Hacky: exclude NaNs in the rendering.
                        radiance += L;
                    }
                }
                img(x, y) = radiance / Real(spp);
            }
        }
        reporter.update(1);
    }, Vector2i(num_tiles_x, num_tiles_y));
    reporter.done();
    return img;
}


Image3 render(const Scene &scene) {
    if (scene.options.integrator == Integrator::Depth ||
            scene.options.integrator == Integrator::ShadingNormal ||
            scene.options.integrator == Integrator::MeanCurvature ||
            scene.options.integrator == Integrator::RayDifferential ||
            scene.options.integrator == Integrator::MipmapLevel) {
        return aux_render(scene);
    } else if (scene.options.integrator == Integrator::Path) {
        return path_render(scene);
    } else if (scene.options.integrator == Integrator::VolPath) {
        return vol_path_render(scene);
    } else if (scene.options.integrator == Integrator::ReSTIR) {
        return restir_path_render(scene);
    } else {
        assert(false);
        return Image3();
    }
}
