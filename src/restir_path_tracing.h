#pragma once

#include "scene.h"
#include "debug.h"

enum Unbiased {
    NONE = 0,
    NAIVE = 1,
    MIS = 2,
    HMIS = 3
};

Real eval_target_pdf(const Scene& scene,
                   const PathVertex& vertex, 
                   const Vector3& dir_view,
                   const PointAndNormal light_sample,
                   const int light_id,
                   bool visibitliy,
                   int x, int y) {
    Vector3 dir_light = normalize(light_sample.position - vertex.position);
    // assert(neighbour_vertex.material_id != -1);
    if(vertex.material_id == -1) {
        if(debug(x, y)) {
            std::cout << "invalid mat" << std::endl;
        }
        return 0;
    }
    const Material &neighbour_mat = scene.materials[vertex.material_id];

    Real G = max(-dot(dir_light, light_sample.normal), Real(0)) /
                distance_squared(light_sample.position, vertex.position);
    const Light &light = scene.lights[light_id];
    Spectrum L = emission(light, -dir_light, Real(0), light_sample, scene);
    // Vector3 dir_view = -ray.dir;
    //? need to change dir_view?
    Spectrum f = eval(neighbour_mat, dir_view, dir_light, vertex, scene.texture_pool);
    Real target_pdf = luminance(L * G * f);

    if(debug(x, y)) {
        std::cout << "evaling target_pdf ..." << target_pdf << std::endl;
        std::cout << "G " << G << std::endl;
        std::cout << "f " << f << std::endl;
        std::cout << "light pos " << light_sample.position << std::endl;
        std::cout << "vertex pos " << vertex.position << std::endl;
        std::cout << normalize(dir_view) << std::endl;
        std::cout << normalize(vertex.geometric_normal) << std::endl;
        std::cout << dot(vertex.geometric_normal, dir_view) << std::endl;
        std::cout << dot(vertex.geometric_normal, dir_light) << std::endl;
        std::cout << "L " << L << std::endl;

    }
    //todo: visibility check
    if(visibitliy) {
        Ray shadow_ray{vertex.position, dir_light, 
                get_shadow_epsilon(scene),
                (1 - get_shadow_epsilon(scene)) *
                    distance(light_sample.position, vertex.position)};

        if (occluded(scene, shadow_ray)) {
            if(debug(x, y)) {
                std::cout << "blocked" << std::endl;
            }
            return 0.;
        }
    }

    return target_pdf;
}

void regular_ris(const Scene& scene, Reservoir& rsv, pcg32_state& rng, const PathVertex& vertex, const Vector3& dir_view, int x, int y) {
    if(vertex.material_id == -1) return;
    for (int r = 0; r < scene.options.ris_samples; r++) {
        // First, we sample a point on the light source.
        // We do this by first picking a light source, then pick a point on it.
        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        // if(debug(x,y)) {
        //     std::cout << "light uv " << light_uv << std::endl;
        //     std::cout << "light w " << light_w << std::endl;
        //     std::cout << "shape w " << shape_w << std::endl;
        // }
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light = sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);
        // if(debug(x, y))
        //     std::cout << "(loop) light position" << point_on_light.position << std::endl;
        
        Real source_pdf = light_pmf(scene, light_id) *
            pdf_point_on_light(light, point_on_light, vertex.position, scene);

        //* target_pdf = Le * f * G
        //* luminance turn spectrum into real
        // Vector3 dir_light = normalize(point_on_light.position - vertex.position);
        // Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
        //             distance_squared(point_on_light.position, vertex.position);
        // Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
        // Spectrum f = eval(scene.materials[vertex.material_id], dir_view, dir_light, vertex, scene.texture_pool);
        // Real target_pdf = luminance(L * G * f);
        Real target_pdf = eval_target_pdf(scene, vertex, dir_view, point_on_light, light_id, false, x, y);
        

        Real contri_W = 1.0 / source_pdf;
        assert(source_pdf > 0);
        // Real unbiased_reuse_m = 1.0 / scene.options.ris_samples;
        // Real w = target_pdf * contri_W * unbiased_reuse_m;
        Real w = target_pdf * contri_W;

        assert(w == w);
        if(debug(x, y)) {
            std::cout << "in ris: " << r << std::endl;
            std::cout << "contri_w " << contri_W << " w " << w  << " target pdf " << target_pdf << std::endl;
        }

        // Real reservior_r = next_pcg32_real<Real>(first_bounce_rng);
        Real reservior_r = next_pcg32_real<Real>(rng);
        if(debug(x,y)) {
            update_reservoir_debug(rsv, light_id, point_on_light, w, reservior_r);
        } else {
            update_reservoir(rsv, light_id, point_on_light, w, reservior_r);
        }
    }
    rsv.ref_vertex = vertex;
    rsv.prev_dir_view = dir_view;
    return;
}

Real combine_reservoirs(const Scene &scene, 
                        const std::vector<Reservoir>& reservoirs,
                        const PathVertex& vertex, 
                        pcg32_state &rng,
                        const Vector3& dir_view,
                        Reservoir& new_reservoir,
                        int x, int y) {
    Real new_reservoir_M = 0;
    for(int r = 0; r < size(reservoirs); r++) {
        //* target_pdf = Le * f * G
        //* luminance turn spectrum into real
        //* pˆq(ri.y)
        
        // if(reservoirs[r].M == 0) {
        //     continue;
        // }

        Real target_pdf = eval_target_pdf(scene, vertex, dir_view, reservoirs[r].y, reservoirs[r].light_id, true, x, y);
        
        //todo: normalized
        //* pˆq(r .y) · r .W · r .M
        Real w = target_pdf * reservoirs[r].W * reservoirs[r].M;
        
        // Real m = 1.0 / reservoirs[r].M;
        // Real w = m * target_pdf * reservoirs[r].W;

        if(debug(x, y)) {
            std::cout << "combine reservior " << r << std::endl;
            std::cout << "vertex " << reservoirs[r].ref_vertex.position << std::endl;
            std::cout << "sample y " << reservoirs[r].y.position << std::endl;
            std::cout << "w " << w << std::endl;
            std::cout << "ri.W " << reservoirs[r].W  << std::endl;
            std::cout << "ri.M " << reservoirs[r].M  << std::endl;
            std::cout << "target " << target_pdf  << std::endl;
            if(update_reservoir_debug(new_reservoir, reservoirs[r].light_id, reservoirs[r].y, w, next_pcg32_real<Real>(rng))) {
                new_reservoir.ref_vertex = reservoirs[r].ref_vertex;
                new_reservoir.prev_dir_view = reservoirs[r].prev_dir_view;
            }
        } else {
            if(update_reservoir(new_reservoir, reservoirs[r].light_id, reservoirs[r].y, w, next_pcg32_real<Real>(rng))) {
                new_reservoir.ref_vertex = reservoirs[r].ref_vertex;
                new_reservoir.prev_dir_view = reservoirs[r].prev_dir_view;
            }
        }

        new_reservoir_M += reservoirs[r].M;
    }
    if(debug(x, y)) {
        std::cout << "after combining ..." << std::endl;
        std::cout << "chose sample " << new_reservoir.y.position << std::endl;
    }
    return new_reservoir_M;
}

/// restir path tracing (DI)
Spectrum restir_path_tracing(const Scene &scene,
                      int x, int y, /* pixel coordinates */
                      pcg32_state &rng, 
                      Reservoir& rsv, 
                      bool reuse,
                      std::vector<Reservoir> &reservoirs) {
    
    int w = scene.camera.width, h = scene.camera.height;
    Real x_r = next_pcg32_real<Real>(rng);
    Real y_r = next_pcg32_real<Real>(rng);

    Vector2 screen_pos((x + x_r) / w,
                       (y + y_r) / h);
    // if(debug(x,y)) {
    //     std::cout << "=========== new round ==============" << std::endl;
    //     std::cout << "x r" << x_r << std::endl;
    //     std::cout << "y r" << y_r << std::endl;
    // }
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (!vertex_) {
        // Hit background. Account for the environment map if needed.
        if (has_envmap(scene)) {
            const Light &envmap = get_envmap(scene);
            return emission(envmap,
                            -ray.dir, // pointing outwards from light
                            ray_diff.spread,
                            PointAndNormal{}, // dummy parameter for envmap
                            scene);
        }
        return make_zero_spectrum();
    }
    PathVertex vertex = *vertex_;

    Spectrum radiance = make_zero_spectrum();
    // A path's contribution is 
// C(v) = W(v0, v1) * G(v0, v1) * f(v0, v1, v2) * 
    //                    G(v1, v2) * f(v1, v2, v3) * 
    //                  ........
    //                  * G(v_{n-1}, v_n) * L(v_{n-1}, v_n)
    // where v is the path vertices, W is the sensor response
    // G is the geometry term, f is the BSDF, L is the emission
    //
    // "sample_primary" importance samples both W and G,
    // and we assume it always has weight 1.

    // current_path_throughput stores the ratio between
    // 1) the path contribution from v0 up to v_{i} (the BSDF f(v_{i-1}, v_i, v_{i+1}) is not included), 
    // where i is where the PathVertex "vertex" lies on, and
    // 2) the probability density for computing the path v from v0 up to v_i,
    // so that we can compute the Monte Carlo estimates C/p. 
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    // eta_scale stores the scale introduced by Snell-Descartes law to the BSDF (eta^2).
    // We use the same Russian roulette strategy as Mitsuba/pbrt-v3
    // and tracking eta_scale and removing it from the
    // path contribution is crucial for many bounces of refraction.
    Real eta_scale = Real(1);

    // We hit a light immediately. 
    // This path has only two vertices and has contribution
    // C = W(v0, v1) * G(v0, v1) * L(v0, v1)
    if (is_light(scene.shapes[vertex.shape_id])) {
        radiance += current_path_throughput *
            emission(vertex, -ray.dir, scene);
    }

    // We iteratively sum up path contributions from paths with different number of vertices
    // If max_depth == -1, we rely on Russian roulette for path termination.
    int max_depth = scene.options.max_depth;
    bool first_bounce = true;
    for (int num_vertices = 3; max_depth == -1 || num_vertices <= max_depth + 1; num_vertices++) {
        // We are at v_i, and all the path contribution on and before has been accounted for.
        // Now we need to somehow generate v_{i+1} to account for paths with more vertices.
        // In path tracing, we generate two vertices:
        // 1) we sample a point on the light source (often called "Next Event Estimation")
        // 2) we randomly trace a ray from the surface point at v_i and hope we hit something.
        //
        // The first importance samples L(v_i, v_{i+1}), and the second
        // importance samples f(v_{i-1}, v_i, v_{i+1}) * G(v_i, v_{i+1})
        //
        // We then combine the two sampling strategies to estimate the contribution using weighted average.
        // Say the contribution of the first sampling is C1 (with probability density p1), 
        // and the contribution of the second sampling is C2 (with probability density p2,
        // then we compute the estimate as w1*C1/p1 + w2*C2/p2.
        //
        // Assuming the vertices for C1 is v^1, and v^2 for C2,
        // Eric Veach showed that it is a good idea setting 
        // w1 = p_1(v^1)^k / (p_1(v^1)^k + p_2(v^1)^k)
        // w2 = p_2(v^2)^k / (p_1(v^2)^k + p_2(v^2)^k),
        // where k is some scalar real number, and p_a(v^b) is the probability density of generating
        // vertices v^b using sampling method "a".
        // We will set k=2 as suggested by Eric Veach.

        // Finally, we set our "next vertex" in the loop to the v_{i+1} generated
        // by the second sampling, and update current_path_throughput using
        // our hemisphere sampling.

        // Let's implement this!
        assert(vertex.material_id != -1);
        const Material &mat = scene.materials[vertex.material_id];

        //* get the light sample        
        //* if it's the first bounce, use ris
        //* otherwise, use normal light sampling
        int light_id = 0;
        PointAndNormal point_on_light;
        Real unbiased_reuse_m = 0.;
        
        if (!first_bounce) {
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            // if(debug(x,y)) {
            //     std::cout << "bounce" << num_vertices << std::endl;
            //     std::cout << "light uv " << light_uv << std::endl;
            //     std::cout << "light w " << light_w << std::endl;
            //     std::cout << "shape w " << shape_w << std::endl;
            // }
            light_id = sample_light(scene, light_w);
            point_on_light = sample_point_on_light(scene.lights[light_id], vertex.position, light_uv, shape_w, scene);
        } else {
            //*=============
            //* RIS
            //*=============
            if(!reuse) {
                regular_ris(scene, rsv, rng, vertex, -ray.dir, x, y);
                assert(rsv.M > 0);
                unbiased_reuse_m = 1.0 / rsv.M;
                light_id = rsv.light_id;
                point_on_light = rsv.y;
                // if(debug(x, y)) {
                //     std::cout << "first pass store rsv as" << std::endl;
                //     std::cout << "ref_vertex " << rsv.ref_vertex.position << std::endl;
                //     std::cout << "prev_dir_view " << rsv.prev_dir_view << std::endl;
                //     std::cout << "M " << rsv.M << std::endl;
                // }
            } else {
                switch(scene.options.unbiased) {
                case static_cast<int>(Unbiased::NONE): {
                    regular_ris(scene, rsv, rng, vertex, -ray.dir, x, y);


                    //* reuse with bias(Algorithm (4))
                    Reservoir new_reservoir = init_reservoir();
                    reservoirs.push_back(rsv);
                    if(debug(x, y)) {
                        std::cout << "reuse reserveroir number " << size(reservoirs) << std::endl;
                    }
                    new_reservoir.M = combine_reservoirs(scene, reservoirs, vertex, rng, -ray.dir, new_reservoir, x, y);
                    //todo: if M == 0 w_sum == 0, no selected reservoir, the y would be zero
                    //todo: just give up this DI...jump to bsdf
                    light_id = new_reservoir.light_id;
                    point_on_light = new_reservoir.y;
                    rsv = new_reservoir; // reuse
                    unbiased_reuse_m = 1.0 / rsv.M;
                    if(debug(x, y)) {
                        std::cout << "rsv.w_sum " << rsv.w_sum << std::endl;
                        std::cout << "rsv.M " << rsv.M << std::endl;
                        std::cout << "rsv.m " << unbiased_reuse_m << std::endl;
                        for (int r = 0; r < size(reservoirs); r++) {
                            std::cout << "weight of rsv " << r << ":" << reservoirs[r].M / rsv.M << std::endl;
                        }
                    }
                    break;
                }
                case static_cast<int>(Unbiased::NAIVE): {
                    //* unbias reuse (Algorithm (6))
                    regular_ris(scene, rsv, rng, vertex, -ray.dir, x, y);
                    reservoirs.push_back(rsv);
                    Reservoir new_reservoir = init_reservoir();
                    new_reservoir.M = combine_reservoirs(scene, reservoirs, vertex, rng, -ray.dir, new_reservoir, x, y);
                    light_id = new_reservoir.light_id;
                    point_on_light = new_reservoir.y;
                    rsv = new_reservoir; // reuse 
                    
                    // light_id = rsv.light_id;
                    // point_on_light = rsv.y;
                    // unbiased_reuse_m = 1. / rsv.M;

                    //* Z add up all the effective samples number
                    Real Z = 0;
                    for(int r = 0; r < size(reservoirs); r++) {
                        //* target_pdf = Le * f * G
                        //* luminance turn spectrum into real
                        //* pˆqi(r.y)
                        
                        if (eval_target_pdf(scene, reservoirs[r].ref_vertex, reservoirs[r].prev_dir_view, new_reservoir.y, new_reservoir.light_id, true, x, y) > 0) {
                            Z += reservoirs[r].M;
                            if(debug(x, y)) {
                                std::cout << "r.M " << reservoirs[r].M << std::endl;
                                std::cout << "Z " << Z << std::endl;
                            }
                        }
                    }
                    if(debug(x, y)) {
                        std::cout << "new ris weight" << std::endl;
                        for(int r = 0; r < size(reservoirs); r++) {
                            std::cout << "reservoir " << r << " reservoirs[r].M / Z: " << reservoirs[r].M / Z << std::endl;
                        }
                    } 
                    
                    if (Z > 0) {
                        unbiased_reuse_m = 1. / Z;
                    } else {
                        unbiased_reuse_m = 0;
                    }
                    
                    if(debug(x, y)) {
                        std::cout << "rsv.w_sum " << rsv.w_sum << std::endl;
                        std::cout << "rsv.M " << rsv.M << std::endl;
                        std::cout << "rsv.m " << unbiased_reuse_m << std::endl;
                    }
                    break;
                }
                case static_cast<int>(Unbiased::MIS): {
                    //* unbias reuse mis version (supplement Algorithm (1))
                    Reservoir new_reservoir = init_reservoir();
                    Real new_reservoir_M = 0;
                    int selected_reservoir_id = 0;
                    
                    for(int r = 0; r < size(reservoirs); r++) {
                        // if (reservoirs[r].W == 0) {
                        //     //* some reservior is "empty"
                        //     //* their W is not updated, which means its G is 0 
                        //     if(debug(x, y)) {
                        //         std::cout << "reserveroir " << r << "is useless" << std::endl;
                        //     }
                        //     continue;
                        // }
                        //* target_pdf = Le * f * G
                        //* luminance turn spectrum into real
                        //* pˆq(ri.y)
                        Vector3 dir_light = normalize(reservoirs[r].y.position - vertex.position);
                        Real G = max(-dot(dir_light, reservoirs[r].y.normal), Real(0)) /
                                    distance_squared(reservoirs[r].y.position, vertex.position);
                        const Light &light = scene.lights[reservoirs[r].light_id];
                        Spectrum L = emission(light, -dir_light, Real(0), reservoirs[r].y, scene);
                        Vector3 dir_view = -ray.dir;
                        Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
                        Real target_pdf = luminance(L * G * f);

                        // pˆq(r .y) · r .W · r .M
                        Real w = target_pdf * reservoirs[r].W;
                        if(debug(x, y)) {
                            std::cout << "reservior " << r << std::endl;
                            std::cout << "vertex " << reservoirs[r].ref_vertex.position << std::endl;
                            std::cout << "sample y " << reservoirs[r].y.position << std::endl;
                            std::cout << "w " << w << std::endl;
                            std::cout << "ri.W " << reservoirs[r].W  << std::endl;
                            std::cout << "ri.M " << reservoirs[r].M  << std::endl;
                            std::cout << "target " << target_pdf  << std::endl;
                            std::cout << "G " << G  << std::endl;
                            std::cout << "f " << f  << std::endl;
                            std::cout << "L " << L  << std::endl;
                            if(update_reservoir_debug(new_reservoir, reservoirs[r].light_id, reservoirs[r].y, w, next_pcg32_real<Real>(rng))) {
                                new_reservoir.ref_vertex = reservoirs[r].ref_vertex;
                                new_reservoir.prev_dir_view = reservoirs[r].prev_dir_view;
                                selected_reservoir_id = r;
                            }
                        } else {
                            if(update_reservoir(new_reservoir, reservoirs[r].light_id, reservoirs[r].y, w, next_pcg32_real<Real>(rng))) {
                                new_reservoir.ref_vertex = reservoirs[r].ref_vertex;
                                new_reservoir.prev_dir_view = reservoirs[r].prev_dir_view;
                                selected_reservoir_id = r;
                            }
                        }
                        new_reservoir_M += reservoirs[r].M;
                    }
                    if(debug(x,y)) {
                        std::cout << "after combined" << std::endl;
                        std::cout << "sample y " << new_reservoir.y.position << std::endl;
                    }

                    //* p_sum add up all the effective \hat{p}i(r.y) at neighbors
                    Real p_sum = 0;
                    Real selected_pdf = 0;
                    for(int r = 0; r < size(reservoirs); r++) {
                        //* target_pdf = Le * f * G
                        //* luminance turn spectrum into real
                        //* pˆqi(r.y)
                        PathVertex neighbour_vertex = reservoirs[r].ref_vertex;
                        PointAndNormal light_sample = new_reservoir.y;
                        Vector3 dir_light = normalize(light_sample.position - neighbour_vertex.position);
                        assert(neighbour_vertex.material_id != -1);
                        const Material &neighbour_mat = scene.materials[neighbour_vertex.material_id];

                        // assert(distance_squared(light_sample.position, neighbour_vertex.position) > 0);
                        Real dist = distance_squared(light_sample.position, neighbour_vertex.position);
                        if (dist == 0) {
                            continue;
                        }
                        Real G = max(-dot(dir_light, light_sample.normal), Real(0)) / dist;
                        const Light &light = scene.lights[new_reservoir.light_id];
                        Spectrum L = emission(light, -dir_light, Real(0), light_sample, scene);
                        // Vector3 dir_view = -ray.dir;
                        Vector3 dir_view = reservoirs[r].prev_dir_view;
                        //? need to change dir_view?
                        Spectrum f = eval(neighbour_mat, dir_view, dir_light, neighbour_vertex, scene.texture_pool);
                        Real target_pdf = luminance(L * G * f);

                        if(r == selected_reservoir_id) {
                            selected_pdf = target_pdf;
                        }

                        if(debug(x, y)) {
                            std::cout << "compare target pdf with 0" << std::endl;
                            std::cout << "target_pdf " << target_pdf << std::endl;
                            std::cout << "G " << G << std::endl;
                            std::cout << "f " << f << std::endl;
                            std::cout << "light pos " << light_sample.position << std::endl;
                            std::cout << "vertex pos " << neighbour_vertex.position << std::endl;
                            std::cout << normalize(dir_view) << std::endl;
                            std::cout << normalize(neighbour_vertex.geometric_normal) << std::endl;
                            std::cout << dot(neighbour_vertex.geometric_normal, dir_view) << std::endl;
                            std::cout << dot(neighbour_vertex.geometric_normal, dir_light) << std::endl;
                            std::cout << "L " << L << std::endl;

                        }
                        //todo: visibility check
                        Ray shadow_ray{neighbour_vertex.position, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(light_sample.position, neighbour_vertex.position)};

                        if (target_pdf > 0 && !occluded(scene, shadow_ray)) {
                            p_sum += target_pdf;
                            if(debug(x, y)) {
                                std::cout << "r.M " << reservoirs[r].M << std::endl;
                                std::cout << "add target_pdf " << target_pdf << " of the " << r << " reversor" << std::endl;

                                std::cout << "p_sum " << p_sum << std::endl;
                            }
                        }
                    }

                    light_id = new_reservoir.light_id;
                    point_on_light = new_reservoir.y;
                    new_reservoir.M = new_reservoir_M;
                    // assert(p_sum > 0);
                    if(p_sum == 0) {
                        //* the p_sum is zero, 
                        unbiased_reuse_m = 0;
                        if(debug(x, y)) {
                            std::cout << "empty m!" << std::endl;

                        }
                    } else {
                        unbiased_reuse_m = selected_pdf / p_sum;
                    }
                    rsv = new_reservoir; // reuse 
                    
                    if(debug(x, y)) {
                        std::cout << "rsv.w_sum " << rsv.w_sum << std::endl;
                        std::cout << "rsv.M " << rsv.M << std::endl;
                        std::cout << "rsv.m " << unbiased_reuse_m << std::endl;
                        std::cout << "select pdf " << selected_pdf << std::endl;
                        std::cout << "p_sum " << p_sum << std::endl;
                    }
                    break;
                }
                case static_cast<int>(Unbiased::HMIS): {
                    //* unbias reuse mis version (supplement Algorithm (1))
                    Reservoir new_reservoir = init_reservoir();
                    Real new_reservoir_M = 0;
                    int selected_reservoir_id = 0;
                    std::vector<Reservoir> H_reservoirs;
                    Reservoir current_pixel_reservoir = reservoirs.back();
                    for(int r = 0; r < size(reservoirs); r++) {
                        //* heuristics trick
                        //* |dot(ni, nj)| < cos(25 degree)
                        //* (|zi - zj|)/|zj| > 0.1
                        if (fabs(dot(reservoirs[r].ref_vertex.geometric_normal, 
                            current_pixel_reservoir.ref_vertex.geometric_normal)) < 0.9 && 
                            (fabs(reservoirs[r].ref_vertex.position.z - current_pixel_reservoir.ref_vertex.position.z) / fabs(current_pixel_reservoir.ref_vertex.position.z)) > 0.1) {
                            if(debug(x, y)) {
                                std::cout << "reserveroir " << r << "is useless" << std::endl;
                            }
                            continue;
                        }
                        H_reservoirs.push_back(reservoirs[r]);
                        //* target_pdf = Le * f * G
                        //* luminance turn spectrum into real
                        //* pˆq(ri.y)
                        Vector3 dir_light = normalize(reservoirs[r].y.position - vertex.position);
                        Real G = max(-dot(dir_light, reservoirs[r].y.normal), Real(0)) /
                                    distance_squared(reservoirs[r].y.position, vertex.position);
                        const Light &light = scene.lights[reservoirs[r].light_id];
                        Spectrum L = emission(light, -dir_light, Real(0), reservoirs[r].y, scene);
                        Vector3 dir_view = -ray.dir;
                        Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
                        Real target_pdf = luminance(L * G * f);

                        // pˆq(r .y) · r .W · r .M
                        Real w = target_pdf * reservoirs[r].W;
                        if(debug(x, y)) {
                            std::cout << "reservior " << r << std::endl;
                            std::cout << "vertex " << reservoirs[r].ref_vertex.position << std::endl;
                            std::cout << "sample y " << reservoirs[r].y.position << std::endl;
                            std::cout << "w " << w << std::endl;
                            std::cout << "ri.W " << reservoirs[r].W  << std::endl;
                            std::cout << "ri.M " << reservoirs[r].M  << std::endl;
                            std::cout << "target " << target_pdf  << std::endl;
                            std::cout << "G " << G  << std::endl;
                            std::cout << "f " << f  << std::endl;
                            std::cout << "L " << L  << std::endl;
                            if(update_reservoir_debug(new_reservoir, reservoirs[r].light_id, reservoirs[r].y, w, next_pcg32_real<Real>(rng))) {
                                new_reservoir.ref_vertex = reservoirs[r].ref_vertex;
                                new_reservoir.prev_dir_view = reservoirs[r].prev_dir_view;
                                selected_reservoir_id = size(H_reservoirs) - 1;
                            }
                        } else {
                            if(update_reservoir(new_reservoir, reservoirs[r].light_id, reservoirs[r].y, w, next_pcg32_real<Real>(rng))) {
                                new_reservoir.ref_vertex = reservoirs[r].ref_vertex;
                                new_reservoir.prev_dir_view = reservoirs[r].prev_dir_view;
                                selected_reservoir_id = size(H_reservoirs) - 1;
                            }
                        }
                        new_reservoir_M += reservoirs[r].M;
                    }
                    if(debug(x,y)) {
                        std::cout << "after combined" << std::endl;
                        std::cout << "sample y " << new_reservoir.y.position << std::endl;
                    }

                    //* p_sum add up all the effective \hat{p}i(r.y) at neighbors
                    Real p_sum = 0;
                    Real selected_pdf = 0;
                    for(int r = 0; r < size(H_reservoirs); r++) {
                        //* target_pdf = Le * f * G
                        //* luminance turn spectrum into real
                        //* pˆqi(r.y)
                        PathVertex neighbour_vertex = H_reservoirs[r].ref_vertex;
                        PointAndNormal light_sample = new_reservoir.y;
                        Vector3 dir_light = normalize(light_sample.position - neighbour_vertex.position);
                        assert(neighbour_vertex.material_id != -1);
                        const Material &neighbour_mat = scene.materials[neighbour_vertex.material_id];

                        //?? why it cannot choose even one sample..
                        Real dist = distance_squared(light_sample.position, neighbour_vertex.position);
                        if (dist == 0) {
                            continue;
                        }
                        Real G = max(-dot(dir_light, light_sample.normal), Real(0)) /
                                    distance_squared(light_sample.position, neighbour_vertex.position);
                        const Light &light = scene.lights[new_reservoir.light_id];
                        Spectrum L = emission(light, -dir_light, Real(0), light_sample, scene);
                        // Vector3 dir_view = -ray.dir;
                        Vector3 dir_view = H_reservoirs[r].prev_dir_view;
                        //? need to change dir_view?
                        Spectrum f = eval(neighbour_mat, dir_view, dir_light, neighbour_vertex, scene.texture_pool);
                        Real target_pdf = luminance(L * G * f);

                        if(r == selected_reservoir_id) {
                            selected_pdf = target_pdf;
                        }
                        if(debug(x, y)) {
                            std::cout << "compare target pdf with 0" << std::endl;
                            std::cout << "target_pdf " << target_pdf << std::endl;
                            std::cout << "G " << G << std::endl;
                            std::cout << "f " << f << std::endl;
                            std::cout << "light pos " << light_sample.position << std::endl;
                            std::cout << "vertex pos " << neighbour_vertex.position << std::endl;
                            std::cout << normalize(dir_view) << std::endl;
                            std::cout << normalize(neighbour_vertex.geometric_normal) << std::endl;
                            std::cout << dot(neighbour_vertex.geometric_normal, dir_view) << std::endl;
                            std::cout << dot(neighbour_vertex.geometric_normal, dir_light) << std::endl;
                            std::cout << "L " << L << std::endl;

                        }
                        //todo: visibility check
                        Ray shadow_ray{neighbour_vertex.position, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(light_sample.position, neighbour_vertex.position)};

                        if (target_pdf > 0 && !occluded(scene, shadow_ray)) {
                            p_sum += target_pdf;
                            if(debug(x, y)) {
                                std::cout << "r.M " << H_reservoirs[r].M << std::endl;
                                std::cout << "add target_pdf " << target_pdf << " of the " << r << " reversor" << std::endl;

                                std::cout << "p_sum " << p_sum << std::endl;
                            }
                        }
                    }

                    light_id = new_reservoir.light_id;
                    point_on_light = new_reservoir.y;
                    new_reservoir.M = new_reservoir_M;
                    // assert(p_sum > 0);
                    if(p_sum == 0) {
                        //* the p_sum is zero, 
                        unbiased_reuse_m = 0;
                        if(debug(x, y)) {
                            std::cout << "empty m!" << std::endl;

                        }
                    } else {
                        unbiased_reuse_m = selected_pdf / p_sum;
                    }
                    rsv = new_reservoir; // reuse 
                    
                    if(debug(x, y)) {
                        std::cout << "rsv.w_sum " << rsv.w_sum << std::endl;
                        std::cout << "rsv.M " << rsv.M << std::endl;
                        std::cout << "rsv.m " << unbiased_reuse_m << std::endl;
                        std::cout << "select pdf " << selected_pdf << std::endl;
                        std::cout << "p_sum " << p_sum << std::endl;
                    }
                    break;
                }
                default:
                    //* unmatch biased method
                    assert(true);
                }
            }
        }
        
        const Light &light = scene.lights[light_id];
        // if(debug(x, y)) {
        //     std::cout << "bounce " << num_vertices << std::endl;
        //     std::cout << "select light " << light_id << " in total " << size(scene.lights) << " lights " << std::endl;
        //     std::cout << "light position" << point_on_light.position << std::endl;
        //     std::cout << "ref vertex" << vertex.position << std::endl; 
        // }


        // Next, we compute w1*C1/p1. We store C1/p1 in C1.
        Spectrum C1 = make_zero_spectrum();
        Real w1 = 0;

        // Remember "current_path_throughput" already stores all the path contribution on and before v_i.
        // So we only need to compute G(v_{i}, v_{i+1}) * f(v_{i-1}, v_{i}, v_{i+1}) * L(v_{i}, v_{i+1})
        {
            // Let's first deal with C1 = G * f * L.
            // Let's first compute G.
            Real G = 0;
            Vector3 dir_light;
            // The geometry term is different between directional light sources and
            // others. Currently we only have environment maps as directional light sources.
            if (!is_envmap(light)) {
                dir_light = normalize(point_on_light.position - vertex.position);
                // If the point on light is occluded, G is 0. So we need to test for occlusion.
                // To avoid self intersection, we need to set the tnear of the ray
                // to a small "epsilon". We set the epsilon to be a small constant times the
                // scale of the scene, which we can obtain through the get_shadow_epsilon() function.
                Ray shadow_ray{vertex.position, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(point_on_light.position, vertex.position)};
                if (!occluded(scene, shadow_ray)) {
                    // geometry term is cosine at v_{i+1} divided by distance squared
                    // this can be derived by the infinitesimal area of a surface projected on
                    // a unit sphere -- it's the Jacobian between the area measure and the solid angle
                    // measure.
                    G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                        distance_squared(point_on_light.position, vertex.position);
                    if (G == 0 && debug(x,y)) {
                        std::cout << "G is 0" << std::endl;
                    }
                } else {
                    if(debug(x, y)) {
                        std::cout << "occluded" << std::endl;
                    }
                }
            } else {
                // The direction from envmap towards the point is stored in
                // point_on_light.normal.
                dir_light = -point_on_light.normal;
                // If the point on light is occluded, G is 0. So we need to test for occlusion.
                // To avoid self intersection, we need to set the tnear of the ray
                // to a small "epsilon" which we define as c_shadow_epsilon as a global constant.
                Ray shadow_ray{vertex.position, dir_light, 
                               get_shadow_epsilon(scene),
                               infinity<Real>() /* envmaps are infinitely far away */};
                if (!occluded(scene, shadow_ray)) {
                    // We integrate envmaps using the solid angle measure,
                    // so the geometry term is 1.
                    G = 1;
                }
            }

            // Before we proceed, we first compute the probability density p1(v1)
            // The probability density for light sampling to sample our point is
            // just the probability of sampling a light times the probability of sampling a point
            Real p1 = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, vertex.position, scene);

            // if(debug(x, y)) {
            //         std::cout << "G: " << G << std::endl;
            //         std::cout << "p1: " << p1 << std::endl;
            // }
            // We don't need to continue the computation if G is 0.
            // Also sometimes there can be some numerical issue such that we generate
            // a light path with probability zero
            if (G > 0 && p1 > 0) {
                // Let's compute f (BSDF) next.
                Vector3 dir_view = -ray.dir;
                assert(vertex.material_id >= 0);
                Spectrum f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);

                // Evaluate the emission
                // We set the footprint to zero since it is not fully clear how
                // to set it in this case.
                // One way is to use a roughness based heuristics, but we have multi-layered BRDFs.
                // See "Real-time Shading with Filtered Importance Sampling" from Colbert et al.
                // for the roughness based heuristics.
                Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);

                // C1 is just a product of all of them!
                C1 = G * f * L;
                if(debug(x,y) && luminance(C1) < 1e6) {
                    std::cout << "C1 too small " << C1 << std::endl;
                    std::cout << "reuse? " << reuse << std::endl;
                    std::cout << "f: " << f << std::endl;
                    std::cout << "L: " << L << std::endl;
                    std::cout << "G: " << G << std::endl;
                }
            
                // Next let's compute w1

                // Remember that we want to set
                // w1 = p_1(v^1)^2 / (p_1(v^1)^2 + p_2(v^1)^2)
                // Notice that all of the probability density share the same path prefix and those cancel out.
                // Therefore we only need to account for the generation of the vertex v_{i+1}.

                // The probability density for our hemispherical sampling to sample 
                Real p2 = pdf_sample_bsdf(
                    mat, dir_view, dir_light, vertex, scene.texture_pool);
                // !!!! IMPORTANT !!!!
                // In general, p1 and p2 now live in different spaces!!
                // our BSDF API outputs a probability density in the solid angle measure
                // while our light probability density is in the area measure.
                // We need to make sure that they are in the same space.
                // This can be done by accounting for the Jacobian of the transformation
                // between the two measures.
                // In general, I recommend to transform everything to area measure 
                // (except for directional lights) since it fits to the path-space math better.
                // Converting a solid angle measure to an area measure is just a
                // multiplication of the geometry term G (let solid angle be dS, area be dA,
                // we have dA/dS = G).
                p2 *= G;

                w1 = (p1*p1) / (p1*p1 + p2*p2); //* keep using p1 to calculate mis 

                // if(debug(x, y)) {
                //     std::cout << "p1" << p1 << std::endl;
                //     std::cout << "p2" << p2 << std::endl;
                //     std::cout << "w1" << w1 << std::endl;
                // }
                
                //* replace 1 / pdf with W
                if (!first_bounce) {
                    C1 /= p1;
                    // if(debug(x, y)) {
                    //     std::cout << "C1: " << C1 << std::endl;
                    // }
                } else {
                    if(rsv.M == 0) {
                        //* empty reservior
                        C1 = make_const_spectrum(0);
                        if(debug(x, y)) {
                            std::cout << "rsv.M is zero" << std::endl;
                        }
                    } else {
                        // if(debug(x, y)) {
                        //     std::cout << "G " << G <<  " F " << f << " Le " << L << std::endl;
                        //     std::cout << "C1 " << C1 <<  " luminace " << luminance(C1) << std::endl;
                        //     std::cout << "M " << rsv.M << " m " << m << " w_sum " << rsv.w_sum << std::endl;
                        // }
                        if (luminance(C1) == 0) {
                            //* if target_pdf is 0, set W to 0
                            rsv.W = 0;
                        } else {
                            // rsv.W = (1. / luminance(C1)) * rsv.w_sum;
                            rsv.W = unbiased_reuse_m * (1. / luminance(C1)) * rsv.w_sum;



                            if(debug(x, y)) {
                                std::cout << "pass " << num_vertices - 3 << std::endl;
                                std::cout << "w_sum " << rsv.w_sum << std::endl;
                                std::cout << "m " << unbiased_reuse_m << std::endl;
                                std::cout << "w_sum * m " << rsv.w_sum * unbiased_reuse_m << std::endl;
                                std::cout << "target pdf inverse " << 1. / luminance(C1) << std::endl;
                                std::cout << "W " << rsv.W << std::endl;
                            }
                        }
                        C1 *= rsv.W;
                    }
                }
                // if (C1.x != C1.x) {
                //     std::cout << "NAN: " << x << ' ' << y << std::endl;
                // }
                // assert(C1.x == C1.x);
                if(debug(x, y)) {
                    std::cout << "C1: " << C1 << std::endl;
                }
            }
        }
        if (first_bounce & ignore_first_bounce()) {
            radiance += make_const_spectrum(0);
            if(debug(x, y))
                std::cout << "first bounce radiance" << radiance << std::endl;
        } else {
            radiance += current_path_throughput * C1 * w1;
            // if(debug(x, y) && first_bounce) {
            //     std::cout << "reuse radiance " <<  radiance << std::endl;
            // }
            if(debug(x, y)) {
                std::cout << "bounce: " << num_vertices << std::endl;
                std::cout << " cpt " << current_path_throughput << std::endl;
                std::cout << " C1 " << C1 << std::endl;
                std::cout << " w1 " << w1 << std::endl;
            

                std::cout << num_vertices << " radiance " << radiance << std::endl;
            }
        }
        first_bounce = false;

        // Let's do the hemispherical sampling next.
        Vector3 dir_view = -ray.dir;
        Vector2 bsdf_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
        std::optional<BSDFSampleRecord> bsdf_sample_ =
            sample_bsdf(mat,
                        dir_view,
                        vertex,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w);
        if (!bsdf_sample_) {
            // BSDF sampling failed. Abort the loop.
            break;
        }
        const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
        Vector3 dir_bsdf = bsdf_sample.dir_out;
        assert(dir_bsdf.x == dir_bsdf.x);
        // Update ray differentials & eta_scale
        if (bsdf_sample.eta == 0) {
            ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
        } else {
            ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
            eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
        }

        // Trace a ray towards bsdf_dir. Note that again we have
        // to have an "epsilon" tnear to prevent self intersection.
        Ray bsdf_ray{vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>()};
        std::optional<PathVertex> bsdf_vertex = intersect(scene, bsdf_ray);

        // if(debug(x, y)) {
        //     std::cout << "bsdf_ray dir" << bsdf_ray.dir << std::endl;
        //     std::cout << "bsdf_ray org" << bsdf_ray.org << std::endl;
        //     std::cout << "sample vertex " << vertex.position << std::endl;
        //     std::cout << "sample dir_view " << dir_view << std::endl;
        //     std::cout << "sample uv " << bsdf_rnd_param_uv << std::endl;
        //     std::cout << "sample w " << bsdf_rnd_param_w << std::endl;
        //     std::cout << "bsdf vertex " << (*bsdf_vertex).position << std::endl;

        // }

        // To update current_path_throughput
        // we need to multiply G(v_{i}, v_{i+1}) * f(v_{i-1}, v_{i}, v_{i+1}) to it
        // and divide it with the pdf for getting v_{i+1} using hemisphere sampling.
        Real G;
        if (bsdf_vertex) {
            G = fabs(dot(dir_bsdf, bsdf_vertex->geometric_normal)) /
                distance_squared(bsdf_vertex->position, vertex.position);
        } else {
            // We hit nothing, set G to 1 to account for the environment map contribution.
            G = 1;
        }

        Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
        Real p2 = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
        if (p2 <= 0) {
            // Numerical issue -- we generated some invalid rays.
            break;
        }

        // Remember to convert p2 to area measure!
        p2 *= G;
        // note that G cancels out in the division f/p, but we still need
        // G later for the calculation of w2.

        // Now we want to check whether dir_bsdf hit a light source, and
        // account for the light contribution (C2 & w2 & p2).
        // There are two possibilities: either we hit an emissive surface,
        // or we hit an environment map.
        // We will handle them separately.
        if (bsdf_vertex && is_light(scene.shapes[bsdf_vertex->shape_id])) {
            // G & f are already computed.
            Spectrum L = emission(*bsdf_vertex, -dir_bsdf, scene);
            Spectrum C2 = G * f * L;
            // Next let's compute p1(v2): the probability of the light source sampling
            // directly drawing the point corresponds to bsdf_dir.
            int light_id = get_area_light_id(scene.shapes[bsdf_vertex->shape_id]);
            assert(light_id >= 0);
            const Light &light = scene.lights[light_id];
            PointAndNormal light_point{bsdf_vertex->position, bsdf_vertex->geometric_normal};
            Real p1 = light_pmf(scene, light_id) *
                pdf_point_on_light(light, light_point, vertex.position, scene);
            Real w2 = (p2*p2) / (p1*p1 + p2*p2);

            C2 /= p2;
            radiance += current_path_throughput * C2 * w2;
        } else if (!bsdf_vertex && has_envmap(scene)) {
            // G & f are already computed.
            const Light &light = get_envmap(scene);
            Spectrum L = emission(light,
                                  -dir_bsdf, // pointing outwards from light
                                  ray_diff.spread,
                                  PointAndNormal{}, // dummy parameter for envmap
                                  scene);
            Spectrum C2 = G * f * L;
            // Next let's compute p1(v2): the probability of the light source sampling
            // directly drawing the direction bsdf_dir.
            PointAndNormal light_point{Vector3{0, 0, 0}, -dir_bsdf}; // pointing outwards from light
            Real p1 = light_pmf(scene, scene.envmap_light_id) *
                      pdf_point_on_light(light, light_point, vertex.position, scene);
            Real w2 = (p2*p2) / (p1*p1 + p2*p2);

            C2 /= p2;
            radiance += current_path_throughput * C2 * w2;
            if (debug(x, y)) {
                std::cout << "bounce: " << num_vertices << std::endl;
                std::cout << " cpt " << current_path_throughput << std::endl;
                std::cout << " C2 " << C1 << std::endl;
                std::cout << " w2 " << w1 << std::endl;
                std::cout << num_vertices << " radiance " << radiance << std::endl;
            }
        }

        if (!bsdf_vertex) {
            // Hit nothing -- can't continue tracing.
            if (debug(x, y)) {
                std::cout << "cant hit anymore" << num_vertices << std::endl;
            }
            break;
        }

        // Update rays/intersection/current_path_throughput/current_pdf
        // Russian roulette heuristics
        Real rr_prob = 1;
        if (num_vertices - 1 >= scene.options.rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                if (debug(x, y)) {
                    std::cout << "terminated" << num_vertices << std::endl;
                }
                break;
            }
        }

        ray = bsdf_ray;
        vertex = *bsdf_vertex;
        current_path_throughput = current_path_throughput * (G * f) / (p2 * rr_prob);
        // if(current_path_throughput.x > 1 && current_path_throughput.y > 1 && current_path_throughput.z > 1) {
        //         std::cout << p2 << ' ' << f << ' ' << current_path_throughput << std::endl;
        // }
    }

    if (debug(x, y)) {
        std::cout << "==================" << radiance << std::endl;
    }

    return radiance;
}
