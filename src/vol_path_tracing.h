#pragma once
#include <scene.h>
#include <pcg.h>
#include <phase_function.h>

#include "debug.h"
#define INF infinity<Real>()

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    //disable ray differentials
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (!vertex_) {
        return make_zero_spectrum();
    }
    Spectrum radiance = make_zero_spectrum();
    PathVertex vertex = *vertex_;
    Spectrum sigma_a = get_sigma_a(scene.media[vertex.exterior_medium_id], vertex.position);
    Real t = length(vertex.position - ray.org);
    Spectrum transmittance = exp(- sigma_a * t);
    Spectrum Le = make_zero_spectrum();
    if(is_light(scene.shapes[vertex.shape_id])) {
        Le = emission(vertex, -ray.dir, scene);
    }
    radiance = Le * transmittance;
    return radiance;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    //disable ray differentials
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    
    Real t_hit = infinity<Real>();
    //* homogenous, sigma_a & sigma_t can be derived from any points in the medium...
    Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], make_zero_spectrum());
    Spectrum sigma_s = get_sigma_s(scene.media[scene.camera.medium_id], make_zero_spectrum());
    Spectrum sigma_t = sigma_a + sigma_s;
    PhaseFunction phase_function = get_phase_function(scene.media[scene.camera.medium_id]);

    if (vertex_) {
        t_hit = distance((*vertex_).position, ray.org);
    }
    
    Real u = next_pcg32_real<Real>(rng);
    Spectrum t = - log(1 - u) / sigma_t;
    if (t.x < t_hit) {
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Vector3 p = ray.org + t * ray.dir;

        //equation 7
        Spectrum L_s1_estimate = make_zero_spectrum();
        Real L_s1_pdf = 1.;

        Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light &light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, p, light_uv, shape_w, scene);

        Vector3 dir_light = normalize(point_on_light.position - p);
        Ray shadow_ray{p, dir_light, get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                        distance(point_on_light.position, p)};
        if (!occluded(scene, shadow_ray)) {
            Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                        distance_squared(point_on_light.position, p);
        
            L_s1_pdf = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, p, scene);
            Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
            L_s1_estimate = L * exp(-sigma_t * distance(p, point_on_light.position)) * eval(phase_function, -ray.dir, dir_light) * G;
        }
        // if(!vertex_) std::cout << L_s1_estimate << std::endl;
        return transmittance / trans_pdf * sigma_s * L_s1_estimate / L_s1_pdf;
    } 

    //* if ray doesn't hit an object, t_hit would be infinity and never reach here
    Spectrum trans_pdf = exp(-sigma_t * t_hit);
    Spectrum transmittance = exp(-sigma_t * t_hit);
    Spectrum Le = make_zero_spectrum();
    if(is_light(scene.shapes[(*vertex_).shape_id])) {
        Le = emission(*vertex_, -ray.dir, scene);
    }
    return Le * transmittance / trans_pdf;
}

inline int update_medium_id(const PathVertex vertex, const Ray ray, const int current_medium_id) {
    if(vertex.interior_medium_id != vertex.exterior_medium_id) {
        if (dot(ray.dir, vertex.geometric_normal) > 0) {
            return vertex.exterior_medium_id;
        } 
        return vertex.interior_medium_id;
    }

    return current_medium_id;
}
// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    //disable ray differentials
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    //?? when will current_medium not exist (except for configuration error)
    // bool current_medium_exist = true; // assume the configuration is correct...
    // Medium current_medium = scene.media[scene.camera.medium_id];
    int current_medium_id = scene.camera.medium_id;

    Spectrum radiance = make_zero_spectrum();
    Spectrum current_path_throughput = make_const_spectrum(1);
    int bounces = 0;
    
    // bool debug = false;
    while (true) {
        bool scatter = false;
        Real t_hit = infinity<Real>();
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        if(vertex_) {
            t_hit = distance((*vertex_).position, ray.org);
            // if(debug) {
            //     assert(true);
            // }
        }
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);
        
        //* homogenous, sigma_a & sigma_t can be derived from any points in the medium...
        if (current_medium_id != -1) {
            Medium current_medium = scene.media[current_medium_id];
            Spectrum sigma_a = get_sigma_a(current_medium, make_zero_spectrum());
            Spectrum sigma_s = get_sigma_s(current_medium, make_zero_spectrum());
            Spectrum sigma_t = sigma_a + sigma_s;

            Real u = next_pcg32_real<Real>(rng);
            Spectrum t = - log(1 - u) / sigma_t;

            if(t.x < t_hit) {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                // if(t_hit == infinity<Real>() && !debug) {
                    // std::cout << "====================" << std::endl;
                    // debug = true; 
                // }
                ray.org = ray.org + t * ray.dir;
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                ray.org = (*vertex_).position;
            }
            // ray.org = ray.org + t * ray.dir;
        }
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && vertex_) {
            //* reach a surface, include emission 
            //? medium still cannot emit?
            if(is_light(scene.shapes[(*vertex_).shape_id])) {
                radiance += current_path_throughput * emission(*vertex_, -ray.dir, scene);
                // if (debug) {
                //     std::cout << "radiance " << radiance << std::endl;
                // }
            }
        }

        //* add one more vertex in path
        if (!scatter && vertex_) {
            if ((*vertex_).material_id == -1) {
                //* index-matching interface, skip it
                current_medium_id = update_medium_id(*vertex_, ray, current_medium_id);
                bounces++;
                continue;
            }
        }

        //?? go through one more time to get a brighter image??
        if (bounces >= scene.options.max_depth  - 1 && scene.options.max_depth != -1) {
            break;
        }

        if (scatter && (current_medium_id != -1)) {
            Medium current_medium = scene.media[current_medium_id];
            Spectrum sigma_s = get_sigma_s(current_medium, make_zero_spectrum());
            PhaseFunction phase_function = get_phase_function(current_medium);

            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, phase_rnd_param_uv);
            Vector3f next_dir;
            if(next_dir_) {
                next_dir = *next_dir_;
                current_path_throughput *= eval(phase_function, -ray.dir, next_dir) / pdf_sample_phase(phase_function,  -ray.dir, next_dir) * sigma_s;
                ray.dir = next_dir;
                // if (debug) {
                    // std::cout << "current_path_throughput " << current_path_throughput << std::endl;
                // }
            } else {
                // fail to sample next_dir
                assert(true);
                break;
            }
        } else {
            //hit a surface
            break;
        }

            // Russian roulette heuristics
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;

    }
    // if (debug) {
        // std::cout << "return radiance " << radiance << std::endl;
    // }
    return radiance;
}

Spectrum get_sigma_t(const Scene& scene, int current_medium_id) {
    Medium current_medium = scene.media[current_medium_id];
    Spectrum sigma_a = get_sigma_a(current_medium, make_zero_spectrum());
    Spectrum sigma_s = get_sigma_s(current_medium, make_zero_spectrum());
    return sigma_a + sigma_s;
}

Spectrum next_event_estimation(const Scene& scene, Vector3 p, const Ray& ray, const int current_medium_id, 
                            pcg32_state& rng, int& bounces, RayDifferential& raydiff, const int x, const int y) {
    Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light &light = scene.lights[light_id];
    //? update last position that issue next event estimation
    PointAndNormal p_prime =
        sample_point_on_light(light, p, light_uv, shape_w, scene);
    Vector3 ori_p = p;
    Vector3 dir_light = normalize(p_prime.position - ori_p);
    
    //* compute transmittance to light, skip through index-matching shapes
    Spectrum T_light = make_const_spectrum(1);
    int shadow_medium_id = current_medium_id;
    int shadow_bounces = 0;
    Spectrum p_trans_dir = make_const_spectrum(1); // for multiple importance sampling
    while(true) {
        Vector3 dir_light = normalize(p_prime.position - p);
        Ray shadow_ray{p, dir_light, get_shadow_epsilon(scene),
                        (1 - get_shadow_epsilon(scene)) *
                        distance(p_prime.position, p)};
        Real next_t = distance(p, p_prime.position);
        if(debug(x, y)) {
            std::cout << "p" << p << std::endl;
            std::cout << "prime" << p_prime.position << std::endl;
            std::cout << "next_t =====1======" << next_t << std::endl;
        }
        std::optional<PathVertex> isect = intersect(scene, shadow_ray, raydiff);
        if (isect) {
            if(debug(x, y) && is_light(scene.shapes[((*isect).shape_id)])) {
                std::cout << "hit light " << std::endl;
            }
            next_t = distance(p, (*isect).position);
            if(debug(x, y)) {
                std::cout << "isect" << (*isect).position << std::endl;
                std::cout << "shape id" << (*isect).shape_id << std::endl;
                std::cout << "next_t =====2======" << next_t << std::endl;
            }
        }

        // if(!is_light(scene.shapes[(*isect).shape_id])) {
        //     std::cout << "material" << (*isect).material_id << std::endl;
        // }

        if(shadow_medium_id != -1) {
            // account for transmittance to next_t
            Spectrum sigma_t = get_sigma_t(scene, shadow_medium_id);
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
            if(debug(x, y)) {
                std::cout << "next_t " << next_t << std::endl;
                std::cout << "sigma_t " << sigma_t << std::endl;
                std::cout << "exp " << exp(-sigma_t * next_t) << std::endl;
                std::cout << "T_light " << T_light << std::endl;
                std::cout << "p_trans_dir" << p_trans_dir << std::endl;

            }
        }

        if(!isect) {
            break;
        } else {
           

            if((*isect).material_id >= 0) {
                // we're block
                return make_zero_spectrum();
            }
                // otherwise, we're on index-matching surface
            shadow_bounces++;
            if(scene.options.max_depth != -1 && (bounces + shadow_bounces + 1) >= scene.options.max_depth) {
                //reach the max no. of vertices, return 0;
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium_id(*isect, shadow_ray, shadow_medium_id);
            p = p + next_t * shadow_ray.dir;
        }
    }
    //?? need check current_medium_id
    if (T_light.x > 0 && (current_medium_id != -1)) { // still homogeneous
        Real G = max(-dot(dir_light, p_prime.normal), Real(0)) /
            distance_squared(p_prime.position, ori_p);
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
        PhaseFunction phase_function = get_phase_function(scene.media[current_medium_id]);
        Spectrum f = eval(phase_function, -ray.dir, dir_light);
        Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, ori_p, scene);
        Spectrum contrib = T_light * G * f * L / pdf_nee;
        
        // multiple importance sampling
        Spectrum pdf_phase = pdf_sample_phase(phase_function, -ray.dir, dir_light) * G * p_trans_dir;
        // Spectrum w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        Spectrum w = make_const_spectrum(0);
        if(debug(x, y)) {
            std::cout << "contrib " << w * contrib << std::endl;
            std::cout << "G " << G << std::endl;
            std::cout << "f " << f << std::endl;
            std::cout << "L " << L << std::endl;
            std::cout << "pdf_nee " << pdf_nee << std::endl;
        }
        return w * contrib;
    }         
    return make_zero_spectrum(); 
}
// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting

Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    if(debug(x, y))
        std::cout << "========= new round =============" << std::endl;
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    //disable ray differentials
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    int current_medium_id = scene.camera.medium_id;

    Spectrum radiance = make_zero_spectrum();
    Spectrum current_path_throughput = fromRGB(Vector3{1, 1, 1});
    int bounces = 0;
    // the pdf of the latest phase function sampling
    Real dir_pdf = 0; 
    // the last position p that can issue a next event estimation
    Vector3 nee_p_cache; 
    // the product PDF of transmittance sampling going through several index-matching surfaces
    // from the last phase funcion sampling
    Spectrum multi_trans_pdf = make_const_spectrum(1); 
    // a flag to indicate that the light path has never scattered
    bool never_scatter = true;
    
    while (true) {
        if(debug(x, y))
            std::cout << "========= new for loop =============" << std::endl;
        bool scatter = false;
        Real t_hit = infinity<Real>();
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        if(vertex_) {
            if(debug(x, y))
                std::cout << "========= ray hit something =============" << std::endl;
            t_hit = distance((*vertex_).position, ray.org);
        }
        Spectrum transmittance = make_const_spectrum(1);
        Spectrum trans_pdf = make_const_spectrum(1);
        
        if (current_medium_id == -1) {
            if(debug(x,y))
                std::cout << "========= into vacumn =============" << std::endl;
        }
        if (current_medium_id != -1) {
            Medium current_medium = scene.media[current_medium_id];
            Spectrum sigma_a = get_sigma_a(current_medium, make_zero_spectrum());
            Spectrum sigma_s = get_sigma_s(current_medium, make_zero_spectrum());
            Spectrum sigma_t = sigma_a + sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Spectrum t = - log(1 - u) / sigma_t;

            if(t.x < t_hit) {
                if(debug(x, y))
                    std::cout << "========= begin to scatter =============" << std::endl;

                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                // if(debug(x,y)) {
                //     std::cout << "ray.org" << ray.org << std::endl;
                // }
                ray.org = ray.org + t * ray.dir;
                // if(debug(x,y)) {
                //     std::cout << "t" << t << std::endl;
                //     std::cout << "before hit" << ray.org << std::endl;
                // }
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                // ray.org = ray.org + t_hit * ray.dir;
                ray.org = (*vertex_).position;
                if(debug(x,y)) {
                    std::cout << "======== hit surface =========" << std::endl;
                }
            }
            multi_trans_pdf *= trans_pdf;
            if(debug(x, y)) {
                std::cout << "multiple trans pdf " << multi_trans_pdf << std::endl;
            }
            // ray.org = ray.org + t * ray.dir;
        }
        auto origin_current_path_throughput = current_path_throughput;
        current_path_throughput *= (transmittance / trans_pdf);
        if(debug(x,y)) {
            std::cout << "current_path_throughput " << current_path_throughput << std::endl;
            std::cout << "original current_path_throughput " << origin_current_path_throughput << std::endl;
            std::cout << "transmittance " << transmittance << std::endl;
            std::cout << "trans_pdf " << trans_pdf << std::endl;
        }

        if (!scatter && vertex_) {
            //* reach a surface, include emission 
            //? medium still cannot emit?
            if(is_light(scene.shapes[(*vertex_).shape_id])) {
                if (never_scatter) {
                    radiance += current_path_throughput * emission(*vertex_, -ray.dir, scene);
                    if(debug(x,y)) {
                        std::cout << "======== never scatter before the hit ========= " << radiance << std::endl;
                    }
                } else {
                    //TODO: account for multiple sampling
                    // account for the next event estimation
                    PointAndNormal light_point;
                    light_point.position = (*vertex_).position;
                    light_point.normal = (*vertex_).geometric_normal;
                    int light_id = get_area_light_id(scene.shapes[(*vertex_).shape_id]);
                    const Light &light = scene.lights[light_id];
                    //?? pdf light
                    Real pdf_nee = pdf_point_on_light(light, light_point, nee_p_cache, scene) ;

                    // Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, light_point, nee_p_cache, scene) ;
                    // Vector3 dir_light = normalize(light_point.position - nee_p_cache);
                    Real G = max(-dot(ray.dir, light_point.normal), Real(0)) /
                        distance_squared(light_point.position, nee_p_cache);
                    Spectrum dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    // phase sampling weight
                    // Spectrum w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    Spectrum w = make_const_spectrum(1.);

                    radiance += current_path_throughput * emission(*vertex_, -ray.dir, scene) * w;
                    if(debug(x, y)) {
                        std::cout << "======== radiance in phase sampling =========" << radiance << std::endl;
                        std::cout << "======== weight =========" << w << std::endl;
                        std::cout << "======== pdf_phase =========" << dir_pdf_ << std::endl;
                        std::cout << "======== pdf_nee =========" << pdf_nee << std::endl;
                        std::cout << "======== dir_pdf =========" << dir_pdf << std::endl;
                        std::cout << "======== multi_trans_pdf =========" << multi_trans_pdf << std::endl;
                        std::cout << "======== G =========" << G << std::endl;
                        std::cout << "======== LE =========" << emission(*vertex_, -ray.dir, scene) << std::endl;
                    }
                }
            }
        }

        //* add one more vertex in path
        if (!scatter && vertex_) {
            if(debug(x, y) && (*vertex_).material_id != -1) {
                std::cout << "non index matching interface" << std::endl;
            }
            if ((*vertex_).material_id == -1) {
                //* index-matching interface, skip it
                auto origin_medium_id = current_medium_id;
                current_medium_id = update_medium_id(*vertex_, ray, current_medium_id);
                bounces++;
                if(debug(x, y)) {
                    std::cout << (*vertex_).shape_id << " index matching interface" << std::endl;
                    std::cout << (*vertex_).shape_id << " prev medium: " << origin_medium_id << std::endl;
                    std::cout << (*vertex_).shape_id << " inter medium: " << (*vertex_).interior_medium_id << std::endl;
                    std::cout << (*vertex_).shape_id << " exter medium: " << (*vertex_).exterior_medium_id << std::endl;
                    std::cout << (*vertex_).shape_id << " cos: " << dot(ray.dir, (*vertex_).geometric_normal) << std::endl;
                    std::cout << (*vertex_).shape_id << " isect: " << (*vertex_).position << std::endl;


                    std::cout << (*vertex_).shape_id << " next medium: " << current_medium_id << std::endl;

                }
                //!!!!!!!!!!!!update the ray.org after skipping the index-matching medium
                //! it will affect the next event estimation contribution (the p is not correct)
                ray.org = (*vertex_).position;
                continue;
            }
        }

        if (bounces >= scene.options.max_depth  - 1 && scene.options.max_depth != -1) {
            break;
        }
        
        if (scatter && (current_medium_id != -1)) {
            Medium current_medium = scene.media[current_medium_id];
            Spectrum sigma_s = get_sigma_s(current_medium, make_zero_spectrum());
            PhaseFunction phase_function = get_phase_function(current_medium);
            //TODO: combine with next estimation estimation
            //*=========================================================
            //* first do the next estimation estimation
            //*=========================================================
            //??? multiple the rsult with the transmittance and sigma_s
            Vector3 p = ray.org; // current position
            nee_p_cache = p; 
            Spectrum nee = next_event_estimation(scene, p, ray, current_medium_id, rng, bounces, ray_diff, x, y);
            radiance += current_path_throughput * nee * sigma_s;
            if(debug(x, y)) { 
                std::cout << "nee" << nee << std::endl;
                std::cout << "curr path through" << current_path_throughput << std::endl;
                std::cout << "radiance in nee" << radiance << std::endl;
            }

            //*=========================================================
            //* then do the phase sampling
            //*=========================================================
            never_scatter = false;
            Vector2 phase_rnd_param_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, phase_rnd_param_uv);
            Vector3f next_dir;
            if(next_dir_) {
                next_dir = *next_dir_;
                //?? update the dir_pdf 
                dir_pdf = pdf_sample_phase(phase_function,  -ray.dir, next_dir);
                current_path_throughput *= eval(phase_function, -ray.dir, next_dir) / dir_pdf * sigma_s;
                ray.dir = next_dir;
                if(debug(x,y)) {
                    std::cout << "after phase sample" << std::endl;
                    std::cout << "curr path through" << current_path_throughput << std::endl;
                    std::cout << "eval phase function" << eval(phase_function, -ray.dir, next_dir) << std::endl;
                    std::cout << "dir_pdf" << dir_pdf << std::endl;
                    std::cout << "reset the multi_trans_pdf!!" << std::endl;
                    

                }
                //!!!!!!!!!!!!!!!! reset the multi_trans_pdf after doing phase sampling
                //! otherwise the multi_trans_pdf will getting larger and larger (the trans_pdf is larger than 1 if sigma_s is too big)
                //! it only needs to record the product of trans_pdf from the last phase sample 
                //! in other words, when phase sample happen (also next event estimation), reset the multi_trans_pdf
                multi_trans_pdf = make_const_spectrum(1);
            } else {
                // fail to sample next_dir
                assert(true);
                break;
            }
        } else {
            //hit a surface
            break;
        }

            // Russian roulette heuristics
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max(current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // Terminate the path
                break;
            } else {
                current_path_throughput /= rr_prob;
            }
        }
        bounces++;

    }
    if (debug(x,y)) {
        std::cout << "return radiance " << radiance << std::endl;
    }
    // if(radiance.x > 0.2 && x < 255 && y < 255 && x > 200 && y > 200) {
    //     std::cout << x << ' ' << y << std::endl;
    // }
    // if(debug(x, y)) return make_const_spectrum(1000);
    // else {
    //     return make_const_spectrum(0);
    // }
    return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
