#include "../microfacet.h"
#include <iostream>

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    //! pow clamps in schlick?
    //* fetch parameters
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    //* calculate weight for each lobes
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic); 
    Real w_sheen = (1 - metallic) * sheen ;
    Real w_metal = 1 - specular_transmission * (1 - metallic);
    Real w_clearcoat = 0.25 * clearcoat ;
    Real w_glass = (1 - metallic) * specular_transmission;

    //* some operations on vectors
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    Vector3 half_vector;
    if(reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = normalize(dir_in + dir_out * eta);
    }

    //* make sure half_dot_n is >= 0
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    //!
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);

    //* ============================== 
    //* calculate f_glass
    //* ============================== 
    Real F_glass = fresnel_dielectric(h_dot_in, eta);
    Real D_glass = anisotropic_GGX(roughness, anisotropic, to_local(frame, half_vector));
    Real G_glass = smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_in)) * 
                    smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_out));

    Spectrum f_glass;
    if (reflect) {
        f_glass = base_color * F_glass * D_glass * G_glass / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        assert(fabs(dot(frame.n, dir_in)) * pow(h_dot_in + eta * h_dot_out, 2));
        f_glass = sqrt(base_color) * (1 - F_glass) * D_glass * G_glass * fabs(h_dot_in * h_dot_out) / 
           (fabs(dot(frame.n, dir_in)) * pow(h_dot_in + eta * h_dot_out, 2));
    }
    
    //* ============================== 
    //* calculate f_diffuse
    //* ============================== 
    double f_d_90 = 1. / 2. + 2. * roughness * pow(abs(h_dot_out), 2);
    double f_d_in = 1. + (f_d_90 - 1.) * pow(1. - abs(n_dot_in), 5);
    double f_d_out = 1. + (f_d_90 - 1.) * pow(1. - abs(n_dot_out), 5);
    Spectrum f_base_diffuse = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) / c_PI * 
             f_d_in * f_d_out * abs(n_dot_out);

    // subsurface
    double f_ss_90 = roughness * pow(abs(h_dot_out), 2);
    double f_ss_in = 1. + (f_ss_90 - 1.) * pow(1. - abs(n_dot_in), 5);
    double f_ss_out = 1. + (f_ss_90 - 1.) * pow(1. - abs(n_dot_out), 5);
    Spectrum f_subsurface = 1.25 * eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) / c_PI * 
             (f_ss_in * f_ss_out * (1. / (abs(n_dot_in) + abs(n_dot_out)) - 0.5) + 0.5) * abs(n_dot_out);

    // combine
    Spectrum f_diffuse = (1. - subsurface) * f_base_diffuse + subsurface * f_subsurface;

    //* ============================== 
    //* calculate f_metal
    //* ============================== 
    Real R0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    Spectrum c_tint = luminance(base_color)  > 0 ? base_color / luminance(base_color) : make_const_spectrum(1.);
    Spectrum k_s = (1 - specular_tint) + specular_tint * c_tint;
    Spectrum base_color_0 = specular * R0 * (1 - metallic) * k_s + metallic * base_color;
    
    Spectrum F_metal = base_color_0 + (1. - base_color_0) * pow(1. - h_dot_out, 5);
    Real D_metal = anisotropic_GGX(roughness, anisotropic, to_local(frame, half_vector));
    Real G_metal = smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_in)) * 
                    smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_out));
    
    Spectrum f_metal = F_metal * D_metal * G_metal / (4 * abs(n_dot_in));

    //* ============================== 
    //* calculate f_sheen
    //* ==============================
    Real L = luminance(base_color);
    Spectrum tint = make_const_spectrum(1.);
    if (L > 0) tint = base_color / luminance(base_color);
    Spectrum C = (1 - sheen_tint) + tint * sheen_tint;

    Spectrum f_sheen = C * pow((1 - fabs(dot(half_vector, dir_out))), 5) * fabs(dot(frame.n, dir_out));

    //* ============================== 
    //* calculate f_clearcoat
    //* ==============================
    Real ior = 1.5;
    Real F_clearcoat = R0 + (1 - R0) * pow(1 - abs(h_dot_out), 5);
    Real D_clearcoat = clearcoat_GGX(clearcoat_gloss, to_local(frame, half_vector));
    Real G_clearcoat = smith_masking_gtr4(clearcoat_gloss, to_local(frame, dir_in)) *
                    smith_masking_gtr4(clearcoat_gloss,  to_local(frame, dir_out));

    Spectrum f_clearcoat = make_const_spectrum(F_clearcoat * D_clearcoat * G_clearcoat / (4 * abs(n_dot_in)));

    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    if(inside) {
        f_diffuse = make_const_spectrum(0);
        f_metal = make_const_spectrum(0);
        f_sheen = make_const_spectrum(0);
        f_clearcoat = make_const_spectrum(0);
    }

    if(reflect) {
        return  w_diffuse * f_diffuse + 
                w_metal * f_metal + 
                w_sheen * f_sheen + 
                w_clearcoat * f_clearcoat + 
                w_glass * f_glass;
    } else {
        return w_glass * f_glass;
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    //* fetch parameters
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    //* calculate and normalize weight for each lobes
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic); 
    Real w_metal = 1 - specular_transmission * (1 - metallic);
    Real w_clearcoat = 0.25 * clearcoat ;
    Real w_glass = (1 - metallic) * specular_transmission;
    Real factor = w_diffuse + w_metal + w_clearcoat + w_glass;
    w_diffuse /= factor; 
    w_metal /= factor; 
    w_clearcoat /= factor; 
    w_glass /= factor; 

    //* some operations on vectors
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    Vector3 half_vector;
    if(reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        half_vector = normalize(dir_in + dir_out * eta);
    }

    //* make sure half_dot_n is >= 0
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    //!
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);

    //* ============================== 
    //* calculate pdf_glass
    //* ============================== 
    Real F_glass = fresnel_dielectric(h_dot_in, eta);
    Real D_glass = anisotropic_GGX(roughness, anisotropic, to_local(frame, half_vector));
    Real G_in_glass = smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_in));

    Real pdf_glass;
    if (reflect) {
        pdf_glass = (F_glass * D_glass * G_in_glass) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        //? anything need to change here?
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        pdf_glass = (1 - F_glass) * D_glass * G_in_glass * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }

    //* ============================== 
    //* calculate pdf_diffuse
    //* ============================== 
    Real pdf_diffuse = fmax(dot(frame.n, dir_out), Real(0)) / c_PI;

    //* ============================== 
    //* calculate pdf_metal
    //* ============================== 
    Real D_metal = anisotropic_GGX(roughness, anisotropic, to_local(frame, half_vector));
    Real G_in_metal = smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_in));
    Real pdf_metal = D_metal * G_in_metal / (4 * abs(n_dot_in));
    
    //* ============================== 
    //* calculate pdf_clearcoat
    //* ==============================
    Real D_clearcoat = clearcoat_GGX(clearcoat_gloss, to_local(frame, half_vector));
    Real pdf_clearcoat = D_clearcoat * abs(dot(frame.n, half_vector)) / (4 * abs(h_dot_out));
        
    assert(factor > 0);
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    if(inside) {
        w_diffuse = 0;
        w_metal = 0;
        w_glass = 1;
        w_clearcoat = 0;
    }

    Real pdf;
    if(reflect) {
        pdf  = w_diffuse * pdf_diffuse + 
               w_metal * pdf_metal + 
               w_glass * pdf_glass +
               w_clearcoat * pdf_clearcoat;
    } else {
        pdf = pdf_glass;
    }

    assert(pdf >= 0);
    return pdf;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    //* fetch parameters
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    //* calculate and normalize weight for each lobes
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic); 
    Real w_metal = 1 - specular_transmission * (1 - metallic);
    Real w_clearcoat = 0.25 * clearcoat ;
    Real w_glass = (1 - metallic) * specular_transmission;
    Real factor = w_diffuse + w_metal + w_clearcoat + w_glass;
    w_diffuse /= factor; 
    w_metal /= factor; 
    w_clearcoat /= factor; 
    w_glass /= factor; 

    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    if(inside) {
        //* ============================== 
        //*  sample glass
        //* ============================== 
        Real alpha = roughness * roughness;
        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

        Vector3 half_vector = to_world(frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_in = dot(half_vector, dir_in);
        Real F = fresnel_dielectric(h_dot_in, eta);

        if (rnd_param_w <= F) {
            // Reflection
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            // set eta to 0 since we are not transmitting
            return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        } else {
            // Refraction
            // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            // (note that our eta is eta2 / eta1, and l = -dir_in)
            Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
            if (h_dot_out_sq <= 0) {
                // Total internal reflection
                // This shouldn't really happen, as F will be 1 in this case.
                return {};
            }
            // flip half_vector if needed
            if (h_dot_in < 0) {
                half_vector = -half_vector;
            }
            //? anychanges
            Real h_dot_out= sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    }

    if (rnd_param_w <= w_diffuse) {
        //reuse random number
        //* ============================== 
        //*  sample diffuse
        //* ============================== 
        return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            Real(0) /* eta */, roughness/* roughness */};
    } {
        if (rnd_param_w - w_diffuse <= w_metal) {
            //* ============================== 
            //*  sample metal
            //* ============================== 
            Real alpha = roughness * roughness;
            Vector3 dir_in_local = to_local(frame, dir_in);
            Vector3 local_micro_normal = sample_visible_normals(dir_in_local, alpha, rnd_param_uv);
            Vector3 half_vector = to_world(frame, local_micro_normal);
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{
                reflected,
                Real(0) /* eta */, roughness /* roughness */
            };
        } else {
            if (rnd_param_w - w_diffuse - w_metal <= w_clearcoat) {
                //* ============================== 
                //*  sample clearcoat
                //* ============================== 
                Real alpha_g = (1. - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
                Real h_elevation = acos(sqrt((1 - pow(alpha_g * alpha_g, 1 - rnd_param_uv.x))/(1 - alpha_g * alpha_g)));
                Real h_azimuth = 2 * c_PI * rnd_param_uv.y;
                Real h_local_x = sin(h_elevation) * cos(h_azimuth);
                Real h_local_y = sin(h_elevation) * sin(h_azimuth);
                Real h_local_z = cos(h_elevation);
                Vector3 half_vector = to_world(frame, Vector3(h_local_x, h_local_y, h_local_z));
                Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
                return BSDFSampleRecord{
                    reflected,
                    Real(0) /* eta */, Real(0) /* roughness */
                };
            } else {
                Real new_rnd = (rnd_param_w  - w_diffuse - w_metal - w_clearcoat) / w_glass;
                //* ============================== 
                //*  sample glass
                //* ============================== 
                Real alpha = roughness * roughness;
                Vector3 local_dir_in = to_local(frame, dir_in);
                Vector3 local_micro_normal =
                    sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

                Vector3 half_vector = to_world(frame, local_micro_normal);
                // Flip half-vector if it's below surface
                if (dot(half_vector, frame.n) < 0) {
                    half_vector = -half_vector;
                }
                Real h_dot_in = dot(half_vector, dir_in);
                Real F = fresnel_dielectric(h_dot_in, eta);
                if (new_rnd <= F) {
                    // Reflection
                    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
                    // set eta to 0 since we are not transmitting
                    return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
                } else {
                    // Refraction
                    // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
                    // (note that our eta is eta2 / eta1, and l = -dir_in)
                    Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
                    if (h_dot_out_sq <= 0) {
                        // Total internal reflection
                        // This shouldn't really happen, as F will be 1 in this case.
                        return {};
                    }
                    // flip half_vector if needed
                    if (h_dot_in < 0) {
                        half_vector = -half_vector;
                    }
                    //? anychanges
                    Real h_dot_out= sqrt(h_dot_out_sq);
                    Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
                    return BSDFSampleRecord{refracted, eta, roughness};
                }
            }
        }
    }
    assert(true);
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}