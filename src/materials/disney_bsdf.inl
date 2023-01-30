#include "../microfacet.h"
#include <iostream>

Spectrum eval_metal(const DisneyMetal &bsdf, const Vector3 &dir_in, const Vector3 &dir_out,
                    const PathVertex &vertex, const TexturePool &texture_pool, const TransportDirection &dir,
                    Real eta, Real specular_tint, Real metallic, Real specular) {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Vector3 half = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_half = dot(frame.n, half);
    //? do diffuse also need this?
    if (n_dot_out <= 0 || n_dot_half <= 0) {
        return make_zero_spectrum();
    } 
  
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color_ = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real R0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
    Spectrum c_tint = luminance(base_color_)  > 0 ? base_color_ / luminance(base_color_) : make_const_spectrum(1.);
    Spectrum k_s = (1 - specular_tint) + specular_tint * c_tint;
    Spectrum base_color = specular * R0 * (1 - metallic) * k_s + metallic * base_color_;
    
    // fresnel reflection term
    Real half_dot_out = abs(dot(half, dir_out));
    Spectrum fresnel = base_color + (1. - base_color) * pow(1. - half_dot_out, 5);

    //normal distribution funciton term
    Real ndf = anisotropic_GGX(roughness, anisotropic, to_local(frame, half));

    // geometry term
    Real geometry = smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_in)) * 
                    smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_out));
    
    return fresnel * ndf * geometry / (4 * n_dot_in);
}

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    // if(inside) {
    //     return eval(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
    //                       dir_in, dir_out, vertex, texture_pool, dir); 
    // }

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = bsdf.eta;


    Spectrum f_diffuse = eval(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface},
                                dir_in, dir_out, vertex, texture_pool, dir);
                            
    Spectrum f_sheen = eval(DisneySheen{bsdf.base_color, bsdf.sheen_tint},
                            dir_in, dir_out, vertex, texture_pool, dir); 

    Spectrum f_metal = eval_metal(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic},
                            dir_in, dir_out, vertex, texture_pool, dir, eta, specular_tint, metallic, specular);  

    Spectrum f_clear_coat = eval(DisneyClearcoat{bsdf.clearcoat_gloss},
                            dir_in, dir_out, vertex, texture_pool, dir);

    Spectrum f_glass = eval(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, dir_out, vertex, texture_pool, dir);  
    
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    Real w_diffuse = inside ? Real(0) : (1 - specular_transmission) * (1 - metallic);
    Real w_sheen = inside ? Real(0) : (1 - metallic) * sheen;
    Real w_metal = inside ? Real(0) : (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = inside ? Real(0) : 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
                                                            
    if (reflect) {
        return w_diffuse * f_diffuse + w_sheen * f_sheen + w_metal * f_metal + 
               w_clearcoat * f_clear_coat + w_glass * f_glass;
    } else {
        return w_glass * f_glass;
    }                                                   
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    //??? how to use reflect
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real pdf_diffuse = pdf_sample_bsdf(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface},
                          dir_in, dir_out, vertex, texture_pool, dir);

    Real pdf_metal = pdf_sample_bsdf(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic},
                          dir_in, dir_out, vertex, texture_pool, dir);  
    // assert(pdf_metal <= 10000000000 && pdf_metal >= 0.01);

    Real pdf_clear_coat = pdf_sample_bsdf(DisneyClearcoat{bsdf.clearcoat_gloss},
                          dir_in, dir_out, vertex, texture_pool, dir);    
    // assert(pdf_clear_coat <= 10000000000 && pdf_clear_coat >= 0.01);

    Real pdf_glass = pdf_sample_bsdf(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, dir_out, vertex, texture_pool, dir);
    // assert(pdf_glass <= 10000000000 && pdf_glass >= 0.01);

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real normalize = (1 - specular_transmission) * (1 - metallic) +  
                     (1 - specular_transmission * (1 - metallic)) + 
                     0.25 * clearcoat + 
                     (1 - metallic) * specular_transmission;
    assert(normalize > 0);
    
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
   
    Real w_diffuse = inside ? 0 : (1 - specular_transmission) * (1 - metallic) / normalize;
    Real w_metal = inside ? 0 : (1 - specular_transmission * (1 - metallic)) / normalize;
    Real w_glass = inside ? 1 : (1 - metallic) * specular_transmission / normalize;
    Real w_clearcoat = inside ? 0 : 0.25 * clearcoat / normalize;
    
    if(reflect) {
        return w_diffuse * pdf_diffuse + 
            w_metal * pdf_metal + 
            w_glass * pdf_glass +
            w_clearcoat * pdf_clear_coat;
    } else {
        return pdf_glass;
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);


    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    if(inside) {
        return sample_bsdf(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }
    Real normalize = (1 - specular_transmission) * (1 - metallic) + 
                     (1 - specular_transmission * (1 - metallic)) + 
                     (1 - metallic) * specular_transmission + 
                     0.25 * clearcoat;
    assert(normalize > 0);                  
    Real p_diffuse = (1 - specular_transmission) * (1 - metallic) / normalize;
    Real p_metal = (1 - specular_transmission * (1 - metallic)) / normalize;
    Real p_glass = (1 - metallic) * specular_transmission / normalize;
    Real p_clearcoat = 0.25 * clearcoat / normalize;

    if (rnd_param_w <= p_diffuse) {
        //reuse random number
        Real new_rnd = rnd_param_w / p_diffuse;
        return sample_bsdf(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface},
                          dir_in, vertex, texture_pool, rnd_param_uv, new_rnd);
    } {
        if (rnd_param_w - p_diffuse <= p_metal) {
            Real new_rnd = (rnd_param_w  - p_diffuse) / p_metal;
            return sample_bsdf(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic},
                          dir_in, vertex, texture_pool, rnd_param_uv, new_rnd);
        } else {
            if (rnd_param_w - p_diffuse - p_metal <= p_clearcoat) {
                Real new_rnd = (rnd_param_w  - p_diffuse - p_metal) / p_clearcoat;
                return sample_bsdf(DisneyClearcoat{bsdf.clearcoat_gloss},
                          dir_in, vertex, texture_pool, rnd_param_uv, new_rnd);
            } else {
                Real new_rnd = (rnd_param_w  - p_diffuse - p_metal - p_clearcoat) / p_glass;
                return sample_bsdf(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, vertex, texture_pool, rnd_param_uv, new_rnd);
            }
        }
    }
    assert(true);
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}