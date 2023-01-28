#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    // bool reflect = dot(vertex.geometric_normal, dir_in) *
    //                dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    if(inside) {
        return eval(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, dir_out, vertex, texture_pool, dir); 
    }

    Spectrum f_diffuse = eval(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface},
                                dir_in, dir_out, vertex, texture_pool, dir);
                            
    Spectrum f_sheen = eval(DisneySheen{bsdf.base_color, bsdf.sheen_tint},
                            dir_in, dir_out, vertex, texture_pool, dir); 

    Spectrum f_metal = eval(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic},
                            dir_in, dir_out, vertex, texture_pool, dir);  

    Spectrum f_clear_coat = eval(DisneyClearcoat{bsdf.clearcoat_gloss},
                            dir_in, dir_out, vertex, texture_pool, dir);

    Spectrum f_glass = eval(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, dir_out, vertex, texture_pool, dir);  

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
                                                            
    return (1 - specular_transmission) * (1 - metallic) * f_diffuse + 
           (1 - metallic) * sheen * f_sheen + 
           (1 - specular_transmission * (1 - metallic)) * f_metal + 
           0.25 * clearcoat * f_clear_coat + 
           (1 - metallic) * specular_transmission * f_glass; 
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    //??? how to use reflect
    // bool reflect = dot(vertex.geometric_normal, dir_in) *
    //                dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    bool inside = dot(vertex.geometric_normal, dir_in) <= 0;
    if(inside) {
        return pdf_sample_bsdf(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                            dir_in, dir_out, vertex, texture_pool, dir); 
    }

    Real pdf_diffuse = pdf_sample_bsdf(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface},
                          dir_in, dir_out, vertex, texture_pool, dir);

    Real pdf_metal = pdf_sample_bsdf(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic},
                          dir_in, dir_out, vertex, texture_pool, dir);  

    Real pdf_clear_coat = pdf_sample_bsdf(DisneyClearcoat{bsdf.clearcoat_gloss},
                          dir_in, dir_out, vertex, texture_pool, dir);    

    Real pdf_glass = pdf_sample_bsdf(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, dir_out, vertex, texture_pool, dir);  

    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real normalize = (1 - specular_transmission) * (1 - metallic) +  
                     (1 - specular_transmission * (1 - metallic)) + 
                     0.25 * clearcoat + 
                     (1 - metallic) * specular_transmission;

    return (1 - specular_transmission) * (1 - metallic) / normalize * pdf_diffuse + 
        (1 - specular_transmission * (1 - metallic)) / normalize * pdf_metal + 
        0.25 * clearcoat / normalize * pdf_clear_coat + 
        (1 - metallic) / normalize * specular_transmission * pdf_glass; 
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

    Real p_diffuse = (1 - specular_transmission) * (1 - metallic) / normalize;
    Real p_metal = (1 - specular_transmission * (1 - metallic)) / normalize;
    Real p_glass = (1 - metallic) * specular_transmission / normalize;
    Real p_clearcoat = 0.25 * clearcoat / normalize;
    //reuse random number
    if (rnd_param_w <= p_diffuse) {
        // in metallic: clearcoat & metal
        return sample_bsdf(DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface},
                          dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    } {
        Real new_rnd = (rnd_param_w  - p_diffuse) / (1 - p_diffuse);
        if (new_rnd <= p_metal) {
            return sample_bsdf(DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic},
                          dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
        } else {
            Real new_rnd = (rnd_param_w  - p_diffuse -  p_metal) / (1 - p_diffuse - p_metal);
            if (new_rnd <= p_glass) {
            return sample_bsdf(DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta},
                          dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
            } else {
                return sample_bsdf(DisneyClearcoat{bsdf.clearcoat_gloss},
                          dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
            }
        }
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
