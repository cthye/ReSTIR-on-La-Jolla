#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    
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

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // VNDF sampling
    Vector3 half = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_half = dot(frame.n, half);
    //? do diffuse also need this?
    if (n_dot_out <= 0 || n_dot_half <= 0) {
        return 0;
    } 

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real G = smith_masking_gtr3(roughness, anisotropic, to_local(frame, dir_in));
    Real D = anisotropic_GGX(roughness, anisotropic, to_local(frame, half));
    return D * G / (4 * n_dot_in);
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // VNDF sampling

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector2 alpha = aniso_alpha(roughness, anisotropic);
    Vector3 dir_in_local = to_local(frame, dir_in);
    Vector3 local_micro_normal = sample_visible_normals_with_aniso(dir_in_local, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
