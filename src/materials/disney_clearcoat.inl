#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real n_dot_in = fmax(dot(frame.n, dir_in), Real(0));
    Real n_dot_out = fmax(dot(frame.n, dir_out), Real(0));
    Real n_dot_half = fmax(dot(frame.n, half), Real(0));
    //? do diffuse also need this?
    if (n_dot_out <= 0 || n_dot_half <= 0) {
        return make_zero_spectrum();
    } 
    
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    
    Real ior = 1.5;
    Real R0 = (ior - 1) * (ior - 1) / ((ior + 1) * (ior + 1));
    Real half_dot_out = fmax(dot(half, dir_out), Real(0));
    Real fresnel = R0 + (1 - R0) * pow(1 - half_dot_out, 5);


    Real ndf = clearcoat_GGX(clearcoat_gloss, to_local(frame, half));
    
    Real geometry = smith_masking_gtr4(clearcoat_gloss, to_local(frame, dir_in)) *
                    smith_masking_gtr4(clearcoat_gloss,  to_local(frame, dir_out));

    return make_const_spectrum(fresnel * ndf * geometry / (4 * n_dot_in));
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Vector3 half = normalize(dir_in + dir_out);
    Real n_dot_out = fmax(dot(frame.n, dir_out), Real(0));
    Real n_dot_half = fmax(dot(frame.n, half), Real(0));
    //? do diffuse also need this?
    if (n_dot_out <= 0 || n_dot_half <= 0) {
        return 0;
    } 

    Real half_dot_out = fmax(dot(half, dir_out), Real(0));
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real D = clearcoat_GGX(clearcoat_gloss, to_local(frame, half));

    return D * n_dot_half / (4 * abs(half_dot_out));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
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
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1. - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001;
    Real h_elevation = acos(sqrt(1 - pow(alpha_g * alpha_g, 1 - rnd_param_uv.x))/(1 - alpha_g * alpha_g));
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
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
