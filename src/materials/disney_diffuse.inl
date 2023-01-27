Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Vector3 half = normalize((dir_in + dir_out));
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    double n_dot_in = abs(dot(frame.n, dir_in));
    double n_dot_out = abs(dot(frame.n, dir_out));
    double half_dot_out = abs(dot(half, dir_out));

    // base diffuse
    double f_d_90 = 1. / 2. + 2. * roughness * pow(half_dot_out, 2);
    double f_d_in = 1. + (f_d_90 - 1.) * pow(1. - n_dot_in, 5);
    double f_d_out = 1. + (f_d_90 - 1.) * pow(1. - n_dot_out, 5);
    Spectrum f_base_diffuse = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) / c_PI * 
             f_d_in * f_d_out * abs(n_dot_out);

    // subsurface
    double f_ss_90 = roughness * pow(half_dot_out, 2);
    double f_ss_in = 1. + (f_ss_90 - 1.) * pow(1. - n_dot_in, 5);
    double f_ss_out = 1. + (f_ss_90 - 1.) * pow(1. - n_dot_out, 5);
    Spectrum f_subsurface = 1.25 * eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) / c_PI * 
             (f_ss_in * f_ss_out * (1. / (n_dot_in + n_dot_out) - 0.5) + 0.5) * n_dot_out;

    // combine
    double subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    return (1. - subsurface) * f_base_diffuse + subsurface * f_subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // importance sample the cosine hemisphere domain
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // importance sample the cosine hemisphere domain
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness/* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
