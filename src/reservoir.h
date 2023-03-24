#pragma once

#include "point_and_normal.h"
#include "intersection.h"

struct Reservoir {
    Real M;
    Real w_sum;
    Real W;
    PointAndNormal y;
    int light_id;
    PathVertex ref_vertex;
    Vector3 prev_dir_view;
};

inline Reservoir init_reservoir() {
    Reservoir r;
    r.M = 0;
    r.w_sum = 0;
    r.W = 0;
    r.y = PointAndNormal{
        Vector3{0, 0, 0},
        Vector3{0, 0, 1},
    };
    r.light_id = 0;
    r.prev_dir_view = Vector3{0, 0, 0};
    r.ref_vertex = PathVertex{};
    r.ref_vertex.position = Vector3{0, 0, 0};
    r.ref_vertex.geometric_normal = Vector3{0, 0, 1};

    return r;
}

inline bool update_reservoir(Reservoir &rsv, int light_id, PointAndNormal x, Real w, const Real &pdf_w) {
    rsv.w_sum += w;
    rsv.M += 1;
    if (rsv.w_sum == 0) {
        return false;
    }
    if (pdf_w < (w / rsv.w_sum)) {
        rsv.y = x;
        rsv.light_id = light_id;
        return true;
    }
    return false;
}

inline bool update_reservoir_debug(Reservoir &rsv, int light_id, PointAndNormal x, Real w, const Real &pdf_w) {
    std::cout << "before updating reservoir" << pdf_w << std::endl;
    std::cout << "w_sum: " << rsv.w_sum << std::endl;
    std::cout << "adding: " << w << std::endl;
    rsv.w_sum += w;
    rsv.M += 1;
    std::cout << "updating reservoir ...." << std::endl;
    std::cout << "pdf_w" << pdf_w << std::endl;
    std::cout << "rsv.w_sum" << rsv.w_sum << std::endl;
    std::cout << "rsv.M " << rsv.M << std::endl;


    if (rsv.w_sum == 0) {
        return false;
    }
    std::cout << "w / rsv.w_sum" << w / rsv.w_sum << std::endl;

    if (pdf_w < (w / rsv.w_sum)) {
        std::cout << "choosen!" << pdf_w << std::endl;
        rsv.y = x;
        rsv.light_id = light_id;
        return true;
    }
    return false;
}