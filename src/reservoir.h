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
};

inline Reservoir init_reservoir() {
    Reservoir r;
    r.M = 0;
    r.w_sum = 0;
    r.y = PointAndNormal{
        Vector3{0, 0, 0},
        Vector3{0, 0, 1},
    };
    r.light_id = 0;
    return r;
}

inline void update_reservoir(Reservoir &rsv, int light_id, PointAndNormal x, Real w, const Real &pdf_w) {
    rsv.w_sum += w;
    rsv.M += 1;
    if (rsv.w_sum == 0) {
        return;
    }
    if (pdf_w < (w / rsv.w_sum)) {
        rsv.y = x;
        rsv.light_id = light_id;
    }
}