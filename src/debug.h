#pragma once
#include <iostream>
#include "lajolla.h"
#include "vector.h"

inline bool debug(int x, int y) {
    // return false;
    // return x == 255 && y == 255; // middle
    // return x == 205 && y == 201;
    // return x == 204 && y == 203;
    // return x == 249 && y == 243;
    // return x == 323 && y == 229;
    // return x == 146 && y == 246;
    // return x == 216 && y == 306;
    return x == 167 && y == 149;
}

inline bool ignore_first_bounce() {
    return true;
}