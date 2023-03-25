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
    // return x == 167 && y == 149;
    // return x == 392 && y == 468;
    // return x == 228 && y == 490; 
    // return x == 239 && y == 482;
    // return x == 202 && y == 335;
    // return x == 211 && y == 218;
    // return x == 264 && y == 333;
    // return x == 330 && y == 99;
    // return x == 149 && y == 443;
    // return x == 118 && y == 465;

    // return x == 213 && y == 442;
    return x == 184 && y == 455;
}

inline bool ignore_first_bounce() {
    return false;
}