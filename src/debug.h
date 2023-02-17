#pragma once
#include <iostream>
#include "lajolla.h"
#include "vector.h"

bool debug(int x, int y) {
    return false;
    // return x == 255 && y == 255; // middle
    // return x == 205 && y == 201;
    // return x == 204 && y == 203;
    // return x == 249 && y == 243;
}

void debug_print(Real val) {
    std::cout << val << std::endl;
}

void debug_print(std::string msg) {
    std::cout << msg << std::endl;
}

void debug_print(std::string msg, Real val) {
    std::cout << msg << ": "<< val << std::endl;
}

void debug_print(std::string msg, Vector3 val) {
    std::cout << msg << ": "<< val << std::endl;
}