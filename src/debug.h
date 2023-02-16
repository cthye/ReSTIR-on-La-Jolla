#pragma once
#include <iostream>
#include "lajolla.h"
#include "vector.h"

bool debug(int x, int y) {
    return false;
    // return x == 255 && y == 255; // middle
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