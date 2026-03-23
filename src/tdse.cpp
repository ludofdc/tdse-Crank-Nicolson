// tdse.cpp
// Main file - TDSE Crank-Nicolson solver

#include "params.hpp"
#include <iostream>

int main() {

    Params p = read_params("input/input.txt");
    print_params(p);

    return 0;
}