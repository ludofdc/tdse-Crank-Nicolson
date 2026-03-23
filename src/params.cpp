// params.cpp
// Implementation of parameter reading and printing functions
// Atomic units: hbar=1, m=1

#include "params.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

Params read_params(const std::string& filename) {

    Params p;
    std::ifstream file(filename);

    if (!file.is_open())
        throw std::runtime_error("Cannot open input file: " + filename);

    std::string line, key;
    double value;

    while (std::getline(file, line)) {

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        if (!(iss >> key >> value)) continue;

        if      (key == "L")      p.L      = value;
        else if (key == "dx")     p.dx     = value;
        else if (key == "dt")     p.dt     = value;
        else if (key == "T_tot")  p.T_tot  = value;
        else if (key == "x0")     p.x0     = value;
        else if (key == "sigma")  p.sigma  = value;
        else if (key == "k0")     p.k0     = value;
        else if (key == "x_bar")  p.x_bar  = value;
        else if (key == "w_bar")  p.w_bar  = value;
        else if (key == "V0")     p.V0     = value;
    }

    // Compute derived quantities
    p.N = static_cast<int>(std::round(p.L / p.dx));
    p.M = static_cast<int>(std::round(p.T_tot / p.dt));

    return p;
}

void print_params(const Params& p) {

    std::cout << "========================================\n";
    std::cout << "  TDSE Crank-Nicolson solver\n";
    std::cout << "  Atomic units: hbar=1, m=1\n";
    std::cout << "========================================\n";

    std::cout << "\n  Domain and grid\n";
    std::cout << "    L      = " << p.L     << " bohr\n";
    std::cout << "    dx     = " << p.dx    << " bohr\n";
    std::cout << "    dt     = " << p.dt    << " hbar/Eh\n";
    std::cout << "    T_tot  = " << p.T_tot << " hbar/Eh\n";
    std::cout << "    N      = " << p.N     << " points\n";
    std::cout << "    M      = " << p.M     << " steps\n";

    std::cout << "\n  Gaussian wave packet\n";
    std::cout << "    x0     = " << p.x0    << " bohr\n";
    std::cout << "    sigma  = " << p.sigma << " bohr\n";
    std::cout << "    k0     = " << p.k0    << " bohr^-1\n";

    std::cout << "\n  Potential barrier\n";
    std::cout << "    x_bar  = " << p.x_bar << " bohr\n";
    std::cout << "    w_bar  = " << p.w_bar << " bohr\n";
    std::cout << "    V0     = " << p.V0    << " Eh\n";

    std::cout << "========================================\n\n";
}