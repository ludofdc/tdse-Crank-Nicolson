#include "physics.hpp"
#include <vector>
#include <cmath>
#include <complex>
#include <numeric>


std::vector<double> build_grid(const Params& p) {
    double dx = p.dx;
    std::vector<double> x(p.N);
    double x_min = -p.L/2;

    for(int i = 0; i < p.N; i++) {
        x[i] = x_min + dx*i;
    }

    return x;
}

std::vector<double> set_potential(const std::vector<double>& x, 
                                  const Params& p) {
    std::vector<double> v(p.N);

    for (int i = 0; i < p.N; i++) {
        if (x[i] > p.x_bar && x[i] < (p.x_bar + p.w_bar)) 
            v[i] = p.V0;
        else 
            v[i] = 0.0;
    }

    return v;
}

std::vector<std::complex<double>> init_gauss_psi(const std::vector<double>& x,
                                                  const Params& p) {
    using cd = std::complex<double>;
    const cd I(0.0, 1.0);

    std::vector<cd> psi(p.N);

    // Build un-normalised Gaussian wave packet:
    //   psi(x) = exp( -(x-x0)^2 / (4*sigma^2) ) * exp(i*k0*x)
    // The Gaussian envelope localises the packet around x0 with width sigma.
    // The plane wave exp(i*k0*x) gives the packet its mean momentum hbar*k0.
    for (int i = 0; i < p.N; i++) {
        double envelope = std::exp(-std::pow(x[i] - p.x0, 2) / (4.0 * p.sigma * p.sigma));
        cd    phase     = std::exp(I * p.k0 * x[i]);
        psi[i] = envelope * phase;
    }

    // Normalise so that integral |psi|^2 dx = 1  (trapezoidal rule, uniform dx)
    double norm2 = 0.0;
    for (int i = 0; i < p.N; i++)
        norm2 += std::norm(psi[i]);   // std::norm returns |z|^2
    norm2 *= p.dx;

    double inv_sqrt_norm = 1.0 / std::sqrt(norm2);
    for (int i = 0; i < p.N; i++)
        psi[i] *= inv_sqrt_norm;

    return psi;
}

// ---------------------------------------------------------------------------
// Crank-Nicolson matrices
//
// The discretised TDSE in atomic units (hbar=1, m=1) is:
//
//   i * (psi^{n+1} - psi^n) / dt = H * (psi^{n+1} + psi^n) / 2
//
// Rearranging:
//
//   A * psi^{n+1} = B * psi^n
//
// where
//   A = I + (i*dt/2) * H_disc
//   B = I - (i*dt/2) * H_disc
//
// The discretised kinetic operator is the standard second-order finite
// difference: H_kin[i] = -1/(2*dx^2) * (psi[i-1] - 2*psi[i] + psi[i+1])
//
// Defining  r = dt / (4*dx^2)  the matrix elements become:
//
//   A: main diagonal      = 1 + i*(dt/2)*( 1/dx^2 + V[i] )
//      off-diagonals      = -i*r
//
//   B: main diagonal      = 1 - i*(dt/2)*( 1/dx^2 + V[i] )
//      off-diagonals      = +i*r
//
// Boundary conditions: Dirichlet  psi[0] = psi[N-1] = 0 (absorbing walls).
// The first and last rows are therefore trivial (identity rows with rhs=0).
// ---------------------------------------------------------------------------

TridiagMatrix build_diagA(const std::vector<double>& V, const Params& p) {
    using cd = std::complex<double>;
    const cd I(0.0, 1.0);

    int N = p.N;
    TridiagMatrix A;
    A.lower.resize(N, 0.0);
    A.main .resize(N, 0.0);
    A.upper.resize(N, 0.0);

    double r = p.dt / (4.0 * p.dx * p.dx);   // off-diagonal coefficient

    // Boundary rows: psi = 0 enforced by identity
    A.main[0]     = 1.0;
    A.main[N - 1] = 1.0;

    for (int i = 1; i < N - 1; i++) {
        double diag_kin = 1.0 / (p.dx * p.dx);          // kinetic part of diagonal
        A.main [i] = 1.0 + I * (p.dt / 2.0) * (diag_kin + V[i]);
        A.lower[i] = -I * r;
        A.upper[i] = -I * r;
    }

    return A;
}

TridiagMatrix build_diagB(const std::vector<double>& V, const Params& p) {
    using cd = std::complex<double>;
    const cd I(0.0, 1.0);

    int N = p.N;
    TridiagMatrix B;
    B.lower.resize(N, 0.0);
    B.main .resize(N, 0.0);
    B.upper.resize(N, 0.0);

    double r = p.dt / (4.0 * p.dx * p.dx);

    // Boundary rows: rhs = 0 (psi stays zero at walls)
    B.main[0]     = 1.0;
    B.main[N - 1] = 1.0;

    for (int i = 1; i < N - 1; i++) {
        double diag_kin = 1.0 / (p.dx * p.dx);
        B.main [i] = 1.0 - I * (p.dt / 2.0) * (diag_kin + V[i]);
        B.lower[i] = +I * r;
        B.upper[i] = +I * r;
    }

    return B;
}
