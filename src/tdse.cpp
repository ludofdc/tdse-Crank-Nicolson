// tdse.cpp
// Main driver - TDSE Crank-Nicolson solver
// Atomic units: hbar = 1, m = 1
//
// Algorithm outline
// -----------------
// At each time step we must solve the linear system
//
//   A * psi^{n+1} = rhs
//
// where  rhs = B * psi^n  and A, B are complex tridiagonal matrices built
// once and reused.  The tridiagonal system is solved with the Thomas
// algorithm (a.k.a. tridiagonal matrix algorithm, TDMA), which is the
// Gaussian-elimination specialisation for tridiagonals: O(N) in time and
// memory.  This is equivalent to what LAPACK's zgtsv does, implemented here
// explicitly in complex arithmetic to avoid the external dependency.
//
// Output
// ------
// At every snapshot step the solver writes a CSV file
//   output/psi_XXXXXX.csv
// with columns:  x , |psi|^2 , Re(psi) , Im(psi) , V(x)
// A summary file  output/norm.csv  records  t , norm  at every snapshot.

#include "params.hpp"
#include "physics.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>

// ---------------------------------------------------------------------------
// Thomas algorithm for a complex tridiagonal system  M * x = rhs
//
// The system has N equations.  The matrix M is given by three vectors:
//   lo[i]   lower diagonal  (lo[0] unused)
//   md[i]   main  diagonal
//   up[i]   upper diagonal  (up[N-1] unused)
//
// The solution overwrites 'rhs' in place.
// The algorithm is numerically stable for diagonally dominant matrices,
// which the Crank-Nicolson matrix always satisfies.
// ---------------------------------------------------------------------------
static void thomas_solve(const std::vector<std::complex<double>>& lo,
                               std::vector<std::complex<double>>  md,   // copy: modified in place
                         const std::vector<std::complex<double>>& up,
                               std::vector<std::complex<double>>& rhs)
{
    int N = static_cast<int>(rhs.size());

    // ----- Forward sweep: eliminate lower diagonal -----
    // For each row i starting from 1, subtract a multiple of row i-1
    // so that the lower diagonal entry becomes zero.
    for (int i = 1; i < N; i++) {
        std::complex<double> factor = lo[i] / md[i - 1];  // elimination factor
        md [i]   -= factor * up[i - 1];                   // update main diagonal
        rhs[i]   -= factor * rhs[i - 1];                  // update right-hand side
    }

    // ----- Back substitution -----
    rhs[N - 1] /= md[N - 1];
    for (int i = N - 2; i >= 0; i--)
        rhs[i] = (rhs[i] - up[i] * rhs[i + 1]) / md[i];
}

// ---------------------------------------------------------------------------
// Apply the tridiagonal matrix M to the vector psi: result = M * psi
// Used to compute rhs = B * psi^n at each step.
// ---------------------------------------------------------------------------
static std::vector<std::complex<double>>
mat_vec(const TridiagMatrix& M, const std::vector<std::complex<double>>& psi)
{
    int N = static_cast<int>(psi.size());
    std::vector<std::complex<double>> result(N, 0.0);

    result[0] = M.main[0] * psi[0] + M.upper[0] * psi[1];
    for (int i = 1; i < N - 1; i++)
        result[i] = M.lower[i] * psi[i - 1]
                  + M.main [i] * psi[i]
                  + M.upper[i] * psi[i + 1];
    result[N - 1] = M.lower[N - 1] * psi[N - 2] + M.main[N - 1] * psi[N - 1];

    return result;
}

// ---------------------------------------------------------------------------
// Compute norm^2 = integral |psi|^2 dx  (rectangle rule, uniform dx)
// ---------------------------------------------------------------------------
static double compute_norm2(const std::vector<std::complex<double>>& psi, double dx)
{
    double s = 0.0;
    for (const auto& z : psi)
        s += std::norm(z);   // std::norm(z) = |z|^2
    return s * dx;
}

// ---------------------------------------------------------------------------
// Write one snapshot CSV:  x , |psi|^2 , Re(psi) , Im(psi) , V(x)
// ---------------------------------------------------------------------------
static void write_snapshot(const std::string& filename,
                            const std::vector<double>&                x,
                            const std::vector<std::complex<double>>&  psi,
                            const std::vector<double>&                V)
{
    std::ofstream f(filename);
    f << "x,prob,re_psi,im_psi,V\n";
    f << std::scientific << std::setprecision(10);
    for (int i = 0; i < static_cast<int>(x.size()); i++) {
        f << x[i]            << ","
          << std::norm(psi[i]) << ","
          << psi[i].real()   << ","
          << psi[i].imag()   << ","
          << V[i]            << "\n";
    }
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {

    // ---- 1. Read parameters and print summary ----
    Params p = read_params("input/input.txt");
    print_params(p);

    // ---- 2. Build spatial grid and potential ----
    std::vector<double> x = build_grid(p);
    std::vector<double> V = set_potential(x, p);

    // ---- 3. Initialise wave packet ----
    std::vector<std::complex<double>> psi = init_gauss_psi(x, p);

    // ---- 4. Build Crank-Nicolson matrices (built once, reused every step) ----
    TridiagMatrix A = build_diagA(V, p);   // LHS:  A * psi^{n+1} = rhs
    TridiagMatrix B = build_diagB(V, p);   // RHS:  rhs = B * psi^n

    // ---- 5. Decide snapshot frequency ----
    // We want at most 200 snapshots spread uniformly over M steps.
    int snap_every = std::max(1, p.M / 200);

    // Open norm log file
    std::ofstream norm_log("output/norm.csv");
    norm_log << "t,norm2\n";
    norm_log << std::scientific << std::setprecision(10);

    // Write initial snapshot (step 0)
    {
        std::ostringstream fname;
        fname << "output/psi_" << std::setfill('0') << std::setw(6) << 0 << ".csv";
        write_snapshot(fname.str(), x, psi, V);

        double n2 = compute_norm2(psi, p.dx);
        norm_log << 0.0 << "," << n2 << "\n";
        std::cout << "step " << std::setw(6) << 0
                  << "  t = " << std::fixed << std::setprecision(4) << 0.0
                  << "  norm^2 = " << std::scientific << std::setprecision(8) << n2 << "\n";
    }

    // ---- 6. Time evolution loop ----
    int snap_index = 1;   // snapshot counter (for filename numbering)

    for (int n = 1; n <= p.M; n++) {

        // (a) Compute rhs = B * psi^n
        std::vector<std::complex<double>> rhs = mat_vec(B, psi);

        // (b) Solve A * psi^{n+1} = rhs  with Thomas algorithm.
        //     We pass copies of A's diagonals because Thomas modifies them.
        std::vector<std::complex<double>> lo_copy = A.lower;
        std::vector<std::complex<double>> md_copy = A.main;
        std::vector<std::complex<double>> up_copy = A.upper;
        thomas_solve(lo_copy, md_copy, up_copy, rhs);
        // rhs now holds psi^{n+1}

        psi = rhs;

        // (c) Enforce Dirichlet boundary conditions (should already be ~0)
        psi[0]       = 0.0;
        psi[p.N - 1] = 0.0;

        // (d) Write snapshot and norm if it is a snapshot step
        if (n % snap_every == 0 || n == p.M) {
            double t  = n * p.dt;
            double n2 = compute_norm2(psi, p.dx);

            std::ostringstream fname;
            fname << "output/psi_" << std::setfill('0') << std::setw(6) << snap_index << ".csv";
            write_snapshot(fname.str(), x, psi, V);
            snap_index++;

            norm_log << t << "," << n2 << "\n";
            std::cout << "step " << std::setw(6) << n
                      << "  t = " << std::fixed << std::setprecision(4) << t
                      << "  norm^2 = " << std::scientific << std::setprecision(8) << n2 << "\n";
        }
    }

    std::cout << "\nDone. Snapshots written to output/\n";
    return 0;
}
