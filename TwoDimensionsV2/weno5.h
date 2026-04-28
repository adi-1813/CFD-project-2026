#pragma once
#include <vector>
#include "geometry.h"   // needs ConsvState, Flux

using namespace std;

// ── WENO epsilon to prevent division by zero in smoothness weights ───────────
static constexpr double WENO_EPS = 1.0e-6;

// ============================================================================
// WENO5 nonlinear weight functions  (Liu, Osher & Chan 1994; Jiang & Shu 1996)
// ============================================================================
// Both functions reconstruct a numerical flux at the right edge x_{i+1/2}
// from five consecutive split-flux point values {v1,...,v5}.
// Interior mapping: v1=f[i-2], v2=f[i-1], v3=f[i], v4=f[i+1], v5=f[i+2].

// Left-biased (upwind) reconstruction — used for F^+ (positive wave speeds).
// Ideal weights d = {1/10, 6/10, 3/10}.
double wcL5(double v1, double v2, double v3, double v4, double v5);

// Right-biased (downwind) reconstruction — used for F^- (negative wave speeds).
// Ideal weights d = {3/10, 6/10, 1/10}  (mirror of wcL5).
double wcR5(double v1, double v2, double v3, double v4, double v5);

// ============================================================================
// Scalar WENO5 sweep over a 1-D split-flux array
// ============================================================================
// Both functions operate on a single scalar component of the split flux
// (e.g. fp1[], fm1[], ...) and write the reconstructed interface flux into hL/hR.
//
// Ghost values required by the 5-point stencil at the boundaries:
//   gLm1 = f[-1],  gLm2 = f[-2]   (two left ghosts)
//   gRp0 = f[N],   gRp1 = f[N+1]  (two right ghosts)
//
// wenoL_scalar writes h^+[i] for all i = 0..N-1.
// wenoR_scalar writes h^-[i] for all i = 0..N-1.
// These are later differenced as:
//   dF^+/dx|_i ≈ (h^+[i] - h^+[i-1]) / dx   (backward)
//   dF^-/dx|_i ≈ (h^-[i+1] - h^-[i]) / dx   (forward)

void wenoL_scalar(int N, const vector<double>& f,
                  vector<double>& hL,
                  double gLm1, double gLm2,   // f[-1], f[-2]
                  double gRp0, double gRp1);  // f[N],  f[N+1]

void wenoR_scalar(int N, const vector<double>& f,
                  vector<double>& hR,
                  double gLm1, double gLm2,   // f[-1], f[-2]
                  double gRp0, double gRp1);  // f[N],  f[N+1]


// ── 1D overload (original) ───────────────────────────────────────────────────
void compute_rhs(const vector<ConsvState>& U, int N, double dx, double gamma,
                 vector<ConsvState>& rhs_out);

// ============================================================================
// Full spatial residual  L(U) = -dF/dx
// ============================================================================
// Computes the FDM-WENO5 spatial residual using Lax-Friedrichs flux splitting.
// Ghost states for wall BCs are built internally from U[1], U[2] (left wall)
// and U[N-2], U[N-3] (right wall) with rhou sign-flipped.
// Output: rhs_out[i] for all i = 0..N-1.
// ============================================================================
// Full x-direction spatial residual  Lx(U) = -dF/dx   [2D overload]
// ============================================================================
// Processes every row j independently.  Left/right wall ghosts flip rhou.
// rhs_out must be pre-sized to [Ny][Nx].
void compute_rhs(const vector<vector<ConsvState>>& U,
                 const Grid2D& g, double gamma,
                 vector<vector<ConsvState>>& rhs_out);

// ============================================================================
// Full y-direction spatial residual  Ly(U) = -dG/dy   [2D only]
// ============================================================================
// Mirrors compute_rhs (x-sweep) but sweeps each column i=0..Nx-1 using the
// y-direction flux G(U) from compute_flux_y().
// Ghost states at top/bottom walls reflect rhov (zero normal velocity).
// rhs_out is ADDED to (+=), not overwritten — caller accumulates Lx+Ly.
void compute_rhs_y(const vector<vector<ConsvState>>& U,
                   const Grid2D& g, double gamma,
                   vector<vector<ConsvState>>& rhs_out);