#pragma once
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

using namespace std;

// ============================================================================
// Conservative state vector
// ============================================================================
// NOTE: rhov added for 2D. In the current 1D x-sweep solver it is carried
// but not yet evolved. Wire up the y-sweep in compute_rhs when ready.
struct ConsvState {
    double rho;   // U1: density
    double rhou;  // U2: x-momentum  (= rho * u)
    double rhov;  // U3: y-momentum  (= rho * v)
    double rhoE;  // U4: total energy (= rho * E)
};

// ============================================================================
// Physical flux vectors
// ============================================================================
// F(U): x-direction flux.   G(U): y-direction flux (compute_flux_y below).
//
// NOTE: The F3/G2 cross-momentum terms below use the average form
// (rhou*v + rhov*u)/2.  The correct 2D x-flux for the rhov equation is
// simply rho*u*v = (rhou/rho)*rhov.  Fix when wiring up the full 2D residual.
struct Flux {
    double F1;  // rho * u
    double F2;  // rho*u^2 + p
    double F3;  // rho*u*v  (cross term)
    double F4;  // u*(rho*E + p)
};

// ============================================================================
// Thermodynamic helpers
// ============================================================================

static inline double pressure(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double v = U.rhov / U.rho;
    double E = U.rhoE / U.rho;
    return (gamma - 1.0) * U.rho * (E - 0.5*(u*u) - 0.5*(v*v));
}

static inline double sound_speed(const ConsvState& U, double gamma)
{
    return sqrt(gamma * pressure(U, gamma) / U.rho);
}

// max(|u|,|v|) + a: correct per-direction CFL bound for 2D Euler.
// Do NOT use sqrt(u^2+v^2)+a — see geometry.h notes.
static inline double max_wave_speed(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double v = U.rhov / U.rho;
    double a = sound_speed(U, gamma);
    return max(fabs(u), fabs(v)) + a;
}

// x-direction physical flux F(U)
static inline Flux compute_flux(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double v = U.rhov / U.rho;
    double p = pressure(U, gamma);
    return { U.rhou,
             U.rhou * u + p,
             (U.rhov * u + U.rhou * v) / 2.0,
             u * (U.rhoE + p) };
}

// y-direction physical flux G(U) — used when y-sweep is added
static inline Flux compute_flux_y(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double v = U.rhov / U.rho;
    double p = pressure(U, gamma);
    return { U.rhov,
             (U.rhou * v + U.rhov * u) / 2.0,
             U.rhov * v + p,
             v * (U.rhoE + p) };
}

// ============================================================================
// Grid — 2D uniform Cartesian
// ============================================================================
// U[j][i] is the state at node (x[i], y[j]).
// For 1D usage: set Ny=1, y_min=0, y_max=0.  dy will be 0 and the y-arrays
// have a single entry y[0]=0.  The 2D loops degenerate to a single j=0 row.
//
// Grid1D is aliased to Grid2D for backward compatibility with 1D main.cpp.

struct Grid2D {
    int    Nx, Ny;
    double x_min, x_max, dx;
    double y_min, y_max, dy;
    vector<double> x;   // x[i] = x_min + i*dx,  i = 0..Nx-1
    vector<double> y;   // y[j] = y_min + j*dy,  j = 0..Ny-1

    Grid2D(int Nx_, double x_min_, double x_max_,
           int Ny_ = 1, double y_min_ = 0.0, double y_max_ = 0.0)
        : Nx(Nx_), Ny(Ny_), x_min(x_min_), x_max(x_max_),
          y_min(y_min_), y_max(y_max_)
    {
        dx = (x_max - x_min) / (Nx - 1);
        dy = (Ny > 1) ? (y_max - y_min) / (Ny - 1) : 0.0;
        x.resize(Nx);  for (int i = 0; i < Nx; i++) x[i] = x_min + i*dx;
        y.resize(Ny);  for (int j = 0; j < Ny; j++) y[j] = y_min + j*dy;
    }
};

using Grid1D = Grid2D;   // backward-compatibility alias

// ============================================================================
// WallCurve — an immersed curved wall defined by  y = f(x)
// ============================================================================
// The wall occupies x in [x0, x1] and its height at any x is given by
// the user-supplied function f(x).  The wall is sampled at n_samples points
// on construction so that height_at(), is_inside(), and normal_at() are all
// O(log n) or O(1) at query time without calling f() again.
//
// Convention: the FLUID side is ABOVE the curve  (y > f(x)).
// If the fluid is below the curve, negate f or flip the normal in normal_at().
//
// ── How to define a wall ────────────────────────────────────────────────────
//
// Pass any callable as wall_func.  A lambda is the most concise syntax:
//
//   // Flat bottom wall at y = 0 (the shock-tube floor used below):
//   WallCurve bottom(-0.5, 0.5, 400, [](double x){ return 0.0; });
//
//   // DUMMY — sinusoidal bump of amplitude A and wavelength equal to the
//   // tube length:
//   // double A = 0.05;
//   // WallCurve bump(-0.5, 0.5, 400,
//   //     [A](double x){ return A * sin(M_PI * (x + 0.5)); });
//   // The normal tilts ±arctan(A*pi) from vertical at the flanks and is
//   // perfectly vertical at the crest/trough.
//
//   // DUMMY — straight wedge at angle theta (oblique-shock test case):
//   // double theta = 15.0 * M_PI / 180.0;
//   // WallCurve wedge(-0.5, 0.0, 400,
//   //     [theta](double x){ return tan(theta) * (x + 0.5); });
//   // The wall starts at (x=-0.5, y=0) and rises to the right.
//
//   // DUMMY — circular arc bump (radius R, centred below the domain):
//   // double R = 2.0, xc = 0.0, yc = -R + 0.1;
//   // WallCurve arc(-0.3, 0.3, 400,
//   //     [R,xc,yc](double x){ return yc + sqrt(R*R - (x-xc)*(x-xc)); });
//
// The WallCurve only covers [x0,x1].  Outside this range is_inside() always
// returns false, so you can stack multiple WallCurve objects to describe a
// piecewise wall.

struct WallCurve {
    double x0, x1;              // x-extent of this wall segment
    int    n;                   // number of sample points
    vector<double> sx, sy;      // sampled (x, y=f(x)) pairs, uniform in x
    function<double(double)> f; // original wall function (kept for reference)

    WallCurve(double x0_, double x1_, int n_samples,
              function<double(double)> wall_func)
        : x0(x0_), x1(x1_), n(n_samples), f(wall_func)
    {
        sx.resize(n);
        sy.resize(n);
        double h = (x1 - x0) / (n - 1);
        for (int k = 0; k < n; k++) {
            sx[k] = x0 + k * h;
            sy[k] = f(sx[k]);
        }
    }

    // Wall height y = f(x) at any x, by linear interpolation on stored samples.
    // Clamps to the endpoint value outside [x0, x1].
    double height_at(double x) const
    {
        if (x <= sx[0])   return sy[0];
        if (x >= sx[n-1]) return sy[n-1];
        // Binary search for the bracketing interval [lo, lo+1]
        int lo = 0, hi = n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            (sx[mid] <= x ? lo : hi) = mid;
        }
        double t = (x - sx[lo]) / (sx[lo+1] - sx[lo]);
        return sy[lo] + t * (sy[lo+1] - sy[lo]);
    }

    // Returns true if the point (px, py) is inside the solid (below the wall).
    // Points outside [x0, x1] are never inside.
    bool is_inside(double px, double py) const
    {
        if (px < x0 || px > x1) return false;
        return py < height_at(px);
    }
};

// ============================================================================
// normal_at — outward unit normal to a WallCurve at a given x
// ============================================================================
// The wall tangent at x is  t = (1, f'(x))  (unnormalised).
// The outward normal (pointing into the fluid, away from the solid) is obtained
// by rotating t by -90°:  n_raw = (f'(x), -1)  ... wait, let's be careful:
//
//   t = (1, f')    →   rotate +90° (CCW) →   n_raw = (-f', 1)
//
// This points upward (into the fluid above the wall) for any slope.
// The unit normal is  n = n_raw / |n_raw|.
//
// f'(x) is estimated by a central finite difference on the stored sample
// array; forward/backward differences are used at the endpoints.
//
// Returns: {nx, ny} as a std::pair, guaranteed unit length.
pair<double,double> normal_at(const WallCurve& wall, double x);

// ============================================================================
// bilinear_interpolate — interpolate ConsvState at arbitrary (px, py)
// ============================================================================
// Given a 2D state array U[j][i] on grid g, finds the four surrounding nodes:
//
//   (i,   j  )  (i+1, j  )
//   (i,   j+1)  (i+1, j+1)
//
// and returns:
//
//   U(px,py) = (1-tx)*(1-ty)*U[j][i]     +  tx*(1-ty)*U[j][i+1]
//            + (1-tx)*   ty *U[j+1][i]   +  tx*   ty *U[j+1][i+1]
//
// with  tx = (px - x[i]) / dx,  ty = (py - y[j]) / dy.
//
// Index clamping is applied so that query points outside the grid use the
// nearest interior cell — the caller does not need to bounds-check.
ConsvState bilinear_interpolate(const vector<vector<ConsvState>>& U,
                                const Grid2D& g,
                                double px, double py);

// ============================================================================
// apply_immersed_wall_bc — ghost-node reflection BC for a curved wall
// ============================================================================
// Iterates over every grid node (i,j).  For nodes that lie INSIDE the wall
// (wall.is_inside(x[i], y[j]) == true), the following steps are applied:
//
//   Step 1 — Wall-normal at the ghost node's x-position:
//               n = normal_at(wall, x[i])      (unit outward normal)
//
//   Step 2 — Signed distance from the ghost node G=(x[i], y[j]) to the wall
//            surface, measured along n.  For a bottom wall:
//               d = y[j] - wall.height_at(x[i])    (negative: G is below wall)
//            More precisely, the foot F on the wall surface is:
//               F = G - d * n   (the nearest point on the wall to G,
//                                approximated along the normal direction)
//
//   Step 3 — Mirror point M of G reflected through F:
//               M = G - 2*d * n   =  G + 2*|d| * n
//            M is on the fluid side (above the wall).
//
//   Step 4 — Bilinear-interpolate the flow state at M:
//               U_m = bilinear_interpolate(U, g, M.x, M.y)
//
//   Step 5 — Reflect the momentum across the wall normal to enforce
//            no-penetration (u · n = 0 at the wall):
//               u_m  = U_m.rhou / U_m.rho
//               v_m  = U_m.rhov / U_m.rho
//               un   = u_m*nx + v_m*ny         (normal component)
//               u_ghost = (u_m - 2*un*nx)      (flip normal, keep tangential)
//               v_ghost = (v_m - 2*un*ny)
//
//   Step 6 — Set the ghost cell state:
//               rho_ghost  = rho_m             (density: same as mirror)
//               rhou_ghost = rho_m * u_ghost
//               rhov_ghost = rho_m * v_ghost
//               rhoE_ghost = rhoE_m            (total energy: same as mirror)
//
// This is the standard ghost-fluid / image-point reflection used in
// immersed boundary methods for inviscid compressible flow.
void apply_immersed_wall_bc(vector<vector<ConsvState>>& U,
                             const Grid2D& g,
                             const WallCurve& wall,
                             double gamma);

// ============================================================================
// Axis-aligned (flat) wall BCs — shock tube and simple box
// ============================================================================
// 2D version: zeros rhou at left/right walls (i=0, i=Nx-1) for all j,
//             zeros rhov at bottom/top walls (j=0, j=Ny-1) for all i.
void apply_wall_bc(vector<vector<ConsvState>>& U, const Grid2D& g);

// 1D overload kept for backward compatibility with single-row solver.
void apply_wall_bc(vector<ConsvState>& U, int N);

// ============================================================================
// Initial conditions — Sod shock tube
// ============================================================================
// Discontinuity at x=0.
//   Left  (x<0): rho=1.0,   u=0, v=0, p=1.0
//   Right (x>0): rho=0.125, u=0, v=0, p=0.1
void set_initial_conditions(vector<vector<ConsvState>>& U,
                             const Grid2D& g, double gamma);

// 1D overload
void set_initial_conditions(vector<ConsvState>& U,
                             const Grid2D& g, double gamma);

// ============================================================================
// Initial conditions — supersonic flow over a wedge
void set_ic_wedge(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma, double Mach, double theta);

void apply_inflow_bc(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma, double Mach, double theta);

void apply_outflow_bc(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma);

