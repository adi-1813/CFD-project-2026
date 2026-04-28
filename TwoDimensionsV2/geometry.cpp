#include "geometry.h"
#include <stdexcept>

// ============================================================================
// normal_at — outward unit normal to a WallCurve at position x
// ============================================================================
// Strategy:
//   Find the sample interval [lo, lo+1] containing x, then estimate f'(x)
//   with a central difference on the stored sample array.  At the endpoints
//   fall back to a one-sided (forward/backward) difference.
//
//   Tangent  t = (1, f')
//   Raw normal (rotate t by +90° CCW) = (-f', 1)   — points upward into fluid
//   Unit normal n = (-f', 1) / sqrt(f'^2 + 1)
pair<double,double> normal_at(const WallCurve& wall, double x)
{
    // Clamp x to the wall extent
    x = max(wall.sx[0], min(wall.sx[wall.n-1], x));

    // Locate the bracketing interval by binary search
    int lo = 0, hi = wall.n - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        (wall.sx[mid] <= x ? lo : hi) = mid;
    }

    // Estimate f'(x) by finite difference on stored samples
    double dfdx;
    if (lo == 0)
    {
        // Forward difference at left endpoint
        double h = wall.sx[1] - wall.sx[0];
        dfdx = (wall.sy[1] - wall.sy[0]) / h;
    }
    else if (lo == wall.n - 2 && x >= wall.sx[wall.n-1] - 1e-14)
    {
        // Backward difference at right endpoint
        int k = wall.n - 1;
        double h = wall.sx[k] - wall.sx[k-1];
        dfdx = (wall.sy[k] - wall.sy[k-1]) / h;
    }
    else
    {
        // Central difference using the neighbours of the bracket
        // Use lo-1 and lo+2 if they exist, otherwise fall back to lo and lo+1
        int left  = max(0,          lo - 1);
        int right = min(wall.n - 1, lo + 2);
        double h  = wall.sx[right] - wall.sx[left];
        dfdx = (wall.sy[right] - wall.sy[left]) / h;
    }

    // Outward unit normal: rotate tangent (1, dfdx) by +90° CCW → (-dfdx, 1)
    double len = sqrt(dfdx*dfdx + 1.0);
    return { -dfdx / len, 1.0 / len };
}

// ============================================================================
// bilinear_interpolate
// ============================================================================
// Finds the cell (i,j) such that x[i] <= px < x[i+1] and y[j] <= py < y[j+1],
// clamps to the valid range, and performs the standard bilinear blend.
//
// Index layout (standard row-major, y varies on outer index):
//   U[j  ][i  ]   U[j  ][i+1]
//   U[j+1][i  ]   U[j+1][i+1]
ConsvState bilinear_interpolate(const vector<vector<ConsvState>>& U,
                                const Grid2D& g,
                                double px, double py)
{
    // ── x index ──────────────────────────────────────────────────────────────
    // Fractional index along x:  fi = (px - x_min) / dx
    // Integer part gives the left cell edge, fractional part gives tx.
    int i = (int)((px - g.x_min) / g.dx);
    i = max(0, min(g.Nx - 2, i));          // clamp so i+1 is always valid
    double tx = (px - g.x[i]) / g.dx;
    tx = max(0.0, min(1.0, tx));

    // ── y index ──────────────────────────────────────────────────────────────
    int j = (g.Ny > 1) ? (int)((py - g.y_min) / g.dy) : 0;
    j = max(0, min(g.Ny - 2, j));
    double ty = (g.Ny > 1) ? (py - g.y[j]) / g.dy : 0.0;
    ty = max(0.0, min(1.0, ty));

    // ── Bilinear blend ────────────────────────────────────────────────────────
    // w00=(1-tx)(1-ty), w10=tx(1-ty), w01=(1-tx)ty, w11=tx*ty
    double w00 = (1.0 - tx) * (1.0 - ty);
    double w10 =        tx  * (1.0 - ty);
    double w01 = (1.0 - tx) *        ty;
    double w11 =        tx  *        ty;

    // For 1D grids (Ny==1), j+1 would be out of range — guard with min().
    int jp1 = min(j + 1, g.Ny - 1);

    ConsvState result;
    result.rho  = w00*U[j][i].rho   + w10*U[j][i+1].rho
                + w01*U[jp1][i].rho  + w11*U[jp1][i+1].rho;
    result.rhou = w00*U[j][i].rhou  + w10*U[j][i+1].rhou
                + w01*U[jp1][i].rhou + w11*U[jp1][i+1].rhou;
    result.rhov = w00*U[j][i].rhov  + w10*U[j][i+1].rhov
                + w01*U[jp1][i].rhov + w11*U[jp1][i+1].rhov;
    result.rhoE = w00*U[j][i].rhoE  + w10*U[j][i+1].rhoE
                + w01*U[jp1][i].rhoE + w11*U[jp1][i+1].rhoE;
    return result;
}

// ============================================================================
// apply_immersed_wall_bc
// ============================================================================
void apply_immersed_wall_bc(vector<vector<ConsvState>>& U,
                             const Grid2D& g,
                             const WallCurve& wall,
                             double gamma)
{
    #pragma omp parallel for collapse(2) schedule(dynamic, 16)
    for (int j = 0; j < g.Ny; j++)
    {
        for (int i = 0; i < g.Nx; i++)
        {
            double px = g.x[i];
            double py = g.y[j];

            // Only act on nodes that lie inside the solid wall
            if (!wall.is_inside(px, py)) continue;

            // ── Step 1: Outward unit normal at this x position ───────────────
            auto [nx, ny] = normal_at(wall, px);

            // ── Step 2: Signed distance from ghost node G to the wall surface
            // The wall surface at x=px has height  yw = wall.height_at(px).
            // The ghost node is at py < yw (inside the wall).
            // We approximate the signed distance along the normal as:
            //   d = (py - yw)   [negative, since the ghost is below the wall]
            // A more accurate formulation would project G onto the wall curve
            // by minimising |G - F(s)|^2, but for small d this is equivalent.
            double yw = wall.height_at(px);
            double d  = (py - yw)*ny;   // d < 0  (ghost is inside/below wall)
            //double dx = 0.0;
            //double dy = py - yw;
            //double d  = dx*nx + dy*ny;
            // ── Step 3: Mirror point M = G - 2*d*n ──────────────────────────
            // Since d < 0, -2*d > 0, so M is displaced in the +n direction,
            // i.e., into the fluid above the wall.
            double mx = px - 2.0 * d * nx;
            double my = py - 2.0 * d * ny;

            // Clamp M to the grid domain so bilinear_interpolate doesn't
            // extrapolate outside the grid.
            mx = max(g.x_min, min(g.x_max, mx));
            my = max(g.y_min, min(g.y_max, my));

            // ── Step 4: Interpolate fluid state at M ─────────────────────────
            ConsvState Um = bilinear_interpolate(U, g, mx, my);

            // ── Step 5: Reflect velocity across wall normal ──────────────────
            // u_m = velocity at mirror point
            double um = Um.rhou / Um.rho;
            double vm = Um.rhov / Um.rho;

            // Normal component of u_m  (scalar)
            double un = um * nx + vm * ny;

            // Reflected velocity: u_ghost = u_m - 2*(u_m · n)*n
            // This zeroes the normal component and keeps the tangential part.
            double u_ghost = um - 2.0 * un * nx;
            double v_ghost = vm - 2.0 * un * ny;

            // ── Step 6: Write ghost cell state ───────────────────────────────
            // Density and total energy are copied from the mirror (no jump).
            // Momentum is set from the reflected velocity.
            U[j][i].rho  = Um.rho;
            U[j][i].rhou = Um.rho * u_ghost;
            U[j][i].rhov = Um.rho * v_ghost;
            U[j][i].rhoE = Um.rhoE;
        }
    }
}

// ============================================================================
// apply_wall_bc — axis-aligned flat walls
// ============================================================================

// 2D version: left/right walls zero rhou; top/bottom walls zero rhov.
void apply_wall_bc(vector<vector<ConsvState>>& U, const Grid2D& g)
{
   // Left wall  (i=0):   zero x-momentum for all rows
   // Right wall (i=Nx-1): zero x-momentum for all rows
   // for (int j = 0; j < g.Ny; j++) {
   //     U[j][0].rhou      = 0.0;
   //     U[j][g.Nx-1].rhou = 0.0;
   // }

    // Bottom wall (j=0):    zero y-momentum for all columns
    // Top wall    (j=Ny-1): zero y-momentum for all columns
    if (g.Ny > 1) {
        for (int i = 0; i < g.Nx; i++) {
            U[0][i].rhov      = 0.0;
            U[g.Ny-1][i].rhov = 0.0;
        }
    }
}

// 1D overload: only left/right walls, only rhou.
void apply_wall_bc(vector<ConsvState>& U, int N)
{
    U[0].rhou     = 0.0;
    U[N-1].rhou   = 0.0;
}

// ============================================================================
// set_initial_conditions — Sod shock tube
// ============================================================================

void set_initial_conditions(vector<vector<ConsvState>>& U,
                             const Grid2D& g, double gamma)
{
    for (int j = 0; j < g.Ny; j++) {
        for (int i = 0; i < g.Nx; i++) {
            double rho_i, u_i, p_i;
            if (g.x[i] < 0.0) { rho_i = 1.0;   u_i = 0.0; p_i = 1.0; }
            else               { rho_i = 0.125; u_i = 0.0; p_i = 0.1; }
            double E_i = p_i / ((gamma - 1.0) * rho_i) + 0.5 * u_i * u_i;
            // rhov = 0: no transverse velocity in the initial Sod problem
            U[j][i] = { rho_i, rho_i * u_i, 0.0, rho_i * E_i };
        }
    }
}

void set_ic_wedge(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma,
                    double Mach, double theta)
{
    double rho0=1.0;
    double p0=1.0;
    double a0=sqrt(gamma*p0/rho0);
    double u0=Mach*a0;
    double v0=0.0;
    double E0=p0/((gamma-1.0)*rho0)+0.5*(u0*u0+v0*v0);
    for (int j = 0; j < g.Ny; j++) {
        for (int i = 0; i < g.Nx; i++) {
            U[j][i] = { rho0, rho0 * u0, rho0 * v0, rho0 * E0 };
        }
    }
}

void apply_inflow_bc(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma,
                    double Mach, double theta)
{
    double rho0=1.0;
    double p0=1.0;
    double a0=sqrt(gamma*p0/rho0);
    double u0=Mach*a0;
    double v0=0.0;
    double E0=p0/((gamma-1.0)*rho0)+0.5*(u0*u0+v0*v0);

    for (int j = 0; j < g.Ny; j++) {
        U[j][0] = { rho0, rho0 * u0, rho0 * v0, rho0 * E0 };
    }
}   

void apply_outflow_bc(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma)
{
    int Nx=g.Nx;
    for (int j = 0; j < g.Ny; j++) {
        //U[j][Nx-1] = U[j][Nx-2];
        U[j][Nx-1] = U[j][Nx-2];  // zero-gradient outflow (copy last interior state)
        U[j][Nx]=U[j][Nx-2]; 
       
         // Extra ghost cells for WENO5 stencil at the right boundary

    }
}

// 1D overload (single row, flat 1D vector)
void set_initial_conditions(vector<ConsvState>& U,
                             const Grid2D& g, double gamma)
{
    for (int i = 0; i < g.Nx; i++) {
        double rho_i, u_i, p_i;
        if (g.x[i] < 0.0) { rho_i = 1.0;   u_i = 0.0; p_i = 1.0; }
        else               { rho_i = 0.125; u_i = 0.0; p_i = 0.1; }
        double E_i = p_i / ((gamma - 1.0) * rho_i) + 0.5 * u_i * u_i;
        U[i] = { rho_i, rho_i * u_i, 0.0, rho_i * E_i };
    }
}
