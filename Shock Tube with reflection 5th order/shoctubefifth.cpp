/*
 * ============================================================================
 * 1D Euler Shock Tube Solver
 * ============================================================================
 * Governing equations: 1D Euler equations in CONSERVATIVE form
 *
 *   dU/dt + dF/dx = 0
 *
 * where the conservative state vector U and flux vector F are:
 *
 *   U = [ rho       ]      F = [ rho*u           ]
 *       [ rho*u     ]          [ rho*u^2 + p      ]
 *       [ rho*E     ]          [ u*(rho*E + p)    ]
 *
 * and rho = density, u = velocity, p = pressure, E = total specific energy.
 * The equation of state closes the system: p = (gamma-1)*rho*(E - 0.5*u^2)
 *
 * ── Spatial Discretisation: FDM-WENO5 with Lax-Friedrichs flux splitting ───
 *
 * This solver uses the FINITE DIFFERENCE METHOD (FDM) formulation of WENO,
 * NOT the finite volume (FVM) formulation.  U_i are POINT VALUES at node x_i.
 * WENO reconstructs a numerical flux function h(x) such that dh/dx = dF/dx
 * to 5th-order accuracy at each node.  The spatial derivative is:
 *
 *   dF/dx|_i  ≈  (h^+[i] - h^+[i-1]) / dx     <- backward diff on F^+
 *               + (h^-[i+1] - h^-[i]) / dx     <- forward  diff on F^-
 *
 * WENO5 uses a 5-point stencil {i-2,..,i+2} with THREE candidate sub-stencils,
 * giving 5th-order accuracy in smooth regions.
 *
 * Lax-Friedrichs flux splitting decomposes F into F^+ and F^- using the
 * global maximum wave speed alpha:
 *
 *   F^+(U) = 0.5 * (F(U) + alpha * U)    <- all "speeds" >= 0
 *   F^-(U) = 0.5 * (F(U) - alpha * U)    <- all "speeds" <= 0
 *
 * F^+ is reconstructed with the left-biased (upwind)   WENO stencil wcL5.
 * F^- is reconstructed with the right-biased (downwind) WENO stencil wcR5.
 *
 * ── Ghost Points for WENO5 ──────────────────────────────────────────────────
 * WENO5 needs a 5-point stencil, so boundary nodes require TWO ghost points
 * outside each end of the domain (compared to one for WENO3).
 *
 * Ghost states are constructed by reflecting interior nodes with the normal
 * velocity sign-flipped (no-penetration / solid wall):
 *
 *   Left wall ghosts:   g_{-1}  mirrors node 1  (rhou sign-flipped)
 *                       g_{-2}  mirrors node 2  (rhou sign-flipped)
 *
 *   Right wall ghosts:  g_{N}   mirrors node N-2 (rhou sign-flipped)
 *                       g_{N+1} mirrors node N-3 (rhou sign-flipped)
 *
 * ── Boundary Conditions ─────────────────────────────────────────────────────
 * Both ends are solid walls (reflection / no-penetration).
 * After each RK stage, rhou at the two wall nodes is explicitly zeroed.
 *
 * ── Time Integration: SSP-RK3 ───────────────────────────────────────────────
 * 3rd-order Strong-Stability-Preserving Runge-Kutta (Shu & Osher, 1988):
 *
 *   U^(1)   = U^n + dt * L(U^n)                        ; enforce wall BC
 *   U^(2)   = 3/4*U^n + 1/4*(U^(1) + dt*L(U^(1)))     ; enforce wall BC
 *   U^(n+1) = 1/3*U^n + 2/3*(U^(2) + dt*L(U^(2)))     ; enforce wall BC
 *
 * where L(U) = -dF/dx is the FDM-WENO5 spatial residual.
 * ============================================================================
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <windows.h>
#pragma comment(lib, "Shell32.lib")

using namespace std;

// ── Conservative state vector: U = [rho, rho*u, rho*E] ──────────────────────
struct ConsvState {
    double rho;   // U1: density
    double rhou;  // U2: momentum  (= rho * u)
    double rhoE;  // U3: total energy (= rho * E)
};

// ── Physical flux vector: F(U) = [F1, F2, F3] ───────────────────────────────
struct Flux {
    double F1;  // F1 = rho * u
    double F2;  // F2 = rho*u^2 + p
    double F3;  // F3 = u * (rho*E + p)
};

// ── WENO epsilon to prevent division by zero in smoothness weights ───────────
static constexpr double WENO_EPS = 1.0e-6;

// ============================================================================
// Primitive / thermodynamic helpers
// ============================================================================

static inline double pressure(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double E = U.rhoE  / U.rho;
    return (gamma - 1.0) * U.rho * (E - 0.5 * u * u);
}

static inline Flux compute_flux(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double p = pressure(U, gamma);
    return { U.rhou,
             U.rhou * u + p,
             u * (U.rhoE + p) };
}

static inline double sound_speed(const ConsvState& U, double gamma)
{
    return sqrt(gamma * pressure(U, gamma) / U.rho);
}

static inline double max_wave_speed(const ConsvState& U, double gamma)
{
    return fabs(U.rhou / U.rho) + sound_speed(U, gamma);
}

// ============================================================================
// WENO5 nonlinear weight functions  (Liu, Osher & Chan 1994; Jiang & Shu 1996)
// ============================================================================
/*
 * Both functions take five consecutive split-flux values {v1, v2, v3, v4, v5}
 * and return the WENO5 reconstructed value at the right edge of cell i,
 * i.e., at x_{i+1/2}.
 *
 * The mapping to grid nodes for an interior node i is:
 *   v1=f[i-2], v2=f[i-1], v3=f[i], v4=f[i+1], v5=f[i+2]
 *
 * ── wcL5 — left-biased (upwind), used for F^+ ───────────────────────────────
 *   Ideal weights: d = {1/10, 6/10, 3/10}  (sum = 1)
 *
 *   Three candidate 3rd-order polynomials evaluated at x_{i+1/2}:
 *     q0 =  1/3*v1 - 7/6*v2 + 11/6*v3          (stencil {i-2,i-1,i})
 *     q1 = -1/6*v2 + 5/6*v3 +  1/3*v4          (stencil {i-1, i, i+1})
 *     q2 =  1/3*v3 + 5/6*v4 -  1/6*v5          (stencil { i, i+1,i+2})
 *
 *   Smoothness indicators (Jiang-Shu):
 *     beta0 = 13/12*(v1-2v2+v3)^2 + 1/4*(v1-4v2+3v3)^2
 *     beta1 = 13/12*(v2-2v3+v4)^2 + 1/4*(v2-v4)^2
 *     beta2 = 13/12*(v3-2v4+v5)^2 + 1/4*(3v3-4v4+v5)^2
 *
 * ── wcR5 — right-biased (downwind), used for F^- ────────────────────────────
 *   Ideal weights: d = {3/10, 6/10, 1/10}  (mirror of wcL5)
 *
 *   Three candidate 3rd-order polynomials evaluated at x_{i+1/2}:
 *     q0 = -1/6*v1 + 5/6*v2 +  1/3*v3          (stencil {i-2,i-1,i})
 *     q1 =  1/3*v2 + 5/6*v3 -  1/6*v4          (stencil {i-1, i, i+1})
 *     q2 = 11/6*v3 - 7/6*v4 +  1/3*v5          (stencil { i, i+1,i+2})
 *
 *   Same smoothness indicators as wcL5.
 *
 * The nonlinear weights for both:
 *   alpha_k = d_k / (eps + beta_k)^2,   w_k = alpha_k / sum(alpha)
 */
static double wcL5(double v1, double v2, double v3, double v4, double v5)
{
    // Smoothness indicators
    double b0 = (13.0/12.0)*pow(v1 - 2.0*v2 + v3, 2) + 0.25*pow(v1 - 4.0*v2 + 3.0*v3, 2);
    double b1 = (13.0/12.0)*pow(v2 - 2.0*v3 + v4, 2) + 0.25*pow(v2 - v4, 2);
    double b2 = (13.0/12.0)*pow(v3 - 2.0*v4 + v5, 2) + 0.25*pow(3.0*v3 - 4.0*v4 + v5, 2);

    // Ideal weights (left-biased / upwind)
    double a0 = 0.1 / pow(WENO_EPS + b0, 2);
    double a1 = 0.6 / pow(WENO_EPS + b1, 2);
    double a2 = 0.3 / pow(WENO_EPS + b2, 2);
    double asum = a0 + a1 + a2;
    double w0 = a0/asum, w1 = a1/asum, w2 = a2/asum;

    // Candidate stencil reconstructions at x_{i+1/2}
    double q0 =  v1/3.0 - 7.0/6.0*v2 + 11.0/6.0*v3;
    double q1 = -v2/6.0 + 5.0/6.0*v3 +       v4/3.0;
    double q2 =  v3/3.0 + 5.0/6.0*v4 -       v5/6.0;

    return w0*q0 + w1*q1 + w2*q2;
}

static double wcR5(double v1, double v2, double v3, double v4, double v5)
{
    // Same smoothness indicators (symmetric)
    double b0 = (13.0/12.0)*pow(v1 - 2.0*v2 + v3, 2) + 0.25*pow(v1 - 4.0*v2 + 3.0*v3, 2);
    double b1 = (13.0/12.0)*pow(v2 - 2.0*v3 + v4, 2) + 0.25*pow(v2 - v4, 2);
    double b2 = (13.0/12.0)*pow(v3 - 2.0*v4 + v5, 2) + 0.25*pow(3.0*v3 - 4.0*v4 + v5, 2);

    // Ideal weights (right-biased / downwind) — mirror of wcL5
    double a0 = 0.3 / pow(WENO_EPS + b0, 2);
    double a1 = 0.6 / pow(WENO_EPS + b1, 2);
    double a2 = 0.1 / pow(WENO_EPS + b2, 2);
    double asum = a0 + a1 + a2;
    double w0 = a0/asum, w1 = a1/asum, w2 = a2/asum;

    // Candidate stencil reconstructions at x_{i+1/2}
    double q0 = -v1/6.0 + 5.0/6.0*v2 +       v3/3.0;
    double q1 =  v2/3.0 + 5.0/6.0*v3 -       v4/6.0;
    double q2 = 11.0/6.0*v3 - 7.0/6.0*v4 + v5/3.0;

    return w0*q0 + w1*q1 + w2*q2;
}

// ============================================================================
// wenoL_scalar / wenoR_scalar: reconstruct a split-flux array at all nodes
// ============================================================================
/*
 * WENO5 requires a 5-point stencil {i-2,i-1,i,i+1,i+2}.
 * At the boundaries the stencil reaches outside the domain, so we need
 * two ghost values on each side.
 *
 * ghost_Lm1 = f_{-1},  ghost_Lm2 = f_{-2}   (left-side ghosts)
 * ghost_Rp0 = f_{N},   ghost_Rp1 = f_{N+1}   (right-side ghosts)
 *
 * wenoL_scalar reconstructs h^+[i] = wcL5(...) at all nodes i = 0..N-1.
 * wenoR_scalar reconstructs h^-[i] = wcR5(...) at all nodes i = 0..N-1.
 *
 * The residual uses:
 *   dF^+/dx|_i ≈ (h^+[i] - h^+[i-1]) / dx     (backward difference)
 *   dF^-/dx|_i ≈ (h^-[i+1] - h^-[i]) / dx     (forward  difference)
 *
 * For the backward diff at i=0  we need h^+[-1], and for the forward diff
 * at i=N-1 we need h^-[N].  These are computed in compute_rhs as "outer
 * ghost" reconstructions (hL_ext and hR_ext below).
 */
static void wenoL_scalar(int N, const vector<double>& f,
                          vector<double>& hL,
                          double gLm1, double gLm2,   // f[-1], f[-2]
                          double gRp0, double gRp1)   // f[N],  f[N+1]
{
    // i = 0: needs f[-2], f[-1], f[0], f[1], f[2]
    hL[0] = wcL5(gLm2, gLm1, f[0], f[1], f[2]);

    // i = 1: needs f[-1], f[0], f[1], f[2], f[3]
    hL[1] = wcL5(gLm1, f[0], f[1], f[2], f[3]);

    // Interior
    for (int i = 2; i < N - 2; i++)
        hL[i] = wcL5(f[i-2], f[i-1], f[i], f[i+1], f[i+2]);

    // i = N-2: needs f[N-4], f[N-3], f[N-2], f[N-1], f[N]
    hL[N-2] = wcL5(f[N-4], f[N-3], f[N-2], f[N-1], gRp0);

    // i = N-1: needs f[N-3], f[N-2], f[N-1], f[N], f[N+1]
    hL[N-1] = wcL5(f[N-3], f[N-2], f[N-1], gRp0, gRp1);
}

static void wenoR_scalar(int N, const vector<double>& f,
                          vector<double>& hR,
                          double gLm1, double gLm2,   // f[-1], f[-2]
                          double gRp0, double gRp1)   // f[N],  f[N+1]
{
    // i = 0: needs f[-2], f[-1], f[0], f[1], f[2]
    hR[0] = wcR5(gLm2, gLm1, f[0], f[1], f[2]);

    // i = 1: needs f[-1], f[0], f[1], f[2], f[3]
    hR[1] = wcR5(gLm1, f[0], f[1], f[2], f[3]);

    // Interior
    for (int i = 2; i < N - 2; i++)
        hR[i] = wcR5(f[i-2], f[i-1], f[i], f[i+1], f[i+2]);

    // i = N-2: needs f[N-4], f[N-3], f[N-2], f[N-1], f[N]
    hR[N-2] = wcR5(f[N-4], f[N-3], f[N-2], f[N-1], gRp0);

    // i = N-1: needs f[N-3], f[N-2], f[N-1], f[N], f[N+1]
    hR[N-1] = wcR5(f[N-3], f[N-2], f[N-1], gRp0, gRp1);
}

// ============================================================================
// Spatial residual  L(U) = -dF/dx  via FDM-WENO5 + Lax-Friedrichs splitting
// ============================================================================
/*
 * Steps:
 * 1. Compute global alpha (max wave speed across all nodes).
 * 2. Build two ghost states on each side via the reflection/wall BC.
 * 3. Compute split fluxes F^+/F^- at interior nodes and all four ghost nodes.
 * 4. WENO5 reconstruct hL^+ and hR^- passing the four ghost split-flux values.
 * 5. Compute "outer ghost" reconstructions hL[-1] and hR[N] needed by the
 *    backward/forward difference at the boundary nodes.
 * 6. Assemble residual over ALL nodes 0..N-1.
 */
static void compute_rhs(const vector<ConsvState>& U, int N, double dx, double gamma,
                         vector<ConsvState>& rhs_out)
{
    // ── Step 1: Global Lax-Friedrichs splitting parameter ───────────────────
    double alpha = 0.0;
    for (int i = 0; i < N; i++)
        alpha = max(alpha, max_wave_speed(U[i], gamma));

    // ── Step 2: Ghost states — reflection / no-penetration BC ───────────────
    //
    // Two ghost nodes are needed on each side for the WENO5 5-point stencil.
    //
    // Left wall: ghost_{-1} mirrors node 1, ghost_{-2} mirrors node 2.
    // Using nodes 1 and 2 (not 0) gives proper odd symmetry in rhou so the
    // flux derivative sees u = 0 at the wall.
    ConsvState gL1, gL2;   // g_{-1} and g_{-2}
    gL1.rho  =  U[1].rho;   gL1.rhou = -U[1].rhou;   gL1.rhoE =  U[1].rhoE;
    gL2.rho  =  U[2].rho;   gL2.rhou = -U[2].rhou;   gL2.rhoE =  U[2].rhoE;

    // Right wall: ghost_{N} mirrors node N-2, ghost_{N+1} mirrors node N-3.
    ConsvState gR0, gR1;   // g_{N} and g_{N+1}
    gR0.rho  =  U[N-2].rho; gR0.rhou = -U[N-2].rhou; gR0.rhoE =  U[N-2].rhoE;
    gR1.rho  =  U[N-3].rho; gR1.rhou = -U[N-3].rhou; gR1.rhoE =  U[N-3].rhoE;

    // Physical fluxes at ghost nodes
    Flux FgL1 = compute_flux(gL1, gamma);
    Flux FgL2 = compute_flux(gL2, gamma);
    Flux FgR0 = compute_flux(gR0, gamma);
    Flux FgR1 = compute_flux(gR1, gamma);

    // Split fluxes at ghost nodes: F^± = 0.5*(F ± alpha*U)
    // Component 1 (mass), ghost left -1 and -2
    double fp1_gL1 = 0.5*(FgL1.F1 + alpha*gL1.rho);   double fm1_gL1 = 0.5*(FgL1.F1 - alpha*gL1.rho);
    double fp1_gL2 = 0.5*(FgL2.F1 + alpha*gL2.rho);   double fm1_gL2 = 0.5*(FgL2.F1 - alpha*gL2.rho);
    // Component 2 (momentum), ghost left
    double fp2_gL1 = 0.5*(FgL1.F2 + alpha*gL1.rhou);  double fm2_gL1 = 0.5*(FgL1.F2 - alpha*gL1.rhou);
    double fp2_gL2 = 0.5*(FgL2.F2 + alpha*gL2.rhou);  double fm2_gL2 = 0.5*(FgL2.F2 - alpha*gL2.rhou);
    // Component 3 (energy), ghost left
    double fp3_gL1 = 0.5*(FgL1.F3 + alpha*gL1.rhoE);  double fm3_gL1 = 0.5*(FgL1.F3 - alpha*gL1.rhoE);
    double fp3_gL2 = 0.5*(FgL2.F3 + alpha*gL2.rhoE);  double fm3_gL2 = 0.5*(FgL2.F3 - alpha*gL2.rhoE);

    // Component 1, ghost right +0 and +1
    double fp1_gR0 = 0.5*(FgR0.F1 + alpha*gR0.rho);   double fm1_gR0 = 0.5*(FgR0.F1 - alpha*gR0.rho);
    double fp1_gR1 = 0.5*(FgR1.F1 + alpha*gR1.rho);   double fm1_gR1 = 0.5*(FgR1.F1 - alpha*gR1.rho);
    // Component 2, ghost right
    double fp2_gR0 = 0.5*(FgR0.F2 + alpha*gR0.rhou);  double fm2_gR0 = 0.5*(FgR0.F2 - alpha*gR0.rhou);
    double fp2_gR1 = 0.5*(FgR1.F2 + alpha*gR1.rhou);  double fm2_gR1 = 0.5*(FgR1.F2 - alpha*gR1.rhou);
    // Component 3, ghost right
    double fp3_gR0 = 0.5*(FgR0.F3 + alpha*gR0.rhoE);  double fm3_gR0 = 0.5*(FgR0.F3 - alpha*gR0.rhoE);
    double fp3_gR1 = 0.5*(FgR1.F3 + alpha*gR1.rhoE);  double fm3_gR1 = 0.5*(FgR1.F3 - alpha*gR1.rhoE);

    // ── Step 3: Split flux arrays at all interior nodes ──────────────────────
    vector<double> fp1(N), fm1(N);
    vector<double> fp2(N), fm2(N);
    vector<double> fp3(N), fm3(N);

    for (int i = 0; i < N; i++) {
        Flux F = compute_flux(U[i], gamma);
        fp1[i] = 0.5*(F.F1 + alpha*U[i].rho);   fm1[i] = 0.5*(F.F1 - alpha*U[i].rho);
        fp2[i] = 0.5*(F.F2 + alpha*U[i].rhou);  fm2[i] = 0.5*(F.F2 - alpha*U[i].rhou);
        fp3[i] = 0.5*(F.F3 + alpha*U[i].rhoE);  fm3[i] = 0.5*(F.F3 - alpha*U[i].rhoE);
    }

    // ── Step 4: WENO5 reconstruction with two reflection ghost values each side
    vector<double> hL1(N,0), hR1(N,0);
    vector<double> hL2(N,0), hR2(N,0);
    vector<double> hL3(N,0), hR3(N,0);

    wenoL_scalar(N, fp1, hL1, fp1_gL1, fp1_gL2, fp1_gR0, fp1_gR1);
    wenoR_scalar(N, fm1, hR1, fm1_gL1, fm1_gL2, fm1_gR0, fm1_gR1);
    wenoL_scalar(N, fp2, hL2, fp2_gL1, fp2_gL2, fp2_gR0, fp2_gR1);
    wenoR_scalar(N, fm2, hR2, fm2_gL1, fm2_gL2, fm2_gR0, fm2_gR1);
    wenoL_scalar(N, fp3, hL3, fp3_gL1, fp3_gL2, fp3_gR0, fp3_gR1);
    wenoR_scalar(N, fm3, hR3, fm3_gL1, fm3_gL2, fm3_gR0, fm3_gR1);

    // ── Step 5: "Outer ghost" reconstructions needed at boundary residuals ───
    //
    // Backward diff at i=0  needs h^+[-1]:
    //   wcL5 centered at node -1 → stencil {-3,-2,-1,0,1}
    //   We use double-ghost values in place of nodes -3 and -2:
    //   approx: v1=gL2, v2=gL2, v3=gL1, v4=f[0], v5=f[1]
    //   (gL2 plays the role of both f[-3] and f[-2]; sufficient for stability)
    double hL1_ext = wcL5(fp1_gL2, fp1_gL2, fp1_gL1, fp1[0], fp1[1]);
    double hL2_ext = wcL5(fp2_gL2, fp2_gL2, fp2_gL1, fp2[0], fp2[1]);
    double hL3_ext = wcL5(fp3_gL2, fp3_gL2, fp3_gL1, fp3[0], fp3[1]);

    // Forward diff at i=N-1 needs h^-[N]:
    //   wcR5 centered at node N → stencil {N-2,N-1,N,N+1,N+2}
    //   approx: v1=f[N-2], v2=f[N-1], v3=gR0, v4=gR1, v5=gR1
    double hR1_ext = wcR5(fm1[N-2], fm1[N-1], fm1_gR0, fm1_gR1, fm1_gR1);
    double hR2_ext = wcR5(fm2[N-2], fm2[N-1], fm2_gR0, fm2_gR1, fm2_gR1);
    double hR3_ext = wcR5(fm3[N-2], fm3[N-1], fm3_gR0, fm3_gR1, fm3_gR1);

    // ── Step 6: Residual over ALL nodes including boundaries ─────────────────

    // Left boundary node (i = 0)
    rhs_out[0].rho  = -((hL1[0]-hL1_ext) + (hR1[1]-hR1[0])) / dx;
    rhs_out[0].rhou = -((hL2[0]-hL2_ext) + (hR2[1]-hR2[0])) / dx;
    rhs_out[0].rhoE = -((hL3[0]-hL3_ext) + (hR3[1]-hR3[0])) / dx;

    // Interior nodes
    for (int i = 1; i < N - 1; i++) {
        rhs_out[i].rho  = -((hL1[i]-hL1[i-1]) + (hR1[i+1]-hR1[i])) / dx;
        rhs_out[i].rhou = -((hL2[i]-hL2[i-1]) + (hR2[i+1]-hR2[i])) / dx;
        rhs_out[i].rhoE = -((hL3[i]-hL3[i-1]) + (hR3[i+1]-hR3[i])) / dx;
    }

    // Right boundary node (i = N-1)
    rhs_out[N-1].rho  = -((hL1[N-1]-hL1[N-2]) + (hR1_ext-hR1[N-1])) / dx;
    rhs_out[N-1].rhou = -((hL2[N-1]-hL2[N-2]) + (hR2_ext-hR2[N-1])) / dx;
    rhs_out[N-1].rhoE = -((hL3[N-1]-hL3[N-2]) + (hR3_ext-hR3[N-1])) / dx;
}

// ── Enforce zero normal velocity at both walls ────────────────────────────────
static inline void apply_wall_bc(vector<ConsvState>& U, int N)
{
    U[0].rhou   = 0.0;
    U[N-1].rhou = 0.0;
}

// ============================================================================
// 3rd-order SSP / TVD Runge-Kutta time step  (Shu & Osher 1988)
// ============================================================================
static double ssp_rk3_step(vector<ConsvState>& U, int N, double dx, double gamma,
                            double CFL, double time, double final_time)
{
    double max_speed = 0.0;
    for (int i = 0; i < N; i++)
        max_speed = max(max_speed, max_wave_speed(U[i], gamma));
    double dt = CFL * dx / max_speed;
    if (time + dt > final_time) dt = final_time - time;

    vector<ConsvState> L0(N), L1(N), L2(N);
    vector<ConsvState> U1(N), U2(N);

    // Stage 1: U^(1) = U^n + dt * L(U^n)
    compute_rhs(U, N, dx, gamma, L0);
    for (int i = 0; i < N; i++) {
        U1[i].rho  = U[i].rho  + dt*L0[i].rho;
        U1[i].rhou = U[i].rhou + dt*L0[i].rhou;
        U1[i].rhoE = U[i].rhoE + dt*L0[i].rhoE;
    }
    apply_wall_bc(U1, N);

    // Stage 2: U^(2) = 3/4*U^n + 1/4*(U^(1) + dt*L(U^(1)))
    compute_rhs(U1, N, dx, gamma, L1);
    for (int i = 0; i < N; i++) {
        U2[i].rho  = 0.75*U[i].rho  + 0.25*(U1[i].rho  + dt*L1[i].rho);
        U2[i].rhou = 0.75*U[i].rhou + 0.25*(U1[i].rhou + dt*L1[i].rhou);
        U2[i].rhoE = 0.75*U[i].rhoE + 0.25*(U1[i].rhoE + dt*L1[i].rhoE);
    }
    apply_wall_bc(U2, N);

    // Stage 3: U^(n+1) = 1/3*U^n + 2/3*(U^(2) + dt*L(U^(2)))
    compute_rhs(U2, N, dx, gamma, L2);
    for (int i = 0; i < N; i++) {
        U[i].rho  = (1.0/3.0)*U[i].rho  + (2.0/3.0)*(U2[i].rho  + dt*L2[i].rho);
        U[i].rhou = (1.0/3.0)*U[i].rhou + (2.0/3.0)*(U2[i].rhou + dt*L2[i].rhou);
        U[i].rhoE = (1.0/3.0)*U[i].rhoE + (2.0/3.0)*(U2[i].rhoE + dt*L2[i].rhoE);
    }
    apply_wall_bc(U, N);

    return dt;
}

// ============================================================================
// HTML live-plot output
// ============================================================================
void write_html(const vector<double>&     x_arr,
                const vector<ConsvState>& U,
                double time, int step,
                double gamma, double Rgas,
                bool finished)
{
    int N = (int)x_arr.size();

    FILE* f = fopen("shock_tube_tmp.html", "w");
    if (!f) { cerr << "ERROR: Cannot write HTML\n"; return; }

    vector<double> u_arr(N), T_arr(N);
    for (int i = 0; i < N; i++) {
        double u_i = U[i].rhou / U[i].rho;
        double p_i = pressure(U[i], gamma);
        u_arr[i]   = u_i;
        T_arr[i]   = p_i / (U[i].rho * Rgas);
    }

    fprintf(f, "<!DOCTYPE html><html><head>\n");
    fprintf(f, "<style>\n");
    fprintf(f, "  body{background:#1a1a2e;color:#eee;font-family:monospace;margin:16px;}\n");
    fprintf(f, "  h2{color:#4e9af1;margin-bottom:6px;}\n");
    fprintf(f, "  .status{font-size:13px;margin-bottom:10px;color:%s;}\n",
            finished ? "#4ef1a4" : "#f1a44e");
    fprintf(f, "  canvas{display:block;background:#0f1923;border:1px solid #2a2a4a;margin-bottom:10px;}\n");
    fprintf(f, "</style></head><body>\n");
    fprintf(f, "<h2>1D Shock Tube (Wall BC) &nbsp;|&nbsp; t = %.5f &nbsp;|&nbsp; step = %d</h2>\n", time, step);
    fprintf(f, "<div class='status'>%s</div>\n",
            finished ? "&#10003; Simulation complete." : "&#9654; Running...");
    fprintf(f, "<canvas id='c1' width='900' height='180'></canvas>\n");
    fprintf(f, "<canvas id='c2' width='900' height='180'></canvas>\n");
    fprintf(f, "<canvas id='c3' width='900' height='180'></canvas>\n");

    fprintf(f, "<script>\n");
    fprintf(f, "var xd=["); for(int i=0;i<N;i++) fprintf(f,"%.5f%s",x_arr[i],  i<N-1?",":""); fprintf(f,"];\n");
    fprintf(f, "var d1=["); for(int i=0;i<N;i++) fprintf(f,"%.5f%s",U[i].rho,  i<N-1?",":""); fprintf(f,"];\n");
    fprintf(f, "var d2=["); for(int i=0;i<N;i++) fprintf(f,"%.5f%s",u_arr[i],  i<N-1?",":""); fprintf(f,"];\n");
    fprintf(f, "var d3=["); for(int i=0;i<N;i++) fprintf(f,"%.5f%s",T_arr[i],  i<N-1?",":""); fprintf(f,"];\n");

    fprintf(f, "function drawChart(id,data,label,color,ymin,ymax){\n");
    fprintf(f, "  var cv=document.getElementById(id),ctx=cv.getContext('2d');\n");
    fprintf(f, "  var W=cv.width,H=cv.height;\n");
    fprintf(f, "  var pl=60,pr=15,pt=25,pb=35;\n");
    fprintf(f, "  var pw=W-pl-pr,ph=H-pt-pb;\n");
    fprintf(f, "  ctx.fillStyle='#0f1923';ctx.fillRect(0,0,W,H);\n");
    fprintf(f, "  if(ymin===undefined){ymin=data.reduce((a,b)=>Math.min(a,b),Infinity);ymax=data.reduce((a,b)=>Math.max(a,b),-Infinity);}\n");
    fprintf(f, "  var yp=(ymax-ymin)*0.08||0.05;ymin-=yp;ymax+=yp;\n");
    fprintf(f, "  var xmin=xd[0],xmax=xd[xd.length-1];\n");
    fprintf(f, "  function sx(x){return pl+(x-xmin)/(xmax-xmin)*pw;}\n");
    fprintf(f, "  function sy(y){return pt+ph-(y-ymin)/(ymax-ymin)*ph;}\n");
    fprintf(f, "  ctx.strokeStyle='#1e2d3d';ctx.lineWidth=1;\n");
    fprintf(f, "  for(var i=0;i<=5;i++){\n");
    fprintf(f, "    var gy=ymin+i*(ymax-ymin)/5,gx=xmin+i*(xmax-xmin)/5;\n");
    fprintf(f, "    ctx.beginPath();ctx.moveTo(pl,sy(gy));ctx.lineTo(pl+pw,sy(gy));ctx.stroke();\n");
    fprintf(f, "    ctx.beginPath();ctx.moveTo(sx(gx),pt);ctx.lineTo(sx(gx),pt+ph);ctx.stroke();\n");
    fprintf(f, "  }\n");
    fprintf(f, "  ctx.strokeStyle='#4a5a6a';ctx.lineWidth=1.5;ctx.strokeRect(pl,pt,pw,ph);\n");
    fprintf(f, "  ctx.fillStyle='#8899aa';ctx.font='11px monospace';ctx.textAlign='center';\n");
    fprintf(f, "  for(var i=0;i<=5;i++){var gx=xmin+i*(xmax-xmin)/5;ctx.fillText(gx.toFixed(2),sx(gx),pt+ph+18);}\n");
    fprintf(f, "  ctx.textAlign='right';\n");
    fprintf(f, "  for(var i=0;i<=5;i++){var gy=ymin+i*(ymax-ymin)/5;ctx.fillText(gy.toFixed(3),pl-5,sy(gy)+4);}\n");
    fprintf(f, "  ctx.fillStyle='#aabbcc';ctx.font='bold 12px monospace';ctx.textAlign='center';\n");
    fprintf(f, "  ctx.fillText('x',pl+pw/2,H-4);\n");
    fprintf(f, "  ctx.save();ctx.translate(13,pt+ph/2);ctx.rotate(-Math.PI/2);ctx.fillText(label,0,0);ctx.restore();\n");
    fprintf(f, "  ctx.strokeStyle=color;ctx.lineWidth=2;ctx.beginPath();\n");
    fprintf(f, "  for(var i=0;i<xd.length;i++){if(i===0)ctx.moveTo(sx(xd[i]),sy(data[i]));else ctx.lineTo(sx(xd[i]),sy(data[i]));}\n");
    fprintf(f, "  ctx.stroke();\n");
    fprintf(f, "  ctx.fillStyle=color;ctx.font='bold 13px monospace';ctx.textAlign='left';\n");
    fprintf(f, "  ctx.fillText(label,pl+6,pt+16);\n");
    fprintf(f, "}\n");

    fprintf(f, "drawChart('c1',d1,'Density',    '#4e9af1',0,1.1);\n");
    fprintf(f, "drawChart('c2',d2,'Velocity',   '#f1a44e');\n");
    fprintf(f, "drawChart('c3',d3,'Temperature','#4ef1a4');\n");

    if (!finished)
        fprintf(f, "setTimeout(function(){location.reload();}, 80);\n");

    fprintf(f, "</script></body></html>\n");
    fclose(f);

    remove("shock_tube.html");
    rename("shock_tube_tmp.html", "shock_tube.html");
}

// ============================================================================
// Main
// ============================================================================
int main()
{
    const int    N          = 400;
    const double x_min      = -0.5;
    const double x_max      =  0.5;
    const double dx         = (x_max - x_min) / (N - 1);
    const double gamma      = 1.4;
    const double Rgas       = 1.0;
    const double CFL        = 0.5;
    const double final_time = 0.6;

    vector<double> x_arr(N);
    for (int i = 0; i < N; i++) x_arr[i] = x_min + i * dx;

    vector<ConsvState> U(N);

    // Initial conditions: Sod shock tube
    for (int i = 0; i < N; i++) {
        double rho_i, u_i, p_i;
        if (x_arr[i] < 0.0) { rho_i = 1.0;   u_i = 0.0; p_i = 1.0; }
        else                 { rho_i = 0.125; u_i = 0.0; p_i = 0.1; }
        double E_i = p_i / ((gamma - 1.0) * rho_i) + 0.5 * u_i * u_i;
        U[i] = { rho_i, rho_i * u_i, rho_i * E_i };
    }
    apply_wall_bc(U, N);

    double time = 0.0;
    int    step = 0;

    write_html(x_arr, U, time, step, gamma, Rgas, false);
    char full_path[MAX_PATH];
    GetFullPathNameA("shock_tube.html", MAX_PATH, full_path, nullptr);
    ShellExecuteA(nullptr, "open", full_path, nullptr, nullptr, SW_SHOWNORMAL);
    Sleep(2000);
    cout << "Browser opened.  Scheme: FDM-WENO5 + SSP-RK3 + Wall (reflection) BC" << endl;

    while (time < final_time)
    {
        if (isnan(U[N/2].rho) || isinf(U[N/2].rho)) {
            cout << "DIVERGED at step " << step << " t=" << time << endl;
            break;
        }

        double dt = ssp_rk3_step(U, N, dx, gamma, CFL, time, final_time);
        time += dt;
        step++;

        if (step % 5 == 0 || time >= final_time)
            write_html(x_arr, U, time, step, gamma, Rgas, time >= final_time);
    }

    cout << "Finished at t = " << time << " (" << step << " steps)" << endl;
    return 0;
}