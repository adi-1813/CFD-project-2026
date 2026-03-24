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
 * ── Spatial Discretisation: FDM-WENO3 with Lax-Friedrichs flux splitting ───
 *
 * This solver uses the FINITE DIFFERENCE METHOD (FDM) formulation of WENO,
 * NOT the finite volume (FVM) formulation.  The distinction matters:
 *
 *   FDM-WENO: U_i are POINT VALUES at node x_i.  WENO reconstructs a
 *             numerical flux function h(x) such that dh/dx = dF/dx to
 *             3rd-order accuracy at each node.  The spatial derivative is:
 *
 *               dF/dx|_i  ≈  (h^+[i] - h^+[i-1]) / dx     <- backward diff
 *                           + (h^-[i+1] - h^-[i]) / dx     <- forward  diff
 *
 *   FVM-WENO: U_i are CELL AVERAGES.  WENO reconstructs interface *states*
 *             U^L and U^R, a numerical flux F_hat is formed at each interface,
 *             and the residual is (F_hat_{i+1/2} - F_hat_{i-1/2}) / dx.
 *
 * Because the Euler flux F(U) is nonlinear, the two approaches are NOT
 * equivalent.
 *
 * WENO3 uses a 3-point stencil {i-1, i, i+1} with TWO candidate sub-stencils
 * (versus WENO5's 5-point stencil with three).  This gives 3rd-order accuracy
 * in smooth regions at lower stencil cost.
 *
 * Lax-Friedrichs flux splitting decomposes F into a positive part F^+ and
 * a negative part F^- using the global maximum wave speed alpha:
 *
 *   F^+(U) = 0.5 * (F(U) + alpha * U)    <- all "speeds" >= 0
 *   F^-(U) = 0.5 * (F(U) - alpha * U)    <- all "speeds" <= 0
 *
 * F^+ is reconstructed with the left-biased (upwind) WENO stencil wcL3 and
 * differenced backward.  F^- is reconstructed with the right-biased (downwind)
 * WENO stencil wcR3 and differenced forward.
 *
 * ── Boundary Conditions: Reflection / No-Penetration (solid wall) ───────────
 * Both ends of the tube are treated as solid walls.  The WENO3 stencil at the
 * boundary nodes needs one ghost point outside the domain.  The ghost state is
 * constructed by mirroring the first interior node with the normal (x) velocity
 * sign-flipped, enforcing zero normal velocity at the wall:
 *
 *   Ghost state at left  wall (node -1):  rho_{-1}    =  rho_1
 *                                         (rho*u)_{-1} = -(rho*u)_1   <- flip
 *                                         (rho*E)_{-1} =  (rho*E)_1
 *
 *   Ghost state at right wall (node N):   rho_N        =  rho_{N-2}
 *                                         (rho*u)_N    = -(rho*u)_{N-2} <- flip
 *                                         (rho*E)_N    =  (rho*E)_{N-2}
 *
 * The ghost state is converted into split fluxes F^+/F^- and passed directly
 * into wenoL_scalar / wenoR_scalar as explicit ghost values, so the residual
 * loop covers ALL nodes (i = 0..N-1) with no residual copy.
 * After each RK stage, rhou at the two boundary nodes is explicitly zeroed to
 * prevent velocity drift from accumulating over many time steps.
 *
 * ── Time Integration: SSP-RK3 ───────────────────────────────────────────────
 * 3rd-order Strong-Stability-Preserving Runge-Kutta (Shu & Osher, 1988):
 *
 *   U^(1)   = U^n + dt * L(U^n)                        ; enforce wall BC
 *   U^(2)   = 3/4*U^n + 1/4*(U^(1) + dt*L(U^(1)))     ; enforce wall BC
 *   U^(n+1) = 1/3*U^n + 2/3*(U^(2) + dt*L(U^(2)))     ; enforce wall BC
 *
 * where L(U) = -dF/dx is the FDM-WENO spatial residual.
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

// Pressure from conservative state via ideal-gas EOS:  p = (gamma-1)*rho*(E - u^2/2)
static inline double pressure(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double E = U.rhoE  / U.rho;
    return (gamma - 1.0) * U.rho * (E - 0.5 * u * u);
}

// Physical flux F(U)
static inline Flux compute_flux(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double p = pressure(U, gamma);
    return { U.rhou,
             U.rhou * u + p,
             u * (U.rhoE + p) };
}

// Sound speed: a = sqrt(gamma * p / rho)
static inline double sound_speed(const ConsvState& U, double gamma)
{
    return sqrt(gamma * pressure(U, gamma) / U.rho);
}

// Spectral radius of Euler Jacobian = |u| + a  (max characteristic speed)
static inline double max_wave_speed(const ConsvState& U, double gamma)
{
    return fabs(U.rhou / U.rho) + sound_speed(U, gamma);
}

// ============================================================================
// WENO3 nonlinear weight functions
// ============================================================================
/*
 * WENO3 uses a 3-point stencil {v1, v2, v3} with TWO candidate sub-stencils,
 * giving 3rd-order accuracy in smooth regions.
 *
 * Sub-stencil 0: {v1, v2}  (left/upwind pair)
 * Sub-stencil 1: {v2, v3}  (right/downwind pair)
 *
 * Smoothness indicators (each a squared first difference):
 *   beta_0 = (v2 - v1)^2
 *   beta_1 = (v3 - v2)^2
 *
 * Nonlinear weights:  alpha_k = d_k / (eps + beta_k)^2,  w_k = alpha_k/sum
 *
 * wcL3 — left-biased (upwind):  ideal weights d = {1/3, 2/3}
 *         Candidate polynomials evaluated at right edge i+1/2:
 *           q0 = -v1/2 + 3*v2/2          (stencil {v1,v2} extrapolated right)
 *           q1 =  v2/2 +   v3/2          (stencil {v2,v3} interpolated)
 *
 * wcR3 — right-biased (downwind): ideal weights d = {2/3, 1/3}
 *         Candidate polynomials evaluated at right edge i+1/2:
 *           q0 =  v1/2 +   v2/2          (stencil {v1,v2} interpolated)
 *           q1 =  3*v2/2 - v3/2          (stencil {v2,v3} extrapolated left)
 *
 * Note: v1=f[i-1], v2=f[i], v3=f[i+1] for the interior reconstruction.
 */
static double wcL3(double v1, double v2, double v3)
{
    double b0 = pow(v2 - v1, 2);
    double b1 = pow(v3 - v2, 2);

    double a0 = (1.0/3.0) / pow(WENO_EPS + b0, 2);
    double a1 = (2.0/3.0) / pow(WENO_EPS + b1, 2);
    double asum = a0 + a1;
    double w0 = a0/asum,  w1 = a1/asum;

    double q0 = -0.5*v1 + 1.5*v2;
    double q1 =  0.5*v2 + 0.5*v3;

    return w0*q0 + w1*q1;
}

static double wcR3(double v1, double v2, double v3)
{
    double b0 = pow(v2 - v1, 2);
    double b1 = pow(v3 - v2, 2);

    double a0 = (2.0/3.0) / pow(WENO_EPS + b0, 2);
    double a1 = (1.0/3.0) / pow(WENO_EPS + b1, 2);
    double asum = a0 + a1;
    double w0 = a0/asum,  w1 = a1/asum;

    double q0 =  0.5*v1 + 0.5*v2;
    double q1 =  1.5*v2 - 0.5*v3;

    return w0*q0 + w1*q1;
}

// ============================================================================
// wenoL_scalar / wenoR_scalar: reconstruct a split-flux array at all nodes
// ============================================================================
/*
 * ghost_L = f_{-1}: ghost value just outside the left boundary.
 * ghost_R = f_{N}:  ghost value just outside the right boundary.
 *
 * For transmissive outflow:  ghost = copy of boundary node.
 * For reflection / wall BC:  ghost = split flux of the mirror state (rhou
 *   sign-flipped), computed in compute_rhs and passed in here.
 *
 * Both functions now cover ALL nodes 0..N-1, so the residual loop can run
 * over the full domain with no special-case boundary copies.
 */
static void wenoL_scalar(int N, const vector<double>& f, vector<double>& hL,
                          double ghost_L, double ghost_R)
{
    // i = 0: left ghost supplied by caller (reflection or transmissive)
    hL[0] = wcL3(ghost_L, f[0], f[1]);

    for (int i = 1; i < N - 1; i++)
        hL[i] = wcL3(f[i-1], f[i], f[i+1]);

    // i = N-1: right ghost supplied by caller
    hL[N-1] = wcL3(f[N-2], f[N-1], ghost_R);
}

static void wenoR_scalar(int N, const vector<double>& f, vector<double>& hR,
                          double ghost_L, double ghost_R)
{
    // i = 0: left ghost supplied by caller
    hR[0] = wcR3(ghost_L, f[0], f[1]);

    for (int i = 1; i < N - 1; i++)
        hR[i] = wcR3(f[i-1], f[i], f[i+1]);

    // i = N-1: right ghost supplied by caller
    hR[N-1] = wcR3(f[N-2], f[N-1], ghost_R);
}

// ============================================================================
// Spatial residual  L(U) = -dF/dx  via FDM-WENO3 + Lax-Friedrichs splitting
// ============================================================================
/*
 * Steps:
 * 1. Compute global alpha.
 * 2. Build ghost states for both walls (rhou sign-flipped mirror).
 * 3. Compute split fluxes F^+/F^- at ghost nodes and all interior nodes.
 * 4. WENO3 reconstruct hL^+ and hR^- passing ghost split-flux values.
 * 5. Compute the "outer ghost" reconstructions hL[-1] and hR[N] needed at
 *    the very boundary of the backward/forward difference stencils.
 * 6. Assemble residual over ALL nodes 0..N-1.
 */
static void compute_rhs(const vector<ConsvState>& U, int N, double dx, double gamma,
                         vector<ConsvState>& rhs_out)
{
    // ── Step 1: Global Lax-Friedrichs splitting parameter ───────────────────
    double alpha = 0.0;
    for (int i = 0; i < N; i++)
        alpha = max(alpha, max_wave_speed(U[i], gamma));

    // ── Step 2: Ghost states for reflection / no-penetration BC ─────────────
    //
    // Left ghost (node -1) mirrors node 1: rho and rhoE even, rhou odd.
    // Using node 1 (not node 0) gives a proper odd-symmetric difference
    // across the wall so that the flux derivative sees u = 0 at x = x_L.
    ConsvState ghost_L;
    ghost_L.rho  =  U[1].rho;
    ghost_L.rhou = -U[1].rhou;   // <- sign flip: reflection condition
    ghost_L.rhoE =  U[1].rhoE;

    // Right ghost (node N) mirrors node N-2.
    ConsvState ghost_R;
    ghost_R.rho  =  U[N-2].rho;
    ghost_R.rhou = -U[N-2].rhou; // <- sign flip: reflection condition
    ghost_R.rhoE =  U[N-2].rhoE;

    // Physical fluxes at ghost nodes
    Flux F_gL = compute_flux(ghost_L, gamma);
    Flux F_gR = compute_flux(ghost_R, gamma);

    // Split fluxes at ghost nodes:  F^± = 0.5*(F ± alpha*U)
    double fp1_gL = 0.5*(F_gL.F1 + alpha*ghost_L.rho);
    double fm1_gL = 0.5*(F_gL.F1 - alpha*ghost_L.rho);
    double fp2_gL = 0.5*(F_gL.F2 + alpha*ghost_L.rhou);
    double fm2_gL = 0.5*(F_gL.F2 - alpha*ghost_L.rhou);
    double fp3_gL = 0.5*(F_gL.F3 + alpha*ghost_L.rhoE);
    double fm3_gL = 0.5*(F_gL.F3 - alpha*ghost_L.rhoE);

    double fp1_gR = 0.5*(F_gR.F1 + alpha*ghost_R.rho);
    double fm1_gR = 0.5*(F_gR.F1 - alpha*ghost_R.rho);
    double fp2_gR = 0.5*(F_gR.F2 + alpha*ghost_R.rhou);
    double fm2_gR = 0.5*(F_gR.F2 - alpha*ghost_R.rhou);
    double fp3_gR = 0.5*(F_gR.F3 + alpha*ghost_R.rhoE);
    double fm3_gR = 0.5*(F_gR.F3 - alpha*ghost_R.rhoE);

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

    // ── Step 4: WENO3 reconstruction with reflection ghost values ───────────
    vector<double> hL1(N,0), hR1(N,0);
    vector<double> hL2(N,0), hR2(N,0);
    vector<double> hL3(N,0), hR3(N,0);

    wenoL_scalar(N, fp1, hL1, fp1_gL, fp1_gR);
    wenoR_scalar(N, fm1, hR1, fm1_gL, fm1_gR);
    wenoL_scalar(N, fp2, hL2, fp2_gL, fp2_gR);
    wenoR_scalar(N, fm2, hR2, fm2_gL, fm2_gR);
    wenoL_scalar(N, fp3, hL3, fp3_gL, fp3_gR);
    wenoR_scalar(N, fm3, hR3, fm3_gL, fm3_gR);

    // ── Step 5: "Outer ghost" reconstructions needed at boundary residuals ───
    //
    // The backward diff at i=0  needs hL^+[-1]:
    //   i.e., wcL3 evaluated one node further left than the domain,
    //   using the ghost node as the central value.
    double hL1_ext = wcL3(fp1_gL, fp1_gL, fp1[0]);
    double hL2_ext = wcL3(fp2_gL, fp2_gL, fp2[0]);
    double hL3_ext = wcL3(fp3_gL, fp3_gL, fp3[0]);

    // The forward diff at i=N-1 needs hR^-[N]:
    //   i.e., wcR3 evaluated one node further right than the domain,
    //   using the ghost node as the central value.
    double hR1_ext = wcR3(fm1[N-1], fm1_gR, fm1_gR);
    double hR2_ext = wcR3(fm2[N-1], fm2_gR, fm2_gR);
    double hR3_ext = wcR3(fm3[N-1], fm3_gR, fm3_gR);

    // ── Step 6: Residual over ALL nodes including boundaries ─────────────────

    // Left boundary node (i = 0): backward diff uses hL_ext as hL[-1]
    rhs_out[0].rho  = -((hL1[0]-hL1_ext) + (hR1[1]-hR1[0])) / dx;
    rhs_out[0].rhou = -((hL2[0]-hL2_ext) + (hR2[1]-hR2[0])) / dx;
    rhs_out[0].rhoE = -((hL3[0]-hL3_ext) + (hR3[1]-hR3[0])) / dx;

    // Interior nodes: standard stencil
    for (int i = 1; i < N - 1; i++) {
        rhs_out[i].rho  = -((hL1[i]-hL1[i-1]) + (hR1[i+1]-hR1[i])) / dx;
        rhs_out[i].rhou = -((hL2[i]-hL2[i-1]) + (hR2[i+1]-hR2[i])) / dx;
        rhs_out[i].rhoE = -((hL3[i]-hL3[i-1]) + (hR3[i+1]-hR3[i])) / dx;
    }

    // Right boundary node (i = N-1): forward diff uses hR_ext as hR[N]
    rhs_out[N-1].rho  = -((hL1[N-1]-hL1[N-2]) + (hR1_ext-hR1[N-1])) / dx;
    rhs_out[N-1].rhou = -((hL2[N-1]-hL2[N-2]) + (hR2_ext-hR2[N-1])) / dx;
    rhs_out[N-1].rhoE = -((hL3[N-1]-hL3[N-2]) + (hR3_ext-hR3[N-1])) / dx;
}

// ── Enforce zero normal velocity at both walls ────────────────────────────────
// Called after every RK stage to prevent rhou drift accumulating over time.
// Only rhou is touched; rho and rhoE evolve freely from the residual.
static inline void apply_wall_bc(vector<ConsvState>& U, int N)
{
    U[0].rhou   = 0.0;   // left wall:  u = 0
    U[N-1].rhou = 0.0;   // right wall: u = 0
}

// ============================================================================
// 3rd-order SSP / TVD Runge-Kutta time step  (Shu & Osher 1988)
// ============================================================================
static double ssp_rk3_step(vector<ConsvState>& U, int N, double dx, double gamma,
                            double CFL, double time, double final_time)
{
    // CFL-based adaptive time step
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
    const double final_time = 0.6;   // run longer to observe wall reflections

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
    apply_wall_bc(U, N);   // u already 0, but enforce explicitly

    double time = 0.0;
    int    step = 0;

    write_html(x_arr, U, time, step, gamma, Rgas, false);
    char full_path[MAX_PATH];
    GetFullPathNameA("shock_tube.html", MAX_PATH, full_path, nullptr);
    ShellExecuteA(nullptr, "open", full_path, nullptr, nullptr, SW_SHOWNORMAL);
    Sleep(2000);
    cout << "Browser opened.  Scheme: FDM-WENO3 + SSP-RK3 + Wall (reflection) BC" << endl;

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