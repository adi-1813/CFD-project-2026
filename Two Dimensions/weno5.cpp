#include "weno5.h"
#include "geometry.h"
#include <cmath>
#include <algorithm>

using namespace std;

// ============================================================================
// wcL5 — left-biased WENO5 reconstruction at x_{i+1/2}
// ============================================================================
// Three candidate 3rd-order polynomials (sub-stencils):
//   q0 =  1/3*v1 - 7/6*v2 + 11/6*v3    {i-2, i-1, i}
//   q1 = -1/6*v2 + 5/6*v3 +  1/3*v4    {i-1,  i, i+1}
//   q2 =  1/3*v3 + 5/6*v4 -  1/6*v5    { i, i+1, i+2}
//
// Jiang-Shu smoothness indicators:
//   beta0 = 13/12*(v1-2v2+v3)^2 + 1/4*(v1-4v2+3v3)^2
//   beta1 = 13/12*(v2-2v3+v4)^2 + 1/4*(v2-v4)^2
//   beta2 = 13/12*(v3-2v4+v5)^2 + 1/4*(3v3-4v4+v5)^2
double wcL5(double v1, double v2, double v3, double v4, double v5)
{
    double b0 = (13.0/12.0)*pow(v1 - 2.0*v2 + v3, 2) + 0.25*pow(v1 - 4.0*v2 + 3.0*v3, 2);
    double b1 = (13.0/12.0)*pow(v2 - 2.0*v3 + v4, 2) + 0.25*pow(v2 - v4, 2);
    double b2 = (13.0/12.0)*pow(v3 - 2.0*v4 + v5, 2) + 0.25*pow(3.0*v3 - 4.0*v4 + v5, 2);

    double a0 = 0.1 / pow(WENO_EPS + b0, 2);
    double a1 = 0.6 / pow(WENO_EPS + b1, 2);
    double a2 = 0.3 / pow(WENO_EPS + b2, 2);
    double asum = a0 + a1 + a2;
    double w0 = a0/asum, w1 = a1/asum, w2 = a2/asum;

    double q0 =  v1/3.0 - 7.0/6.0*v2 + 11.0/6.0*v3;
    double q1 = -v2/6.0 + 5.0/6.0*v3 +       v4/3.0;
    double q2 =  v3/3.0 + 5.0/6.0*v4 -       v5/6.0;

    return w0*q0 + w1*q1 + w2*q2;
}

// ============================================================================
// wcR5 — right-biased WENO5 reconstruction at x_{i+1/2}
// ============================================================================
// Three candidate 3rd-order polynomials (sub-stencils):
//   q0 = -1/6*v1 + 5/6*v2 +  1/3*v3    {i-2, i-1, i}
//   q1 =  1/3*v2 + 5/6*v3 -  1/6*v4    {i-1,  i, i+1}
//   q2 = 11/6*v3 - 7/6*v4 +  1/3*v5    { i, i+1, i+2}
//
// Same Jiang-Shu smoothness indicators as wcL5 (symmetric in beta).
// Ideal weights are the mirror of wcL5: d = {3/10, 6/10, 1/10}.
double wcR5(double v1, double v2, double v3, double v4, double v5)
{
    double b0 = (13.0/12.0)*pow(v1 - 2.0*v2 + v3, 2) + 0.25*pow(v1 - 4.0*v2 + 3.0*v3, 2);
    double b1 = (13.0/12.0)*pow(v2 - 2.0*v3 + v4, 2) + 0.25*pow(v2 - v4, 2);
    double b2 = (13.0/12.0)*pow(v3 - 2.0*v4 + v5, 2) + 0.25*pow(3.0*v3 - 4.0*v4 + v5, 2);

    double a0 = 0.3 / pow(WENO_EPS + b0, 2);
    double a1 = 0.6 / pow(WENO_EPS + b1, 2);
    double a2 = 0.1 / pow(WENO_EPS + b2, 2);
    double asum = a0 + a1 + a2;
    double w0 = a0/asum, w1 = a1/asum, w2 = a2/asum;

    double q0 = -v1/6.0 + 5.0/6.0*v2 +        v3/3.0;
    double q1 =  v2/3.0 + 5.0/6.0*v3 -        v4/6.0;
    double q2 = 11.0/6.0*v3 - 7.0/6.0*v4 + v5/3.0;

    return w0*q0 + w1*q1 + w2*q2;
}

// ============================================================================
// wenoL_scalar — reconstruct h^+[i] at all nodes from a F^+ component array
// ============================================================================
// Boundary nodes i=0,1 and i=N-2,N-1 use ghost values because the 5-point
// stencil extends outside the domain there.
void wenoL_scalar(int N, const vector<double>& f,
                  vector<double>& hL,
                  double gLm1, double gLm2,
                  double gRp0, double gRp1)
{
    // i=0: stencil needs f[-2], f[-1], f[0], f[1], f[2]
    hL[0] = wcL5(gLm2, gLm1, f[0], f[1], f[2]);

    // i=1: stencil needs f[-1], f[0], f[1], f[2], f[3]
    hL[1] = wcL5(gLm1, f[0], f[1], f[2], f[3]);

    // Interior: full 5-point stencil available
    for (int i = 2; i < N - 2; i++)
        hL[i] = wcL5(f[i-2], f[i-1], f[i], f[i+1], f[i+2]);

    // i=N-2: stencil needs f[N-4]..f[N-1], f[N]
    hL[N-2] = wcL5(f[N-4], f[N-3], f[N-2], f[N-1], gRp0);

    // i=N-1: stencil needs f[N-3]..f[N-1], f[N], f[N+1]
    hL[N-1] = wcL5(f[N-3], f[N-2], f[N-1], gRp0, gRp1);
}

// ============================================================================
// wenoR_scalar — reconstruct h^-[i] at all nodes from a F^- component array
// ============================================================================
void wenoR_scalar(int N, const vector<double>& f,
                  vector<double>& hR,
                  double gLm1, double gLm2,
                  double gRp0, double gRp1)
{
    hR[0] = wcR5(gLm2, gLm1, f[0], f[1], f[2]);
    hR[1] = wcR5(gLm1, f[0], f[1], f[2], f[3]);

    for (int i = 2; i < N - 2; i++)
        hR[i] = wcR5(f[i-2], f[i-1], f[i], f[i+1], f[i+2]);

    hR[N-2] = wcR5(f[N-4], f[N-3], f[N-2], f[N-1], gRp0);
    hR[N-1] = wcR5(f[N-3], f[N-2], f[N-1], gRp0, gRp1);
}


// ============================================================================
// compute_rhs  [1D] — x-sweep on a flat vector, original single-row solver
// ============================================================================
void compute_rhs(const vector<ConsvState>& U, int N, double dx, double gamma,
                 vector<ConsvState>& rhs_out)
{
    double alpha = 0.0;
    for (int i = 0; i < N; i++)
        alpha = max(alpha, max_wave_speed(U[i], gamma));

    ConsvState gL1{ U[1].rho, -U[1].rhou,  U[1].rhov,  U[1].rhoE };
    ConsvState gL2{ U[2].rho, -U[2].rhou,  U[2].rhov,  U[2].rhoE };
    ConsvState gR0{ U[N-2].rho, -U[N-2].rhou, U[N-2].rhov, U[N-2].rhoE };
    ConsvState gR1{ U[N-3].rho, -U[N-3].rhou, U[N-3].rhov, U[N-3].rhoE };

    Flux FgL1 = compute_flux(gL1, gamma);
    Flux FgL2 = compute_flux(gL2, gamma);
    Flux FgR0 = compute_flux(gR0, gamma);
    Flux FgR1 = compute_flux(gR1, gamma);

    vector<double> fp1(N),fm1(N), fp2(N),fm2(N), fp3(N),fm3(N), fp4(N),fm4(N);
    for (int i = 0; i < N; i++) {
        Flux F = compute_flux(U[i], gamma);
        fp1[i] = 0.5*(F.F1 + alpha*U[i].rho);   fm1[i] = 0.5*(F.F1 - alpha*U[i].rho);
        fp2[i] = 0.5*(F.F2 + alpha*U[i].rhou);  fm2[i] = 0.5*(F.F2 - alpha*U[i].rhou);
        fp3[i] = 0.5*(F.F3 + alpha*U[i].rhov);  fm3[i] = 0.5*(F.F3 - alpha*U[i].rhov);
        fp4[i] = 0.5*(F.F4 + alpha*U[i].rhoE);  fm4[i] = 0.5*(F.F4 - alpha*U[i].rhoE);
    }

    auto sp = [&](const ConsvState& g, const Flux& Fg,
                  double& p1, double& m1, double& p2, double& m2,
                  double& p3, double& m3, double& p4, double& m4) {
        p1=0.5*(Fg.F1+alpha*g.rho);  m1=0.5*(Fg.F1-alpha*g.rho);
        p2=0.5*(Fg.F2+alpha*g.rhou); m2=0.5*(Fg.F2-alpha*g.rhou);
        p3=0.5*(Fg.F3+alpha*g.rhov); m3=0.5*(Fg.F3-alpha*g.rhov);
        p4=0.5*(Fg.F4+alpha*g.rhoE); m4=0.5*(Fg.F4-alpha*g.rhoE);
    };

    double fp1_gL1,fm1_gL1,fp2_gL1,fm2_gL1,fp3_gL1,fm3_gL1,fp4_gL1,fm4_gL1;
    double fp1_gL2,fm1_gL2,fp2_gL2,fm2_gL2,fp3_gL2,fm3_gL2,fp4_gL2,fm4_gL2;
    double fp1_gR0,fm1_gR0,fp2_gR0,fm2_gR0,fp3_gR0,fm3_gR0,fp4_gR0,fm4_gR0;
    double fp1_gR1,fm1_gR1,fp2_gR1,fm2_gR1,fp3_gR1,fm3_gR1,fp4_gR1,fm4_gR1;

    sp(gL1,FgL1, fp1_gL1,fm1_gL1, fp2_gL1,fm2_gL1, fp3_gL1,fm3_gL1, fp4_gL1,fm4_gL1);
    sp(gL2,FgL2, fp1_gL2,fm1_gL2, fp2_gL2,fm2_gL2, fp3_gL2,fm3_gL2, fp4_gL2,fm4_gL2);
    sp(gR0,FgR0, fp1_gR0,fm1_gR0, fp2_gR0,fm2_gR0, fp3_gR0,fm3_gR0, fp4_gR0,fm4_gR0);
    sp(gR1,FgR1, fp1_gR1,fm1_gR1, fp2_gR1,fm2_gR1, fp3_gR1,fm3_gR1, fp4_gR1,fm4_gR1);

    vector<double> hL1(N,0),hR1(N,0), hL2(N,0),hR2(N,0);
    vector<double> hL3(N,0),hR3(N,0), hL4(N,0),hR4(N,0);

    wenoL_scalar(N,fp1,hL1, fp1_gL1,fp1_gL2, fp1_gR0,fp1_gR1);
    wenoR_scalar(N,fm1,hR1, fm1_gL1,fm1_gL2, fm1_gR0,fm1_gR1);
    wenoL_scalar(N,fp2,hL2, fp2_gL1,fp2_gL2, fp2_gR0,fp2_gR1);
    wenoR_scalar(N,fm2,hR2, fm2_gL1,fm2_gL2, fm2_gR0,fm2_gR1);
    wenoL_scalar(N,fp3,hL3, fp3_gL1,fp3_gL2, fp3_gR0,fp3_gR1);
    wenoR_scalar(N,fm3,hR3, fm3_gL1,fm3_gL2, fm3_gR0,fm3_gR1);
    wenoL_scalar(N,fp4,hL4, fp4_gL1,fp4_gL2, fp4_gR0,fp4_gR1);
    wenoR_scalar(N,fm4,hR4, fm4_gL1,fm4_gL2, fm4_gR0,fm4_gR1);

    double hL1_ext = wcL5(fp1_gL2,fp1_gL2,fp1_gL1,fp1[0],fp1[1]);
    double hL2_ext = wcL5(fp2_gL2,fp2_gL2,fp2_gL1,fp2[0],fp2[1]);
    double hL3_ext = wcL5(fp3_gL2,fp3_gL2,fp3_gL1,fp3[0],fp3[1]);
    double hL4_ext = wcL5(fp4_gL2,fp4_gL2,fp4_gL1,fp4[0],fp4[1]);

    double hR1_ext = wcR5(fm1[N-2],fm1[N-1],fm1_gR0,fm1_gR1,fm1_gR1);
    double hR2_ext = wcR5(fm2[N-2],fm2[N-1],fm2_gR0,fm2_gR1,fm2_gR1);
    double hR3_ext = wcR5(fm3[N-2],fm3[N-1],fm3_gR0,fm3_gR1,fm3_gR1);
    double hR4_ext = wcR5(fm4[N-2],fm4[N-1],fm4_gR0,fm4_gR1,fm4_gR1);

    rhs_out[0].rho  = -((hL1[0]-hL1_ext) + (hR1[1]-hR1[0])) / dx;
    rhs_out[0].rhou = -((hL2[0]-hL2_ext) + (hR2[1]-hR2[0])) / dx;
    rhs_out[0].rhov = -((hL3[0]-hL3_ext) + (hR3[1]-hR3[0])) / dx;
    rhs_out[0].rhoE = -((hL4[0]-hL4_ext) + (hR4[1]-hR4[0])) / dx;

    for (int i = 1; i < N-1; i++) {
        rhs_out[i].rho  = -((hL1[i]-hL1[i-1]) + (hR1[i+1]-hR1[i])) / dx;
        rhs_out[i].rhou = -((hL2[i]-hL2[i-1]) + (hR2[i+1]-hR2[i])) / dx;
        rhs_out[i].rhov = -((hL3[i]-hL3[i-1]) + (hR3[i+1]-hR3[i])) / dx;
        rhs_out[i].rhoE = -((hL4[i]-hL4[i-1]) + (hR4[i+1]-hR4[i])) / dx;
    }

    rhs_out[N-1].rho  = -((hL1[N-1]-hL1[N-2]) + (hR1_ext-hR1[N-1])) / dx;
    rhs_out[N-1].rhou = -((hL2[N-1]-hL2[N-2]) + (hR2_ext-hR2[N-1])) / dx;
    rhs_out[N-1].rhov = -((hL3[N-1]-hL3[N-2]) + (hR3_ext-hR3[N-1])) / dx;
    rhs_out[N-1].rhoE = -((hL4[N-1]-hL4[N-2]) + (hR4_ext-hR4[N-1])) / dx;
}



// ============================================================================
// compute_rhs  [2D overload] — x-sweep, one row at a time
// ============================================================================
void compute_rhs(const vector<vector<ConsvState>>& U,
                 const Grid2D& g, double gamma,
                 vector<vector<ConsvState>>& rhs_out)
{
    int Nx = g.Nx, Ny = g.Ny;

    // Global LF alpha: max spectral radius over the entire 2D domain
    double alpha = 0.0;
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
            alpha = max(alpha, max_wave_speed(U[j][i], gamma));

    for (int j = 0; j < Ny; j++)
    {
        const vector<ConsvState>& row = U[j];
        vector<ConsvState>&       rhs = rhs_out[j];

        // ── Ghost states: left/right wall — flip rhou ─────────────────────
        //ConsvState gL1{ row[1].rho, -row[1].rhou, row[1].rhov, row[1].rhoE };
        //ConsvState gL2{ row[2].rho, -row[2].rhou, row[2].rhov, row[2].rhoE };
       // ConsvState gR0{ row[Nx-2].rho, -row[Nx-2].rhou, row[Nx-2].rhov, row[Nx-2].rhoE };
       //ConsvState gR1{ row[Nx-3].rho, -row[Nx-3].rhou, row[Nx-3].rhov, row[Nx-3].rhoE };
        
        ConsvState gL1=U[j][1];
        ConsvState gL2=U[j][2];
        ConsvState gR0=U[j][Nx-1];
        ConsvState gR1=U[j][Nx-1];

        Flux FgL1 = compute_flux(gL1, gamma);
        Flux FgL2 = compute_flux(gL2, gamma);
        Flux FgR0 = compute_flux(gR0, gamma);
        Flux FgR1 = compute_flux(gR1, gamma);

        // ── Split flux arrays ─────────────────────────────────────────────
        vector<double> fp1(Nx), fm1(Nx);
        vector<double> fp2(Nx), fm2(Nx);
        vector<double> fp3(Nx), fm3(Nx);
        vector<double> fp4(Nx), fm4(Nx);

        for (int i = 0; i < Nx; i++) {
            Flux F = compute_flux(row[i], gamma);
            fp1[i] = 0.5*(F.F1 + alpha*row[i].rho);   fm1[i] = 0.5*(F.F1 - alpha*row[i].rho);
            fp2[i] = 0.5*(F.F2 + alpha*row[i].rhou);  fm2[i] = 0.5*(F.F2 - alpha*row[i].rhou);
            fp3[i] = 0.5*(F.F3 + alpha*row[i].rhov);  fm3[i] = 0.5*(F.F3 - alpha*row[i].rhov);
            fp4[i] = 0.5*(F.F4 + alpha*row[i].rhoE);  fm4[i] = 0.5*(F.F4 - alpha*row[i].rhoE);
        }

        // Ghost split fluxes
        auto gfs = [&](const ConsvState& g, const Flux& Fg,
                        double& fp1g, double& fm1g,
                        double& fp2g, double& fm2g,
                        double& fp3g, double& fm3g,
                        double& fp4g, double& fm4g) {
            fp1g = 0.5*(Fg.F1 + alpha*g.rho);   fm1g = 0.5*(Fg.F1 - alpha*g.rho);
            fp2g = 0.5*(Fg.F2 + alpha*g.rhou);  fm2g = 0.5*(Fg.F2 - alpha*g.rhou);
            fp3g = 0.5*(Fg.F3 + alpha*g.rhov);  fm3g = 0.5*(Fg.F3 - alpha*g.rhov);
            fp4g = 0.5*(Fg.F4 + alpha*g.rhoE);  fm4g = 0.5*(Fg.F4 - alpha*g.rhoE);
        };

        double fp1_gL1,fm1_gL1,fp2_gL1,fm2_gL1,fp3_gL1,fm3_gL1,fp4_gL1,fm4_gL1;
        double fp1_gL2,fm1_gL2,fp2_gL2,fm2_gL2,fp3_gL2,fm3_gL2,fp4_gL2,fm4_gL2;
        double fp1_gR0,fm1_gR0,fp2_gR0,fm2_gR0,fp3_gR0,fm3_gR0,fp4_gR0,fm4_gR0;
        double fp1_gR1,fm1_gR1,fp2_gR1,fm2_gR1,fp3_gR1,fm3_gR1,fp4_gR1,fm4_gR1;

        gfs(gL1,FgL1, fp1_gL1,fm1_gL1, fp2_gL1,fm2_gL1, fp3_gL1,fm3_gL1, fp4_gL1,fm4_gL1);
        gfs(gL2,FgL2, fp1_gL2,fm1_gL2, fp2_gL2,fm2_gL2, fp3_gL2,fm3_gL2, fp4_gL2,fm4_gL2);
        gfs(gR0,FgR0, fp1_gR0,fm1_gR0, fp2_gR0,fm2_gR0, fp3_gR0,fm3_gR0, fp4_gR0,fm4_gR0);
        gfs(gR1,FgR1, fp1_gR1,fm1_gR1, fp2_gR1,fm2_gR1, fp3_gR1,fm3_gR1, fp4_gR1,fm4_gR1);

        // ── WENO5 reconstruction ──────────────────────────────────────────
        vector<double> hL1(Nx,0),hR1(Nx,0), hL2(Nx,0),hR2(Nx,0);
        vector<double> hL3(Nx,0),hR3(Nx,0), hL4(Nx,0),hR4(Nx,0);

        wenoL_scalar(Nx,fp1,hL1, fp1_gL1,fp1_gL2, fp1_gR0,fp1_gR1);
        wenoR_scalar(Nx,fm1,hR1, fm1_gL1,fm1_gL2, fm1_gR0,fm1_gR1);
        wenoL_scalar(Nx,fp2,hL2, fp2_gL1,fp2_gL2, fp2_gR0,fp2_gR1);
        wenoR_scalar(Nx,fm2,hR2, fm2_gL1,fm2_gL2, fm2_gR0,fm2_gR1);
        wenoL_scalar(Nx,fp3,hL3, fp3_gL1,fp3_gL2, fp3_gR0,fp3_gR1);
        wenoR_scalar(Nx,fm3,hR3, fm3_gL1,fm3_gL2, fm3_gR0,fm3_gR1);
        wenoL_scalar(Nx,fp4,hL4, fp4_gL1,fp4_gL2, fp4_gR0,fp4_gR1);
        wenoR_scalar(Nx,fm4,hR4, fm4_gL1,fm4_gL2, fm4_gR0,fm4_gR1);

        // ── Outer ghost reconstructions for boundary residuals ────────────
        double hL1_ext = wcL5(fp1_gL2,fp1_gL2,fp1_gL1,fp1[0],fp1[1]);
        double hL2_ext = wcL5(fp2_gL2,fp2_gL2,fp2_gL1,fp2[0],fp2[1]);
        double hL3_ext = wcL5(fp3_gL2,fp3_gL2,fp3_gL1,fp3[0],fp3[1]);
        double hL4_ext = wcL5(fp4_gL2,fp4_gL2,fp4_gL1,fp4[0],fp4[1]);

        double hR1_ext = wcR5(fm1[Nx-2],fm1[Nx-1],fm1_gR0,fm1_gR1,fm1_gR1);
        double hR2_ext = wcR5(fm2[Nx-2],fm2[Nx-1],fm2_gR0,fm2_gR1,fm2_gR1);
        double hR3_ext = wcR5(fm3[Nx-2],fm3[Nx-1],fm3_gR0,fm3_gR1,fm3_gR1);
        double hR4_ext = wcR5(fm4[Nx-2],fm4[Nx-1],fm4_gR0,fm4_gR1,fm4_gR1);

        double dx = g.dx;

        // ── Assemble x-residual (overwrite rhs) ───────────────────────────
        rhs[0].rho  = -((hL1[0]-hL1_ext) + (hR1[1]-hR1[0])) / dx;
        rhs[0].rhou = -((hL2[0]-hL2_ext) + (hR2[1]-hR2[0])) / dx;
        rhs[0].rhov = -((hL3[0]-hL3_ext) + (hR3[1]-hR3[0])) / dx;
        rhs[0].rhoE = -((hL4[0]-hL4_ext) + (hR4[1]-hR4[0])) / dx;

        for (int i = 1; i < Nx-1; i++) {
            rhs[i].rho  = -((hL1[i]-hL1[i-1]) + (hR1[i+1]-hR1[i])) / dx;
            rhs[i].rhou = -((hL2[i]-hL2[i-1]) + (hR2[i+1]-hR2[i])) / dx;
            rhs[i].rhov = -((hL3[i]-hL3[i-1]) + (hR3[i+1]-hR3[i])) / dx;
            rhs[i].rhoE = -((hL4[i]-hL4[i-1]) + (hR4[i+1]-hR4[i])) / dx;
        }

        rhs[Nx-1].rho  = -((hL1[Nx-1]-hL1[Nx-2]) + (hR1_ext-hR1[Nx-1])) / dx;
        rhs[Nx-1].rhou = -((hL2[Nx-1]-hL2[Nx-2]) + (hR2_ext-hR2[Nx-1])) / dx;
        rhs[Nx-1].rhov = -((hL3[Nx-1]-hL3[Nx-2]) + (hR3_ext-hR3[Nx-1])) / dx;
        rhs[Nx-1].rhoE = -((hL4[Nx-1]-hL4[Nx-2]) + (hR4_ext-hR4[Nx-1])) / dx;
    }
}

// ============================================================================
// compute_rhs_y — y-sweep: Ly(U) = -dG/dy, accumulated into rhs_out
// ============================================================================
// Operates column-by-column (fixed i, varying j).
// Uses compute_flux_y() which returns G(U) = {rhov, cross, rhov*v+p, v*(rhoE+p)}.
// Ghost states at bottom (j=-1,-2) and top (j=Ny, Ny+1) walls flip rhov.
void compute_rhs_y(const vector<vector<ConsvState>>& U,
                   const Grid2D& g, double gamma,
                   vector<vector<ConsvState>>& rhs_out)
{
    int Nx = g.Nx, Ny = g.Ny;
    double dy = g.dy;

    // Global LF alpha for y-direction (max |v| + a over domain)
    double alpha = 0.0;
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            double v = U[j][i].rhov / U[j][i].rho;
            double a = sound_speed(U[j][i], gamma);
            alpha = max(alpha, fabs(v) + a);
        }

    for (int i = 0; i < Nx; i++)
    {
        // ── Extract column i into a flat array for the WENO sweep ─────────
        vector<ConsvState> col(Ny);
        for (int j = 0; j < Ny; j++) col[j] = U[j][i];

        // ── Ghost states: bottom/top walls — flip rhov ────────────────────
        // Bottom wall ghosts mirror rows 1 and 2 with rhov negated
        ConsvState gB1{ col[1].rho,  col[1].rhou, -col[1].rhov,  col[1].rhoE };
        ConsvState gB2{ col[2].rho,  col[2].rhou, -col[2].rhov,  col[2].rhoE };
        // Top wall ghosts mirror rows Ny-2 and Ny-3 with rhov negated
        ConsvState gT0{ col[Ny-2].rho, col[Ny-2].rhou, -col[Ny-2].rhov, col[Ny-2].rhoE };
        ConsvState gT1{ col[Ny-3].rho, col[Ny-3].rhou, -col[Ny-3].rhov, col[Ny-3].rhoE };

        Flux FgB1 = compute_flux_y(gB1, gamma);
        Flux FgB2 = compute_flux_y(gB2, gamma);
        Flux FgT0 = compute_flux_y(gT0, gamma);
        Flux FgT1 = compute_flux_y(gT1, gamma);

        // ── Split flux arrays (using G flux components) ───────────────────
        // G = {G1=rhov, G2=cross, G3=rhov*v+p, G4=v*(rhoE+p)}
        vector<double> fp1(Ny),fm1(Ny), fp2(Ny),fm2(Ny);
        vector<double> fp3(Ny),fm3(Ny), fp4(Ny),fm4(Ny);

        for (int j = 0; j < Ny; j++) {
            Flux G = compute_flux_y(col[j], gamma);
            fp1[j] = 0.5*(G.F1 + alpha*col[j].rho);   fm1[j] = 0.5*(G.F1 - alpha*col[j].rho);
            fp2[j] = 0.5*(G.F2 + alpha*col[j].rhou);  fm2[j] = 0.5*(G.F2 - alpha*col[j].rhou);
            fp3[j] = 0.5*(G.F3 + alpha*col[j].rhov);  fm3[j] = 0.5*(G.F3 - alpha*col[j].rhov);
            fp4[j] = 0.5*(G.F4 + alpha*col[j].rhoE);  fm4[j] = 0.5*(G.F4 - alpha*col[j].rhoE);
        }

        // Ghost split fluxes (bottom = "left", top = "right" in the sweep)
        auto sf = [&](const ConsvState& g, const Flux& Fg,
                      double& p1, double& m1, double& p2, double& m2,
                      double& p3, double& m3, double& p4, double& m4) {
            p1 = 0.5*(Fg.F1+alpha*g.rho);   m1 = 0.5*(Fg.F1-alpha*g.rho);
            p2 = 0.5*(Fg.F2+alpha*g.rhou);  m2 = 0.5*(Fg.F2-alpha*g.rhou);
            p3 = 0.5*(Fg.F3+alpha*g.rhov);  m3 = 0.5*(Fg.F3-alpha*g.rhov);
            p4 = 0.5*(Fg.F4+alpha*g.rhoE);  m4 = 0.5*(Fg.F4-alpha*g.rhoE);
        };

        double fp1_gB1,fm1_gB1,fp2_gB1,fm2_gB1,fp3_gB1,fm3_gB1,fp4_gB1,fm4_gB1;
        double fp1_gB2,fm1_gB2,fp2_gB2,fm2_gB2,fp3_gB2,fm3_gB2,fp4_gB2,fm4_gB2;
        double fp1_gT0,fm1_gT0,fp2_gT0,fm2_gT0,fp3_gT0,fm3_gT0,fp4_gT0,fm4_gT0;
        double fp1_gT1,fm1_gT1,fp2_gT1,fm2_gT1,fp3_gT1,fm3_gT1,fp4_gT1,fm4_gT1;

        sf(gB1,FgB1, fp1_gB1,fm1_gB1, fp2_gB1,fm2_gB1, fp3_gB1,fm3_gB1, fp4_gB1,fm4_gB1);
        sf(gB2,FgB2, fp1_gB2,fm1_gB2, fp2_gB2,fm2_gB2, fp3_gB2,fm3_gB2, fp4_gB2,fm4_gB2);
        sf(gT0,FgT0, fp1_gT0,fm1_gT0, fp2_gT0,fm2_gT0, fp3_gT0,fm3_gT0, fp4_gT0,fm4_gT0);
        sf(gT1,FgT1, fp1_gT1,fm1_gT1, fp2_gT1,fm2_gT1, fp3_gT1,fm3_gT1, fp4_gT1,fm4_gT1);

        // ── WENO5 reconstruction along the column ─────────────────────────
        vector<double> hL1(Ny,0),hR1(Ny,0), hL2(Ny,0),hR2(Ny,0);
        vector<double> hL3(Ny,0),hR3(Ny,0), hL4(Ny,0),hR4(Ny,0);

        // bottom ghost = "left ghost" (gLm1=gB1, gLm2=gB2)
        // top ghost    = "right ghost" (gRp0=gT0, gRp1=gT1)
        wenoL_scalar(Ny,fp1,hL1, fp1_gB1,fp1_gB2, fp1_gT0,fp1_gT1);
        wenoR_scalar(Ny,fm1,hR1, fm1_gB1,fm1_gB2, fm1_gT0,fm1_gT1);
        wenoL_scalar(Ny,fp2,hL2, fp2_gB1,fp2_gB2, fp2_gT0,fp2_gT1);
        wenoR_scalar(Ny,fm2,hR2, fm2_gB1,fm2_gB2, fm2_gT0,fm2_gT1);
        wenoL_scalar(Ny,fp3,hL3, fp3_gB1,fp3_gB2, fp3_gT0,fp3_gT1);
        wenoR_scalar(Ny,fm3,hR3, fm3_gB1,fm3_gB2, fm3_gT0,fm3_gT1);
        wenoL_scalar(Ny,fp4,hL4, fp4_gB1,fp4_gB2, fp4_gT0,fp4_gT1);
        wenoR_scalar(Ny,fm4,hR4, fm4_gB1,fm4_gB2, fm4_gT0,fm4_gT1);

        // ── Outer ghost reconstructions for boundary residuals ─────────────
        double hL1_ext = wcL5(fp1_gB2,fp1_gB2,fp1_gB1,fp1[0],fp1[1]);
        double hL2_ext = wcL5(fp2_gB2,fp2_gB2,fp2_gB1,fp2[0],fp2[1]);
        double hL3_ext = wcL5(fp3_gB2,fp3_gB2,fp3_gB1,fp3[0],fp3[1]);
        double hL4_ext = wcL5(fp4_gB2,fp4_gB2,fp4_gB1,fp4[0],fp4[1]);

        double hR1_ext = wcR5(fm1[Ny-2],fm1[Ny-1],fm1_gT0,fm1_gT1,fm1_gT1);
        double hR2_ext = wcR5(fm2[Ny-2],fm2[Ny-1],fm2_gT0,fm2_gT1,fm2_gT1);
        double hR3_ext = wcR5(fm3[Ny-2],fm3[Ny-1],fm3_gT0,fm3_gT1,fm3_gT1);
        double hR4_ext = wcR5(fm4[Ny-2],fm4[Ny-1],fm4_gT0,fm4_gT1,fm4_gT1);

        // ── Accumulate y-residual into rhs_out ────────────────────────────
        rhs_out[0][i].rho  += -((hL1[0]-hL1_ext) + (hR1[1]-hR1[0])) / dy;
        rhs_out[0][i].rhou += -((hL2[0]-hL2_ext) + (hR2[1]-hR2[0])) / dy;
        rhs_out[0][i].rhov += -((hL3[0]-hL3_ext) + (hR3[1]-hR3[0])) / dy;
        rhs_out[0][i].rhoE += -((hL4[0]-hL4_ext) + (hR4[1]-hR4[0])) / dy;

        for (int j = 1; j < Ny-1; j++) {
            rhs_out[j][i].rho  += -((hL1[j]-hL1[j-1]) + (hR1[j+1]-hR1[j])) / dy;
            rhs_out[j][i].rhou += -((hL2[j]-hL2[j-1]) + (hR2[j+1]-hR2[j])) / dy;
            rhs_out[j][i].rhov += -((hL3[j]-hL3[j-1]) + (hR3[j+1]-hR3[j])) / dy;
            rhs_out[j][i].rhoE += -((hL4[j]-hL4[j-1]) + (hR4[j+1]-hR4[j])) / dy;
        }

        rhs_out[Ny-1][i].rho  += -((hL1[Ny-1]-hL1[Ny-2]) + (hR1_ext-hR1[Ny-1])) / dy;
        rhs_out[Ny-1][i].rhou += -((hL2[Ny-1]-hL2[Ny-2]) + (hR2_ext-hR2[Ny-1])) / dy;
        rhs_out[Ny-1][i].rhov += -((hL3[Ny-1]-hL3[Ny-2]) + (hR3_ext-hR3[Ny-1])) / dy;
        rhs_out[Ny-1][i].rhoE += -((hL4[Ny-1]-hL4[Ny-2]) + (hR4_ext-hR4[Ny-1])) / dy;
    }
}