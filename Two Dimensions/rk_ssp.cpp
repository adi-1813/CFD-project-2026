#include "rk_ssp.h"
#include "weno5.h"
#include "geometry.h"
#include <algorithm>

using namespace std;

// ============================================================================
// ssp_rk3_step  [1D overload — unchanged]
// ============================================================================
double ssp_rk3_step(vector<ConsvState>& U, int N, double dx, double gamma,
                    double CFL, double time, double final_time)
{
    double max_speed = 0.0;
    for (int i = 0; i < N; i++)
        max_speed = max(max_speed, max_wave_speed(U[i], gamma));
    double dt = CFL * dx / max_speed;
    if (time + dt > final_time) dt = final_time - time;

    vector<ConsvState> L0(N), L1(N), L2(N);
    vector<ConsvState> U1(N), U2(N);

    compute_rhs(U, N, dx, gamma, L0);
    for (int i = 0; i < N; i++) {
        U1[i].rho  = U[i].rho  + dt*L0[i].rho;
        U1[i].rhou = U[i].rhou + dt*L0[i].rhou;
        U1[i].rhov = U[i].rhov + dt*L0[i].rhov;
        U1[i].rhoE = U[i].rhoE + dt*L0[i].rhoE;
    }
    apply_wall_bc(U1, N);

    compute_rhs(U1, N, dx, gamma, L1);
    for (int i = 0; i < N; i++) {
        U2[i].rho  = 0.75*U[i].rho  + 0.25*(U1[i].rho  + dt*L1[i].rho);
        U2[i].rhou = 0.75*U[i].rhou + 0.25*(U1[i].rhou + dt*L1[i].rhou);
        U2[i].rhov = 0.75*U[i].rhov + 0.25*(U1[i].rhov + dt*L1[i].rhov);
        U2[i].rhoE = 0.75*U[i].rhoE + 0.25*(U1[i].rhoE + dt*L1[i].rhoE);
    }
    apply_wall_bc(U2, N);

    compute_rhs(U2, N, dx, gamma, L2);
    for (int i = 0; i < N; i++) {
        U[i].rho  = (1.0/3.0)*U[i].rho  + (2.0/3.0)*(U2[i].rho  + dt*L2[i].rho);
        U[i].rhou = (1.0/3.0)*U[i].rhou + (2.0/3.0)*(U2[i].rhou + dt*L2[i].rhou);
        U[i].rhov = (1.0/3.0)*U[i].rhov + (2.0/3.0)*(U2[i].rhov + dt*L2[i].rhov);
        U[i].rhoE = (1.0/3.0)*U[i].rhoE + (2.0/3.0)*(U2[i].rhoE + dt*L2[i].rhoE);
    }
    apply_wall_bc(U, N);

    return dt;
}

// ============================================================================
// ssp_rk3_step  [2D overload]
// ============================================================================
double ssp_rk3_step(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma,
                    double CFL, double time, double final_time, double Mach, double theta, const WallCurve& wedge)
{
    int Nx = g.Nx, Ny = g.Ny;

    double max_speed = 0.0;
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
            max_speed = max(max_speed, max_wave_speed(U[j][i], gamma));
    double dt = CFL * min(g.dx, g.dy) / max_speed;
    if (time + dt > final_time) dt = final_time - time;

    auto make2D = [&]{ return vector<vector<ConsvState>>(Ny, vector<ConsvState>(Nx)); };

    vector<vector<ConsvState>> U1 = make2D();
    vector<vector<ConsvState>> U2 = make2D();
    vector<vector<ConsvState>> L(Ny, vector<ConsvState>(Nx, {0,0,0,0}));

    auto compute_L = [&](const vector<vector<ConsvState>>& Uin,
                          vector<vector<ConsvState>>& Lout) {
        for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
                Lout[j][i] = {0,0,0,0};
        compute_rhs  (Uin, g, gamma, Lout);
        compute_rhs_y(Uin, g, gamma, Lout);
    };

    // Stage 1
    compute_L(U, L);
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            U1[j][i].rho  = U[j][i].rho  + dt*L[j][i].rho;
            U1[j][i].rhou = U[j][i].rhou + dt*L[j][i].rhou;
            U1[j][i].rhov = U[j][i].rhov + dt*L[j][i].rhov;
            U1[j][i].rhoE = U[j][i].rhoE + dt*L[j][i].rhoE;
        }
    apply_inflow_bc(U, g, gamma, Mach, theta);
    apply_outflow_bc(U, g,gamma);
    apply_wall_bc(U, g);
    apply_immersed_wall_bc(U, g, wedge, gamma);
    // Stage 2
    compute_L(U1, L);
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            U2[j][i].rho  = 0.75*U[j][i].rho  + 0.25*(U1[j][i].rho  + dt*L[j][i].rho);
            U2[j][i].rhou = 0.75*U[j][i].rhou + 0.25*(U1[j][i].rhou + dt*L[j][i].rhou);
            U2[j][i].rhov = 0.75*U[j][i].rhov + 0.25*(U1[j][i].rhov + dt*L[j][i].rhov);
            U2[j][i].rhoE = 0.75*U[j][i].rhoE + 0.25*(U1[j][i].rhoE + dt*L[j][i].rhoE);
        }
    apply_inflow_bc(U, g, gamma, Mach, theta);
    apply_outflow_bc(U, g,gamma);
    apply_wall_bc(U, g);
    apply_immersed_wall_bc(U, g, wedge, gamma);
    // Stage 3
    compute_L(U2, L);
    for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++) {
            U[j][i].rho  = (1.0/3.0)*U[j][i].rho  + (2.0/3.0)*(U2[j][i].rho  + dt*L[j][i].rho);
            U[j][i].rhou = (1.0/3.0)*U[j][i].rhou + (2.0/3.0)*(U2[j][i].rhou + dt*L[j][i].rhou);
            U[j][i].rhov = (1.0/3.0)*U[j][i].rhov + (2.0/3.0)*(U2[j][i].rhov + dt*L[j][i].rhov);
            U[j][i].rhoE = (1.0/3.0)*U[j][i].rhoE + (2.0/3.0)*(U2[j][i].rhoE + dt*L[j][i].rhoE);
        }
    apply_inflow_bc(U, g, gamma, Mach, theta);
    apply_outflow_bc(U, g,gamma);
    apply_wall_bc(U, g);
    apply_immersed_wall_bc(U, g, wedge, gamma);
    return dt;
}