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
    double rhou;  // U2: momentum
    double rhoE;  // U3: total energy
};

// ── Physical flux vector: F(U) = [F_mass, F_momentum, F_energy] ─────────────
struct Flux {
    double F_mass;      // F1 = rho*u
    double F_momentum;  // F2 = rho*u^2 + p
    double F_energy;    // F3 = u*(rho*E + p)
};

// Compute primitive variables and flux from a conservative state
static inline double pressure(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double E = U.rhoE  / U.rho;
    return (gamma - 1.0) * U.rho * (E - 0.5*u*u);
}

static inline Flux compute_flux(const ConsvState& U, double gamma)
{
    double u = U.rhou / U.rho;
    double p = pressure(U, gamma);
    return { U.rhou,
             U.rhou*u + p,
             u*(U.rhoE + p) };
}

static inline double sound_speed(const ConsvState& U, double gamma)
{
    return sqrt(gamma * pressure(U, gamma) / U.rho);
}

// ── HTML output ──────────────────────────────────────────────────────────────
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
    fprintf(f, "<h2>1D Shock Tube &nbsp;|&nbsp; t = %.5f &nbsp;|&nbsp; step = %d</h2>\n", time, step);
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

void writeData(const vector<double>& x_arr,
               const vector<ConsvState>& U,
               double gamma,
               double Rgas)
{
    ofstream file("shocktube.dat");

    file << "# x rho u T\n";

    for (int i = 0; i < x_arr.size(); i++)
    {
        double rho = U[i].rho;
        double u   = U[i].rhou / U[i].rho;
        double p   = pressure(U[i], gamma);
        double T   = p / (rho * Rgas);

        file << x_arr[i] << " "
             << rho      << " "
             << u        << " "
             << T        << "\n";
    }

    file.close();
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main()
{
    const int    N          = 400;
    const double x_min      = -0.5;
    const double x_max      =  0.5;
    const double dx         = (x_max - x_min) / (N - 1);
    const double gamma      = 1.4;
    const double Rgas       = 1.0;
    const double CFL        = 0.99;
    const double final_time = 0.2;

    // ── Grid ─────────────────────────────────────────────────────────────────
    vector<double> x_arr(N);
    for (int i = 0; i < N; i++) x_arr[i] = x_min + i*dx;

    // ── Conservative state vector U and its update buffer ────────────────────
    vector<ConsvState> U(N), U_new(N);

    // ── Initial conditions (Sod shock tube) ──────────────────────────────────
    for (int i = 0; i < N; i++) {
        double rho_i, u_i, p_i;
        if (x_arr[i] < 0.0) { rho_i = 1.0;   u_i = 0.0; p_i = 1.0; }
        else                 { rho_i = 0.125; u_i = 0.0; p_i = 0.1; }
        double E_i  = p_i / ((gamma - 1.0) * rho_i) + 0.5*u_i*u_i;
        U[i] = { rho_i, rho_i*u_i, rho_i*E_i };
    }

    double time = 0.0;
    int    step = 0;

    write_html(x_arr, U, time, step, gamma, Rgas, false);

    char full_path[MAX_PATH];
    GetFullPathNameA("shock_tube.html", MAX_PATH, full_path, nullptr);
    ShellExecuteA(nullptr, "open", full_path, nullptr, nullptr, SW_SHOWNORMAL);
    Sleep(3000);
    cout << "Browser opened. Simulating..." << endl;

    // ── Time integration (Lax-Friedrichs with scalar viscosity) ──────────────
    while (time < final_time)
    {
        if (isnan(U[N/2].rho) || isinf(U[N/2].rho)) {
            cout << "DIVERGED at step " << step << " t=" << time << endl;
            break;
        }

        // CFL time step
        double max_speed = 0.0;
        for (int i = 0; i < N; i++)
            max_speed = max(max_speed, fabs(U[i].rhou / U[i].rho) + sound_speed(U[i], gamma));
        double dt = CFL * dx / max_speed;
        if (time + dt > final_time) dt = final_time - time;

        // Precompute flux vector F(U) at every node
        vector<Flux> F(N);
        for (int i = 0; i < N; i++)
            F[i] = compute_flux(U[i], gamma);

        // Interior update:  U_new = U - dt/dx * (F_{i+1} - F_{i-1})/2  +  dt * visc
        for (int i = 1; i < N - 1; i++) {
            double alpha = fabs(U[i].rhou / U[i].rho) + sound_speed(U[i], gamma);

            // Central-difference flux divergence
            double dF_mass_dx      = (F[i+1].F_mass      - F[i-1].F_mass)      / (2.0*dx);
            double dF_momentum_dx  = (F[i+1].F_momentum  - F[i-1].F_momentum)  / (2.0*dx);
            double dF_energy_dx    = (F[i+1].F_energy     - F[i-1].F_energy)    / (2.0*dx);

            // Lax-Friedrichs scalar viscosity (stabilisation): alpha * d²U/dx²
            double visc_mass      = alpha * (U[i+1].rho  - 2.0*U[i].rho  + U[i-1].rho)  / (dx*dx);
            double visc_momentum  = alpha * (U[i+1].rhou - 2.0*U[i].rhou + U[i-1].rhou) / (dx*dx);
            double visc_energy    = alpha * (U[i+1].rhoE - 2.0*U[i].rhoE + U[i-1].rhoE) / (dx*dx);

            U_new[i].rho  = U[i].rho  - dt*dF_mass_dx     + dt*0.5*visc_mass;
            U_new[i].rhou = U[i].rhou - dt*dF_momentum_dx + dt*0.5*visc_momentum;
            U_new[i].rhoE = U[i].rhoE - dt*dF_energy_dx   + dt*0.5*visc_energy;
        }

        // Boundary conditions (zero-gradient / transmissive)
        U_new[0]     = U[0];
        U_new[N - 1] = U[N - 1];

        U = U_new;
        time += dt;
        step++;

        if (step % 3 == 0 || time >= final_time)
            write_html(x_arr, U, time, step, gamma, Rgas, time >= final_time);
    }
    
    writeData(x_arr, U, gamma, Rgas);
    cout << "Finished at t = " << time << " (" << step << " steps)" << endl;
    return 0;
}