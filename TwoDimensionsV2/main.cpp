#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <windows.h>
#pragma comment(lib, "Shell32.lib")
#include<omp.h>
#include "geometry.h"
#include "weno5.h"
#include "rk_ssp.h"

using namespace std;

// ============================================================================
// HTML live-plot: 2D density heatmap + 1D slice plots
// ============================================================================
// The heatmap uses a Canvas ImageData pixel buffer (one pixel per grid node).
// Three line charts below it show density, x-velocity, and pressure along
// the horizontal centreline (j = Ny/2).
//
// The page auto-reloads every 80 ms while the simulation is running.
void write_html(const Grid2D& g,
                const vector<vector<ConsvState>>& U,
                double time, int step, double gamma, bool finished)
{
    int Nx = g.Nx, Ny = g.Ny;
    int jmid = Ny / 2;

    FILE* f = fopen("shock_tube_tmp.html", "w");
    if (!f) { cerr << "ERROR: Cannot write HTML\n"; return; }

    // Gather primitive variables for the centreline plots
    vector<double> rho_mid(Nx), u_mid(Nx), p_mid(Nx);
    for (int i = 0; i < Nx; i++) {
        rho_mid[i] = U[jmid][i].rho;
        u_mid[i]   = U[jmid][i].rhou / U[jmid][i].rho;
        p_mid[i]   = pressure(U[jmid][i], gamma);
    }

    // Flatten the 2D density field for the heatmap
    // We pass it as a flat JSON array, row-major (j=0 is bottom)
    fprintf(f, "<!DOCTYPE html><html><head>\n");
    fprintf(f, "<style>\n");
    fprintf(f, "  body{background:#1a1a2e;color:#eee;font-family:monospace;margin:16px;}\n");
    fprintf(f, "  h2{color:#4e9af1;margin-bottom:6px;}\n");
    fprintf(f, "  .status{font-size:13px;margin-bottom:10px;color:%s;}\n",
            finished ? "#4ef1a4" : "#f1a44e");
    fprintf(f, "  canvas{display:block;background:#0f1923;border:1px solid #2a2a4a;margin-bottom:10px;}\n");
    fprintf(f, "</style></head><body>\n");
    fprintf(f, "<h2>2D Shock Tube — WENO5+SSP-RK3 "
               "&nbsp;|&nbsp; t = %.5f &nbsp;|&nbsp; step = %d "
               "&nbsp;|&nbsp; %dx%d</h2>\n", time, step, Nx, Ny);
    fprintf(f, "<div class='status'>%s</div>\n",
            finished ? "&#10003; Simulation complete." : "&#9654; Running...");

    // Heatmap canvas (scaled up for visibility: each cell = 2px)
    int scale = max(1, 600 / Nx);
    fprintf(f, "<canvas id='hmap' width='%d' height='%d'></canvas>\n",
            Nx*scale, Ny*scale);
    fprintf(f, "<canvas id='c1' width='900' height='160'></canvas>\n");
    fprintf(f, "<canvas id='c2' width='900' height='160'></canvas>\n");
    fprintf(f, "<canvas id='c3' width='900' height='160'></canvas>\n");

    fprintf(f, "<script>\n");

    // Pass grid x-coords for line plots
    fprintf(f, "var xd=[");
    for (int i=0;i<Nx;i++) fprintf(f,"%.5f%s",g.x[i],i<Nx-1?",":"");
    fprintf(f,"];\n");

    // Centreline data
    fprintf(f, "var drho=[");
    for (int i=0;i<Nx;i++) fprintf(f,"%.5f%s",rho_mid[i],i<Nx-1?",":"");
    fprintf(f,"];\n");
    fprintf(f, "var du=[");
    for (int i=0;i<Nx;i++) fprintf(f,"%.5f%s",u_mid[i],i<Nx-1?",":"");
    fprintf(f,"];\n");
    fprintf(f, "var dp=[");
    for (int i=0;i<Nx;i++) fprintf(f,"%.5f%s",p_mid[i],i<Nx-1?",":"");
    fprintf(f,"];\n");

    // Full 2D density for heatmap
    fprintf(f, "var rho2d=[");
    for (int j=0;j<Ny;j++)
        for (int i=0;i<Nx;i++)
            fprintf(f,"%.4f%s",U[j][i].rho,(j==Ny-1&&i==Nx-1)?"":",");
    fprintf(f,"];\n");
    fprintf(f, "var NX=%d,NY=%d,SC=%d;\n", Nx, Ny, scale);

    // --- Heatmap renderer (blue→white plasma-style colourmap) ---
    fprintf(f, "function drawHmap(){\n");
    fprintf(f, "  var cv=document.getElementById('hmap');\n");
    fprintf(f, "  var ctx=cv.getContext('2d');\n");
    fprintf(f, "  var img=ctx.createImageData(NX*SC,NY*SC);\n");
    // Find global min/max density for normalisation
    fprintf(f, "  var mn=rho2d[0],mx=rho2d[0];\n");
    fprintf(f, "  for(var k=1;k<rho2d.length;k++){mn=Math.min(mn,rho2d[k]);mx=Math.max(mx,rho2d[k]);}\n");
    fprintf(f, "  var rng=mx-mn||1;\n");
    // Jet-like colourmap: t in [0,1] -> RGB
    fprintf(f, "  function jet(t){\n");
    fprintf(f, "    var r=Math.min(1,Math.max(0, 1.5-Math.abs(4*t-3)));\n");
    fprintf(f, "    var g=Math.min(1,Math.max(0, 1.5-Math.abs(4*t-2)));\n");
    fprintf(f, "    var b=Math.min(1,Math.max(0, 1.5-Math.abs(4*t-1)));\n");
    fprintf(f, "    return [Math.round(r*255),Math.round(g*255),Math.round(b*255)];\n");
    fprintf(f, "  }\n");
    // Fill pixels (flip j so j=0=bottom)
    fprintf(f, "  for(var j=0;j<NY;j++){\n");
    fprintf(f, "    for(var i=0;i<NX;i++){\n");
    fprintf(f, "      var t=(rho2d[(NY-1-j)*NX+i]-mn)/rng;\n");
    fprintf(f, "      var c=jet(t);\n");
    fprintf(f, "      for(var sy=0;sy<SC;sy++) for(var sx=0;sx<SC;sx++){\n");
    fprintf(f, "        var idx=4*((j*SC+sy)*NX*SC+(i*SC+sx));\n");
    fprintf(f, "        img.data[idx]=c[0];img.data[idx+1]=c[1];\n");
    fprintf(f, "        img.data[idx+2]=c[2];img.data[idx+3]=255;\n");
    fprintf(f, "      }\n");
    fprintf(f, "    }\n");
    fprintf(f, "  }\n");
    fprintf(f, "  ctx.putImageData(img,0,0);\n");
    // Axis labels on heatmap
    fprintf(f, "  ctx.fillStyle='#aabbcc';ctx.font='12px monospace';\n");
    fprintf(f, "  ctx.fillText('Density (2D)',4,14);\n");
    fprintf(f, "  ctx.fillStyle='#4ef1a4';ctx.fillText('t=%.4f',NX*SC-70,14);\n", time);
    fprintf(f, "}\n");
    fprintf(f, "drawHmap();\n");

    // --- Line chart renderer (same as 1D version) ---
    fprintf(f, "function drawChart(id,data,label,color){\n");
    fprintf(f, "  var cv=document.getElementById(id),ctx=cv.getContext('2d');\n");
    fprintf(f, "  var W=cv.width,H=cv.height,pl=60,pr=15,pt=20,pb=30;\n");
    fprintf(f, "  var pw=W-pl-pr,ph=H-pt-pb;\n");
    fprintf(f, "  ctx.fillStyle='#0f1923';ctx.fillRect(0,0,W,H);\n");
    fprintf(f, "  var ymin=data.reduce((a,b)=>Math.min(a,b),Infinity);\n");
    fprintf(f, "  var ymax=data.reduce((a,b)=>Math.max(a,b),-Infinity);\n");
    fprintf(f, "  var yp=(ymax-ymin)*0.08||0.05;ymin-=yp;ymax+=yp;\n");
    fprintf(f, "  var xmin=xd[0],xmax=xd[xd.length-1];\n");
    fprintf(f, "  function sx(x){return pl+(x-xmin)/(xmax-xmin)*pw;}\n");
    fprintf(f, "  function sy(y){return pt+ph-(y-ymin)/(ymax-ymin)*ph;}\n");
    fprintf(f, "  ctx.strokeStyle='#1e2d3d';ctx.lineWidth=1;\n");
    fprintf(f, "  for(var i=0;i<=5;i++){\n");
    fprintf(f, "    var gy=ymin+i*(ymax-ymin)/5;\n");
    fprintf(f, "    ctx.beginPath();ctx.moveTo(pl,sy(gy));ctx.lineTo(pl+pw,sy(gy));ctx.stroke();\n");
    fprintf(f, "  }\n");
    fprintf(f, "  ctx.strokeStyle='#4a5a6a';ctx.lineWidth=1.5;ctx.strokeRect(pl,pt,pw,ph);\n");
    fprintf(f, "  ctx.fillStyle='#8899aa';ctx.font='10px monospace';\n");
    fprintf(f, "  ctx.textAlign='center';\n");
    fprintf(f, "  for(var i=0;i<=5;i++){var gx=xmin+i*(xmax-xmin)/5;ctx.fillText(gx.toFixed(2),sx(gx),pt+ph+16);}\n");
    fprintf(f, "  ctx.textAlign='right';\n");
    fprintf(f, "  for(var i=0;i<=5;i++){var gy=ymin+i*(ymax-ymin)/5;ctx.fillText(gy.toFixed(3),pl-4,sy(gy)+4);}\n");
    fprintf(f, "  ctx.strokeStyle=color;ctx.lineWidth=2;ctx.beginPath();\n");
    fprintf(f, "  for(var i=0;i<xd.length;i++){\n");
    fprintf(f, "    if(i===0)ctx.moveTo(sx(xd[i]),sy(data[i]));\n");
    fprintf(f, "    else ctx.lineTo(sx(xd[i]),sy(data[i]));\n");
    fprintf(f, "  }\n");
    fprintf(f, "  ctx.stroke();\n");
    fprintf(f, "  ctx.fillStyle=color;ctx.font='bold 12px monospace';ctx.textAlign='left';\n");
    fprintf(f, "  ctx.fillText(label+' (y=mid)',pl+4,pt+14);\n");
    fprintf(f, "}\n");

    fprintf(f, "drawChart('c1',drho,'Density',   '#4e9af1');\n");
    fprintf(f, "drawChart('c2',du,  'Velocity-x','#f1a44e');\n");
    fprintf(f, "drawChart('c3',dp,  'Pressure',  '#4ef1a4');\n");

    if (!finished)
        fprintf(f, "setTimeout(function(){location.reload();},80);\n");

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
    // ── Parameters ───────────────────────────────────────────────────────────
    const int    Nx         = 400;      // x grid points
    const int    Ny         = 400;      // y grid points (square domain)
    const double x_min      = 0.0;
    const double x_max      = 3.0;
    const double y_min      = 0.0;
    const double y_max      = 2.0;
    const double gamma      = 1.4;
    const double Rgas       = 1.0;
    const double CFL        = 0.4;     // slightly conservative for 2D
    const double final_time = 5.0;
    const double M_PI       = 3.14159265358979323846;
    // For wedge flow
    const double Mach=2.0;
    const double theta=15.0*M_PI/180.0;  // wedge angle
    omp_set_num_threads(4);

    // ── Grid and state ────────────────────────────────────────────────────────
    Grid2D g(Nx, x_min, x_max, Ny, y_min, y_max);
    vector<vector<ConsvState>> U(Ny, vector<ConsvState>(Nx));

    WallCurve wedge(x_min+0.5, x_max, Nx-1, [&](double x){ return tan(theta)*(x - (x_min+0.5)); });

    set_ic_wedge(U, g, gamma, Mach, theta);
    apply_inflow_bc(U, g, gamma, Mach, theta);
    apply_outflow_bc(U, g,gamma);
    apply_wall_bc(U, g);
    apply_immersed_wall_bc(U, g, wedge, gamma);

    // ── Open browser ──────────────────────────────────────────────────────────
    double time = 0.0;
    int    step = 0;
    write_html(g, U, time, step, gamma, false);
    char full_path[MAX_PATH];
    GetFullPathNameA("shock_tube.html", MAX_PATH, full_path, nullptr);
    ShellExecuteA(nullptr, "open", full_path, nullptr, nullptr, SW_SHOWNORMAL);
    Sleep(2000);
    cout << "Browser opened.  Scheme: FDM-WENO5 + SSP-RK3 + Wall BC  ["
         << Nx << "x" << Ny << "]" << endl;

    auto make2D = [&]{
    return vector<vector<ConsvState>>(Ny, vector<ConsvState>(Nx, {0,0,0,0}));
};
vector<vector<ConsvState>> U1 = make2D();
vector<vector<ConsvState>> U2 = make2D();
vector<vector<ConsvState>> L  = make2D();

    // ── Solver loop ───────────────────────────────────────────────────────────
    while (time < final_time)
    {
        if (isnan(U[Ny/2][Nx/2].rho) || isinf(U[Ny/2][Nx/2].rho)) {
            cout << "DIVERGED at step " << step << " t=" << time << endl;
            break;
        }


        double dt = ssp_rk3_step(U,U1,U2,L,g, gamma, CFL, time, final_time, Mach, theta, wedge);
        time += dt;
        step++;

        if (step % 5 == 0 || time >= final_time)
            write_html(g, U, time, step, gamma, time >= final_time);
    }

    cout << "Finished at t = " << time << " (" << step << " steps)" << endl;
    return 0;
}