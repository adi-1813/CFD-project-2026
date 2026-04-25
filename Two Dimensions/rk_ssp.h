#pragma once
#include <vector>
#include "geometry.h"

using namespace std;

// 1D overload
double ssp_rk3_step(vector<ConsvState>& U, int N, double dx, double gamma,
                    double CFL, double time, double final_time);

// 2D overload
double ssp_rk3_step(vector<vector<ConsvState>>& U,
                    const Grid2D& g, double gamma,
                    double CFL, double time, double final_time, double Mach, double theta, const WallCurve& wedge);