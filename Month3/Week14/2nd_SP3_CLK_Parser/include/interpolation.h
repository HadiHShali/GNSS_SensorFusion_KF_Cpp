#pragma once
#include <vector>
#include <string>
#include "sp3_reader.h"

// 10th-order Lagrange interpolation for one coordinate
inline double lagrangeInterp(const std::vector<double>& times,
    const std::vector<double>& values,
    double t)
{
    double result = 0.0;
    int n = static_cast<int>(times.size());
    for (int i = 0; i < n; i++) {
        double term = values[i];
        for (int j = 0; j < n; j++) {
            if (j != i) term *= (t - times[j]) / (times[i] - times[j]);
        }
        result += term;
    }
    return result;
}

// Linear interpolation for clocks -- deliberately NOT Lagrange 
inline double linearClockInterp(double t1, double clk1,
    double t2, double clk2, double t)
{
    return clk1 + (clk2 - clk1) * (t - t1) / (t2 - t1);
}

// Convenience wrapper: find the 10 nearest SP3 epochs to target_t for
// a given satellite, run lagrangeInterp() on each of X/Y/Z.
// Returns false if satellite not found or insufficient surrounding epochs.

//
//Given a satellite's PRN and a target time, it finds that satellite's 10 nearest SP3 epochs,
// then uses Lagrange interpolation to compute the satellite's exact ECEF
// position at that target time — writing the result into out_pos[3].

inline bool interpolateSatPosition(const std::vector<Sp3Epoch>& epochs,
    const std::string& prn,
    double target_t,
    double out_pos[3])
{

    // Phase 1 — Collect This Satellite's Data Across All Epochs
    std::vector<double> times;
    std::vector<double> xs, ys, zs;
    for (const auto& ep : epochs) {
        auto it = ep.sat_pos.find(prn);
        if (it != ep.sat_pos.end()) {
            times.push_back(ep.t_gps);
            xs.push_back(it->second[0]);
            ys.push_back(it->second[1]);
            zs.push_back(it->second[2]);
        }
    }
    if (times.size() < 10) return false;   // not enough data for order-10

    // Phase 1: Find the 10 epochs closest to target_t (5 before, 5 after ideally)
    // Simplified: find nearest index, take a centered window of 10
    int nearest = 0;
    double best_diff = 1e18;
    for (size_t i = 0; i < times.size(); i++) {
        double diff = std::abs(times[i] - target_t);
        if (diff < best_diff) { best_diff = diff; nearest = static_cast<int>(i); }
    }
    int start = nearest - 5;
    if (start < 0) start = 0;
    if (start + 10 > static_cast<int>(times.size()))
        start = static_cast<int>(times.size()) - 10;

    std::vector<double> wt(times.begin() + start, times.begin() + start + 10);
    std::vector<double> wx(xs.begin() + start, xs.begin() + start + 10);
    std::vector<double> wy(ys.begin() + start, ys.begin() + start + 10);
    std::vector<double> wz(zs.begin() + start, zs.begin() + start + 10);

    out_pos[0] = lagrangeInterp(wt, wx, target_t);
    out_pos[1] = lagrangeInterp(wt, wy, target_t);
    out_pos[2] = lagrangeInterp(wt, wz, target_t);
    return true;
}



