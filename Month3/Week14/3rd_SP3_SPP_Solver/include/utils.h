#pragma once
#include <vector>
#include <string>
#include "gps_ls_solver.h"

std::vector<Observation> loadObservations(const std::string& filename);
void ecefToLla(double X, double Y, double Z, double& lat, double& lon, double& h);
void llaToEcef(double lat_deg, double lon_deg, double h_m, double& x, double& y, double& z);

// ─────────────────────────────────────────────────────────────
// Dual-frequency ionosphere-free combination (Week 13)
// Formula from IS-GPS-200: P_IF = K1*P1 - K2*P2
// Verified current practice (2024-2025 literature) -- unchanged since 1991.
// ─────────────────────────────────────────────────────────────
constexpr double GPS_F1 = 1575.42e6;   // Hz, L1 carrier frequency
constexpr double GPS_F2 = 1227.60e6;   // Hz, L2 carrier frequency
constexpr double GPS_K1 = (GPS_F1 * GPS_F1) / (GPS_F1 * GPS_F1 - GPS_F2 * GPS_F2);  // ~2.546
constexpr double GPS_K2 = (GPS_F2 * GPS_F2) / (GPS_F1 * GPS_F1 - GPS_F2 * GPS_F2);  // ~1.546

// Computes the ionosphere-free pseudorange from L1 (P1) and L2 (P2)
// raw pseudoranges, both in meters. Eliminates first-order ionospheric
// delay algebraically -- do NOT apply Klobuchar correction on top of
// the result (would double-correct a term already physically removed).
inline double ionoFreePseudorange(double P1, double P2) {
    return GPS_K1 * P1 - GPS_K2 * P2;
}