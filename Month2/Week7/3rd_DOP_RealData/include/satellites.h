#pragma once
#include <vector>
#include <string>

struct GpsEphemeris {
    int prn;
    int year, month, day, hour, minute;
    double second;
    double clk_bias, clk_drift, clk_drift_rate;
    double IODE, Crs, delta_n, M0;
    double Cuc, e, Cus, sqrt_a;
    double toe, Cic, Omega0, Cis;
    double i0, Crc, omega, Omega_dot;
    double i_dot, L2_codes, gps_week, L2P_flag;
    double sv_accuracy, sv_health, TGD, IODC;
    double trans_time, fit_interval;
};

struct SatPosition { double X, Y, Z; };

// Public API
std::vector<GpsEphemeris> parseRinexNav(const std::string& filename);

SatPosition computeSatPosECEF(const GpsEphemeris& eph, double t_gps);

const GpsEphemeris* findBestEphemeris(
    const std::vector<GpsEphemeris>& ephs, int prn, double t_gps);
	
double computeElevation(double rec_x, double rec_y, double rec_z,
                        double sat_x, double sat_y, double sat_z);