// real_data_position.cpp - Compute receiver position from real RINEX data
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "gps_ls_solver.h"        // Week 6 solver
#include "utils.h"                 // ecefToLla
#include "satellites.h"            // parseRinexNav, computeSatPosECEF, atmospheric models
#include "rinex_obs.h"             // RINEX observation parser

const double C_LIGHT = 299792458.0;
const double OMEGA_EARTH = 7.2921151467e-5;
using namespace std;

int main() {
    cout << "=== REAL GPS POSITIONING WITH ATMOSPHERIC CORRECTIONS ===" << endl << endl;

    // ─────────────────────────────────────────────────────────────
    // STEP 0: Parse RINEX NAV (broadcast ephemerides)
    // ─────────────────────────────────────────────────────────────
    string nav_file = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    vector<GpsEphemeris> ephs = parseRinexNav(nav_file);
    cout << "Parsed " << ephs.size() << " ephemerides" << endl;

    // ─────────────────────────────────────────────────────────────
    // STEP 1: Parse Klobuchar iono parameters
    // ─────────────────────────────────────────────────────────────
    KlobucharParams kp = parseKlobucharFromNav(nav_file);
    if (kp.valid) {
        cout << "Loaded Klobuchar Iono parameters" << endl;
    } else {
        cout << "WARNING: no Klobuchar params - iono correction disabled" << endl;
    }

    // ─────────────────────────────────────────────────────────────
    // STEP 2: Parse RINEX OBS (pseudoranges) at noon UTC
    // ─────────────────────────────────────────────────────────────
    string obs_file = "../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx";
    ObsEpoch epoch = parseRinexObsOneEpoch(obs_file, 43200.0);
    cout << "Epoch: " << epoch.year << "-" << epoch.month << "-" << epoch.day
         << " " << epoch.hour << ":" << epoch.minute << ":" << epoch.second << endl;
    cout << "GPS sec of week: " << epoch.t_gps << endl;
    cout << "Number of GPS observations: " << epoch.gps_records.size() << endl;

    if (epoch.gps_records.size() < 4) {
        cerr << "Need at least 4 satellites!" << endl;
        return 1;
    }

    // ─────────────────────────────────────────────────────────────
    // STEP 3: Build initial observations (sat clock + Earth rotation only)
    // ─────────────────────────────────────────────────────────────
    vector<Observation> obs;
    cout << endl << "Building observations..." << endl;

    for (const auto& rec : epoch.gps_records) {
        const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
        if (!eph) continue;

        SatPosition sp = computeSatPosECEF(*eph, epoch.t_gps);

        // Earth rotation correction (Sagnac)
        double travel_time = rec.pseudorange / C_LIGHT;
        double rotation_angle = OMEGA_EARTH * travel_time;
        double cos_r = cos(rotation_angle);
        double sin_r = sin(rotation_angle);
        double sx_new =  cos_r * sp.X + sin_r * sp.Y;
        double sy_new = -sin_r * sp.X + cos_r * sp.Y;
        sp.X = sx_new;
        sp.Y = sy_new;

        // ─── Full satellite clock correction ────────────────────
		double dt_sat = epoch.t_gps - eph->toe;
		if (dt_sat >  302400.0) dt_sat -= 604800.0;
		if (dt_sat < -302400.0) dt_sat += 604800.0;

		// Polynomial clock model
		double sat_clock_correction = eph->clk_bias
			+ eph->clk_drift * dt_sat
			+ eph->clk_drift_rate * dt_sat * dt_sat;

		// Relativistic correction
		double a_sma = eph->sqrt_a * eph->sqrt_a;
		double n_mean = sqrt(3.986005e14 / (a_sma * a_sma * a_sma)) + eph->delta_n;
		double Mk = eph->M0 + n_mean * dt_sat;
		double Ek = Mk;
		for (int j = 0; j < 10; j++) {
			Ek = Mk + eph->e * sin(Ek);
		}
		const double F_REL = -4.442807633e-10;
		double rel_corr = F_REL * eph->e * eph->sqrt_a * sin(Ek);
		sat_clock_correction += rel_corr;

		// TGD (Timing Group Delay) - for single-frequency L1 users
		sat_clock_correction -= eph->TGD;
		// ────────────────────────────────────────────────────────

        Observation o;
        o.prn = rec.prn;
        o.sat_x = sp.X;
        o.sat_y = sp.Y;
        o.sat_z = sp.Z;
        o.pseudorange = rec.pseudorange + C_LIGHT * sat_clock_correction;
        obs.push_back(o);
    }
    cout << "Built " << obs.size() << " observation tuples" << endl << endl;

    if (obs.size() < 4) { cerr << "Not enough observations!" << endl; return 1; }

    // ─────────────────────────────────────────────────────────────
    // PASS 1a: Dirty rough solve with all satellites (to get initial position)
    // ─────────────────────────────────────────────────────────────
    PositionSolution sol1a = solveGpsPosition(obs, false);
    cout << "PASS 1a (all sats): " << obs.size() << " sats, PDOP = "
         << sol1a.dop.PDOP << ", sigma0 = " << sol1a.sigma0 << " m" << endl;

    // ─────────────────────────────────────────────────────────────
    // Filter by elevation using PASS 1a rough position
    // ─────────────────────────────────────────────────────────────
    vector<Observation> filtered_obs;
    const double ELEV_MASK = 10.0;  // degrees

    for (const auto& o : obs) {
        double elev = computeElevation(sol1a.x, sol1a.y, sol1a.z,
                                        o.sat_x, o.sat_y, o.sat_z);
        if (elev >= ELEV_MASK) {
            filtered_obs.push_back(o);
        } else {
            cout << "  G" << o.prn << " rejected (elev = " << elev << " deg)" << endl;
        }
    }

    if (filtered_obs.size() < 4) {
        cerr << "Not enough satellites after masking!" << endl;
        return 1;
    }

    // ─────────────────────────────────────────────────────────────
    // PASS 1: CLEAN rough solve with masked satellites (no atmo yet)
    // ─────────────────────────────────────────────────────────────
    PositionSolution sol1 = solveGpsPosition(filtered_obs, false);
    cout << endl << "PASS 1 (masked, no atmo): " << filtered_obs.size()
         << " sats, PDOP = " << sol1.dop.PDOP
         << ", sigma0 = " << sol1.sigma0 << " m" << endl;

    // Convert clean rough position to lat/lon/height for atmo models
    double rough_lat, rough_lon, rough_h;
    ecefToLla(sol1.x, sol1.y, sol1.z, rough_lat, rough_lon, rough_h);
    cout << "Rough position lat=" << rough_lat << " lon=" << rough_lon
         << " h=" << rough_h << endl << endl;


	// After computing sol1:
	const double TRUE_X = -2420422.0410;
	const double TRUE_Y = -4737132.4707;
	const double TRUE_Z =  3507827.6034;

	double sol1_err_x = sol1.x - TRUE_X;
	double sol1_err_y = sol1.y - TRUE_Y;
	double sol1_err_z = sol1.z - TRUE_Z;

	cout << "PASS 1 sol vs truth:" << endl;
	cout << "  dX = " << sol1_err_x << " m" << endl;
	cout << "  dY = " << sol1_err_y << " m" << endl;
	cout << "  dZ = " << sol1_err_z << " m" << endl;
	cout << "  clock bias = " << sol1.clock_bias << " m" << endl;
	cout << "  iterations = " << sol1.iterations << endl;
	cout << "  converged  = " << sol1.converged << endl;
	
	cout << "PASS 1 residuals:" << endl;
	for (int i = 0; i < sol1.residuals.size(); i++) {
		cout << "  " << i << ": " << sol1.residuals(i) << " m" << endl;
	}
	
	
	cout << "PASS 1 ECEF: X=" << sol1.x 
     << " Y=" << sol1.y 
     << " Z=" << sol1.z << endl;
	cout << "Truth  ECEF: X=" << TRUE_X 
     << " Y=" << TRUE_Y 
     << " Z=" << TRUE_Z << endl;
    // ─────────────────────────────────────────────────────────────
    // PASS 2: Apply atmospheric corrections to SAME masked satellites
    // ─────────────────────────────────────────────────────────────
    vector<Observation> corrected_obs;
    cout << "Applying atmospheric corrections:" << endl;

    for (const auto& o : filtered_obs) {
        double el = computeElevation(sol1.x, sol1.y, sol1.z,
                                      o.sat_x, o.sat_y, o.sat_z);
        double az = computeAzimuth(sol1.x, sol1.y, sol1.z,
                                    o.sat_x, o.sat_y, o.sat_z);

        double iono  = klobucharIonoDelay(kp, rough_lat, rough_lon,
                                           az, el, epoch.t_gps);
        double tropo = saastamoinenTropoDelay(rough_lat, rough_h, el);

        Observation o_new = o;
        o_new.pseudorange = o.pseudorange - iono - tropo;
        corrected_obs.push_back(o_new);

        cout << "  G" << o.prn << " el=" << el
             << " iono=" << iono << " tropo=" << tropo << endl;
    }

    PositionSolution sol2 = solveGpsPosition(corrected_obs, false);

    // ─────────────────────────────────────────────────────────────
    // Compare PASS 1 vs PASS 2 (residual metrics)
    // ─────────────────────────────────────────────────────────────
    cout << endl << "==============================" << endl;
    cout << "    BEFORE vs AFTER CORRECTIONS   " << endl;
    cout << "==============================" << endl;
    cout << fixed << setprecision(2);
    cout << "PASS 1 (no atmo):   sigma-0 = " << sol1.sigma0
         << " m, PDOP = " << sol1.dop.PDOP << endl;
    cout << "PASS 2 (with atmo): sigma-0 = " << sol2.sigma0
         << " m, PDOP = " << sol2.dop.PDOP << endl;
    cout << "Sigma-0 reduction: " << (sol1.sigma0 - sol2.sigma0)
         << " m (" << (100.0 * (sol1.sigma0 - sol2.sigma0) / sol1.sigma0)
         << "%)" << endl;

    // ─────────────────────────────────────────────────────────────
    // Print final position in LLA
    // ─────────────────────────────────────────────────────────────
    double final_lat, final_lon, final_h;
    ecefToLla(sol2.x, sol2.y, sol2.z, final_lat, final_lon, final_h);
    cout << endl << "=== FINAL POSITION ===" << endl;
    cout << "  Latitude:  " << final_lat << " deg" << endl;
    cout << "  Longitude: " << final_lon << " deg" << endl;
    cout << "  Height:    " << final_h << " m" << endl;

    // ─────────────────────────────────────────────────────────────
    // Compare against true BILL00USA position (from IGS site log)
    // ─────────────────────────────────────────────────────────────

    double err1 = sqrt(pow(sol1.x - TRUE_X, 2)
                     + pow(sol1.y - TRUE_Y, 2)
                     + pow(sol1.z - TRUE_Z, 2));

    double err2 = sqrt(pow(sol2.x - TRUE_X, 2)
                     + pow(sol2.y - TRUE_Y, 2)
                     + pow(sol2.z - TRUE_Z, 2));

    cout << endl << "=== ACCURACY vs TRUE POSITION ===" << endl;
    cout << "PASS 1 (no atmo):   3D error = " << err1 << " m" << endl;
    cout << "PASS 2 (with atmo): 3D error = " << err2 << " m" << endl;
    cout << "Position improvement: " << (err1 - err2) << " m" << endl;

    return 0;
}