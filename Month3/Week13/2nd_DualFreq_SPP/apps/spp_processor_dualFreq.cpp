#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "gps_ls_solver.h"
#include "utils.h"
#include "satellites.h"
#include "rinex_obs.h"
 
using namespace std;
using namespace Eigen;
 
const double C_LIGHT = 299792458.0;
const double OMEGA_EARTH = 7.2921151467e-5;

int main()
{
	cout << "=== FULL SPP PROCESSOR ===" << endl;
    
	// ---------------------------------------------------------------------------------//
	//								One_time Step (before loop)							//
	// ---------------------------------------------------------------------------------//
    string nav_file = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    string obs_file = "../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx";
	
	// Parse Nav file
    vector<GpsEphemeris> ephs = parseRinexNav(nav_file);
	KlobucharParams kp = parseKlobucharFromNav(nav_file);
    cout << "Parsed " << ephs.size() << " ephemerides" << endl;
	
	// Parse All observation epochs
	vector<ObsEpoch> epochs = parseRinexAllEpochs(obs_file);
	cout << "Parsed " << epochs.size() << " epochs" << endl;
	
	// Open trajectoory CSV
	ofstream csv("../../data/trajectory_dualfreq.csv");
	csv << "t_gps,lat_deg,lon_deg,height_m,sigma0_m,pdop,n_sats,converged,n_iono_free,n_single_freq,err_3d" << endl;
	
	// Track statistics
	int successful = 0;
	int failed = 0;
	// NEW: for 3D error summary
	const double TRUE_X = -2420420.5012;
	const double TRUE_Y = -4737131.5508;
	const double TRUE_Z =  3507827.665;
	vector<double> errors_3d;   // store every successful epoch's error
	// ---------------------------------------------------------------------------------//
	//								Loop. Process every Epoch							//
	// ---------------------------------------------------------------------------------//
    
	for (const auto& epoch: epochs)
	{
		// skip epochs with too few satellites
		if (epoch.gps_records.size() < 4){failed++; continue;}
		
		// Step1: Build observations (Sagnac + sat clock + relativity + TGD)
		int n_iono_free_epoch = 0, n_single_freq_epoch = 0;   // NEW: per-epoch counters
		vector<Observation> obs;
		for (const auto& rec: epoch.gps_records)
		{
			const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
			if (!eph) continue;
			
			SatPosition sp = computeSatPosECEF(*eph, epoch.t_gps);
			
			// Choose pseudorange source based on C2W availability
			double base_pr;
			bool is_iono_free=false;
			if (rec.has_c2){
				base_pr = ionoFreePseudorange(rec.pseudorange, rec.pseudorange_c2);
				is_iono_free = true;
				n_iono_free_epoch++;
				
			} else{
				base_pr = rec.pseudorange;
				is_iono_free = false;
				n_single_freq_epoch++;
			}
			
			    // ────────── ADD DEBUG PRINT RIGHT HERE ──────────
			static bool printed_once = false;
			if (!printed_once && is_iono_free) {
				cout << fixed << setprecision(3);   // ADD THIS — show 3 decimal places
				cout << "  DEBUG sat G" << rec.prn
					 << "  P1(C1W)=" << rec.pseudorange
					 << "  P2(C2W)=" << rec.pseudorange_c2
					 << "  P_IF="    << base_pr
					 << "  (P2-P1)=" << (rec.pseudorange_c2 - rec.pseudorange)
					 << endl;
				cout.unsetf(ios::fixed);   // reset formatting for the rest of the program
				printed_once = true;
			}
			// ──────────────────────────────────────────────
			
			//Sangac
			double travel_time = base_pr / C_LIGHT;
			double rot = OMEGA_EARTH * travel_time;
			double cr = cos(rot), sr = sin(rot);
			double sx = cr * sp.X + sr * sp.Y;
			double sy = -sr * sp.X + cr * sp.Y;
			sp.X = sx; sp.Y = sy;
			
			// Full satellite clock correction
			double dt_sat = epoch.t_gps - eph->toe; 
			if (dt_sat >  302400.0) dt_sat -= 604800.0;
            if (dt_sat < -302400.0) dt_sat += 604800.0;
            double clk = eph->clk_bias + eph->clk_drift * dt_sat
                       + eph->clk_drift_rate * dt_sat * dt_sat;
			
			// Relativistic
			double a = eph->sqrt_a * eph->sqrt_a;
            double n = sqrt(3.986005e14 / (a*a*a)) + eph->delta_n;
            double Mk = eph->M0 + n * dt_sat;
            double Ek = Mk;
            for (int j = 0; j < 10; j++) Ek = Mk + eph->e * sin(Ek);
            clk += -4.442807633e-10 * eph->e * eph->sqrt_a * sin(Ek);
			//clk -= eph->TGD;
			// TGD: Only for single frequency (see day 3 research check note)
			if (!is_iono_free){
				clk -= eph->TGD;
			}
            

            Observation o;
            o.prn = rec.prn;
            o.sat_x = sp.X; o.sat_y = sp.Y; o.sat_z = sp.Z;
            o.pseudorange = base_pr + C_LIGHT * clk;
			o.used_iono_free = is_iono_free;   // NEW
            obs.push_back(o);
		}
		if (obs.size() < 4) { failed++; continue; }
		
		
		// Step2: Rough Solve
		        // ── STEP 2: Rough solve ──
        PositionSolution rough = solveGpsPosition(obs, false);
        if (!rough.converged) { failed++; continue; }


// ── STEP 3: Elevation mask ──
        vector<Observation> masked;
        for (const auto& o : obs) {
            double el = computeElevation(rough.x, rough.y, rough.z,
                                          o.sat_x, o.sat_y, o.sat_z);
            if (el >= 10.0) masked.push_back(o);
        }
        if (masked.size() < 4) { failed++; continue; }

		
		
		double rough_lat, rough_lon, rough_h;
        ecefToLla(rough.x, rough.y, rough.z, rough_lat, rough_lon, rough_h);


        // ── STEP 4: Apply atmospheric corrections ──
        vector<Observation> corrected;
        for (const auto& o : masked) {
            double el = computeElevation(rough.x, rough.y, rough.z,
                                          o.sat_x, o.sat_y, o.sat_z);
            double az = computeAzimuth(rough.x, rough.y, rough.z,
                                        o.sat_x, o.sat_y, o.sat_z);
            double tropo = saastamoinenTropoDelay(rough_lat, rough_h, el);
			
			double iono = 0.0;
			if (!o.used_iono_free){ // Skip Klobuchar for already iono free obse
				iono = klobucharIonoDelay(kp, rough_lat, rough_lon, az, el, epoch.t_gps);
			}
			
            Observation o_new = o;
            o_new.pseudorange -= (iono + tropo);
            corrected.push_back(o_new);
        }


        // ── STEP 5: Final solve ──
        PositionSolution sol = solveGpsPosition(corrected, false);
        if (!sol.converged) { failed++; continue; }
        
		// NEW: compute 3D error vs truth for this epoch
		double err_x = sol.x - TRUE_X;
		double err_y = sol.y - TRUE_Y;
		double err_z = sol.z - TRUE_Z;
		double err_3d = sqrt(err_x*err_x + err_y*err_y + err_z*err_z);
		errors_3d.push_back(err_3d);
        // ── STEP 6: Save to CSV ──
        double lat, lon, h;
        ecefToLla(sol.x, sol.y, sol.z, lat, lon, h);
        csv << fixed << setprecision(6)
            << epoch.t_gps << "," << lat << "," << lon << ","
            << h << "," << sol.sigma0 << "," << sol.dop.PDOP << ","
            << corrected.size() << "," << sol.converged << ","
			<< n_iono_free_epoch << "," << n_single_freq_epoch << ","
			<< err_3d << endl;   // NEW columns
        
        successful++;

	}
	
	csv.close();
    
	cout << endl << "=================================" << endl;
	cout << "       PROCESSING COMPLETE" << endl;
	cout << "==========================================" << endl;
	cout << "Total epochs:      " << epochs.size() << endl;
	cout << "Successful fixes:  " << successful << endl;
	cout << "Failed:            " << failed << endl;
	cout << "Success rate:      "
		 << (100.0 * successful / epochs.size()) << "%" << endl;
	cout << "Output: trajectory_dualfreq.csv" << endl;

	// NEW: 3D error statistics
	if (!errors_3d.empty()) {
		double sum = 0.0, sum_sq = 0.0;
		double min_err = errors_3d[0], max_err = errors_3d[0];
		for (double e : errors_3d) {
			sum += e;
			sum_sq += e * e;
			if (e < min_err) min_err = e;
			if (e > max_err) max_err = e;
		}
		double mean_err = sum / errors_3d.size();
		double rms_err  = sqrt(sum_sq / errors_3d.size());

		cout << endl << "=== 3D ACCURACY vs TRUE POSITION ===" << endl;
		cout << fixed << setprecision(3);
		cout << "Mean 3D error:  " << mean_err << " m" << endl;
		cout << "RMS  3D error:  " << rms_err  << " m" << endl;
		cout << "Min  3D error:  " << min_err  << " m" << endl;
		cout << "Max  3D error:  " << max_err  << " m" << endl;
	}
    
    return 0;

}