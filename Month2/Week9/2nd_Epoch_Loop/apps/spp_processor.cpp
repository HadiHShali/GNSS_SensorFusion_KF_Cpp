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
	ofstream csv("../../data/trajectory.csv");
	csv << "t_gps,lat_deg,lon_deg,height_m,sigma0_m,pdop,n_sats,converged" << endl;
	
	// Track statistics
	int successful = 0;
	int failed = 0;
	
	// ---------------------------------------------------------------------------------//
	//								Loop. Process every Epoch							//
	// ---------------------------------------------------------------------------------//
    
	for (const auto& epoch: epochs)
	{
		// skip epochs with too few satellites
		if (epoch.gps_records.size() < 4){failed++; continue;}
		
		// Step1: Build observations (Sagnac + sat clock + relativity + TGD)
		vector<Observation> obs;
		for (const auto& rec: epoch.gps_records)
		{
			const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
			if (!eph) continue;
			
			SatPosition sp = computeSatPosECEF(*eph, epoch.t_gps);
			
			//Sangac
			double travel_time = rec.pseudorange / C_LIGHT;
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
            clk -= eph->TGD;

            Observation o;
            o.prn = rec.prn;
            o.sat_x = sp.X; o.sat_y = sp.Y; o.sat_z = sp.Z;
            o.pseudorange = rec.pseudorange + C_LIGHT * clk;
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
            double iono = klobucharIonoDelay(kp, rough_lat, rough_lon,
                                              az, el, epoch.t_gps);
            double tropo = saastamoinenTropoDelay(rough_lat, rough_h, el);
            Observation o_new = o;
            o_new.pseudorange -= (iono + tropo);
            corrected.push_back(o_new);
        }


        // ── STEP 5: Final solve ──
        PositionSolution sol = solveGpsPosition(corrected, false);
        if (!sol.converged) { failed++; continue; }
        
        // ── STEP 6: Save to CSV ──
        double lat, lon, h;
        ecefToLla(sol.x, sol.y, sol.z, lat, lon, h);
        csv << fixed << setprecision(6)
            << epoch.t_gps << "," << lat << "," << lon << ","
            << h << "," << sol.sigma0 << "," << sol.dop.PDOP << ","
            << corrected.size() << "," << sol.converged << endl;
        
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
    cout << "Output: trajectory.csv" << endl;
    
    return 0;

}