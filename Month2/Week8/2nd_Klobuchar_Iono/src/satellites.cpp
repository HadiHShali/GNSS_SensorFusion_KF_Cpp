#include "satellites.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;

const double MU_EARTH = 3.986005e14;
const double OMEGA_EARTH = 7.2921151467e-5;

// Helper functions (static = private to this file)
static double solveKepler(double M, double e) {
    double E = M;
    for (int i = 0; i < 20; i++) {
        double E_new = M + e * sin(E);
        if (abs(E_new - E) < 1e-12) return E_new;
        E = E_new;
    }
    return E;
}

static string fortranToCpp(const string& s) {
    string out = s;
    for (char& c : out) if (c == 'D' || c == 'd') c = 'E';
    return out;
}

static double parseField(const string& line, size_t start, size_t len = 19) {
    if (start >= line.length()) return 0.0;
    string field = fortranToCpp(line.substr(start, len));
    try { return stod(field); }
    catch (...) { return 0.0; }
}

static void readOrbitLine(const string& line, double& v1, double& v2,
    double& v3, double& v4) {
    v1 = parseField(line, 4);  v2 = parseField(line, 23);
    v3 = parseField(line, 42); v4 = parseField(line, 61);
}

// Then paste your three functions here:
// - parseRinexNav()
// Main Parser
vector<GpsEphemeris> parseRinexNav(const string& filename)
{
	vector<GpsEphemeris> ephemerides;
	ifstream file(filename);
	if (!file.is_open())
	{
		cerr << "Error: Can not Open the" << filename << endl;
		return ephemerides;
	}

	string line;
	bool header_done = false;

	// skip the header
	while (getline(file, line))
	{
		if (line.find("END OF HEADER") != string::npos)
		{
			header_done = true;
			break;
		}
	}
	if (!header_done)
	{
		cerr << "Error: No END OF HEADER found" << endl;
		return ephemerides;
	}


	// Step2 : Read records
	while (getline(file, line))
	{
		//skip empty lines
		if (line.length() < 23) continue;

		// GPS records start with 'G'
		if (line[0] != 'G')
		{
			// skip non-GPS satellite for today (R, E, C. etc.)
			// each line is 8 lines - skip 7 more lines
			for (int i = 0; i < 7; i++)
			{
				if (!getline(file, line)) break;//!getline(file, line) becomes true if: End-of-file (EOF) is reached, or A read error occurs.
			}
			continue;
		}

		// Step3: We have a GPS record. Parse all 8 line.
		GpsEphemeris eph = {}; // zero initialize
		// line 0: PRN, Epoch, Clock parameters
		try
		{
			eph.prn = stoi(line.substr(1, 2));
			eph.year = stoi(line.substr(4, 4));
			eph.month = stoi(line.substr(9, 2));
			eph.day = stoi(line.substr(12, 2));
			eph.hour = stoi(line.substr(15, 2));
			eph.minute = stoi(line.substr(18, 2));
			eph.second = parseField(line, 21, 2);
		}
		catch (...) { continue; }

		eph.clk_bias = parseField(line, 23);
		eph.clk_drift = parseField(line, 42);
		eph.clk_drift_rate = parseField(line, 61);


		// Lines 1-7: broadcast orbit parameters (4 numbers each)
		if (!getline(file, line)) break;
		readOrbitLine(line, eph.IODE, eph.Crs, eph.delta_n, eph.M0);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.Cuc, eph.e, eph.Cus, eph.sqrt_a);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.toe, eph.Cic, eph.Omega0, eph.Cis);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.i0, eph.Crc, eph.omega, eph.Omega_dot);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.i_dot, eph.L2_codes, eph.gps_week, eph.L2P_flag);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.sv_accuracy, eph.sv_health, eph.TGD, eph.IODC);

		if (!getline(file, line)) break;
		double spare1, spare2;
		readOrbitLine(line, eph.trans_time, eph.fit_interval, spare1, spare2);


		// Step 4: Store the complete ephemeris
		ephemerides.push_back(eph);
	}

	file.close();
	return ephemerides;
}

// - computeSatPosECEF()
SatPosition computeSatPosECEF(const GpsEphemeris& eph, double t_gps)
{
	SatPosition pos;
	double a = eph.sqrt_a * eph.sqrt_a;
	double n = sqrt(MU_EARTH / (a * a * a)) + eph.delta_n;
	double tk = t_gps - eph.toe;
	if (tk > 302400.0) tk -= 604800.0;
	if (tk < -302400.0) tk += 604800.0;
	double Mk = eph.M0 + n * tk;
	double Ek = solveKepler(Mk, eph.e);
	double nuk = atan2(sqrt(1 - eph.e * eph.e) * sin(Ek), cos(Ek) - eph.e);
	double phi_k = nuk + eph.omega;
	double sin_2phi = sin(2 * phi_k), cos_2phi = cos(2 * phi_k);
	double u_k = phi_k + eph.Cuc * cos_2phi + eph.Cus * sin_2phi;
	double r_k = a * (1 - eph.e * cos(Ek)) + eph.Crc * cos_2phi + eph.Crs * sin_2phi;
	double i_k = eph.i0 + eph.i_dot * tk + eph.Cic * cos_2phi + eph.Cis * sin_2phi;
	double x_orb = r_k * cos(u_k);
	double y_orb = r_k * sin(u_k);
	double Omega_k = eph.Omega0 + (eph.Omega_dot - OMEGA_EARTH) * tk - OMEGA_EARTH * eph.toe;
	pos.X = x_orb * cos(Omega_k) - y_orb * cos(i_k) * sin(Omega_k);
	pos.Y = x_orb * sin(Omega_k) + y_orb * cos(i_k) * cos(Omega_k);
	pos.Z = y_orb * sin(i_k);
	return pos;
}
// - findBestEphemeris()
// New Today: Pick the best ephemeris for a satellite at time t
// Returns iterators to the ephemeris closest in time to t_gps for this prn
// findBestEphemeris() answers this question:
// "I have hundreds of ephemerides in my vector. Out of all of them, 
// which one belongs to satellite #5 AND has a time-of-ephemeris closest to right now?"

// Why this matters: Remember from Day 3, your BRDC file has 440 ephemerides(about 14 per satellite).
// For each satellite, you have multiple versions recorded throughout the day.You want the one that's
// most accurate for your target time — the one whose toe is closest to t_gps.
// This function searches the entire vector and returns a pointer (address in memory) to the best match.

//const: "Whatever I return, the caller can't modify it"
//GpsEphemeris: The type of thing we're pointing to. 
// *: POINTER — an address (where it is located in memory), not the data itself
// So this function returns a pointer to a GpsEphemeris. NOT a copy.NOT the struct itself. Just its address in memory.
const GpsEphemeris* findBestEphemeris(const vector<GpsEphemeris>& ephs, int prn, double t_gps)
{
	const GpsEphemeris* best = nullptr;      //nullptr means "no pointer yet — empty". We start with nothing found, then update as we discover matches.
	double best_dt = 1e18;                  //will hold the time difference of the best match so far. We want to find the SMALLEST dt, so we start at a deliberately HUGE value.
	for (const auto& e : ephs) {            //read-only reference to each element in the ALL 440 ephemerides. 
		if (e.prn != prn) continue;
		double dt = abs(e.toe - t_gps);
		if (dt > 302400.0) dt = 604800.0 - dt;  // week rollover
		if (dt < best_dt) { best_dt = dt; best = &e; }  //The & operator means "give me the address of". So, '&e' gives the address of e in the memory
	}
	return best;
}

// compute elevation angle (in degrees) from receiver to satellite
double computeElevation(double rec_x, double rec_y, double rec_z, double sat_x, double sat_y, double sat_z)
{
	// receiver lat lon (approximate from ECEF)
	double r = sqrt(rec_x*rec_x + rec_y*rec_y + rec_z*rec_z);
	double lat = atan2(rec_z, sqrt(rec_x*rec_x+rec_y*rec_y));
	double lon = atan2(rec_y, rec_x);

	// satellite vector relative to receiver
	double dx = sat_x - rec_x;
	double dy = sat_y - rec_y;
	double dz = sat_z - rec_z;

	// rotate ECEF to ENU
	double sin_lat = sin(lat);
	double cos_lat = cos(lat);
	double sin_lon = sin(lon);
	double cos_lon = cos(lon);

	double east = -sin_lon * dx + cos_lon * dy;
	double north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
	double up    =  cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;

	// elevation angle
	double horizontal = sqrt(east*east + north*north);
	double elev_rad = atan2(up, horizontal);

	return elev_rad * 180.0 / 3.14159265358979;
}

double computeAzimuth (double rec_x, double rec_y, double rec_z, double sat_x, double sat_y, double sat_z)
{
	// same ENU transform as computeElevation, but return azimuth
	double lat = atan2(rec_z, sqrt(rec_x*rec_x + rec_y*rec_y));
	double lon = atan2(rec_y, rec_x);
	double dx = sat_x - rec_x;
	double dy = sat_y - rec_y;
	double dz = sat_z - rec_z;
	double sin_lat = sin(lat), cos_lat = cos(lat);
	double sin_lon = sin(lon), cos_lon = cos(lon);
	double east  = -sin_lon * dx + cos_lon * dy;
	double north = -sin_lat*cos_lon*dx - sin_lat*sin_lon*dy + cos_lat*dz;
    double az_rad = atan2(east, north);  // atan2 returns -pi..pi
    if (az_rad < 0) az_rad += 2 * 3.14159265358979;
    return az_rad * 180.0 / 3.14159265358979;
}


// Parse Klobuchar paramaters from RINEX Nav file header
KlobucharParams parseKlobucharFromNav(const std::string& filename)
{
	KlobucharParams kp;
	
	ifstream file(filename);
	if (!file.is_open()) return kp;
	
	string line;
	bool found_gpsa = false;
	bool found_gpsb = false;
	
	while (getline(file, line))
	{
		if (line.find("END OF HEADER") !=string::npos) break;
		
		if (line.find("IONOSPHERIC CORR") != string::npos)
		{
			if (line.substr(0, 4) == "GPSA")
			{
				kp.alpha[0] = parseField(line, 5, 12);
				kp.alpha[1] = parseField(line, 17, 12);
				kp.alpha[2] = parseField(line, 29, 12);
				kp.alpha[3] = parseField(line, 41, 12);
				found_gpsa = true;
			}
			else if (line.substr(0, 4)=="GPSB")
			{
				kp.beta[0] = parseField(line, 5, 12);
				kp.beta[1] = parseField(line, 17, 12);
				kp.beta[2] = parseField(line, 29, 12);
				kp.beta[3] = parseField(line, 41, 12);
				found_gpsb = true;
			}
		}
	}
	kp.valid = (found_gpsa && found_gpsb);
	return kp;
} 

// Klobuchar Ionospheric Delay Model (IS-GPS-200)
// Returns delay in meters for L1 frequency
double klobucharIonoDelay (const KlobucharParams& kp, double rec_lat_deg, double rec_lon_deg, double sat_az_deg, double sat_el_deg, double t_gps_sec)
{
	if (!kp.valid) return 0.0; //safe fallback
	    
	const double C_LIGHT = 299792458.0;
    
    // Convert to semi-circles (1 semi-circle = 180 deg = pi rad)
    double phi_u   = rec_lat_deg / 180.0;
    double lambda_u = rec_lon_deg / 180.0;
    double A       = sat_az_deg * M_PI / 180.0;       // radians
    double E       = sat_el_deg / 180.0;              // semi-circles
    
    // STEP 1: Earth-centered angle psi (semi-circles)
    double psi = 0.0137 / (E + 0.11) - 0.022;
	
	// STEP 2: Geodetic latitude of IPP (Ionospheric Pierce Point)
    double phi_i = phi_u + psi * cos(A);
    if (phi_i >  0.416) phi_i =  0.416;
    if (phi_i < -0.416) phi_i = -0.416;
	
	// STEP 3: Geodetic longitude of IPP
    double lambda_i = lambda_u + psi * sin(A) / cos(phi_i * M_PI);

    // STEP 4: Geomagnetic latitude
    double phi_m = phi_i + 0.064 * cos((lambda_i - 1.617) * M_PI);

    // STEP 5: Local time at IPP (seconds)
    double t = 43200.0 * lambda_i + t_gps_sec;
    while (t >= 86400.0) t -= 86400.0;
    while (t < 0.0)      t += 86400.0;

    // STEP 6: Compute amplitude and period from alpha/beta
    double AMP = kp.alpha[0] + kp.alpha[1]*phi_m
               + kp.alpha[2]*phi_m*phi_m + kp.alpha[3]*phi_m*phi_m*phi_m;
    if (AMP < 0.0) AMP = 0.0;
	double PER = kp.beta[0] + kp.beta[1]*phi_m
               + kp.beta[2]*phi_m*phi_m + kp.beta[3]*phi_m*phi_m*phi_m;
    if (PER < 72000.0) PER = 72000.0;

    // STEP 7: Phase
    double x = 2.0 * M_PI * (t - 50400.0) / PER;
    
    // STEP 8: Slant factor F (1 at zenith, ~3 at horizon)
    double F = 1.0 + 16.0 * pow(0.53 - E, 3);
    
    // STEP 9: Compute delay in SECONDS
    double T_iono_sec;
    if (fabs(x) < 1.57) {
        T_iono_sec = F * (5.0e-9 + AMP * (1.0 - x*x/2.0 + x*x*x*x/24.0));
    } else {
        T_iono_sec = F * 5.0e-9;
    }
    
    // Convert delay from SECONDS to METERS
    return T_iono_sec * C_LIGHT;
}

