// ── LOAD CSV ─────────────────────────────────────────────────────────────
#include "utils.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;

const double WGS84_A  = 6378137.0;
const double WGS84_F  = 1.0/298.257223563;
const double WGS84_E2 = 2*WGS84_F - WGS84_F*WGS84_F;
const double PI       = 3.14159265358979;
const double RAD2DEG  = 180.0 / PI;

vector<Observation> loadObservations(const string& filename) 
{
	vector<Observation> obs;
	ifstream file(filename);
	if (!file.is_open()) 
	{
		cerr << "Cannot open " << filename << endl;
		return obs;
	}
	string line;
	getline(file, line);  // skip CSV header
	while (getline(file, line)) 
	{
		if (line.empty()) continue;
		stringstream ss(line);
		string field;
		Observation o;
		try {
			getline(ss, field, ','); o.prn = stoi(field);
			getline(ss, field, ','); o.sat_x = stod(field);
			getline(ss, field, ','); o.sat_y = stod(field);
			getline(ss, field, ','); o.sat_z = stod(field);
			getline(ss, field, ','); o.pseudorange = stod(field);
			obs.push_back(o);
		}
		catch (...) { continue; }
	}
	return obs;
}

void ecefToLla(double X, double Y, double Z,
	double& lat_deg, double& lon_deg, double& h_m) {
	double lon_rad = atan2(Y, X);
	double p = sqrt(X * X + Y * Y);
	double lat_rad = atan2(Z, p * (1 - WGS84_E2));
	// Iterate (geodetic latitude is implicit)
	for (int i = 0; i < 5; i++) {
		double N = WGS84_A / sqrt(1 - WGS84_E2 * sin(lat_rad) * sin(lat_rad));
		h_m = p / cos(lat_rad) - N;
		lat_rad = atan2(Z, p * (1 - WGS84_E2 * N / (N + h_m)));
	}
	lat_deg = lat_rad * RAD2DEG;
	lon_deg = lon_rad * RAD2DEG;
}
