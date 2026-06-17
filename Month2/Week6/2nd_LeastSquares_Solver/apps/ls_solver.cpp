// ====================================================================//
//       ls_solver.cpp - Least Squares GPS Position Solver             //
// ====================================================================//

// Implement Δρ = H·Δx  iteratively untill convergence

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;
// 'using Eigen::MatrixXd;' lets you write MatrixXd instead of Eigen::MatrixXd
// This is just a name shortcut, like 'using namespace std;' for cout
using Eigen::MatrixXd;
using Eigen::VectorXd;

// ── CONSTANTS ────────────────────────────────────────────────────────────
const double C_LIGHT = 299792458.0;          // speed of light (m/s)
const double DEG2RAD = 3.14159265358979 / 180.0;
const double RAD2DEG = 180.0 / 3.14159265358979;

// WGS84 (for ECEF -> LLA conversion at the end)
const double WGS84_A = 6378137.0;
const double WGS84_F = 1.0 / 298.257223563;
const double WGS84_E2 = 2 * WGS84_F - WGS84_F * WGS84_F;

// ── DATA STRUCTURES ──────────────────────────────────────────────────────
struct Observation {
	int prn;
	double sat_x, sat_y, sat_z;   // satellite ECEF (m)
	double pseudorange;            // measured pseudorange (m)
};

// ── LOAD CSV ─────────────────────────────────────────────────────────────
vector<Observation> loadObservations(const string& filename) {
	vector<Observation> obs;
	ifstream file(filename);
	if (!file.is_open()) {
		cerr << "Cannot open " << filename << endl;
		return obs;
	}
	string line;
	getline(file, line);  // skip CSV header
	while (getline(file, line)) {
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



/// ── BUILD H MATRIX (Nx4) ─────────────────────────────────────────────────
// Row i: [ -ux, -uy, -uz, 1 ]  where u = unit vec from receiver to sat i
MatrixXd buildH(const vector<Observation>& obs,
	double rec_x, double rec_y, double rec_z) {
	int N = obs.size();
	MatrixXd H(N, 4);
	for (int i = 0; i < N; i++) {
		double dx = obs[i].sat_x - rec_x;
		double dy = obs[i].sat_y - rec_y;
		double dz = obs[i].sat_z - rec_z;
		double r = sqrt(dx * dx + dy * dy + dz * dz);

		// Note: PARENTHESES, not brackets! Eigen is different from C arrays.
		// H[i][j] does NOT work for Eigen matrices.
		H(i, 0) = -dx / r;
		H(i, 1) = -dy / r;
		H(i, 2) = -dz / r;
		H(i, 3) = 1.0;
	}
	return H;
}


// ── BUILD Δρ VECTOR (Nx1) ────────────────────────────────────────────────
// dr[i] = measured pseudorange - predicted pseudorange at current guess
VectorXd buildDeltaRho(const vector<Observation>& obs,
	double rec_x, double rec_y, double rec_z,
	double rec_dt) {
	int N = obs.size();
	VectorXd dr(N);
	for (int i = 0; i < N; i++) {
		double dx = obs[i].sat_x - rec_x;
		double dy = obs[i].sat_y - rec_y;
		double dz = obs[i].sat_z - rec_z;
		double r = sqrt(dx * dx + dy * dy + dz * dz);
		double rho_predicted = r + C_LIGHT * rec_dt;
		dr(i) = obs[i].pseudorange - rho_predicted;
	}
	return dr;
}


// ── SOLVE LEAST SQUARES Δx = (HᵀH)⁻¹ HᵀΔρ ────────────────────────────────
VectorXd solveLeastSquares(const MatrixXd& H, const VectorXd& dr) {
	// Use ldlt().solve() — numerically stable, avoids explicit inverse
	MatrixXd HtH = H.transpose() * H;
	VectorXd Htdr = H.transpose() * dr;

	// .ldlt() = perform L·D·Lᵀ decomposition (a kind of Cholesky)
    // .solve(...) = solve the linear system using that decomposition
   // Equivalent to (HᵀH)⁻¹ Hᵀ Δρ but FASTER and MORE STABLE
   // This is THE standard idiom in modern C++ numerical code

	return HtH.ldlt().solve(Htdr);
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


// ── MAIN ─────────────────────────────────────────────────────────────────
int main() {
	cout << "=== GPS Least-Squares Position Solver ===" << endl << endl;

	// STEP 1: Load test data
	string filename = "../../data/test_pseudoranges.csv";
	vector<Observation> obs = loadObservations(filename);
	if (obs.size() < 4) {
		cerr << "Need at least 4 satellites! Got " << obs.size() << endl;
		return 1;
	}
	cout << "Loaded " << obs.size() << " pseudorange observations" << endl << endl;

	// STEP 2: Initial guess — center of Earth, no clock bias
	double rec_x = 0.0, rec_y = 0.0, rec_z = 0.0;
	double rec_dt = 0.0;  // receiver clock bias (seconds)

	// STEP 3: Iterate
	cout << "Iterating..." << endl;
	cout << "Iter |    Dposition    |  |Dx|     |  cDt(m)" << endl;
	cout << "------------------------------------------------" << endl;
	cout << fixed << setprecision(3);

	const int MAX_ITER = 20;
	const double TOLERANCE = 0.001;  // 1 mm

	for (int iter = 0; iter < MAX_ITER; iter++) {
		// Build linearized system
		MatrixXd H = buildH(obs, rec_x, rec_y, rec_z);
		VectorXd dr = buildDeltaRho(obs, rec_x, rec_y, rec_z, rec_dt);

		// Solve
		VectorXd dx = solveLeastSquares(H, dr);

		// Update guess
		rec_x += dx(0);
		rec_y += dx(1);
		rec_z += dx(2);
		rec_dt += dx(3) / C_LIGHT;

		// Check convergence (position correction only)
		double dp_mag = sqrt(dx(0) * dx(0) + dx(1) * dx(1) + dx(2) * dx(2));

		cout << setw(4) << iter << " | "
			<< setw(15) << dp_mag << " | "
			<< setw(10) << dp_mag << " | "
			<< setw(10) << dx(3) << endl;

		if (dp_mag < TOLERANCE) {
			cout << "\nConverged after " << iter + 1 << " iterations!" << endl;
			break;
		}
	}

	// STEP 4: Convert to LLA and print results
	double lat_deg, lon_deg, h_m;
	ecefToLla(rec_x, rec_y, rec_z, lat_deg, lon_deg, h_m);

	cout << endl << "=== SOLUTION ===" << endl;
	cout << fixed << setprecision(3);
	cout << "ECEF (m):" << endl;
	cout << "  X = " << setw(15) << rec_x << endl;
	cout << "  Y = " << setw(15) << rec_y << endl;
	cout << "  Z = " << setw(15) << rec_z << endl;
	cout << endl << "LLA:" << endl;
	cout << setprecision(6);
	cout << "  Latitude:  " << lat_deg << " deg" << endl;
	cout << "  Longitude: " << lon_deg << " deg" << endl;
	cout << setprecision(2);
	cout << "  Height:    " << h_m << " m" << endl;
	cout << endl;
	cout << "Receiver clock bias: " << rec_dt * 1e9 << " ns" << endl;
	cout << "  (equivalent distance: " << rec_dt * C_LIGHT << " m)" << endl;

	return 0;
}
