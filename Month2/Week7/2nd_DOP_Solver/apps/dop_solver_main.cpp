// ====================================================//
// // robust_solver_main.cpp — uses the solver library //
// ====================================================//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "gps_ls_solver.h"
#include "utils.h"

using namespace std;

std::vector<Observation> loadObservations(const std::string& filename);

void ecefToLla(double X, double Y, double Z, double& lat, double& lon, double& h);


int main() {
    cout << "=== Robust GPS Least-Squares Solver ===" << endl << endl;

    auto obs = loadObservations("../../data/test_pseudoranges.csv");
    cout << "Loaded " << obs.size() << " observations" << endl;

    // Set elevation for each observation (compute from geometry)
    // For test data with no known receiver, default to 60° (medium weight)
    for (auto& o : obs) o.elevation_deg = 60.0;


    //
    // ── UNWEIGHTED SOLUTION ────────────────────────────────────
    cout << endl << "--- Unweighted LS ---" << endl;
    auto sol1 = solveGpsPosition(obs, false);
    if (!sol1.converged) {
        cerr << "FAILED to converge!" << endl;
        return 1;
    }
    cout << fixed << setprecision(3);
    cout << "  Iterations: " << sol1.iterations << endl;
    cout << "  Position:   (" << sol1.x << ", " << sol1.y << ", " << sol1.z << ")" << endl;
    cout << "  Clock bias: " << sol1.clock_bias * 1e9 << " ns" << endl;
    cout << "  RMS resid:  " << sol1.rms_residual << " m" << endl;
    cout << "  Sigma-0:    " << sol1.sigma0 << " m" << endl;

    // ── WEIGHTED SOLUTION ──────────────────────────────────────
    cout << endl << "--- Weighted LS ---" << endl;
    auto sol2 = solveGpsPosition(obs, true);
    cout << "  Iterations: " << sol2.iterations << endl;
    cout << "  Position:   (" << sol2.x << ", " << sol2.y << ", " << sol2.z << ")" << endl;
    cout << "  Clock bias: " << sol2.clock_bias * 1e9 << " ns" << endl;
    cout << "  RMS resid:  " << sol2.rms_residual << " m" << endl;
    cout << "  Sigma-0:    " << sol2.sigma0 << " m" << endl;

    // ── PER-SAT RESIDUALS ──────────────────────────────────────
    cout << endl << "Per-satellite residuals (unweighted):" << endl;
    cout << "  PRN     residual (m)" << endl;
    cout << "  -------------------" << endl;
    for (size_t i = 0; i < obs.size(); i++) {
        cout << "  G" << setw(2) << setfill('0') << obs[i].prn << setfill(' ')
            << setw(15) << sol1.residuals(i) << endl;
    }

    // ── CONVERT TO LLA ─────────────────────────────────────────
    double lat, lon, h;
    ecefToLla(sol1.x, sol1.y, sol1.z, lat, lon, h);
    cout << endl << "LLA solution:" << endl;
    cout << setprecision(6);
    cout << "  Latitude:  " << lat << " deg" << endl;
    cout << "  Longitude: " << lon << " deg" << endl;
    cout << setprecision(2);
    cout << "  Height:    " << h << " m" << endl;

	// Print DOP quality metrics
	cout << endl << "=== POSITION QUALITY (DOP) ===" << endl;
	cout << fixed << setprecision(2);
	cout << "  GDOP = " << sol1.dop.GDOP << "  (overall)" << endl;
	cout << "  PDOP = " << sol1.dop.PDOP << "  (3D position)" << endl;
	cout << "  HDOP = " << sol1.dop.HDOP << "  (horizontal)" << endl;
	cout << "  VDOP = " << sol1.dop.VDOP << "  (vertical)" << endl;
	cout << "  TDOP = " << sol1.dop.TDOP << "  (time)" << endl;
	 
	// Quality assessment
	if (sol1.dop.PDOP < 3.0)
		cout << "  Quality: EXCELLENT geometry" << endl;
	else if (sol1.dop.PDOP < 6.0)
		cout << "  Quality: Good geometry" << endl;
	else if (sol1.dop.PDOP < 10.0)
		cout << "  Quality: Marginal geometry" << endl;
	else
		cout << "  Quality: POOR geometry — fix may be unreliable" << endl;

    return 0;

}