#include <iostream>
#include <iomanip>
#include "sp3_reader.h"
#include "clk_reader.h"
#include "interpolation.h"

using namespace std;

int main() {
    cout << "=== SP3/CLK PARSER TEST ===" << endl;

    // Parse SP3 orbit file
    auto sp3_epochs = parseSp3("../../data/orbit.sp3");
    cout << "Parsed " << sp3_epochs.size() << " SP3 epochs" << endl;

    // Parse CLK file
    auto clk_epochs = parseClk("../../data/clock.clk");
    cout << "Parsed " << clk_epochs.size() << " CLK epochs" << endl;

    if (sp3_epochs.empty() || clk_epochs.empty()) {
        cerr << "Parsing failed -- check file paths!" << endl;
        return 1;
    }

    // Sanity check: print first epoch's G06 position (compare to broadcast
    // ephemeris position for the same satellite/time as a rough cross-check)
    string test_prn = "G06";
    if (sp3_epochs[0].sat_pos.count(test_prn)) {
        auto& pos = sp3_epochs[0].sat_pos[test_prn];
        cout << fixed << setprecision(3);
        cout << test_prn << " first epoch ECEF (m): "
            << pos[0] << ", " << pos[1] << ", " << pos[2] << endl;
    }

    if (clk_epochs[0].sat_clock.count(test_prn)) {
        cout << test_prn << " first epoch clock (s): "
            << scientific << clk_epochs[0].sat_clock[test_prn] << endl;
    }

    // Test Lagrange interpolation: interpolate G06's position at a
    // time BETWEEN two SP3 epochs (e.g., 7.5 minutes after epoch 0)
    double target_t = sp3_epochs[0].t_gps + 450.0;  // 7.5 min later
    double interp_pos[3];
    bool ok = interpolateSatPosition(sp3_epochs, test_prn, target_t, interp_pos);
    if (ok) {
        cout << fixed << setprecision(3);
        cout << test_prn << " interpolated ECEF at t=" << target_t
            << ": " << interp_pos[0] << ", " << interp_pos[1]
            << ", " << interp_pos[2] << endl;
    }

    cout << "Test complete." << endl;
    return 0;
}
