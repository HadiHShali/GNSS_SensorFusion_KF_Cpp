#include <clk_reader.h>
#include <fstream>
#include <sstream>
using namespace std;

vector<ClkEpoch> parseClk(const string& filename)
{
	vector<ClkEpoch> epochs;

	ifstream file(filename);

	string line;

	map<double, ClkEpoch> by_time; // group AS lines sharing the same epoch
	// groups CLK file records by their timestamp — so all satellites' clock corrections
	// for the same epoch end up bundled together in one ClkEpoch, even though the raw
	// file lists them as separate AS lines, one per satellite.
	// Key type: double — the epoch's GPS time (t_gps)
	// Value type: ClkEpoch — your struct holding t_gps plus the sat_clock map (PRN → clock bias)

    while (getline(file, line)) {
        if (line.substr(0, 2) != "AS") continue;   // skip AR (receiver) + header lines

        istringstream ss(line);
        string tag, prn;
        int yr, mo, dy, hr, mi, nvals;
        double sec, clk_bias;

        //This line reads 10 space-separated pieces of text off the line, one after another,
        // automatically converting each into the correct variable type (string, int, or double) as it goes
        ss >> tag >> prn >> yr >> mo >> dy >> hr >> mi >> sec >> nvals >> clk_bias;
        //Exa: ↑      ↑      ↑    ↑     ↑
        //     AS    G01    2024  1     7     0     0  0.000000   2   -1.234567890123E-04   1.2E-11 (we don't read the last one)
 
        // clk_bias is ALREADY in seconds -- no unit conversion needed here

        double t_gps = hr * 3600.0 + mi * 60.0 + sec;
        by_time[t_gps].t_gps = t_gps;
        by_time[t_gps].sat_clock[prn] = clk_bias;
    }
    for (auto& [t, ep] : by_time) epochs.push_back(ep);
    return epochs;
}
