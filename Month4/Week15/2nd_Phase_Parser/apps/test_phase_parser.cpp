#include <iostream>
#include <iomanip>
#include "rinex_obs.h"
#include "cycle_slip.h"

using namespace std;

int main()
{
	cout << "=== PHASE PARSER + CYCLE SLIP TEST ===" << endl;
	auto epochs = parseRinexAllEpochs("../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx");
	cout << "Parsed " << epochs.size() << " epochs" << endl;

	// Track G06's LG combination across all epochs, flag slips
	double lg_prev = 0.0;
	bool have_prev = false;
	int slip_count = 0, lli_flag_count = 0;

	for (const auto& ep : epochs)
	{
		for (const auto& rec : ep.gps_records)
		{
			if (rec.prn != 6) continue; // G06 only, for a focused first look
			if (!rec.has_phase_l1 || !rec.has_phase_l2) continue;

			double lg = geometryFree(rec.phase_l1, rec.phase_l2);

			// TEMP DEBUG — add inside the loop, right after computing lg:
			if (rec.prn == 6 && ep.t_gps < 90) {
				cout << "Test the estimated LG:"
					<< "  t=" << ep.t_gps
					<< "  L1=" << fixed << setprecision(3) << rec.phase_l1
					<< "  L2=" << rec.phase_l2
					<< "  LG=" << lg << endl;
			}

			if (have_prev && likelySlip(lg_prev, lg)) {
				slip_count++;
				cout << "  Possible slip at t=" << ep.t_gps
					 << "  LG jump=" << (lg - lg_prev) << " m"
					 << "  LLI_L1=" << rec.lli_l1
					 << "  LLI_L2=" << rec.lli_l2 << endl;
			}
			if (rec.lli_l1 & 1) lli_flag_count++;  // bit 1 = possible slip
			lg_prev = lg;
			have_prev = true;
		}
	}
	cout << "G06: " << slip_count << " LG-detected slips, " << lli_flag_count << " RINEX LLI flags" << endl;
	return 0;
}