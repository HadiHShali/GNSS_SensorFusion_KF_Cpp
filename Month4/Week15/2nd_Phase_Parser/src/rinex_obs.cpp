#include "rinex_obs.h"
#include <fstream>
#include <sstream>
#include <iostream>
 
using namespace std;
 
// Compute GPS time (seconds of GPS week) from calendar date.
// GPS week starts at midnight Saturday/Sunday (UTC).
static double computeGpsSecOfWeek(int year, int month, int day,
                                   int hour, int min, double sec) {
    // Zeller-style algorithm to find day-of-week
    // Adjust for Jan/Feb being treated as months 13/14 of previous year
    int y = year, m = month;
    if (m < 3) { m += 12; y -= 1; }
    int K = y % 100;
    int J = y / 100;
    
    // h: 0 = Saturday, 1 = Sunday, 2 = Monday, ..., 6 = Friday (Zeller convention)
    int h = (day + (13*(m+1))/5 + K + K/4 + J/4 + 5*J) % 7;
    
    // Convert to GPS day-of-week: 0 = Sunday, 1 = Monday, ..., 6 = Saturday
    int gps_dow = (h + 6) % 7;
    
    return gps_dow * 86400.0 + hour * 3600.0 + min * 60.0 + sec;
}
 
ObsEpoch parseRinexObsOneEpoch(const string& filename, double target_t_gps) {
    ObsEpoch epoch = {};
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open OBS file: " << filename << endl;
        return epoch;
    }
    
    string line;
    // Skip header
    while (getline(file, line)) {
        if (line.find("END OF HEADER") != string::npos) break;
    }
    
    // Find target epoch (or first valid one if target = 0)
    bool found = false;
    int n_sats = 0;
    while (getline(file, line)) {
        if (line.empty() || line[0] != '>') continue;
        
        // Parse epoch line: > YYYY MM DD HH MM SS  flag  N
        try {
            epoch.year   = stoi(line.substr(2, 4));
            epoch.month  = stoi(line.substr(7, 2));
            epoch.day    = stoi(line.substr(10, 2));
            epoch.hour   = stoi(line.substr(13, 2));
            epoch.minute = stoi(line.substr(16, 2));
            epoch.second = stod(line.substr(19, 11));
            n_sats       = stoi(line.substr(32, 3));
        } catch (...) { continue; }
        
        epoch.t_gps = computeGpsSecOfWeek(epoch.year, epoch.month, epoch.day,
                                           epoch.hour, epoch.minute, epoch.second);
        
        // Take first epoch if no target, or one closest to target
        if (target_t_gps == 0.0 || epoch.t_gps >= target_t_gps) {
            found = true;
            break;
        } else {
            // skip n_sats lines and continue
            for (int i = 0; i < n_sats; i++) getline(file, line);
        }
    }
    
    if (!found) {
        cerr << "No epoch found" << endl;
        return epoch;
    }
    
    // Read N satellite observation lines
    for (int i = 0; i < n_sats; i++) {
        if (!getline(file, line)) break;
        if (line.empty() || line.length() < 17) continue;
        
        // GPS records only (start with 'G')
        if (line[0] != 'G') continue;
        
        try {
            ObsRecord rec;
            rec.prn = stoi(line.substr(1, 2));
            // First observable is C1C (pseudorange) — columns 3-16
            rec.pseudorange = stod(line.substr(3, 14));
            if (rec.pseudorange > 1.0e6) {  // sanity check
                epoch.gps_records.push_back(rec);
            }
        } catch (...) { continue; }
    }
    
    return epoch;
}

vector<ObsEpoch> parseRinexAllEpochs(const std::string& filename)
{
	vector<ObsEpoch> epochs;
	ifstream file(filename);
	if (!file.is_open())
	{
		cerr << "Cannot open:" << filename << endl;
		return epochs;
	}

	string line;
	vector<string> gps_obs_types;   // NEW: holds the ordered observable list for GPS ('G')

	//skip the header -- NOW also captures SYS / # / OBS TYPES for GPS
	while (getline(file, line))
	{
		if (line.find("END OF HEADER") != string::npos) break;

		if (line.find("SYS / # / OBS TYPES") != string::npos && line[0] == 'G')
		{
			int n_obs = stoi(line.substr(3, 3));
			gps_obs_types.clear();

			int col = 7;
			int parsed = 0;
			string cur_line = line;

			while (parsed < n_obs)
			{
				if (col + 3 > 60) {
					getline(file, cur_line);
					col = 7;
				}
				string code = cur_line.substr(col, 3);
				code.erase(code.find_last_not_of(' ') + 1);
				gps_obs_types.push_back(code);
				col += 4;
				parsed++;
			}
		}
	}

	// resolve dynamic column indices
	// for pseudoranges:
	int idx_C1W = findObsIndex(gps_obs_types, "C1W");
	int idx_C2W = findObsIndex(gps_obs_types, "C2W");
	//for carrier phase:
	int idx_L1C = findObsIndex(gps_obs_types, "L1C");
	int idx_L2W = findObsIndex(gps_obs_types, "L2W");
	
	if (idx_C1W < 0) idx_C1W = findObsIndex(gps_obs_types, "C1C");  // fallback

	cout << "Resolved indices: C1W=" << idx_C1W << " C2W=" << idx_C2W
		 << " L1C=" << idx_L1C << " L2W=" << idx_L2W << endl;


	int col_C1 = (idx_C1W >= 0) ? (3 + 16 * idx_C1W) : -1;   // fallback to old hard-coded 3
	int col_C2 = (idx_C2W >= 0) ? (3 + 16 * idx_C2W) : -1;
	int col_L1 = (idx_L1C >= 0) ? (3 + 16 * idx_L1C) : -1;   // fallback to old hard-coded 3
	int col_L2 = (idx_L2W >= 0) ? (3 + 16 * idx_L2W) : -1;

	int lli_col_L1 = col_L1 + 14;   // 1 character, '0'-'9' or blank
	int lli_col_L2 = col_L2 + 14;

	//parse epoch by epoch
	while (getline(file, line))
	{
		if (line.empty() || line[0] != '>') continue;

		ObsEpoch epochh;
		epochh.year = stoi(line.substr(2, 4));
		epochh.month = stoi(line.substr(7, 2));
		epochh.day = stoi(line.substr(10, 2));
		epochh.hour = stoi(line.substr(13, 2));
		epochh.minute = stoi(line.substr(16, 2));
		epochh.second = stod(line.substr(19, 11));

		int n_satss = stoi(line.substr(32, 3));

		epochh.t_gps = epochh.hour * 3600.0
		            + epochh.minute * 60.0
		            + epochh.second;

		for (int i = 0; i < n_satss; i++)
		{
			getline(file, line);
			if (line.empty() || line[0] != 'G') continue;

			ObsRecord recc;
			recc.prn = stoi(line.substr(1, 2));

			// C1W pseudorange (dynamic column, was hard-coded 3)
			try
			{
				recc.pseudorange = stod(line.substr(col_C1, 14));
			}
			catch (...) { continue; }

			// C2W pseudorange (dynamic column, optional)
			if (col_C2 >= 0 && (int)line.size() >= col_C2 + 14) {
				try {
					string c2_field = line.substr(col_C2, 14);
					c2_field.erase(0, c2_field.find_first_not_of(' '));
					if (!c2_field.empty()) {
						recc.pseudorange_c2 = stod(line.substr(col_C2, 14));
						recc.has_c2 = true;
					}
				}
				catch (...) { recc.has_c2 = false; }
			}

			// New: L1C phase (cycles + LLI)
			if (col_L1 >=0 && (int)line.size() >= col_L1 + 14)
			{
				try 
				{
					string l1_field = line.substr(col_L1, 14);
					l1_field.erase(0, l1_field.find_first_not_of(' '));
					if (!l1_field.empty())
					{
						recc.phase_l1 = stod(line.substr(col_L1, 14));
						recc.has_phase_l1 = true;
						if ((int)line.size() > lli_col_L1)
						{
							char lli_ch = line[lli_col_L1];
							recc.lli_l1 = (lli_ch == ' ') ? 0 : (lli_ch - '0');
						}
					}
					
				} catch (...) {recc.has_phase_l1 = false;}
				
			}
			
			// New: L2W phase (cycles + LLI)
			
			if (col_L2 >=0 && (int)line.size() >= col_L2 + 14)
			{
				try 
				{
					string l2_field = line.substr(col_L2, 14);
					l2_field.erase(0, l2_field.find_first_not_of(' '));
					if (!l2_field.empty())
					{
						recc.phase_l2 = stod(line.substr(col_L2, 14));
						recc.has_phase_l2 = true;
						if ((int)line.size() > lli_col_L2)
						{
							char lli_ch = line[lli_col_L2];
							recc.lli_l2 = (lli_ch == ' ') ? 0 : (lli_ch - '0');
						}
					}
					
				} catch (...) {recc.has_phase_l2 = false;}
				
			}
			epochh.gps_records.push_back(recc);
		}

		if (!epochh.gps_records.empty())
		{
			epochs.push_back(epochh);
		}
	}

	return epochs;
}


//static_cast<int>(i) converts i (which is size_t, unsigned) into a plain int (signed)
//, so it matches the function's declared return type int. This conversion is needed 
// because the function needs to be able to return -1 (a negative number) for "not found," and size_t can't represent negative values at all.
int findObsIndex(const vector<string>& obs_types, const string& code)
{
    for (size_t i = 0; i < obs_types.size(); i++) {
        if (obs_types[i] == code) return static_cast<int>(i);
    }
    return -1;
}