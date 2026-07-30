#pragma once
#include <string>
#include <vector>
 
struct ObsRecord {
    int prn;
    double pseudorange;   // C1W in meters (existing — keep as-is)
    double pseudorange_c2 = 0.0;   // NEW: C2W in meters
    bool has_c2 = false;           // NEW: was C2W present this epoch?
};
 
struct ObsEpoch {
    int year, month, day, hour, minute;
    double second;
    double t_gps;         // seconds of GPS week
    std::vector<ObsRecord> gps_records;   // GPS only
};
 
 
// Column position for observable index i in a RINEX 3 obs line
// (short + trivial -> safe to define directly here)
inline int obsColumn(int index) {
    return 3 + 16 * index;   // 14-char value + 2-char LLI/SNR flags per field
}



// Search a parsed observable-type list for a given code (e.g. "C1W") 
// Returns -1 if not found. (declaration only -- body goes in the .cpp) 

//Given the list of observable codes you parsed from the header (e.g., 
// ["C1C", "L1C", "D1C", "S1C", "C1W", "S1W", "C2W", ...]) 
// and a target string like "C2W", it searches through the list
// and returns the position where that string is found — or -1 if it's not in the list at all.
int findObsIndex(const std::vector<std::string>& obs_types, const std::string& code);


// Parse ONE epoch from RINEX OBS at the given time (or first epoch if t=0)
ObsEpoch parseRinexObsOneEpoch(const std::string& filename, double target_t_gps);

// parse All Epoch from a RINEX Observation file
std::vector<ObsEpoch> parseRinexAllEpochs(const std::string& filename);
