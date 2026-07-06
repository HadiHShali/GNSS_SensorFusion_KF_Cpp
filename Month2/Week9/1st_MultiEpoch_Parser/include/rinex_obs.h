#pragma once
#include <string>
#include <vector>
 
struct ObsRecord {
    int prn;
    double pseudorange;   // C1C in meters
};
 
struct ObsEpoch {
    int year, month, day, hour, minute;
    double second;
    double t_gps;         // seconds of GPS week
    std::vector<ObsRecord> gps_records;   // GPS only
};
 
// Parse ONE epoch from RINEX OBS at the given time (or first epoch if t=0)
ObsEpoch parseRinexObsOneEpoch(const std::string& filename, double target_t_gps);

// parse All Epoch from a RINEX Observation file
std::vector<ObsEpoch> parseRinexAllEpochs(const std::string& filename);
