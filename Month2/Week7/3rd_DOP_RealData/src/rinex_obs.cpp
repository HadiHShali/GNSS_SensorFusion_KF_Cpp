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
