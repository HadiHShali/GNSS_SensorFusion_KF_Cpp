#pragma once
#include <string>
#include <vector>

// One row of the trajectory.csv
struct TrajectoryEpoch
{
	double t_gps; // GPS seconds of week
	double lat_deg; // latitude in degrees
	double lon_deg; // longitude degrees
	double height_m; // ellipsoida height (meters)
	double sigma0_m; // solve sigma-0 (used to weight R)
	double pdop;
	int n_sats;
	int converged;

};

// parse trajectory.csv from week 9 SPP output
std::vector<TrajectoryEpoch> readTrajectoryCsv(const std::string & filename);