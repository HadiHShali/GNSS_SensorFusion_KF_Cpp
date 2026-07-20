#include "trajectory_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

vector<TrajectoryEpoch> readTrajectoryCsv(const string& filename)
{
	vector<TrajectoryEpoch> epochs;

	ifstream file(filename);
	if (!file.is_open())
	{
		cerr << "Cannot Open: " << filename << endl;
		return epochs;
	}

	string line;
	getline(file, line); // skip the header

	// main body of the csv file
	while (getline(file, line))
	{
		if (line.empty()) continue;

		stringstream ss(line);
		string cell;
		TrajectoryEpoch e;

		// Format: t_gps,lat_deg,lon_deg,height_m,sigma0_m,pdop,n_sats,converged
		getline(ss, cell, ',');
		e.t_gps = stod(cell);

		getline(ss, cell, ','); 
		e.lat_deg = stod(cell);
		
		getline(ss, cell, ','); 
		e.lon_deg = stod(cell);

		getline(ss, cell, ','); 
		e.height_m = stod(cell);

		getline(ss, cell, ','); 
		e.sigma0_m = stod(cell);

		getline(ss, cell, ','); 
		e.pdop = stod(cell);

		getline(ss, cell, ','); 
		e.n_sats = stoi(cell);

		getline(ss, cell, ','); 
		e.converged = stoi(cell);

		epochs.push_back(e);
	}
	return epochs;
}

