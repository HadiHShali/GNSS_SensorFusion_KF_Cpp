#include "sp3_reader.h"
#include <fstream>
#include <sstream>

using namespace std;

vector<Sp3Epoch> parseSp3(const string& filename) {
	
	vector<Sp3Epoch> epochs;

	ifstream file(filename);
	string line;
	Sp3Epoch current;
	bool have_epoch = false;

	while (getline(file, line))
	{
		if (line.empty()) continue;

		if (line[0] == '*') //'*' : start of a new epoch
		{
			if (have_epoch) epochs.push_back(current);
			current = Sp3Epoch{};
			int yr = stoi(line.substr(3, 4)), mo = stoi(line.substr(8, 2));
			int dy = stoi(line.substr(11, 2)), hr = stoi(line.substr(14, 2));
			int mi = stoi(line.substr(17, 2));
			double sec = stod(line.substr(20, 10));
			current.t_gps = hr * 3600.0 + mi * 60.0 + sec;  // simplified -- see Day 1 GPS-week note
			have_epoch = true;
		}
		else if (line[0] == 'P')  // Sat Position line: PG01
		{
			string prn = line.substr(1, 3);  // e.g. "G01"
			double x_km = stod(line.substr(4, 14));
			double y_km = stod(line.substr(18, 14));
			double z_km = stod(line.substr(32, 14));
			// Critical : Conver km to meter
			current.sat_pos[prn][0] = x_km * 1000.0;
			current.sat_pos[prn][1] = y_km * 1000.0;
			current.sat_pos[prn][2] = z_km * 1000.0;

		}
	}
	if (have_epoch) epochs.push_back(current);
	
	return epochs;
}