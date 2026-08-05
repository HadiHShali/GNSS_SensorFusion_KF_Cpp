#pragma once
#include <string>
#include <vector>
#include <map>

struct ClkEpoch {
	double t_gps; 
	std::map<std::string, double> sat_clock; // PRN -> clock bias in seconds
};

std::vector<ClkEpoch> parseClk(const std::string& filename);
