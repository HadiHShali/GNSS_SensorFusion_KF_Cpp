#pragma once
#include <string>
#include <vector>
#include <map>

// Explaining the std::map<std::string, double[3]> sat_pos;
// A std::map that looks up a satellite's ECEF position by its
// PRN string (like "G06"), storing that position as a 3-element array [X, Y, Z] in meters.
//
// std::map<KeyType, ValueType> — an associative container: give it a key, get back a value
// Key type : std::string — the PRN, e.g., "G06"
// Value type : double[3] — a fixed - size array of 3 doubles
//
//sat_pos["G06"][0] = 14563123.456;  // X
//sat_pos["G06"][1] = -21234654.321; // Y
//sat_pos["G06"][2] = 3456789.012;   // Z

struct Sp3Epoch {
	double t_gps; // seconds of GPS week
	std::map<std::string, double[3]> sat_pos; // PRN -> [X,Y,Z] in meters
};


std::vector<Sp3Epoch> parseSp3(const std::string& filename);