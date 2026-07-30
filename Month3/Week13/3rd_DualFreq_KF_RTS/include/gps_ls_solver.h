// ==========================================================//
// gps_ls_solver.h - public interface for the GPS LS solver  //
// ==========================================================//
// API : Application Programming Interface
#pragma once
#include <vector>
#include <Eigen/Dense>

struct Observation
{
	int prn; 
	double sat_x, sat_y, sat_z;
	double pseudorange;
	double elevation_deg = 90.0; // optional, default = overhead
};

struct DopValues
{
	double GDOP; // Geometric (all 4 unknowns)
	double PDOP; // Position (X, Y, Z)
	double HDOP; // Horizontal (X, Y)
	double VDOP; // Vertical (Z)
	double TDOP; // Time (Clock Bias)
};

struct PositionSolution
{
	double x, y, z; // ECEF in m
	double clock_bias;  // seconds
	int iterations;
	bool converged;
	double rms_residual; // Root Mean Square of final residuals in m
	double sigma0;  // residual standard deviation
	Eigen::VectorXd residuals;
	DopValues dop;
};

// Main API
// An API is the set of rules for how OTHER code can talk to YOUR code.
// So whenever you hear or read "API" in software contexts, just substitute:
//		"a clean interface for using somebody's code." That's the whole concept.
// The best way to understand API is to think of a restaurant:
//				You(the customer) → MENU → Kitchen(the chef)
//
//									  ↑
//							   This is the API
//You don't need to know HOW the chef cooks the food. 
// You just need the menu to know what you can order. Same with code.
PositionSolution solveGpsPosition(const std::vector<Observation>& obs, bool use_weights = false);

// std::vector<Observati....  
// std::vector: the vector class that's part of the standard library
// std -> standard library
// 
// So next time you see std:: anywhere, just translate it in your head as "this comes from the standard library" 
// — and don't be confused by the double-colon syntax. It's just C++'s way of saying "look in this folder."
// 
// :: -> The :: is called the "scope resolution operator" — it's how C++ asks "which namespace?"
// std::vector<int> nums;     // standard library's vector
// Eigen::Vector3d vec;        // Eigen's vector


// Namespaces — The C++ "Folder System" : Think of namespaces like folders for code:
// C++ Universe/
//├── std / ← standard library
//│   ├── vector
//│   ├── string
//│   ├── cout
//│   └── map
//├── Eigen / ← Eigen library
//│   ├── MatrixXd
//│   ├── VectorXd
//│   └── Vector3d
//└──(your code) / ← your own code
//└── Observation
