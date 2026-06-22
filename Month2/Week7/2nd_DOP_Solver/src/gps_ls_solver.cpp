// =====================================================//
// gps_ls_solver.cpp robust GPS LS solver              //
// =====================================================//

#include "gps_ls_solver.h"
#include <cmath>

using namespace Eigen; // imports EVERYTHING from Eigen folder

const double C_LIGHT = 299792458.0;
const double PI = 3.14159265358979;


// Build H Matrix (N*4)
// static: So whenever you see static at the start of a function in a .cpp file, 
// just translate it to: "this is a private helper for this file only." That's the whole concept.
static MatrixXd buildH(const std::vector<Observation>& obs, double rec_x, double rec_y, double rec_z)
{
	int N = obs.size();
	MatrixXd H(N, 4);
	for (int i = 0; i < N; i++)
	{
		double dx = obs[i].sat_x - rec_x;
		double dy = obs[i].sat_y - rec_y;
		double dz = obs[i].sat_z - rec_z;
		double r = sqrt(dx * dx + dy * dy + dz * dz);
		H(i, 0) = -dx / r;
		H(i, 1) = -dy / r;
		H(i, 2) = -dz / r;
		H(i, 3) = 1.0;
	}
	return H;
}

// Build Δρ vector (Nx1)
static VectorXd buildDeltaRho(const std::vector<Observation>& obs,
	double rec_x, double rec_y, double rec_z,
	double rec_dt) 
{
	int N = obs.size();
	VectorXd dr(N);
	for (int i = 0; i < N; i++) {
		double dx = obs[i].sat_x - rec_x;
		double dy = obs[i].sat_y - rec_y;
		double dz = obs[i].sat_z - rec_z;
		double r = sqrt(dx * dx + dy * dy + dz * dz);
		double rho_predicted = r + C_LIGHT * rec_dt;
		dr(i) = obs[i].pseudorange - rho_predicted;
	}
	return dr;
}

// Build Weight Matrix W (N*N diagonal matrix)
static MatrixXd buildW(const std::vector<Observation>& obs)
{
	int N = obs.size();
	MatrixXd W = MatrixXd::Zero(N, N);
	for (int i = 0; i < N; i++)
	{
		double s = sin(obs[i].elevation_deg * PI / 180.0);
		W(i, i) = s * s; // Weight = sin^2(elevation)
	}
	return W;
}

// compute all five DOP values from the Design matrix H
static DopValues computeDOP(const MatrixXd& H)
{
	DopValues dop = {};

	// Q = (H'H)^-1 is the covariance matrix
	MatrixXd Q = (H.transpose() * H).inverse();

	// pull out diagonal elements: [X, Y, Z, clock bias]
	double qxx = Q(0, 0);
	double qyy = Q(1, 1);
	double qzz = Q(2, 2);
	double qtt = Q(3, 3);
	
	dop.HDOP = sqrt(qxx + qyy);
	dop.VDOP = sqrt(qzz);
	dop.PDOP = sqrt(qxx + qyy + qzz);
	dop.TDOP = sqrt(qtt);
	dop.GDOP = sqrt(qxx + qyy+ qzz+ qtt);

	return dop;
}


// Solver
PositionSolution solveGpsPosition(const std::vector<Observation>& obs, bool use_weights)
{
	PositionSolution sol = {};  // zero init

	// edge case: need at least 4 satellites
	if (obs.size() < 4)
	{
		sol.converged = false;
		return sol;
	}

	// initial guess: center of the Earth
	sol.x = 0.0; sol.y = 0.0; sol.z = 0.0; sol.clock_bias = 0.0;

	const int MAX_ITER = 20;
	const double TOLERANCE = 1e-3; // 1mm

	for (int iter = 0; iter < MAX_ITER; iter++)
	{
		MatrixXd H = buildH(obs, sol.x, sol.y, sol.z);
		VectorXd dr = buildDeltaRho(obs, sol.x, sol.y, sol.z, sol.clock_bias);

		VectorXd dx;
		
		if (use_weights)
		{
			MatrixXd W = buildW(obs);
			MatrixXd HtWH = H.transpose() * W * H;
			VectorXd HtWdr = H.transpose() * W * dr;
			dx = HtWH.ldlt().solve(HtWdr);
		}
		else
		{
			MatrixXd HtWH = H.transpose() * H;
			VectorXd HtWdr = H.transpose() * dr;
			dx = HtWH.ldlt().solve(HtWdr);
		}

		sol.x += dx(0);
		sol.y += dx(1);
		sol.z += dx(2);
		sol.clock_bias += dx(3) / C_LIGHT;
		sol.iterations = iter + 1;

		double dp_mag = sqrt(dx(0) * dx(0)+ dx(1) * dx(1)+ dx(2) * dx(2)); //dp_mag is NOT a pseudorange. It's the size of the position update in this iteration:
		if (dp_mag < TOLERANCE)
		{
			sol.converged = true;
			break;
		}
	}

	// Compute the final residuals and statistics
	sol.residuals = buildDeltaRho(obs, sol.x, sol.y, sol.z, sol.clock_bias);
	int N = obs.size();
	sol.rms_residual = sqrt(sol.residuals.squaredNorm() / N);
	sol.sigma0 = (N > 4) ? sqrt(sol.residuals.squaredNorm() / (N - 4)) : 0.0;
	MatrixXd H_final = buildH(obs, sol.x, sol.y, sol.z);
	sol.dop = computeDOP(H_final);

	return sol;

}

