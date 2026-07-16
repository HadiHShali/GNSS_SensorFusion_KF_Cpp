// test the Kalman Filter for a static point at the center (0, 0, 0)

// ==============  Section1: Header + Setup =======================//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include "kalman_filter_gps.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;

int main()
{
	cout << "=== KF Static Station Test ===" << endl;

	// Truth: Static station at origin
	const double TRUE_X = 0.0, TRUE_Y = 0.0, TRUE_Z = 0.0;
	const double MEAS_NOISE = 10.0; // 10m 1-sigma
	const int N_EPOCHS = 300;
	const double DT = 1.0; 

	// Random Noise Genetrator (fixed seed = reproducible)
	default_random_engine gen(42);
	normal_distribution<double> noise(0.0, MEAS_NOISE);

	// Initialization KF (static station -> small q_accel)

	KalmanFilterGps kf(0.01); // its equivalent to write: KalmanFilterGps kf(0.01, 0.01, 0.04);
	// Whithin the header file, we set the default value. So, we can set the first one and leave the defaults 
	VectorXd x0 = VectorXd::Zero(8);
	MatrixXd P0 = MatrixXd::Identity(8, 8) * 100.0; 
	kf.initialize(x0, P0); 

	// Measurement noise covariance
	Matrix3d R = Matrix3d::Identity() * (MEAS_NOISE * MEAS_NOISE);

	// Output CSV file
	ofstream csv("../../data/kf_static_test.csv");
	csv << "epoch,time,truth_x,meas_x,kf_x,P00" << endl;

	// ==============  Section2: Main Loop =======================// 
	for (int i = 0; i < N_EPOCHS; i++)
	{
		double t = i * DT; 

		// Generate Noisy Measurement of Truth
		Vector3d z;
		z(0) = TRUE_X + noise(gen); 
		z(1) = TRUE_Y + noise(gen);
		z(2) = TRUE_Z + noise(gen);

		// Run KF: Predict then Update (skip predict on first epoch)
		if (i > 0) kf.predict(DT);
		kf.update(z, R);

		VectorXd x_est = kf.getState();
		MatrixXd P_est = kf.getCovariance();

		csv << i << "," << t << "," << TRUE_X << "," << z(0) << "," << x_est(0) << "," << P_est(0, 0) << endl;

		// print first and every 50th epoch
		if (i == 0 || (i + 1) % 50 == 0)
		{
			cout << fixed << setprecision(3);
			cout << "Epoch " << setw(3) << i << " z_x = " << setw(7) << z(0)
				<< " kf_x = " << setw(7) << x_est(0)
				<< " sqrt(P00) = " << setw(6) << sqrt(P_est(0, 0))
				<< endl;

		}
	}

	csv.close();

	// Final Summary
	VectorXd xf = kf.getState();
	MatrixXd Pf = kf.getCovariance();
	cout << endl << "=== FINAL STATE ===" << endl;
	cout << " x = " << xf(0) << " m (truth = 0)" << endl;
	cout << " y = " << xf(1) << " m (truth = 0)" << endl;
	cout << " z = " << xf(2) << " m (truth = 0)" << endl;
	cout << " Final 1-sigma = " << sqrt(Pf(0, 0)) << " m" << endl;
	cout << endl << "CSV written: kf_static_test.csv" << endl;
	return 0;

}