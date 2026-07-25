#pragma once
#include <Eigen/Dense>
#include <vector>

// Kalman filter and RTS Smoother for GPS positioning
//
// State vector (8x1): [x, y, z, vx, vy, vz, dt_bias, dt_drift]
class KalmanFilterGps {
public:
    // Constructor - takes noise tuning parameters
    KalmanFilterGps(double q_accel = 0.01,
                    double q_clk_bias  = 0.01,
                    double q_clk_drift = 0.04);
    
    // Initialize with starting state and covariance
    void initialize(const Eigen::VectorXd& x0,
                    const Eigen::MatrixXd& P0);
    
    // Time update: extrapolate state forward by dt seconds
    void predict(double dt);


	// Measurement update: correct using position measurement and its uncertainty
	// z (3*1) vector : [x_measured, y_measured, z_measured]
	// R (3*3) matrix : [measurement noise covariance]
	void update(const Eigen::Vector3d& z, const Eigen::Matrix3d& R);

	// Accessors
	Eigen::VectorXd getState() const { return x_; }
	// This means calling getState() will not modify the object. So inside the function, 
	// something like this is not allowed: x_(0) = 10;    // ❌ Error
	// but  reading from x_ is perfectly fine:	return x_;     // ✅ OK

	Eigen::MatrixXd getCovariance() const { return P_; }
	bool isInitialized() const { return initialized_; }

private:
	// state + covariance
	Eigen::VectorXd x_;  //8*1 state
	Eigen::MatrixXd P_;  // 8*8 covariance
	bool initialized_ = false;

	// Noise tuning parameters
	double q_accel_;
	double q_clk_bias_;
	double q_clk_drift_;

	// Helpers - construct F and Q for a given dt
	Eigen::MatrixXd buildF(double dt) const;  //Use const after a member function whenever it only reads the object's state and doesn't change it.
	Eigen::MatrixXd buildQ(double dt) const; 

};

// Eigen::VectorXd getState() const     >>>   Returns a copy of the state; does not modify the object.
// const Eigen::VectorXd& getState() const   >>>  Returns a read-only reference to the internal state; avoids copying.
//


// RTS Smoother
// Rauch-Tung-Striebel backward smoother
// Call AFTER running the forward KF and saving all 4 arrays.
// F must be the SAME state transition matrix used at each
// corresponding forward step (pass per-epoch if dt varies).
void rtsSmooth(
    const std::vector<Eigen::VectorXd>& x_filt,
    const std::vector<Eigen::MatrixXd>& P_filt,
    const std::vector<Eigen::VectorXd>& x_pred,   // x_pred[k] used going INTO epoch k
    const std::vector<Eigen::MatrixXd>& P_pred,
    const std::vector<Eigen::MatrixXd>& F_used,    // F used at each step (epoch 0 unused)
    std::vector<Eigen::VectorXd>& x_smooth,        // OUTPUT
    std::vector<Eigen::MatrixXd>& P_smooth);       // OUTPUT

