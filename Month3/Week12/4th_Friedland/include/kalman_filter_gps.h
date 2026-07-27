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
	// something like this is not allowed: x_(0) = 10;   
	// but  reading from x_ is perfectly fine:	return x_;  
	Eigen::MatrixXd getLastGain() const { return K_; }
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

	Eigen::MatrixXd K_;
	
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


// Friedland separate-bias estimator (measurement-only bias, Bb=0, Cb=I)
// Runs IN PARALLEL with the main forward KF loop, using its intermediate
// outputs (L_k, P_pred, z_hat). Estimates a near-constant 3x1 ECEF bias.
class BiasEstimator {
public:
    BiasEstimator();  // initializes bhat=0, SigmaB=large
    
    // Call ONCE per epoch, AFTER kf.predict() but paired with kf.update()
    //   F        : same F used by kf.predict() this step
    //   H        : measurement matrix (3x8 position selector, same as kf.update)
    //   L        : Kalman gain just computed by kf.update() (8x3)
    //   P_pred   : predicted covariance BEFORE this update (8x8)
    //   R        : measurement noise used this step (3x3)
    //   z        : actual measurement (3x1)
    //   z_hat    : predicted measurement = H * x_pred (3x1)
    void step(const Eigen::MatrixXd& F, const Eigen::MatrixXd& H,
              const Eigen::MatrixXd& L, const Eigen::MatrixXd& P_pred,
              const Eigen::Matrix3d& R, const Eigen::Vector3d& z,
              const Eigen::Vector3d& z_hat);
    
    Eigen::Vector3d getBias()       const { return bhat_; }
    Eigen::Matrix3d getBiasCov()    const { return SigmaB_; }
    Eigen::VectorXd getDelta()      const { return delta_; }  // 8x1 state correction
private:
    Eigen::Vector3d bhat_;      // current bias estimate
    Eigen::Matrix3d SigmaB_;    // bias covariance
    Eigen::MatrixXd V_;         // 8x3, V(0)=0
    Eigen::VectorXd delta_;     // 8x1 correction = V * bhat
};
