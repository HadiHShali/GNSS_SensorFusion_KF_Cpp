#include "kalman_filter_gps.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::MatrixXi;


KalmanFilterGps::KalmanFilterGps(double q_accel,
	double q_clk_bias,
	double q_clk_drift)
	: q_accel_(q_accel),
	q_clk_bias_(q_clk_bias),
	q_clk_drift_(q_clk_drift),
	x_(VectorXd::Zero(8)),
	P_(MatrixXd::Identity(8, 8) * 1e6) {  // huge initial uncertainty
}


//The '::' means: The function 'initialize' belongs to the class 'KalmanFilterGps'
void KalmanFilterGps::initialize(const VectorXd& x0, const MatrixXd& P0)
{
	if (x0.size() != 8 || P0.rows() != 8 || P0.cols() != 8)
	{
		throw std::runtime_error("KF initialize: dimensions must be 8");
	}
	x_ = x0;
	P_ = P0;
	initialized_ = true;
}

// --------- buildF() helper ----------------//
MatrixXd KalmanFilterGps::buildF(double dt) const
{
	// 8x8 identity + dt in position-velocity + clock-drift couplings
	MatrixXd F = MatrixXd::Identity(8, 8);
	F(0, 3) = dt;  // x += vx * dt
	F(1, 4) = dt;  // y += vy * dt
	F(2, 5) = dt;  // z += vz * dt
	F(6, 7) = dt;  // dt_bias += dt_drift * dt
	return F;
}

// --------- buildQ() helper ----------------//
MatrixXd KalmanFilterGps::buildQ(double dt) const
{
	double dt2 = dt * dt;
	double dt3 = dt2 * dt;
	double q_a_sq = q_accel_ * q_accel_;

	// 8x8 zero matrix
	MatrixXd Q = MatrixXd::Zero(8, 8);

	// position and velocity blocks (3 axes: x, y, z) 
	for (int i; i < 3; i++)
	{
		Q(i, i) = q_a_sq * dt3 / 3.0; //pos-pos
		Q(i+3, i+3) = q_a_sq * dt; //vel-vel
		Q(i, i+3) = q_a_sq * dt2 / 2.0; //pos-vel cross
		Q(i+3, i) = q_a_sq * dt2 / 2.0; //vel-pos cross
	}
	// clock bias and drift block (states 6 and 7)
	Q(6, 6) = q_clk_bias_ * dt + q_clk_drift_ * dt3 / 3.0; 
	Q(7, 7) = q_clk_drift_ * dt;
	Q(6, 7) = q_clk_drift_ * dt2 / 2.0;
	Q(7, 6) = q_clk_drift_ * dt2 / 2.0;
	return Q;
}

// --------- predict() helper ----------------//
void KalmanFilterGps::predict(double dt)
{
	if (!initialized_)
	{
		throw std::runtime_error("KF predict: not initialized");
	}

	MatrixXd F = buildF(dt);
	MatrixXd Q = buildQ(dt);

	// the tow prediction equations:
	x_ = F * x_;
	P_ = F * P_ * F.transpose() + Q; 
}

// --------- update() helper ----------------//
void KalmanFilterGps::update(const Vector3d& z, const Matrix3d& R)
{
	if (!initialized_)
	{
		throw std::runtime_error("KF update: not initialized");
	}

	// H = Observation Matrix (3*8)
	// Measurement is z = [x, y, z] - first 3 states directly observed
	MatrixXd H = MatrixXd::Zero(3, 8);
	H(0, 0) = 1.0; // z_x = state_x
	H(1, 1) = 1.0; // z_y = state_z
	H(2, 2) = 1.0; // z_z = state_z

	// Innovation (measurement residual)
	Vector3d y = z - H * x_;

	// Innovation Covariance
	Matrix3d S = H * P_ * H.transpose() + R;

	// Kalman gain
	MatrixXd K = P_ * H.transpose() * S.inverse();

	// state + covariance update
	x_ = x_ + K * y;
	P_ = (MatrixXd::Identity(8, 8) - K * H) * P_; 

}

// -------------------RTS Smoother------------------------//
void rtsSmooth(
    const std::vector<VectorXd>& x_filt,
    const std::vector<MatrixXd>& P_filt,
    const std::vector<VectorXd>& x_pred,
    const std::vector<MatrixXd>& P_pred,
    const std::vector<MatrixXd>& F_used,
    std::vector<VectorXd>& x_smooth,
    std::vector<MatrixXd>& P_smooth)
{
    int N = static_cast<int>(x_filt.size());
    x_smooth.resize(N);
    P_smooth.resize(N);
    
    // Last epoch: smoothed = filtered (no future data available)
    x_smooth[N - 1] = x_filt[N - 1];
    P_smooth[N - 1] = P_filt[N - 1];
    
    // Sweep backward: k = N-2 down to 0
    for (int k = N - 2; k >= 0; --k) {
        // F_used[k+1] is the F that was used to go FROM epoch k
        //   TO the prediction x_pred[k+1]
        const MatrixXd& F = F_used[k + 1];
        
        // Smoother gain: C_k = P_filt[k] * F^T * inv(P_pred[k+1])
        MatrixXd Ck = P_filt[k] * F.transpose() * P_pred[k + 1].inverse();
        
        // Smoothed state
        x_smooth[k] = x_filt[k] +
            Ck * (x_smooth[k + 1] - x_pred[k + 1]);
        
        // Smoothed covariance
        P_smooth[k] = P_filt[k] +
            Ck * (P_smooth[k + 1] - P_pred[k + 1]) * Ck.transpose();
    }
}

