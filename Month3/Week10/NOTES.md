# Week 10 — Kalman Filter Foundation
 
## Goal
Build a working Kalman filter class in modern C++17 for GPS positioning,
with a clean public API and standalone verification.
 
## What was built
- 5 canonical Kalman equations documented
- 8-state vector: [x, y, z, vx, vy, vz, dt_bias, dt_drift]
- State transition matrix F (identity + dt couplings)
- Process noise Q (continuous white noise acceleration model)
- KalmanFilterGps class (~200 lines with Eigen)
  - Public: initialize(), predict(), update(), getState(), getCovariance()
  - Private: buildF(dt), buildQ(dt), state x_, covariance P_
- Synthetic static-station test
- Python plot script (kf_static_test.png)
 
## Test results
- 300 epochs of measurements: sigma = 10 m Gaussian noise
- Final position error: ~0.1 m (100x improvement over single measurement)
- Final 1-sigma uncertainty: ~1 m (10x reduction)
- Filter converges within ~50 epochs, then tracks truth cleanly
 
## Folders
- 1st_KF_Theory/       Theory notes
- 2nd_F_Q_Derivation/  F and Q derivation notes
- 3rd_KF_Class/         Class 
- 4th_KF_Class_Test/   Class + test + plot

## Next: Week 11 — integrate KF with real BILL00USA trajectory
