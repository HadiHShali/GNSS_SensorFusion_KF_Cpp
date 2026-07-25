# Week 11 — Real Data Integration
 
## Goal
Integrate the Week 10 Kalman filter with real BILL00USA GPS data
from the Week 9 SPP processor.
 
## What was built
- TrajectoryReader: parses trajectory.csv (5,760 epochs)
- kf_spp_processor: full KF pipeline on real ECEF measurements
- plot_kf_comparison.py: 4-panel raw vs filtered analysis
 
## Debugging log
- Segfault traced to R matrix near-singularity when sigma0 was
  very small; fixed with a measurement noise floor (15m)
- Initial KF was overconfident (1-sigma = 0.57m vs 36.7m actual error)
 
## Tuning experiments
  Run 1: q=0.01, R=raw sigma0        -> 36.7m error, 0.57m sigma
  Run 2: q=0.05, R=2.5x floor15      -> ~34m error,  ~2-4m sigma
  Run 3: q=0.01, R=2.5x floor15      -> 33.5m error, 1.44m sigma  [LOCKED]
 
## Conclusion
Mean error improved 28% (39.0m -> 28.0m). Remaining error is
SYSTEMATIC (ionospheric residual + multipath), not random noise --
confirmed by persistent confidence-vs-error mismatch across all
three tuning configurations. No amount of Q/R tuning removes it.
 
## Next: Month 4 — dual-frequency L1/L2 ionospheric elimination
