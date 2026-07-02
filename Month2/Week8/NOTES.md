\# Week 8 - Atmospheric Corrections

&#x20;

\## Goal

Add ionospheric (Klobuchar) and tropospheric (Saastamoinen) delay

corrections to the SPP solver.



\## What was built

\- Klobuchar iono model (IS-GPS-200)

\- Saastamoinen tropo model (standard atmosphere)

\- Full nav header parser for GPSA/GPSB Klobuchar params

\- Relativistic satellite clock correction

\- TGD (Timing Group Delay) correction

\- Two-pass integration in real-data solver

\- Truth comparison against IGS reference coords

&#x20;

\## Results (BILL00USA, Jan 7 2024, 12:00 UTC, 8 satellites)

\- PASS 1 (no atmo):   sigma0 = 38 m, 3D error = 53 m

\- PASS 2 (with atmo): sigma0 = 38 m, 3D error = 63 m

\- NOTE: single-epoch SPP is inherently noisy at this level



\## Why sigma0 didn't drop dramatically

Single-frequency SPP has fundamental limits with broadcast

products alone. 30-60m single-epoch error is expected.

Kalman filtering (Week 10+) reduces this to \~5-15m.

&#x20;

\## Folders

\- 1st\_Atmospheric\_Theory/  Theory notes

\- 2nd\_Klobuchar\_Iono/       Standalone iono test

\- 3rd\_Saastamoinen/          Standalone tropo test

\- 4th\_Integration/            Full solver with both corrections

&#x20;

\## Next: Week 9 = FULL SPP PROCESSOR (Month 2 MILESTONE)



