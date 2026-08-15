# Week 15 -- Carrier Phase Observable Theory
 
## Goal
Establish carrier phase observable foundations: parsing, theory,
and cycle-slip detection -- prerequisite for RTK/PPP (Weeks 16-19).
 
## What was built
- Carrier phase equation theory (Day 1) -- iono sign-flip, integer
  ambiguity N, wavelength-driven precision (~1000x better than code)
- RINEX parser extended for L1C/L2W/LLI (Day 2)
- Geometry-free (LG) cycle-slip detector, unit bug found+fixed
- Melbourne-Wubbena (MW) independent detector added (Day 3)
- Full 30-satellite scan with proper gap-vs-slip distinction
 
## Key results
- G06 (detailed case): 4 tracking arcs, 2 genuine slips within-arc,
  cross-validated by LG, MW, AND RINEX LLI flags simultaneously
- Constellation-wide: LG/MW agree almost exactly on all 30 satellites
- L2 slips 2-3x more often than L1 across the ENTIRE constellation
  (semi-codeless P(Y) tracking is weaker) -- confirmed, not anecdotal
- G20 flagged as an outlier (12 LG / 12 MW / 12 LLI_L1 / 15 LLI_L2)
  -- worth an elevation cross-check in a future session
 
## Bugs found + fixed (both via numerical hand-verification)
1. LG formula: raw cycle subtraction -> must convert to meters first
2. Gap handling: silently skipped short dropouts where real slips
   live -- fixed with short-gap (<=4 epoch) vs new-arc (>4) split
 
## Next: Week 16 -- Differencing techniques (single/double-difference)
Requires sourcing a nearby CORS/IGS base station for baseline data.

