"""
plot_final_sp3clk_comparison.py

Month 3 FINAL comparison — SP3/CLK precise products + KF + RTS
The capstone plot for the entire Month 3 arc (Weeks 9-14).

Produces a 4-panel figure:
    Top-left:     Horizontal scatter — KF trajectory vs RTS-smoothed vs Truth
    Top-right:    3D error vs time — KF vs RTS
    Bottom-left:  Error distribution histogram — KF vs RTS
    Bottom-right: Position uncertainty vs time — KF vs RTS (log scale)

Run from: Month3/Week14/4th_SP3_KF_RTS/
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────
# CONFIG
# ─────────────────────────────────────────────────────────────
CSV_PATH = "../data/trajectory_dualfreq_SP3CLK_RTS.csv"

TRUE_LAT, TRUE_LON, TRUE_H = 33.578381, -117.064265, 472.0
R_EARTH = 6378137.0
LAT_RAD = np.radians(TRUE_LAT)

KF_COLOR  = "#2E75B6"
RTS_COLOR = "#1E7145"
TRUTH_COLOR = "gold"

# Reference benchmark from Week 13 (dual-freq + KF + RTS, no precise products)
BENCHMARK_MEAN = 13.84
BENCHMARK_RMS  = 13.91


def err_components(lat, lon, h):
    """East / North / Height error components (m) vs BILL00USA truth."""
    east  = np.radians(lon - TRUE_LON) * R_EARTH * np.cos(LAT_RAD)
    north = np.radians(lat - TRUE_LAT) * R_EARTH
    dz    = h - TRUE_H
    return east, north, dz


def err3d(lat, lon, h):
    e, n, dz = err_components(lat, lon, h)
    return np.sqrt(e**2 + n**2 + dz**2)


# ─────────────────────────────────────────────────────────────
# LOAD + VERIFY COLUMNS + PROCESS
# ─────────────────────────────────────────────────────────────
df = pd.read_csv(CSV_PATH)

# Safety check: confirm the actual column names before assuming they match
# kf_lat/kf_lon/kf_h/rts_lat/rts_lon/rts_h -- catches silent mismatches if
# your C++ output used different column names than expected.
print("Columns found in CSV:", df.columns.tolist())
print()

required_cols = ['kf_lat', 'kf_lon', 'kf_h', 'rts_lat', 'rts_lon', 'rts_h']
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise SystemExit(
        f"ERROR: expected columns {missing} not found in CSV.\n"
        f"Update the column names below to match what's actually in your file."
    )

df['t_hours'] = df['t_gps'] / 3600.0

kf_east, kf_north, kf_dz   = err_components(df['kf_lat'],  df['kf_lon'],  df['kf_h'])
rts_east, rts_north, rts_dz = err_components(df['rts_lat'], df['rts_lon'], df['rts_h'])

df['kf_3d']  = np.sqrt(kf_east**2  + kf_north**2  + kf_dz**2)
df['rts_3d'] = np.sqrt(rts_east**2 + rts_north**2 + rts_dz**2)

kf_err  = df['kf_3d']
rts_err = df['rts_3d']
kf_mean, kf_rms   = kf_err.mean(),  np.sqrt((kf_err**2).mean())
rts_mean, rts_rms = rts_err.mean(), np.sqrt((rts_err**2).mean())

print(f"SP3/CLK + KF:        mean={kf_mean:.2f}m   RMS={kf_rms:.2f}m")
print(f"SP3/CLK + KF + RTS:  mean={rts_mean:.2f}m   RMS={rts_rms:.2f}m")
print(f"[Benchmark] Dual-freq + KF + RTS (Wk13): mean={BENCHMARK_MEAN}m  RMS={BENCHMARK_RMS}m")
beat = "BEATS" if rts_mean < BENCHMARK_MEAN else "does not beat"
print(f"-> Final result {beat} the Week 13 benchmark")
print()

# ─────────────────────────────────────────────────────────────
# FIGURE — 4 panels
# ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(15, 11), facecolor='white')

# ---------- TOP-LEFT: Horizontal scatter ----------
ax = axes[0, 0]
ax.plot(kf_east, kf_north, color=KF_COLOR, lw=0.8, alpha=0.6,
        label=f'KF (mean={kf_mean:.1f}m)')
ax.plot(rts_east, rts_north, color=RTS_COLOR, lw=1.8,
        label=f'RTS (mean={rts_mean:.1f}m)')
ax.scatter([0], [0], s=400, c=TRUTH_COLOR, marker='*',
           edgecolors='black', linewidths=1.8, zorder=10, label='Truth')
ax.axhline(0, color='#404040', lw=0.5)
ax.axvline(0, color='#404040', lw=0.5)
ax.set_xlabel('East error (m)', fontsize=11)
ax.set_ylabel('North error (m)', fontsize=11)
ax.set_title('Horizontal Scatter — SP3/CLK + KF vs RTS',
              fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=10, loc='best')
ax.grid(alpha=0.3)
ax.set_aspect('equal')

# ---------- TOP-RIGHT: 3D error vs time ----------
ax = axes[0, 1]
ax.plot(df['t_hours'], df['kf_3d'], color=KF_COLOR, lw=0.8, alpha=0.6, label='KF')
ax.plot(df['t_hours'], df['rts_3d'], color=RTS_COLOR, lw=1.8, label='RTS')
ax.axhline(BENCHMARK_MEAN, color='#C55A11', lw=1.5, linestyle='--', alpha=0.8,
           label=f'Wk13 benchmark ({BENCHMARK_MEAN}m)')
ax.set_xlabel('Time (hours)', fontsize=11)
ax.set_ylabel('3D error (m)', fontsize=11)
ax.set_title('3D Error vs Time', fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
ax.set_xlim(0, 24)

# ---------- BOTTOM-LEFT: Error histogram ----------
ax = axes[1, 0]
bins = np.linspace(0, max(df['kf_3d'].max(), df['rts_3d'].max()), 45)
ax.hist(df['kf_3d'], bins=bins, color=KF_COLOR, alpha=0.5, edgecolor='black',
        label=f'KF (RMS={kf_rms:.1f}m)')
ax.hist(df['rts_3d'], bins=bins, color=RTS_COLOR, alpha=0.6, edgecolor='black',
        label=f'RTS (RMS={rts_rms:.1f}m)')
ax.axvline(BENCHMARK_MEAN, color='#C55A11', lw=1.5, linestyle='--', alpha=0.8)
ax.set_xlabel('3D error (m)', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title('Error Distribution', fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

# ---------- BOTTOM-RIGHT: Uncertainty vs time ----------
ax = axes[1, 1]
if 'kf_sigma' in df.columns and 'rts_sigma' in df.columns:
    ax.plot(df['t_hours'], df['kf_sigma'], color=KF_COLOR, lw=1.3, label='KF 1-sigma')
    ax.plot(df['t_hours'], df['rts_sigma'], color=RTS_COLOR, lw=1.8, label='RTS 1-sigma')
    ax.set_yscale('log')
    ax.legend(fontsize=10)
else:
    ax.text(0.5, 0.5, 'sigma columns not found in CSV\n(check kf_sigma / rts_sigma column names)',
            ha='center', va='center', transform=ax.transAxes, fontsize=10, color='#999')
ax.set_xlabel('Time (hours)', fontsize=11)
ax.set_ylabel('Position 1-sigma (m)', fontsize=11)
ax.set_title('Position Uncertainty vs Time', fontsize=12.5, fontweight='bold', color='#1B3A6B')
ax.grid(alpha=0.3, which='both')
ax.set_xlim(0, 24)

fig.suptitle('Month 3 Finale — BILL00USA: Dual-Freq + SP3/CLK + KF/RTS (24 hours)',
              fontsize=15.5, fontweight='bold', color='#1B3A6B', y=1.01)
plt.tight_layout()
plt.savefig('final_sp3clk_comparison.png', dpi=200, bbox_inches='tight', facecolor='white')
print("\nSaved final_sp3clk_comparison.png")