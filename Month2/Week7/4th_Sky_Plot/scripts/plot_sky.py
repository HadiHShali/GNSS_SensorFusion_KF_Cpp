# plot_sky.py - sky plot for GPS positioning fix
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Load the data
df = pd.read_csv('../data/sky_plot_data.csv')
print(f'loaded {len(df)} satellites')

# setup polar plot
fig, ax = plt.subplots(figsize = (10, 10), facecolor ='white', subplot_kw = dict(projection='polar') )
ax.set_theta_zero_location('N')  # 0 deg at top (North)
ax.set_theta_direction(-1)  # clockwise direction: compass convention
ax.set_rlim(0, 90) # radius limits: 0 to 90 degrees
ax.set_rticks([30, 60, 90]) # radial ticks
ax.set_yticklabels(['60 deg','30 deg','0 deg (Horizon)'], color='#666') # radial tick labels
ax.set_xticks(np.linspace(0, 2*np.pi, 8, endpoint=False))
ax.set_xticklabels(['N','NE','E','SE','S','SW','W','NW'],
                   color='#1B3A6B', fontsize=11, fontweight='bold')

# Elevation mas zone (yellow ring)
ELEV_MASK = 10.0
theta = np.linspace(0, 2*np.pi, 100)
ax.fill_between(theta, 90 - ELEV_MASK, 90,
                color='#FFE699', alpha=0.5, zorder=1)


# Plot the satellites
for _, row in df.iterrows():
    az = np.radians(row['azimuth_deg'])
    r = 90 - row['elevation_deg']  # convert elevation to polar radius

    if row['used']==1:
        ax.scatter(az, r, s=300, c='#1E7145', edgecolors='black', linewidth=2, zorder=5)
        ax.annotate(f"G{int(row['prn']):02d}", (az, r), fontsize=10, color='white', fontweight='bold', ha='center', va='center', zorder=6)

    else:
        ax.scatter(az, r, s=350, c='#C55A11', edgecolors='black', linewidth=2, marker='X', zorder=5)
        ax.annotate(f"G{int(row['prn']):02d}", (az, r-7), fontsize=10, color='#C55A11', fontweight='bold',ha='center', va='center')

# ── TITLE ──
n_used = int(df['used'].sum())
n_total = len(df)
ax.set_title(f'GPS Sky View - BILL00USA Station\n'
    f'Jan 7, 2024  12:00 UTC  |  {n_used} used / {n_total - n_used} rejected',
    fontsize=13, fontweight='bold', color='#1B3A6B', pad=20)

# ── LEGEND ──
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#1E7145',
           markersize=12, label='Used in fix', markeredgecolor='black'),
    Line2D([0],[0], marker='X', color='w', markerfacecolor='#C55A11',
           markersize=14, label=f'Rejected (below {ELEV_MASK} deg)',
           markeredgecolor='black'),
]
ax.legend(handles=legend_elements, loc='upper left',
          bbox_to_anchor=(-0.1, 1.05), fontsize=10, framealpha=0.95)
 
# ── SAVE ──
plt.tight_layout()
plt.savefig('sky_plot_real_data.png', dpi=200, bbox_inches='tight',
            facecolor='white')
plt.savefig('sky_plot_real_data.pdf', bbox_inches='tight',
            facecolor='white')
plt.close()
print('Saved sky_plot_real_data.png and .pdf')
