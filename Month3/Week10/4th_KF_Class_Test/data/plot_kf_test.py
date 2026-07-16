import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
 
df = pd.read_csv('kf_static_test.csv')
 
fig, axes = plt.subplots(2, 1, figsize=(12, 8))
 
ax = axes[0]
ax.axhline(0, color='green', lw=2, label='Truth (0 m)')
ax.scatter(df['time'], df['meas_x'], s=8, color='orange',
           alpha=0.4, label='Noisy measurements')
ax.plot(df['time'], df['kf_x'], color='blue', lw=2, label='KF estimate')
ax.set_xlabel('Time (s)'); ax.set_ylabel('X (m)')
ax.legend(); ax.grid(alpha=0.3)
ax.set_title('KF Smooths Noise Toward Truth')
 
ax = axes[1]
ax.plot(df['time'], np.sqrt(df['P00']), color='purple', lw=2)
ax.set_xlabel('Time (s)'); ax.set_ylabel('sqrt(P00) — 1-sigma (m)')
ax.grid(alpha=0.3)
ax.set_title('Filter Confidence Grows Over Time')
 
plt.tight_layout()
plt.savefig('kf_static_test.png', dpi=200)
print('Saved kf_static_test.png')

