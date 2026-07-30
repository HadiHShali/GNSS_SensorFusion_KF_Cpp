import pandas as pd, numpy as np
 
TRUE_LAT, TRUE_LON, TRUE_H = 33.578381, -117.064265, 472.0
R_EARTH = 6378137.0
lat_rad = np.radians(TRUE_LAT)
 
def err3d(lat, lon, h):
    e = np.radians(lon - TRUE_LON) * R_EARTH * np.cos(lat_rad)
    n = np.radians(lat - TRUE_LAT) * R_EARTH
    dz = h - TRUE_H
    return np.sqrt(e**2 + n**2 + dz**2)
 
df = pd.read_csv('../data/trajectory_rts.csv')   # or whatever your KF/RTS output is named
kf_err  = err3d(df['kf_lat'], df['kf_lon'], df['kf_h'])
rts_err = err3d(df['rts_lat'], df['rts_lon'], df['rts_h'])
 
print(f"Dual-freq raw (Day 2):         mean=48.40m  RMS=52.39m")
print(f"Dual-freq + KF (today):        mean={kf_err.mean():.2f}m  RMS={np.sqrt((kf_err**2).mean()):.2f}m")
print(f"Dual-freq + KF + RTS (today):  mean={rts_err.mean():.2f}m  RMS={np.sqrt((rts_err**2).mean()):.2f}m")
print(f"[Reference] Single-freq + RTS (Week 12): mean=15.56m  RMS=15.58m")

