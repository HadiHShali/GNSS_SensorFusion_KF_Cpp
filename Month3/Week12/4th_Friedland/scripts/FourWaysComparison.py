import pandas as pd, numpy as np

df = pd.read_csv('../data/trajectory_rts_Bias.csv')  # adjust path as needed
TRUE_LAT, TRUE_LON, TRUE_H = 33.578381, -117.064265, 472.0
R_EARTH = 6378137.0
lat_rad = np.radians(TRUE_LAT)

def err3d(lat, lon, h):
    e = np.radians(lon - TRUE_LON) * R_EARTH * np.cos(lat_rad)
    n = np.radians(lat - TRUE_LAT) * R_EARTH
    dz = h - TRUE_H
    return np.sqrt(e**2 + n**2 + dz**2)

variants = {
    'A: KF only':              ('kf_lat','kf_lon','kf_h'),
    'B: KF + Friedland':       ('fc_lat','fc_lon','fc_h'),
    'C: RTS only':              ('rts_lat','rts_lon','rts_h'),
    'D: RTS + Friedland':       ('corrected_lat','corrected_lon','corrected_h'),
}

for name, (lat_c, lon_c, h_c) in variants.items():
    e = err3d(df[lat_c], df[lon_c], df[h_c])
    print(f"{name:24s}  mean={e.mean():6.2f}m   RMS={np.sqrt((e**2).mean()):6.2f}m")