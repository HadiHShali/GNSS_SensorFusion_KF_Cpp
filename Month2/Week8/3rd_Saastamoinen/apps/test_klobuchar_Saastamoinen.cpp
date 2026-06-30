#include <iostream>
#include <iomanip>
#include "satellites.h"

using namespace std;
 
int main() {
    cout << "=== KLOBUCHAR TEST ===" << endl;
    
    // Parse the parameters from your nav file
    string nav_file = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    KlobucharParams kp = parseKlobucharFromNav(nav_file);
    
    if (!kp.valid) {
        cerr << "Failed to parse Klobuchar params!" << endl;
        return 1;
    }
    
    cout << "Alpha: ";
    for (int i = 0; i < 4; i++) cout << kp.alpha[i] << " ";
    cout << endl;
    cout << "Beta:  ";
    for (int i = 0; i < 4; i++) cout << kp.beta[i] << " ";
    cout << endl << endl;
    
    // Test Ionospheric correcton: cases: BILL00USA location, various satellites
    double rec_lat = 33.578;     // BILL00USA
    double rec_lon = -117.064;
    double t_gps   = 43200.0;     // noon UTC
    	
	cout << endl << "===Klobuchar TEST===" << endl;
    cout << fixed << setprecision(2);
    cout << "Iono at zenith (el=90):       "
         << klobucharIonoDelay(kp, rec_lat, rec_lon, 0, 90, t_gps) << " m" << endl;
	cout << "Iono at zenith (el=45):       "
         << klobucharIonoDelay(kp, rec_lat, rec_lon, 0, 45, t_gps) << " m" << endl;
	cout << "Iono at el=30, az=0:           "
         << klobucharIonoDelay(kp, rec_lat, rec_lon, 0, 30, t_gps) << " m" << endl;
    cout << "Iono at el=15, az=0:           "
         << klobucharIonoDelay(kp, rec_lat, rec_lon, 0, 15, t_gps) << " m" << endl;
    cout << "Iono at el=10 (near horizon): "
         << klobucharIonoDelay(kp, rec_lat, rec_lon, 0, 10, t_gps) << " m" << endl;
    
	// Test Tropospheric correcton:
	cout << endl << "===SAASTAMOINEN TEST===" << endl;
	
	double rec_height = 472.0; // meter
	
	cout << "Tropo at el=90: " << saastamoinenTropoDelay (rec_lat, rec_height, 90) << " m" << endl;
	cout << "Tropo at el=45: " << saastamoinenTropoDelay (rec_lat, rec_height, 45) << " m" << endl;
	cout << "Tropo at el=30: " << saastamoinenTropoDelay (rec_lat, rec_height, 30) << " m" << endl;
	cout << "Tropo at el=15: " << saastamoinenTropoDelay (rec_lat, rec_height, 15) << " m" << endl;
	cout << "Tropo at el=10: " << saastamoinenTropoDelay (rec_lat, rec_height, 10) << " m" << endl;
	
    return 0;
}
