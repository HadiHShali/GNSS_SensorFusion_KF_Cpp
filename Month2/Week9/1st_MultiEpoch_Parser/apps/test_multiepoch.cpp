#include <iostream>
#include "rinex_obs.h"

using namespace std;

int main()
{
	string obs_file = "../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx";
	
	cout << "parsing all epochs..." << endl;
	vector<ObsEpoch> epochs = parseRinexAllEpochs(obs_file);
	
	cout << "Total Epochs Parsed: "<< epochs.size() << endl;
	
	if (!epochs.empty()) {
        cout << "First epoch: " << epochs[0].hour << ":" 
             << epochs[0].minute << " - " << epochs[0].gps_records.size()
             << " satellites" << endl;
        cout << "Last epoch:  " << epochs.back().hour << ":"
             << epochs.back().minute << " - " 
             << epochs.back().gps_records.size() << " satellites" << endl;
    }
    
    return 0;
}