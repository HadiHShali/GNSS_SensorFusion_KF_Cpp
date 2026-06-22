#pragma once
#include <vector>
#include <string>
#include "gps_ls_solver.h"

std::vector<Observation> loadObservations(const std::string& filename);

void ecefToLla(double X, double Y, double Z, double& lat, double& lon, double& h);