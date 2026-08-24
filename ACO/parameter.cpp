#include <bits/stdc++.h>

using namespace std;

vector<vector<double>> read_map(string file_path){
    vector<vector<double>> cities;

    ifstream file(file_path);

    string line;

    // Tim NODE_COORD_SECTION
    while (getline(file, line)){
        if (line == "NODE_COORD_SECTION"){
            break;
        }
    }

    // Doc cac dong toa do
    while (getline(file, line)){
        if (line == "EOF"){
            break;
        }

        stringstream ss(line);

        int id;
        double x;
        double y;

        ss >> id >> x >> y;

        cities.push_back({x, y});
    }

    return cities;
}


string PATH_TO_MAP = "Benchmark/kroA100.tsp";

vector<vector<double>> TOWNS = read_map(PATH_TO_MAP);

int n = TOWNS.size();

int ANTS = n / 4;

double EVAPORATION_RATE = 0.1;

double INTENSIFICATION = 5000;

int MAX_ITER_ACO = 200;

double BETA_EVAPORATION_RATE = 0.001;

double ALPHA = 1.0;

double BETA = 5.0;
