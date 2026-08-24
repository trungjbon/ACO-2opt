#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;


// --------------------------------------------------
// Dao nguoc doan tu i den k
// --------------------------------------------------
vector<int> swap_algorithm(vector<int> route, int i, int k){
    vector<int> result = route;

    while (i < k){
        int temp = result[i];
        result[i] = result[k];
        result[k] = temp;

        i++;
        k--;
    }

    return result;
}


// --------------------------------------------------
// Tinh tong do dai cua mot route
// --------------------------------------------------
double cal_score_hc(vector<int> route, vector<vector<double>> distance_matrix){
    double distance = 0;

    int n = route.size();

    for (int i = 0; i < n - 1; i++){
        distance = distance + distance_matrix[route[i]][route[i + 1]];
    }

    // Quay lai thanh pho dau tien
    distance = distance + distance_matrix[route[n - 1]][route[0]];

    return distance;
}


// --------------------------------------------------
// 2-opt plus
// --------------------------------------------------
void two_opt_plus(
    vector<int> best_hc,
    double best_score,
    vector<vector<double>> distance_matrix,
    vector<int>& result_hc,
    double& result_score,
    bool& is_better)
{
    int n = best_hc.size();

    bool improved = true;

    is_better = false;

    while (improved){
        improved = false;

        double best_delta = 0;

        int best_i = -1;
        int best_k = -1;

        for (int i = 1; i < n - 2; i++){
            for (int k = i + 1; k < n - 1; k++){
                int p1 = best_hc[i - 1];
                int p2 = best_hc[i];

                int p3 = best_hc[k];
                int p4 = best_hc[k + 1];

                double delta = (distance_matrix[p1][p2] + distance_matrix[p3][p4])
                    - (distance_matrix[p1][p3] + distance_matrix[p2][p4]);

                if (delta > best_delta){
                    best_delta = delta;

                    best_i = i;
                    best_k = k;

                    improved = true;
                }
            }
        }

        if (improved){
            best_hc = swap_algorithm(best_hc, best_i, best_k);

            best_score -= best_delta;

            is_better = true;
        }
    }

    result_hc = best_hc;
    result_score = best_score;
}
