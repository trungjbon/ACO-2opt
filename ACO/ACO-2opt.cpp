#include <bits/stdc++.h>

using namespace std;

vector<vector<double>> pairwise_distances(vector<vector<double>> points1, vector<vector<double>> points2);
vector<vector<double>> round_matrix(vector<vector<double>> matrix, int decimal);
void fill_diagonal(vector<vector<double>>& matrix, double value);
void two_opt_plus(
    vector<int> best_hc, double best_score, vector<vector<double>> distance_matrix, 
    vector<int>& result_hc, double& result_score, bool& isBetter);

extern vector<vector<double>> TOWNS;
extern int ANTS;
extern double EVAPORATION_RATE;
extern double INTENSIFICATION;
extern int MAX_ITER_ACO;
extern double BETA_EVAPORATION_RATE;
extern double ALPHA;
extern double BETA;


class ACO_2opt{
private:
    int ants;
    double evaporation_rate;
    double Q;
    double pheromone_alpha;
    double heuristic_beta;
    double beta_evaporation_rate;

    vector<vector<double>> distance_matrix;
    vector<vector<double>> pheromone_matrix;
    vector<vector<double>> heuristic_matrix;
    vector<vector<double>> probability_matrix;

    vector<vector<double>> spatial_map;
    int num_nodes;
    vector<int> list_of_nodes;

    vector<double> global_best;
    vector<double> best_series;
    vector<int> best_path;
    double best;

    bool fitted;
    double fit_time;

    bool converged;
    int stopped_at_iteration;

public:
    ACO_2opt(int ants, double evaporation_rate, double intensification, double alpha, double beta, double beta_evaporation_rate){
        this->ants = ants;
        this->evaporation_rate = evaporation_rate;
        this->Q = intensification;
        this->pheromone_alpha = alpha;
        this->heuristic_beta = beta;
        this->beta_evaporation_rate = beta_evaporation_rate;

        fitted = false;
        fit_time = 0;
        converged = false;
        stopped_at_iteration = 0;
    }


    void initialize(vector<vector<double>> spatial_map){
        this->spatial_map = spatial_map;
        num_nodes = spatial_map.size();

        list_of_nodes.clear();
        for (int i = 0; i < num_nodes; i++) list_of_nodes.push_back(i);

        distance_matrix = round_matrix(pairwise_distances(spatial_map, spatial_map), 0);

        pheromone_matrix = vector<vector<double>>(num_nodes, vector<double>(num_nodes, 1));

        heuristic_matrix = vector<vector<double>>(num_nodes, vector<double>(num_nodes));

        for (int i = 0; i < num_nodes; i++){
            for (int j = 0; j < num_nodes; j++){
                if (distance_matrix[i][j] == 0) heuristic_matrix[i][j] = 0;
                else heuristic_matrix[i][j] = 1.0 / distance_matrix[i][j];
            }
        }

        fill_diagonal(pheromone_matrix, 0);
        fill_diagonal(heuristic_matrix, 0);

        update_probabilities();
    }


    void reset_list_of_nodes(){
        list_of_nodes.clear();
        for (int i = 0; i < num_nodes; i++) list_of_nodes.push_back(i);
    }


    void update_probabilities(){
        probability_matrix = vector<vector<double>>(num_nodes, vector<double>(num_nodes));

        for (int i = 0; i < num_nodes; i++){
            for (int j = 0; j < num_nodes; j++){
                probability_matrix[i][j] = pow(pheromone_matrix[i][j], pheromone_alpha) * pow(heuristic_matrix[i][j], heuristic_beta);
            }
        }
    }


    int travel_to_next_node_from(int current_node){
        vector<double> numerator;

        for (int i = 0; i < list_of_nodes.size(); i++){
            numerator.push_back(probability_matrix[current_node][list_of_nodes[i]]);
        }

        double sum_allowed_k = 0;
        for (int i = 0; i < numerator.size(); i++) sum_allowed_k += numerator[i];

        vector<double> probabilities;

        for (int i = 0; i < numerator.size(); i++) probabilities.push_back(numerator[i] / sum_allowed_k);

        double random_number = (double)rand() / RAND_MAX;
        double sum = 0;

        for (int i = 0; i < probabilities.size(); i++){
            sum += probabilities[i];

            if (random_number <= sum) return list_of_nodes[i];
        }

        return list_of_nodes[list_of_nodes.size() - 1];
    }


    void remove_node(int node){
        for (int i = 0; i < list_of_nodes.size(); i++){
            if (list_of_nodes[i] == node){
                list_of_nodes.erase(list_of_nodes.begin() + i);
                return;
            }
        }
    }


    void evaluate(vector<vector<vector<int>>> paths, int iteration, vector<int>& result_path, double& result_score){
        vector<vector<int>> paths_iter = paths[iteration];

        double best_score = 999999999999.0;
        vector<int> best_path;

        for (int i = 0; i < paths_iter.size(); i++){
            vector<int> path = paths_iter[i];
            double score = 0;

            for (int j = 0; j < path.size() - 1; j++) score += distance_matrix[path[j]][path[j + 1]];

            if (score < best_score){
                best_score = score;
                best_path = path;
            }
        }

        result_path = best_path;
        result_score = best_score;
    }


    void evaporate(){
        double evaporation = 1 - evaporation_rate;

        for (int i = 0; i < pheromone_matrix.size(); i++){
            for (int j = 0; j < pheromone_matrix[i].size(); j++){
                pheromone_matrix[i][j] *= evaporation;
            }
        }

        heuristic_beta *= (1 - beta_evaporation_rate);
    }


    void update_pheromone(int iteration, vector<vector<vector<int>>> paths){
        vector<vector<int>> paths_iter = paths[iteration];

        for (int k = 0; k < paths_iter.size(); k++){
            vector<int> path = paths_iter[k];
            double L_k = 0;

            for (int i = 0; i < path.size() - 1; i++) L_k += distance_matrix[path[i]][path[i + 1]];

            double delta = Q / L_k;

            for (int i = 0; i < path.size() - 1; i++){
                pheromone_matrix[path[i]][path[i + 1]] += delta;
                pheromone_matrix[path[i + 1]][path[i]] += delta;
            }
        }
    }


    void intensify(vector<int> best_path, double best_score){
        double delta = Q / best_score;

        for (int i = 0; i < best_path.size() - 1; i++){
            pheromone_matrix[best_path[i]][best_path[i + 1]] += delta;
            pheromone_matrix[best_path[i + 1]][best_path[i]] += delta;
        }
    }


    void fit(vector<vector<double>> spatial_map, int iterations, int conv_crit = 20, bool verbose = true, int decimal = 2){
        clock_t start = clock();

        initialize(spatial_map);

        int num_equal = 0;

        if (verbose) cout << num_nodes << " nodes were given. Beginning ACO Optimization with " << iterations << " iterations..." << endl << endl;

        best_series = vector<double>(iterations, 0);
        global_best = vector<double>(iterations, 0);

        double best_score_so_far = 999999999999.0;

        vector<vector<vector<int>>> paths;

        for (int iteration = 0; iteration < iterations; iteration++){
            clock_t start_iter = clock();

            paths.push_back(vector<vector<int>>());

            for (int ant = 0; ant < ants; ant++){
                vector<int> path(num_nodes + 1);

                int current_node = list_of_nodes[rand() % num_nodes];
                int start_node = current_node;

                for (int node = 0; node < num_nodes; node++){
                    path[node] = current_node;
                    remove_node(current_node);

                    if (list_of_nodes.size() != 0) current_node = travel_to_next_node_from(current_node);
                }

                path[num_nodes] = start_node;

                reset_list_of_nodes();

                paths[iteration].push_back(path);
            }

            vector<int> best_path_now;
            double best_score;

            evaluate(paths, iteration, best_path_now, best_score);

            vector<int> optimized_path;
            double optimized_score;
            bool isBetter;

            two_opt_plus(best_path_now, best_score, distance_matrix, optimized_path, optimized_score, isBetter);

            best_path_now = optimized_path;
            best_score = optimized_score;

            if (isBetter) intensify(best_path_now, best_score);

            if (best_score < best_score_so_far){
                best_score_so_far = best_score;
                best_path = best_path_now;
            }

            if (best_score == best_score_so_far){
                num_equal++;
            } else {
                num_equal = 0;
            }


            best_series[iteration] = best_score;
            global_best[iteration] = best_score_so_far;

            evaporate();
            update_pheromone(iteration, paths);
            update_probabilities();

            if (verbose){
                double time_iter = (double)(clock() - start_iter) / CLOCKS_PER_SEC;

                cout << "Iteration " << iteration << "/" << iterations
                     << " | Score: " << round(best_score * 100) / 100
                     << " | Best so far: " << round(best_score_so_far * 100) / 100
                     << " | " << round(time_iter * 100) / 100 << " s"
                     << " | Evaporation: " << evaporation_rate
                     << " | Alpha: " << pheromone_alpha
                     << " | Beta: " << heuristic_beta << endl;
            }

            if (num_equal == conv_crit){
                converged = true;
                stopped_at_iteration = iteration;

                if (verbose) cout << endl << "Convergence criterion has been met. Stopping...." << endl;

                break;
            }
        }

        if (!converged) stopped_at_iteration = iterations;

        fit_time = round((double)(clock() - start) / CLOCKS_PER_SEC);
        fitted = true;
        best = best_score_so_far;

        if (verbose) cout << endl << "ACO fitted. Runtime: " << fit_time << " second(s). Best score: " << round(best * 100) / 100 << endl;
    }


    void get_result(
        vector<int>& result_path, double& result_best, double& result_fit_time,
        vector<double>& result_best_series, vector<double>& result_global_best)
    {
        result_path = best_path;
        result_best = best;
        result_fit_time = fit_time;
        result_best_series = best_series;
        result_global_best = global_best;
    }


    void plot_best_path(){
        cout << endl << "Best Path:" << endl;

        for (int i = 0; i < best_path.size(); i++){
            int node = best_path[i];

            cout << "(" << spatial_map[node][0] << ", " << spatial_map[node][1] << ")";

            if (i < best_path.size() - 1) cout << " -> ";
        }

        cout << endl;
    }
};


int main(){
    srand(42);

    vector<vector<double>> towns = TOWNS;

    ACO_2opt ACO_optimizer(ANTS, EVAPORATION_RATE, INTENSIFICATION, ALPHA, BETA, BETA_EVAPORATION_RATE);

    ACO_optimizer.fit(towns, MAX_ITER_ACO);

    vector<int> best_path;
    double best;
    double fit_time;
    vector<double> best_series;
    vector<double> global_best;

    ACO_optimizer.get_result(best_path, best, fit_time, best_series, global_best);

    ACO_optimizer.plot_best_path();

    return 0;
}
