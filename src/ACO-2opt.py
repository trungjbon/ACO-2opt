from Local_Optimize import *
import time
from parameter import *
import lib
import random
from typing import Tuple, List
import matplotlib.pyplot as plt

class ACO_2opt:
    def __init__(self, ants: int, evaporation_rate: float, intensification: float,
                 alpha: float, beta: float, beta_evaporation_rate: float) -> None:
        # Parameters
        self.ants = ants
        self.evaporation_rate = evaporation_rate
        self.Q = intensification
        self.pheromone_alpha = alpha
        self.heuristic_beta = beta
        self.beta_evaporation_rate = beta_evaporation_rate

        # Matrices
        self.distance_matrix = None
        self.pheromone_matrix = None
        self.heuristic_matrix = None
        self.probability_matrix = None

        self.spatial_map = None
        self.num_nodes = None
        self.list_of_nodes = None

        # Internal statistics
        self.global_best = None
        self.best_series = None
        self.best_path = None
        self.best = None

        # Optimizer's status
        self.fitted = False
        self.fit_time = None

        # Early stopping 
        self.converged = False
        self.stopped_at_iteration = None


    def __str__(self) -> str:
        string = "Ant Colony Optimizer"
        string += "\n------------------------"
        string += "\nDesigned to solve traveling salesman problem. Optimizes either the minimum or maximum distance traveled."
        string += "\n------------------------"
        string += f"\nNumber of ants: {self.ants}"
        string += f"\nEvaporation rate: {self.evaporation_rate}"
        string += f"\nIntensification factor: {self.Q}"
        string += f"\nAlpha Heuristic: {self.pheromone_alpha}"
        string += f"\nBeta Heuristic: {self.heuristic_beta}"
        string += f"\nBeta Evaporation Rate: {self.beta_evaporation_rate}"
        string += "\n------------------------"

        if (self.fitted):
            string += "\nThis optimizer has been fitted."
        else:
            string += "\nThis optimizer has NOT been fitted."
        return string
    
    
    def __initialize(self, spatial_map: List[List[float]]) -> None:
        for point in spatial_map:
            assert len(point) == 2

        self.spatial_map = spatial_map
        self.num_nodes = len(spatial_map)
        self.list_of_nodes = list(range(0, self.num_nodes, 1))

        self.distance_matrix = lib.round_matrix(lib.pairwise_distances(spatial_map, spatial_map), 0)

        self.pheromone_matrix = [[1 for _ in range(self.num_nodes)] for _ in range(self.num_nodes)]

        self.heuristic_matrix = [[(0.0 if x == 0 else 1.0 / x) for x in row] for row in self.distance_matrix]

        lib.fill_diagonal(self.pheromone_matrix, val=0)
        lib.fill_diagonal(self.heuristic_matrix, val=0)

        self.probability_matrix = [[((self.pheromone_matrix[i][j] ** self.pheromone_alpha) * (self.heuristic_matrix[i][j] ** self.heuristic_beta)) 
            for j in range(len(self.heuristic_matrix[i]))] for i in range(len(self.heuristic_matrix))]

        
    def __reset_list_of_nodes(self) -> None:
        self.list_of_nodes = list(range(0, self.num_nodes, 1))


    def __update_probabilities(self) -> None:
        pheromone = self.pheromone_matrix
        heuristic = self.heuristic_matrix
        alpha = self.pheromone_alpha
        beta = self.heuristic_beta
        
        self.probability_matrix = [[((pheromone[i][j] ** alpha) * (heuristic[i][j] ** beta)) 
            for j in range(len(heuristic[i]))] for i in range(len(heuristic))]


    def __travel_to_next_node_from(self, current_node: int) -> int:
        prob_row = self.probability_matrix[current_node]

        numerator = [prob_row[i] for i in self.list_of_nodes]
        
        sum_allowed_k = sum(numerator)
        probabilities = [x / sum_allowed_k for x in numerator]

        next_node = random.choices(self.list_of_nodes, probabilities)[0]

        return next_node


    def __remove_node(self, node: int) -> None:
        self.list_of_nodes.remove(node)


    def __evaluate(self, iteration: int, paths: List) -> Tuple[List[int], float]:
        paths = paths[iteration]
        distance = self.distance_matrix

        best_score = float('inf')
        best_path = []
        
        for i, path in enumerate(paths):
            score = sum(distance[path[idx]][path[idx + 1]] for idx in range(len(path) - 1))

            if (score < best_score):
                best_score = score
                best_path = paths[i]

        return best_path, best_score
    
  
    def __evaporate(self) -> None:
        evaporation = 1 - self.evaporation_rate

        for i in range (len(self.pheromone_matrix)):
            for j in range (len(self.pheromone_matrix[i])):
                self.pheromone_matrix[i][j] *= evaporation
                
        
        self.heuristic_beta *= (1 - self.beta_evaporation_rate)


    def __update_pheromone(self, iteration : int, paths : List[List[int]]) -> None:
        paths_iter = paths[iteration]
        distance = self.distance_matrix

        for path_kth in paths_iter:
            L_k = sum(distance[path_kth[i]][path_kth[i + 1]] for i in range(len(path_kth) - 1))

            delta = self.Q / L_k
            for i in range(len(path_kth) - 1):
                self.pheromone_matrix[path_kth[i]][path_kth[i + 1]] += delta # (∆τkij(t))
                self.pheromone_matrix[path_kth[i + 1]][path_kth[i]] += delta


    def __intensify(self, best_path: List[int], best: float) -> None:
        delta = self.Q / best
        for i in range(len(best_path) - 1):
            self.pheromone_matrix[best_path[i]][best_path[i + 1]] += delta
            self.pheromone_matrix[best_path[i + 1]][best_path[i]] += delta


    def fit(self, spatial_map: List[List[float]], iterations: int, conv_crit=20, verbose=True, decimal=2) -> None:
        start = time.time()
        self.__initialize(spatial_map=spatial_map)

        num_equal = 0

        if (verbose):
            print(f"{self.num_nodes} nodes were given. Beginning ACO Optimization with {iterations} iterations...\n")

        self.best_series = [0.0] * iterations
        self.global_best = [0.0] * iterations

        best_score_so_far = float('inf')
        paths = []

        for iteration in range(iterations):
            start_iter = time.time()
            paths.append([])

            for ant in range(self.ants):
                paths[iteration].append([])
                path = [None] * (self.num_nodes + 1)

                current_node = self.list_of_nodes[random.randint(0, self.num_nodes - 1)]
                start_node = current_node

                for node in range(self.num_nodes):
                    path[node] = current_node
                    self.__remove_node(node=current_node)
                    if (len(self.list_of_nodes) != 0):
                        current_node = self.__travel_to_next_node_from(current_node=current_node) 


                path[-1] = start_node # go back to start 
                self.__reset_list_of_nodes() 
                paths[iteration][ant] = path


            best_path, best_score = self.__evaluate(iteration=iteration, paths=paths)
            # best_path, best_score, isBetter = local_optimize(best_hc=best_path, best_score=best_score, distance_matrix=self.distance_matrix)
            best_path, best_score, isBetter = two_opt_plus(best_hc=best_path, best_score=best_score, distance_matrix=self.distance_matrix)

            if (isBetter):
                self.__intensify(best_path, best_score)

            if (best_score < best_score_so_far):
                best_score_so_far = best_score
                self.best_path = best_path


            if (best_score == best_score_so_far):
                num_equal += 1
            else:
                num_equal = 0


            self.best_series[iteration] = best_score
            self.global_best[iteration] = best_score_so_far
            self.__evaporate()
            self.__update_pheromone(iteration=iteration, paths=paths)
            self.__update_probabilities()

            if (verbose):
                print(
                    f"Iteration {iteration}/{iterations} | Score: {round(best_score, decimal)} | Best so far: {round(best_score_so_far, decimal)} | {round(time.time() - start_iter, decimal)} s"
                    + f" | Evaporation: {round(self.evaporation_rate, decimal)} | Alpha: {round(self.pheromone_alpha, decimal)} | Beta: {round(self.heuristic_beta, decimal)}"
                )

            if (num_equal == conv_crit):
                self.converged = True
                self.stopped_at_iteration = iteration
                if (verbose): 
                    print("\nConvergence criterion has been met. Stopping....")
                break

        if (not self.converged): 
            self.stopped_at_iteration = iterations

        self.fit_time = round(time.time() - start)
        self.fitted = True
        self.best = self.global_best[-1]

        if (verbose): 
            print(
                f"\nACO fitted. Runtime: {self.fit_time} second(s). Best score: {round(self.best, decimal)}"
        )
        


    def get_result(self) -> Tuple[List[int], float, float, List[float], List[float]]:
        return self.best_path, self.best, self.fit_time, self.best_series, self.global_best
    
    def plot_best_path(self):
        coordinates = [self.spatial_map[i] for i in self.best_path]
        x = [p[0] for p in coordinates]
        y = [p[1] for p in coordinates]
        plt.plot(x, y, 'o-')
        plt.plot(x[0], y[0], 'ro', label='Start/End')
        plt.title(f"Best Path (Score: {self.best})")
        plt.show()


if __name__ == "__main__":
    random.seed(42)
    towns = TOWNS

    ACO_optimizer = ACO_2opt(ants=ANTS, evaporation_rate=EVAPORATION_RATE, intensification=INTENSIFICATION, alpha=ALPHA, beta=BETA, beta_evaporation_rate=BETA_EVAPORATION_RATE)
    ACO_optimizer.fit(towns, iterations=MAX_ITER_ACO)
    best_path, best , fit_time , best_series, global_best = ACO_optimizer.get_result()
    best_path_coor = [towns[i] for i in best_path]
    print("Best path : ", best_path_coor)

    ACO_optimizer.plot_best_path()
    