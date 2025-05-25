from typing import List, Tuple


def swap_algorithm(r, i, k):
    return r[0:i] + r[i:k + 1][::-1] + r[k + 1:len(r)]


def two_opt(best_hc: List[int], best_score: float, distance_matrix: List[List[float]]) -> Tuple[List[int], float, bool]:
    best_hc = list(best_hc[:-1])
    n = len(best_hc)
    dist_matrix = distance_matrix
    improved = False
    
    for swap1 in range(1, n - 2):
        for swap2 in range(swap1 + 1, n):
            new_solution = swap_algorithm(best_hc, swap1, swap2)
            new_distance = cal_score_hc(new_solution, dist_matrix)
            if (new_distance < best_score):
                best_hc = new_solution
                best_score = new_distance
                improved = True

    best_hc.append(best_hc[0])
    return best_hc, best_score, improved

def two_opt_plus(best_hc: List[int], best_score: float, distance_matrix: List[List[float]]) -> Tuple[List[int], float, bool]:
    n = len(best_hc)
    isBetter = False

    improved = True
    while (improved):
        improved = False
        best_delta = 0
        best_i, best_k = -1, -1

        for i in range(1, n - 2):
            for k in range(i + 1, n - 1):
                p1 = best_hc[i - 1]
                p2 = best_hc[i]
                p3 = best_hc[k]
                p4 = best_hc[k + 1] 

                delta = (distance_matrix[p1][p2] + distance_matrix[p3][p4]) - (distance_matrix[p1][p3] + distance_matrix[p2][p4])

                if (delta > best_delta):
                    best_delta = delta
                    best_i, best_k = i, k
                    improved = True

        if (improved):
            best_hc = swap_algorithm(best_hc, best_i, best_k)
            best_score -= best_delta
            isBetter = True        

    return best_hc, best_score, isBetter



def local_optimize(best_hc: List[int], best_score: float, distance_matrix: List[List[float]]) -> Tuple[List[int], float, bool]:
    best_hc = list(best_hc[:-1])
    n = len(best_hc)
    improved = False
    best_found_hc = best_hc
    best_found_score = best_score
        
    for k in range(n - 2):    
        k_ = k + 1
        d_min = distance_matrix[best_found_hc[k]][best_found_hc[k_]]

        for i in range(k + 2, n):
            d = distance_matrix[best_found_hc[k]][best_found_hc[i]]
            if (d < d_min):
                d_min = d
                k_ = i

        if (k_ != k + 1):
            reverse_hc = reverse_oper(solution=best_found_hc, k=k, k_=k_)
            reverse_score = cal_score_hc(new_hc=reverse_hc, distance_matrix=distance_matrix)

            if (reverse_score < best_found_score):
                best_found_hc = reverse_hc
                best_found_score = reverse_score
                improved = True
            
    best_found_hc.append(best_found_hc[0])
    return best_found_hc, best_found_score, improved


def cal_score_hc(new_hc: List[int], distance_matrix: List[List[float]]) -> float:
    dis = sum(distance_matrix[new_hc[i]][new_hc[i + 1]]
              for i in range(len(new_hc) - 1))
    dis += distance_matrix[new_hc[-1]][new_hc[0]]

    return dis
    

def reverse_oper(solution: List[int], k: int, k_: int) -> List[int]:
    new_sol = solution.copy()
    new_sol[k + 1: k_ + 1] = new_sol[k + 1: k_ + 1][::-1]
    return new_sol