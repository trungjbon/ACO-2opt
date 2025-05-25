from typing import List
def read_map(file_path: str) -> List[List[float]]:
    cities = []

    with open(file_path, "r") as file:
        lines = file.readlines()

        start_idx = lines.index("NODE_COORD_SECTION\n") + 1
        end_idx = lines.index("EOF\n")

        for line in lines[start_idx: end_idx]:
            data = line.strip().split()

            cities.append([float(data[1]), float(data[2])])
    
    return cities


PATH_TO_MAP = "C:/Users/trung/DAA/Final_Exam/SOS-ACO-Algorithm-master/Benchmark/kroA100.tsp"
TOWNS = read_map(PATH_TO_MAP)


# Parameter for ACO
n = len(TOWNS)
ANTS = n // 4
EVAPORATION_RATE = 0.1
INTENSIFICATION = 5000
MAX_ITER_ACO = 200  # The number of iterations
BETA_EVAPORATION_RATE = 0.001
ALPHA = 1.0
BETA = 5.0