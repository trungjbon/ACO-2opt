import math

def euclidean_distance(point1: list, point2: list):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(point1, point2)))


def pairwise_distances(points1: list, points2: list):
    return [[euclidean_distance(p1, p2) for p2 in points2] for p1 in points1]


def round_matrix(matrix: list, afterDot: int):
    return [[round(x, afterDot) for x in row] for row in matrix]


def fill_diagonal(matrix: list, val: int):
    for i in range(len(matrix)):
        matrix[i][i] = val

