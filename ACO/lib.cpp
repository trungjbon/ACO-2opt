#include <bits/stdc++.h>

using namespace std;


// Tinh khoang cach Euclidean giua 2 diem
double euclidean_distance(vector<double> point1, vector<double> point2){
    double sum = 0;

    for (int i = 0; i < 2; i++){
        double d = point1[i] - point2[i];
        sum = sum + d * d;
    }

    return sqrt(sum);
}


// Tinh ma tran khoang cach
vector<vector<double>> pairwise_distances(vector<vector<double>> points1, vector<vector<double>> points2){
    int n1 = points1.size();
    int n2 = points2.size();

    vector<vector<double>> result(
        n1, vector<double>(n2, 0)
    );

    for (int i = 0; i < n1; i++){
        for (int j = 0; j < n2; j++){
            result[i][j] = euclidean_distance(
                    points1[i],
                    points2[j]
                );
        }
    }

    return result;
}


// Lam tron cac phan tu trong ma tran
vector<vector<double>> round_matrix(vector<vector<double>> matrix, int number_of_digits){
    double power = pow(10, number_of_digits);

    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            matrix[i][j] = round(matrix[i][j] * power) / power;
        }
    }

    return matrix;
}


// Gan gia tri vao duong cheo chinh
void fill_diagonal(vector<vector<double>>& matrix, double value){
    int n = matrix.size();

    for (int i = 0; i < n; i++){
        matrix[i][i] = value;
    }
}
