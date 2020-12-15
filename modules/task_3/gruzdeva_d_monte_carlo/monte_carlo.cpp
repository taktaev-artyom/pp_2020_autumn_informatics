// Copyright 2020 Gruzdeva Diana

#include "../../../modules/task_3/gruzdeva_d_monte_carlo/monte_carlo.h"
#include <random>
#include <ctime>

double calculateError(int p, int n) {
    double error = p * sqrt(0.125 / n);
    return error;
}

double getSequentialIntegral(double x1, double x2,
          double y1, double y2, double z1, double z2, int n,
          std::function<double(double, double, double)> function, time_t seed) {
    std::mt19937 gen(seed);
    std::normal_distribution<double> normal(0.5, 0.125);
    double result = 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double x = x1 + (x2 - x1) * normal(gen);
        double y = y1 + (y2 - y1) * normal(gen);
        double z = z1 + (z2 - z1) * normal(gen);
        sum += function(x, y, z);
    }
    result = sum * (x2 - x1) * (y2 - y1) * (z2 - z1) / n;
    return result;
}

double getParallelIntegral(double x1, double x2,
          double y1, double y2, double z1, double z2, int n,
          std::function<double(double, double, double)> function, time_t seed) {
    int size, rank;
    double height = (z2 - z1) / n;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double result = 0;
    int local_n = n / size,
        remaining = n - (size - 1) * local_n;
    double local_result = 0.0,
           local_z1 = z1 + rank * local_n * height,
           local_z2 = local_z1 + local_n * height;
    if (rank != size - 1) {
        local_result = getSequentialIntegral(x1, x2, y1, y2, local_z1, local_z2, local_n, function, seed);
    } else {
        local_result = getSequentialIntegral(x1, x2, y1, y2, local_z1, z2, remaining, function, seed);
    }
    MPI_Reduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return result;
}

double callFunction(std::function<double(double, double, double)> function, double x, double y, double z) {
    return function(x, y, z);
}

double linearFunction(double x, double y, double z) {
    return (x + y + z);
}

double polinomFunction(double x, double y, double z) {
    return x * x + y - z * z * z + 6;
}

double compositeFunction(double x, double y, double z) {
    return sin(x * y * z);
}
