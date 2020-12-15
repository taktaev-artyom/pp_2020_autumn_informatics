// Copyright 2020 Gruzdeva Diana
#ifndef MODULES_TASK_3_GRUZDEVA_D_MONTE_CARLO_MONTE_CARLO_H_
#define MODULES_TASK_3_GRUZDEVA_D_MONTE_CARLO_MONTE_CARLO_H_

#include <mpi.h>
#include <random>
#include <ctime>
#include <functional>

double calculateError(int p, int n);

double getSequentialIntegral(double x1, double x2,
          double y1, double y2, double z1, double z2, int n,
          std::function<double(double, double, double)> function, time_t seed);
double getParallelIntegral(double x1, double x2,
        double y1, double y2, double z1, double z2, int n,
        std::function<double(double, double, double)> function, time_t seed);

double callFunction(std::function<double(double)> function, double value);
double linearFunction(double x, double y, double z);
double polinomFunction(double x, double y, double z);
double compositeFunction(double x, double y, double z);

#endif  // MODULES_TASK_3_GRUZDEVA_D_MONTE_CARLO_MONTE_CARLO_H_
