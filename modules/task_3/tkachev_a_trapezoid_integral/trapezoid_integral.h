// Copyright 2020 Tkachev Alexey
#ifndef MODULES_TASK_3_TKACHEV_A_TRAPEZOID_INTEGRAL_TRAPEZOID_INTEGRAL_H
#define MODULES_TASK_3_TKACHEV_A_TRAPEZOID_INTEGRAL_TRAPEZOID_INTEGRAL_H

double integral3D(int integral, int count_processes, double count_data_x, double count_data_y, 
               double count_data_z,  int rank, double dx, double dy, double dz);

double parallelIntegral3D(int integral, double count_data_x, double count_data_y,
                     double count_data_z, double dx, double dy, double dz);

#endif  // MODULES_TASK_3_TKACHEV_A_TRAPEZOID_INTEGRAL_TRAPEZOID_INTEGRAL_H