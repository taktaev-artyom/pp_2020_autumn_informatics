// Copyright 2020 Tkachev Alexey
#ifndef MODULES_TASK_3_TKACHEV_A_TRAPEZOID_INTEGRAL_TRAPEZOID_INTEGRAL_H
#define MODULES_TASK_3_TKACHEV_A_TRAPEZOID_INTEGRAL_TRAPEZOID_INTEGRAL_H

double integral3D(int equation, int count_processes, double count_data_x, double count_data_y, 
               double count_data_z,  int rank, double dx, double dy, double dz);

double parallelIntegral3D(int equation, int delta, double len_x, double len_y, double len_z);

#endif  // MODULES_TASK_3_TKACHEV_A_TRAPEZOID_INTEGRAL_TRAPEZOID_INTEGRAL_H