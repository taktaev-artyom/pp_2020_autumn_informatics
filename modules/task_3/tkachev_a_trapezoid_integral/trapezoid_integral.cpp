// Copyright 2020 Tkachev Alexey
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <vector>
#include <ctime>
#include <iostream>
#include <random>
#include <cassert>
#include <cmath>
#include "../../../../modules/task_3/tkachev_a_trapezoid_integral/trapezoid_integral.h"


double integral3D(int equation, int count_processes, double count_data_x, double count_data_y, 
                double count_data_z, int rank, double dx, double dy, double dz) {

    assert(equation > 0 && count_processes > 0 && rank > 0 && dx > 0 && dy > 0 && dz > 0);
    double x = count_data_x * rank;
    double y = count_data_y * rank;
    double z = count_data_z * rank;
    double OXY = 0;
    double OZX = 0;
    double h = 0;
    double V = 0;

    if (equation == 1) {
        if (rank != (count_processes - 1)) {
            while (x <= count_data_x * (rank + 1) ) {
                OZX = ((x * x + 2 + ((x + dx) * (x + dx) + 2)) / 2 ) * dz;
                h = 2 - x;
                V += OZX  * h;
                x += dx;
                y += dy;
                //std::cout << "V " << V << std::endl;
            } 
        } else {
            while (x <= count_data_x * count_processes) {
                OZX = ((x * x + 2 + ((x + dx) * (x + dx) + 2)) / 2 ) * dz;
                h = y * y + 2;
                V += OZX  * h;
                x += dx;
                y += dy;
                //std::cout << "V " << V << std::endl;
            }  
        }
        
    } else {
        if (equation == 2) {
            if (rank != (count_processes - 1)) {
                while (x <= count_data_x * (rank + 1) ) {
                    OXY = ( ( 4 + 4 ) / 2 ) * dx;
                    h = 5;
                    V += OXY * h;
                    std::cout << "V " << V << std::endl;
                } 
            } else {
                while (x <= count_data_x * count_processes) {
                    OXY = ( ( 4 + 4 ) / 2 ) * dx;
                    h = 5;
                    V += OXY * h;
                    std::cout << "V " << V << std::endl;
                }  
            }
        } else {
            if (equation == 3) {
                double OXY = ( ( 3 + 3 ) / 2 ) * dx;
                double h = sqrt(4 - x * x);
            } else {
                if (equation == 4) {
                    double OXY = ((x * y * y * y + (x + dx) * (y + dy)  *(y + dy) * (y + dy)) /2 ) * dx;
                    double h = 2.718 - 1;
                }
            }
        }
    }
    return V;
}


double parallelIntegral3D(int equation, int count_delta, double len_x, double len_y, double len_z) {

    std::cout << "eq " << equation << std::endl;

    int count_processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    assert(count_processes > 0 && rank >= 0);

    double dx = len_x / count_delta;
    double dy = len_y / count_delta;
    double dz = len_z / count_delta;


    double count_data_x = dx * count_delta / count_processes;
    double count_data_y = dy * count_delta / count_processes;
    double count_data_z = dz * count_delta / count_processes;

    std::cout << "X[" << rank << "] "<< count_data_x << std::endl;
    std::cout << "Y[" << rank << "] "<< count_data_y << std::endl;
    std::cout << "Z[" << rank << "] "<< count_data_z << std::endl;

    double local_result;
    local_result = integral3D(equation, count_processes, count_data_x, count_data_y, count_data_z, rank, dx, dy, dz);

    std::cout << "local[" << rank << "] "<< local_result << std::endl;

    double parralel_integration_result;

    MPI_Reduce(&local_result, &parralel_integration_result, 1, 
    MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    return parralel_integration_result;
}