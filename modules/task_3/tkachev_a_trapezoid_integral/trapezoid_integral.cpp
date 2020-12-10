// Copyright 2020 Tkachev Alexey
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "../../../../modules/task_3/tkachev_a_trapezoid_integral/trapezoid_integral.h"


double integral3D(int integral, int count_processes, double count_data_x, double count_data_y, 
                double count_data_z, int rank, double dx, double dy, double dz) {
    assert(integral > 0);
    assert(count_processes > 0);
    assert(rank >= 0);
    assert(dx >= 0.0);
    assert(dy >= 0.0);
    assert(dz > 0.0);

    double x = count_data_x * rank;
    double y = count_data_y * rank;
    double z = count_data_z * rank;

    double OXY = 0.0;
    double h = 0.0;
    double V = 0.0;

    if (integral == 1) {
        if (rank != (count_processes - 1)) {
            while (x < count_data_x * (rank + 1) ) {
                OXY = ((1 - x + (1 - (dx + x)) ) / 2 ) * dx;
                h = x;
                V += OXY * h;
                x += dx;
            } 
        } else {
            while (x < count_data_x * count_processes) {
                OXY = ((1 - x + (1 - (dx + x)) ) / 2 ) * dx;
                h = x;
                V += OXY * h;
                x += dx;
            }  
        }
    } else {
        if (integral == 2) {
            if (rank != (count_processes - 1)) {
                while (x < count_data_x * (rank + 1) ) {
                    OXY = ( ( 3 + 3 ) / 2 ) * dx;
                    h = 4;
                    V += OXY * h;
                    x += dx;
                } 
            } else {
                while (x < count_data_x * count_processes) {
                    OXY = ( ( 3 + 3 ) / 2 ) * dx;
                    h = 4;
                    V += OXY * h;
                    x += dx;
                }  
            }
        } else {
            if (integral == 3) {
                if (rank != (count_processes - 1)) {
                    while (x < count_data_x * (rank + 1) ) {
                          OXY = ( ( 1 - x + (1 - (dx + x)) ) / 2 ) * dx;
                          h = y * y;
                          V += OXY * h;
                          y += dy;
                          x += dx;
                    } 
                } else {
                    while (x < count_data_x * count_processes) {
                         OXY = ( ( 1 - x + (1 - (dx + x)) ) / 2 ) * dx;
                         h = y * y;
                         V += OXY * h;
                         y += dy;
                         x += dx;
                    } 
                }
            } else {
                if (integral == 4) {
                    // +
                    if (rank != (count_processes - 1)) {
                        while (x < count_data_x * (rank + 1) ) {
                            OXY = ((3 - x + 3 - (dx + x)) / 2 ) * dx;
                            h = 1 - x * x;
                            V += OXY * h;
                            x += dx;
                            }
                        // -
                        x = rank * count_data_x;
                        while (x < count_data_x * (rank + 1) ) {
                            OXY = ((3 + x + 3 + (dx + x)) / 2 ) * dx;
                            h = 1 - x * x;
                            V += OXY * h;
                            x += dx;
                            }
                    } else {
                        while (x < count_data_x * count_processes ) {
                            OXY = ((3 - x + 3 - (dx + x)) / 2 ) * dx;
                            h = 1 - x * x;
                            V += OXY * h;
                            x += dx;
                        }
                        x = rank * count_data_x;
                        while (x < count_data_x * count_processes) {
                            OXY = ((3 + x + 3 + (dx + x)) / 2 ) * dx;
                            h = 1 - x * x;
                            V += OXY * h;
                            x += dx;
                        }
                    }
                } else {
                    if (integral == 5) {
                        if (rank != (count_processes - 1)) {
                            while (x < count_data_x * (rank + 1)) {
                                OXY = ( ( 2 - x + 2 - (dx + x) ) / 2 ) * dx;
                                h = 2 - x;
                                V += OXY * h;
                                x += dx;
                            } 
                        } else {
                            while (x < count_data_x * count_processes) {
                                OXY = ( ( 2 - x + 2 - (dx + x) ) / 2 ) * dx;
                                h = 2 - x;
                                V += OXY * h;
                                x += dx;
                            } 
                        } 
                    } else {
                        if (integral == 6) {
                            if (rank != (count_processes - 1)) {
                                while (z < count_data_z * (rank + 1)) {
                                    OXY = ( ( 6 + 6 ) / 2 ) * dz; // OXY is OZY
                                    h = 1;
                                    V += OXY * h;
                                    z += dz;
                                } 
                            } else {
                                while (z < count_data_z * count_processes) {
                                    OXY = ( ( 6 + 6 ) / 2 ) * dz; // OXY is OZY
                                    h = 1;
                                    V += OXY * h;
                                    z += dz;
                                } 
                            } 
                        } else {
                            if (integral == 7) {
                                if (rank != (count_processes - 1)) {
                                    while (y < count_data_y * (rank + 1)) {
                                        OXY = ( ( y + (y + dy) ) / 2 ) * dy; // OXY is OZY
                                        h = 1 - y;
                                        V += OXY * h;
                                        y += dy;
                                    } 
                                } else {
                                    while (y < count_data_y * count_processes) {
                                        OXY = ( ( y + (y + dy) ) / 2 ) * dy; // OXY is OZY
                                        h = 1 - y;
                                        V += OXY * h;
                                        y += dy;
                                    } 
                                }
                            } else {
                                if (integral == 8) {
                                    if (rank != (count_processes - 1)) {
                                        while (x < count_data_x * (rank + 1)) {
                                            OXY = ( (4 + 4) / 2 ) * dx;
                                            h = 1;
                                            V += OXY * h;
                                            x += dx;
                                        } 
                                    } else {
                                        while (x < count_data_x * count_processes) {
                                            OXY = ( (4 + 4) / 2 ) * dx;
                                            h = 1;
                                            V += OXY * h;
                                            x += dx;
                                        } 
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }   
    return V;
}

double parallelIntegral3D(int integral, double count_data_x, double count_data_y,
                     double count_data_z, double dx, double dy, double dz) {
    int count_processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    assert(count_processes > 0);
    assert(rank >= 0);

    double local_result = integral3D(integral, count_processes, count_data_x,
                                 count_data_y, count_data_z, rank, dx, dy, dz);

    double parralel_integration_result;

    MPI_Reduce(&local_result, &parralel_integration_result, 1, 
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    return parralel_integration_result;
}
