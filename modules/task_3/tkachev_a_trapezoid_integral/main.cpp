// Copyright 2020 Tkachev Alexey
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <cassert>
#include <cmath>
#include "../../../../modules/task_3/tkachev_a_trapezoid_integral/trapezoid_integral.h"


TEST(Trapezoid_Integral_Tests, eq1) {
    // SSS x dxdydz
    // x,y,z = 0
    // 2x + 2y + z - 6 = 0
    int rank, count_processes;
    double len_x, len_y, len_z;
    double parallel_integrate = 0.0;
    double the_answer = 6.75;
    double count_delta = 1000;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    int equation = 1;

    len_x = 2;
    len_y = 2;
    len_z = 6;
    the_answer = 6.75;
    

    double _time1 = MPI_Wtime();
    parallel_integrate = parallelIntegral3D(equation, count_delta, len_x, len_y, len_z);
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        parallel_integrate = parallel_integrate * 2;
        // double not_parallel_integrate = 0.0;
        // printf("PARRALEL TIME: %.5f\n", _time2-_time1);

        // _time1 = MPI_Wtime();
        // not_parallel_integrate = integral3D(equation, count_processes, 0, len_x, len_y, len_z);
        // _time2 = MPI_Wtime();

        // printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        printf("RESULT PARALLEL %.5f\n", parallel_integrate);
        // printf("RESULT NOT PARALLEL %.5f\n", not_parallel_integrate);
        // printf("DIFF = %.5f\n", parallel_integrate - not_parallel_integrate);
        // printf("EPSILON = %.5f\n", the_answer - parallel_integrate);
    }
}

// TEST(Trapezoid_Integral_Tests, eq2) {
//     // SSS x dxdydz
//     // x = [0, 5]
//     // y = [0, 5]
//     // z = [0, 5]
//     int rank, count_processes;
//     double parallel_integrate = 0.0;
//     int count_delta = 1000;

//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

//     int equation = 2;

//     double len_x = 6;
//     double len_y = 4;
//     double len_z = 5;
//     double the_answer = 125000;
 

//     double _time1 = MPI_Wtime();
//     parallel_integrate = parallelIntegral3D(equation, count_delta, len_x, len_y, len_z);
//     double _time2 = MPI_Wtime();

//     if (rank == 0) {
//         parallel_integrate = parallel_integrate;
//         double not_parallel_integrate = 0.0;
//         printf("PARRALEL TIME: %.5f\n", _time2-_time1);

//        //_time1 = MPI_Wtime();
//         //not_parallel_integrate = integral3D(equation, count_processes, 0, len_x, len_y, len_z);
//         //_time2 = MPI_Wtime();

//         //printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
//         printf("RESULT PARALLEL %.5f\n", parallel_integrate);
//        // printf("RESULT NOT PARALLEL %.5f\n", not_parallel_integrate);
//         //printf("DIFF = %.5f\n", parallel_integrate - not_parallel_integrate);
//        // printf("EPSILON = %.5f\n", the_answer - parallel_integrate);
//     }
// }

// TEST(Trapezoid_Integral_Tests, eq3) {
//     // SSS 1 *  dxdydz
//     // x * x + y * y = 4s
//     // x = [0, 2]
//     // y = [0, 3]
//     // z = [0, sqrt(4 - x * x)]
//     int rank, count_processes;
//     double len_x, len_y, len_z;
//     double parallel_integrate = 0.0;
//     double the_answer = -1;
//     double result;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

//     int equation = 3;

//     len_x = 2;
//     len_y = 3;
//     len_z = 2;
//     the_answer = 3 * 3.14;
 
//     double _time1 = MPI_Wtime();
//     parallel_integrate = parallelIntegral3D(equation, len_x, len_y, len_z);
//     double _time2 = MPI_Wtime();

//     if (rank == 0) {
//         parallel_integrate = parallel_integrate;
//         double not_parallel_integrate = 0.0;
//         printf("PARRALEL TIME: %.5f\n", _time2-_time1);

//         _time1 = MPI_Wtime();
//         not_parallel_integrate = integral3D(equation, count_processes, 0, len_x, len_y, len_z);
//         _time2 = MPI_Wtime();

//         printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
//         printf("RESULT PARALLEL %.5f\n", parallel_integrate);
//         printf("RESULT NOT PARALLEL %.5f\n", not_parallel_integrate);
//         printf("DIFF = %.5f\n", parallel_integrate - not_parallel_integrate);
//         printf("EPSILON = %.5f\n", the_answer - parallel_integrate);
//     }
// }

// TEST(Trapezoid_Integral_Tests, eq4) {
//     // ball. R=3
//     // V = 4 * pi / 3 * (R ^ 3)
//     int rank, count_processes;
//     double len_x, len_y, len_z;
//     double parallel_integrate = 0.0;
//     double the_answer = -1;
//     double result;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

//     int equation = 4;

//     len_x = 2;
//     len_y = 3;
//     len_z = 2.718 - 1;
//     the_answer = 10;
 
//     double _time1 = MPI_Wtime();
//     parallel_integrate = parallelIntegral3D(equation, len_x, len_y, len_z);
//     double _time2 = MPI_Wtime();

//     if (rank == 0) {
//         parallel_integrate = parallel_integrate;
//         double not_parallel_integrate = 0.0;
//         printf("PARRALEL TIME: %.5f\n", _time2-_time1);

//         _time1 = MPI_Wtime();
//         not_parallel_integrate = integral3D(equation, count_processes, 0, len_x, len_y, len_z);
//         _time2 = MPI_Wtime();

//         printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
//         printf("RESULT PARALLEL %.5f\n", parallel_integrate);
//         printf("RESULT NOT PARALLEL %.5f\n", not_parallel_integrate);
//         printf("DIFF = %.5f\n", parallel_integrate - not_parallel_integrate);
//         printf("EPSILON = %.5f\n", the_answer - parallel_integrate);
//     }
// }




int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}