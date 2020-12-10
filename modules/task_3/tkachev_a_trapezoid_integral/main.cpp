// Copyright 2020 Tkachev Alexey
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <stdlib.h>
#include <ctime>
#include <iostream>
#include <cassert>
#include "../../../../modules/task_3/tkachev_a_trapezoid_integral/trapezoid_integral.h"


TEST(Trapezoid_Integral_Tests, test_integral_1) {
    // 15 * (y * y + z * z)dxdyzy
    // z = x + y, x + y = 1, x,y,z = 0
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    int integral = 1;

    const double count_delta = 1000;
    const double EPSILON = 0.6;
    const double the_answer = 2;
    
    const double len_x = 1;
    const double len_y = 1;
    const double len_z = 1;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    parallel_integrate *= 15;
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.5f\n", _time2-_time1);
        _time1 = MPI_Wtime();
        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);
        _time2 = MPI_Wtime();
        printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

TEST(Trapezoid_Integral_Tests, test_integral_2) {
    // dxdydz
    // x = [0, 2], y = [0, 3], z = [0, 4]
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    int integral = 2;

    const double count_delta = 1000;
    const double EPSILON = 0.4;
    const double the_answer = 2 * 3 * 4;
    
    const double len_x = 2;
    const double len_y = 3;
    const double len_z = 4;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.6f\n", _time2-_time1);
        _time1 = MPI_Wtime();
        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);
        _time2 = MPI_Wtime();
        printf("NOT PARRALEL TIME %.6f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

TEST(Trapezoid_Integral_Tests, test_integral_3) {
    // dxdydz
    // x,y,z = 0, y = 1 - x, z = y * y
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    int integral = 3;

    const double count_delta = 1000;
    const double EPSILON = 0.01;
    const double the_answer = 0.083333;
    
    const double len_x = 1;
    const double len_y = 1;
    const double len_z = 1;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.5f\n", _time2-_time1);
        _time1 = MPI_Wtime();
        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);

        _time2 = MPI_Wtime();
        printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

TEST(Trapezoid_Integral_Tests, test_integral_4) {
    // dxdydz
    // z = 0, z = 1 - x * x, y = 0, y = 3 - x
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    int integral = 4;

    const double count_delta = 1000;
    const double EPSILON = 0.25;
    const double the_answer = 4;
    
    const double len_x = 1;
    const double len_y = 3;
    const double len_z = 1;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.5f\n", _time2-_time1);
        _time1 = MPI_Wtime();
        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);
        _time2 = MPI_Wtime();
        printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

TEST(Trapezoid_Integral_Tests, test_integral_5) {
    // (y + z)dxdydz
    // x, y = 0, z = 2, x + y + z = 4
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    const double count_delta = 1000;
    const double EPSILON = 0.7;
    const double the_answer = 2;
    
    int integral = 5;

    const double len_x = 2;
    const double len_y = 2;
    const double len_z = 4 - 2;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.5f\n", _time2-_time1);
        _time1 = MPI_Wtime();
        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);
        _time2 = MPI_Wtime();
        printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

TEST(Trapezoid_Integral_Tests, test_integral_6) {
    // dydz ~ 6 dz
    // y,z = 0, y = 6, z = 5
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    const double count_delta = 1000;
    const double EPSILON = 1;
    const double the_answer = 30;
    
    int integral = 6;

    const double len_x = 1;
    const double len_y = 6;
    const double len_z = 5;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.5f\n", _time2-_time1);

        _time1 = MPI_Wtime();

        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);

        _time2 = MPI_Wtime();

        printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

TEST(Trapezoid_Integral_Tests, test_integral_7) {
    // but using OZY instead of OXY
    // 15 * (y * y + z * z)dxdyzy
    // z = x + y, x + y = 1, x,y,z = 0
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    int integral = 7;

    const double count_delta = 1000;
    const double EPSILON = 0.65;
    const double the_answer = 2;
    
    const double len_x = 1;
    const double len_y = 1;
    const double len_z = 1;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    parallel_integrate *= 15;
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.5f\n", _time2-_time1);
        _time1 = MPI_Wtime();
        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);

        _time2 = MPI_Wtime();
        printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

TEST(Trapezoid_Integral_Tests, test_integral_8) {
    // 4 dx ~ dxdy
    // x, y = 0, x = 3, y = 4
    int rank, count_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count_processes);

    int integral = 8;

    const double count_delta = 1000;
    const double EPSILON = 0.25;
    const double the_answer = 3 * 4;
    
    const double len_x = 3;
    const double len_y = 4;
    const double len_z = 1;

    const double dx = len_x / count_delta;
    const double dy = len_y / count_delta;
    const double dz = len_z / count_delta;

    const double count_data_x = dx * count_delta / count_processes;
    const double count_data_y = dy * count_delta / count_processes;
    const double count_data_z = dz * count_delta / count_processes;
    
    double _time1 = MPI_Wtime();
    double parallel_integrate = parallelIntegral3D(integral, count_data_x, count_data_y,
                                                count_data_z, dx, dy, dz);
    double _time2 = MPI_Wtime();

    if (rank == 0) {
        double not_parallel_integrate = 0.0;
        printf("PARRALEL TIME: %.5f\n", _time2-_time1);
        _time1 = MPI_Wtime();
        not_parallel_integrate = integral3D(integral, count_processes, count_data_x, 
                                        count_data_y, count_data_z, 0, dx, dy, dz);
        _time2 = MPI_Wtime();
        printf("NOT PARRALEL TIME %.5f\n", _time2-_time1);
        ASSERT_NEAR(the_answer, parallel_integrate, EPSILON);
    }
}

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
