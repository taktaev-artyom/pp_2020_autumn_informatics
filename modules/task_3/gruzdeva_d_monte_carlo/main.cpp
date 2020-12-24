// Copyright 2020 Gruzdeva Diana
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include "./monte_carlo.h"


TEST(Parallel_Operations_MPI, LINEAR_10_3) {
    int rank;
    int n = 1000;
    double error = calculateError(10, n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start = MPI_Wtime();
    double parallel_result = getParallelIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, linearFunction, time(0));
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << end - start << std::endl;
        start = MPI_Wtime();
        double sequential_result = getSequentialIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, linearFunction, time(0));
        double end = MPI_Wtime();
        std::cout << end - start << std::endl;
        std::cout << sequential_result << std::endl;
        ASSERT_LT(std::fabs(parallel_result - sequential_result), error);
    }
}

TEST(Parallel_Operations_MPI, LINEAR_10_6) {
    int rank;
    int n = 1000000;
    double error = calculateError(10, n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start = MPI_Wtime();
    double parallel_result = getParallelIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, linearFunction, time(0));
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << end - start << std::endl;
        start = MPI_Wtime();
        double sequential_result = getSequentialIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, linearFunction, time(0));
        double end = MPI_Wtime();
        std::cout << end - start << std::endl;
        std::cout << sequential_result << std::endl;
        std::cout << parallel_result << std::endl;
        ASSERT_LT(std::fabs(parallel_result - sequential_result), error);
    }
}

TEST(Parallel_Operations_MPI, POLINOM_10_3) {
    int rank;
    int n = 1000;
    double error = calculateError(300, n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start = MPI_Wtime();
    double parallel_result = getParallelIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, polinomFunction, time(0));
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << end - start << std::endl;
        start = MPI_Wtime();
        double sequential_result = getSequentialIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, polinomFunction, time(0));
        double end = MPI_Wtime();
        std::cout << end - start << std::endl;
        std::cout << sequential_result << std::endl;
        std::cout << parallel_result << std::endl;
        ASSERT_LT(std::fabs(parallel_result - sequential_result), error);
    }
}

TEST(Parallel_Operations_MPI, POLINOM_10_6) {
    int rank;
    int n = 1000000;
    double error = calculateError(300, n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start = MPI_Wtime();
    double parallel_result = getParallelIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, polinomFunction, time(0));
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << end - start << std::endl;
        start = MPI_Wtime();
        double sequential_result = getSequentialIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, polinomFunction, time(0));
        double end = MPI_Wtime();
        std::cout << end - start << std::endl;
        std::cout << sequential_result << std::endl;
        std::cout << parallel_result << std::endl;
        ASSERT_LT(std::fabs(parallel_result - sequential_result), error);
    }
}


TEST(Parallel_Operations_MPI, COMPOSITE_10_3) {
    int rank;
    int n = 1000;
    double error = calculateError(300, n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start = MPI_Wtime();
    double parallel_result = getParallelIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, compositeFunction, time(0));
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << end - start << std::endl;
        start = MPI_Wtime();
        double sequential_result = getSequentialIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, compositeFunction, time(0));
        double end = MPI_Wtime();
        std::cout << end - start << std::endl;
        std::cout << sequential_result << std::endl;
        std::cout << parallel_result << std::endl;
        ASSERT_LT(std::fabs(parallel_result - sequential_result), error);
    }
}

TEST(Parallel_Operations_MPI, COMPOSITE_10_6) {
    int rank;
    int n = 1000000;
    double error = calculateError(300, n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start = MPI_Wtime();
    double parallel_result = getParallelIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, compositeFunction, time(0));
    double end = MPI_Wtime();
    if (rank == 0) {
        std::cout << end - start << std::endl;
        start = MPI_Wtime();
        double sequential_result = getSequentialIntegral(0.0, 1.0, 0.0, 1.0, 0.0, 1.0, n, compositeFunction, time(0));
        double end = MPI_Wtime();
        std::cout << end - start << std::endl;
        std::cout << sequential_result << std::endl;
        std::cout << parallel_result << std::endl;
        ASSERT_LT(std::fabs(parallel_result - sequential_result), error);
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
