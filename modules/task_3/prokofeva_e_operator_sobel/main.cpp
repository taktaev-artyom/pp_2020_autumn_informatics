// Copyright 2020 Prokofeva Elizaveta
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include "random"
#include "./operator_sobel.h"

TEST(test, test1)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 3;
    int cols = 3;
    std::vector<int> test(rows * cols);
    test = getRandomMatrix(rows, cols);
    std::vector<int> paral(rows * cols);
    paral = parallel_sobel(test, rows, cols);
    if (rank== 0)
    {
        std::vector<int> seq(rows * cols);
        seq = calc_sobel(test, rows, cols);
        ASSERT_EQ(seq, paral);
    }
}

TEST(test, test2)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 10;
    int cols = 10;
    std::vector<int> test(rows * cols);
    test = getRandomMatrix(rows, cols);
    std::vector<int> paral(rows * cols);
    paral = parallel_sobel(test, rows, cols);
    if (rank == 0)
    {
        std::vector<int> seq(rows * cols);
        seq = calc_sobel(test, rows, cols);
        ASSERT_EQ(seq, paral);
    }
}
TEST(test, test3)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 1000;
    int cols = 1000;
    std::vector<int> test(rows * cols);
    test = getRandomMatrix(rows, cols);
    std::vector<int> paral(rows * cols);
    double start1 = MPI_Wtime();
    paral = parallel_sobel(test, rows, cols);
    double finish1 = MPI_Wtime();
    if (rank == 0)
    {
        std::vector<int> seq(rows * cols);
        double start2 = MPI_Wtime();
        seq = calc_sobel(test, rows, cols);
        double finish2 = MPI_Wtime();
        ASSERT_EQ(seq, paral);
        std::cout << "Parallel: " << finish1 - start1 << std::endl;
        std::cout << "Sequence: " << finish2 - start2 << std::endl;
    }
}
TEST(test, test4)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 4;
    int cols = 4;
    std::vector<int> test(rows * cols);
    test[0] = 15;
    test[1] = 74;
    test[2] = 123;
    test[3] = 98;
    test[4] = 2;
    test[5] = 10;
    test[6] = 36;
    test[7] = 89;
    test[8] = 55;
    test[9] = 47;
    test[10] = 25;
    test[11] = 81;
    test[12] = 16;
    test[13] = 158;
    test[14] = 200;
    test[15] = 247;
    std::vector<int> paral(rows * cols);
    paral = parallel_sobel(test, rows, cols);
    if (rank == 0)
    {
        std::vector<int> seq(rows * cols);
        seq = calc_sobel(test, rows, cols);
        ASSERT_EQ(seq, paral);
    }
}
TEST(test, test5)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 5;
    int cols = 5;
    std::vector<int> test(rows * cols);
    test = getRandomMatrix(rows, cols);
    std::vector<int> paral(rows * cols);
    paral = parallel_sobel(test, rows, cols);
    if (rank == 0)
    {
        std::vector<int> seq(rows * cols);
        seq = calc_sobel(test, rows, cols);
        ASSERT_EQ(seq, paral);
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
