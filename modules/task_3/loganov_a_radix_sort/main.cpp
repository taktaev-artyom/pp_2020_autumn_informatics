// Copyright 2020 Loganov Andrei
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include "./radix.h"
TEST(TEST_PARALEL_MPI, TEST1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 100;
    std::vector<double> t = getRandomVector(n);
    std::vector<double> lin;
    std::vector<double> Par;
    Par = ParallelSort(t);
    if (rank == 0) {
        lin = seqRadixSort(t);
        ASSERT_EQ(Par, lin);
        }
}
TEST(TEST_PARALEL_MPI, TEST2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 221;
    std::vector<double> t = getRandomVector(n);
    std::vector<double> lin;
    std::vector<double> Par;
    Par = ParallelSort(t);
    if (rank == 0) {
        lin = seqRadixSort(t);
        ASSERT_EQ(Par, lin);
        }
}
TEST(TEST_PARALEL_MPI, TEST3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 500;
    std::vector<double> t = getRandomVector(n);
    std::vector<double> lin;
    std::vector<double> Par;
    Par = ParallelSort(t);
    if (rank == 0) {
        lin = seqRadixSort(t);
        ASSERT_EQ(Par, lin);
        }
}
TEST(TEST_PARALEL_MPI, TEST4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 1036;
    std::vector<double> t = getRandomVector(n);
    std::vector<double> lin;
    std::vector<double> Par;
    Par = ParallelSort(t);
    if (rank == 0) {
        lin = seqRadixSort(t);
        ASSERT_EQ(Par, lin);
        }
}
TEST(TEST_PARALEL_MPI, TEST5) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 1500;
    std::vector<double> t = getRandomVector(n);
    std::vector<double> lin;
    std::vector<double> Par;
    Par = ParallelSort(t);
    if (rank == 0) {
        lin = seqRadixSort(t);
        ASSERT_EQ(Par, lin);
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
