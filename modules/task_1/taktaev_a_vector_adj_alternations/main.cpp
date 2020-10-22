// Copyright 2020 Taktaev Artem
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "../../../modules/task_1/taktaev_a_vector_adj_alternations/vector_adj_alternations.h"

TEST(Parallel_Adj_Alternations_MPI, Test_Size_10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int vec_size = 10;
    std::vector<int> full_vec(vec_size);

    if (rank == 0) full_vec = createRandomVector(vec_size);
    int count_parallel = calculateAdjAlternationsParallel(full_vec, vec_size);
    if (rank == 0) {
        int count_sequential = calculateAdjAlternationsSequential(full_vec, 1, 1);
        ASSERT_EQ(count_parallel, count_sequential);
    }
}

TEST(Parallel_Adj_Alternations_MPI, Test_Wrong_Size_Rand) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int vec_size = 0;

    ASSERT_ANY_THROW(createRandomVector(vec_size));
}

TEST(Parallel_Adj_Alternations_MPI, Test_Wrong_Inc) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int vec_size = 10;
    std::vector<int> full_vec(vec_size);

    int inc = -1;

    if (rank == 0) full_vec = createRandomVector(vec_size);
    ASSERT_ANY_THROW(calculateAdjAlternationsSequential(full_vec, -inc, 1));
}

TEST(Parallel_Adj_Alternations_MPI, Test_Wrong_Start_Index) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int vec_size = 10;
    std::vector<int> full_vec(vec_size);

    int start_index = -1;

    if (rank == 0) full_vec = createRandomVector(vec_size);
    ASSERT_ANY_THROW(calculateAdjAlternationsSequential(full_vec, 1, start_index));
}

TEST(Parallel_Adj_Alternations_MPI, Test_Size_2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int vec_size = 2;
    std::vector<int> full_vec(vec_size);

    if (rank == 0) full_vec = createRandomVector(vec_size);
    int count_parallel = calculateAdjAlternationsParallel(full_vec, vec_size);
    if (rank == 0) {
        int count_sequential = calculateAdjAlternationsSequential(full_vec, 1, 1);
        ASSERT_EQ(count_parallel, count_sequential);
    }
}

TEST(Parallel_Adj_Alternations_MPI, Test_Wrong_Size_Par) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int vec_size = 0;
    std::vector<int> full_vec(vec_size);

    ASSERT_ANY_THROW(calculateAdjAlternationsParallel(full_vec, vec_size));
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
