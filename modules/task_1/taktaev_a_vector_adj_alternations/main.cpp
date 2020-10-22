// Copyright 2020 Taktaev Artem
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include "../../../modules/task_1/taktaev_a_vector_adj_alternations/vector_adj_alternations.h"

TEST(Parallel_Symbol_Count_MPI, Test_10_Same_Symbols_String) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int count_size_string = 10;
    std::vector<int> global_string(count_size_string);
    
    if (rank == 0) global_string = createRandomVector(count_size_string);
    for(int i = 0; i < count_size_string; i++)
        std::cout << global_string[i] << std::endl;

    int global_sum = calculateAdjAlternationsParallel(global_string, count_size_string);

    if (rank == 0) {
        int reference_sum = calculateAdjAlternationsSequential(global_string, 1, 1);
        ASSERT_EQ(reference_sum, global_sum);
    }
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
