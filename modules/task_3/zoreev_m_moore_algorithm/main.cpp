// Copyright 2020 Zoreev Mikhail

#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <mpi.h>

#include "../../../modules/task_3/zoreev_m_moore_algorithm/moore_algorithm.h"

TEST(Moore_Algotithm_MPI, Graph_4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int64_t graph[4 * 4];
    if (rank == 0) {
        randomCompleteGraph(4, graph);
        printGraph(4, graph);
        printPredecessor(4, mooreAlgorithm(4, graph, 0));
    }
    MPI_Bcast(graph, 4*4, MPI_UINT64_T, 0, MPI_COMM_WORLD);
    size_t* result = mooreAlgorithmParallel(4, graph, 0);
    if (rank == 0) {
        printPredecessor(4, result);
    }
    
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}