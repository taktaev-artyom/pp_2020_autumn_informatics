// Copyright 2020 Kirillov Konstantin
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include<iostream>
#include "./dijkstrat_algorithm.h"
TEST(Parallel_Operations_MPI, Test_Size_4) {
    int procRank, procNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    int sizeGraph = 5;
    if (sizeGraph < procNum) {
        procNum = sizeGraph;
    }
    std::vector <int> graph(sizeGraph*sizeGraph);
    std::vector <int> seqDist(sizeGraph);
    std::vector <int> parDist(sizeGraph);
    if (procRank == 0) {
        graph = getRandomGraph(sizeGraph);
        double start_time = MPI_Wtime();
        seqDist = getSequentialDijkstras(graph, 0);
        double end_time = MPI_Wtime();
        std::cout << "Seqential: " << end_time - start_time << std::endl;
        // printGraph(graph);
        printDist(seqDist);
    }
    double start_time = MPI_Wtime();
    parDist = getParallelDijkstras(graph, 0);
    double end_time = MPI_Wtime();
    if (procRank == 0) {
        std::cout << "Parallel: " << end_time - start_time << std::endl;
        printDist(parDist);
        for (int i = 0; i < sizeGraph; i++) {
            ASSERT_EQ(seqDist[i], parDist[i]);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_Size_17) {
    int procRank, procNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    int sizeGraph = 17;
    if (sizeGraph < procNum) {
        procNum = sizeGraph;
    }
    std::vector <int> graph(sizeGraph*sizeGraph);
    std::vector <int> seqDist(sizeGraph);
    std::vector <int> parDist(sizeGraph);
    if (procRank == 0) {
        graph = getRandomGraph(sizeGraph);
        double start_time = MPI_Wtime();
        seqDist = getSequentialDijkstras(graph, 0);
        double end_time = MPI_Wtime();
        std::cout << "Seqential: " << end_time - start_time << std::endl;
        printDist(seqDist);
    }
    double start_time = MPI_Wtime();
    parDist = getParallelDijkstras(graph, 0);
    double end_time = MPI_Wtime();
    if (procRank == 0) {
        std::cout << "Parallel: " << end_time - start_time << std::endl;
        printDist(parDist);
        for (int i = 0; i < sizeGraph; i++) {
            ASSERT_EQ(seqDist[i], parDist[i]);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_Size_11) {
    int procRank, procNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    int sizeGraph = 11;
    if (sizeGraph < procNum) {
        procNum = sizeGraph;
    }
    std::vector <int> graph(sizeGraph*sizeGraph);
    std::vector <int> seqDist(sizeGraph);
    std::vector <int> parDist(sizeGraph);
    if (procRank == 0) {
        graph = getRandomGraph(sizeGraph);
        double start_time = MPI_Wtime();
        seqDist = getSequentialDijkstras(graph, 0);
        double end_time = MPI_Wtime();
        std::cout << "Seqential: " << end_time - start_time << std::endl;
        printDist(seqDist);
    }
    double start_time = MPI_Wtime();
    parDist = getParallelDijkstras(graph, 0);
    double end_time = MPI_Wtime();
    if (procRank == 0) {
        std::cout << "Parallel: " << end_time - start_time << std::endl;
        printDist(parDist);
        for (int i = 0; i < sizeGraph; i++) {
            ASSERT_EQ(seqDist[i], parDist[i]);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_Size_23) {
    int procRank, procNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    int sizeGraph = 23;
    if (sizeGraph < procNum) {
        procNum = sizeGraph;
    }
    std::vector <int> graph(sizeGraph*sizeGraph);
    std::vector <int> seqDist(sizeGraph);
    std::vector <int> parDist(sizeGraph);
    if (procRank == 0) {
        graph = getRandomGraph(sizeGraph);
        double start_time = MPI_Wtime();
        seqDist = getSequentialDijkstras(graph, 0);
        double end_time = MPI_Wtime();
        std::cout << "Seqential: " << end_time - start_time << std::endl;
        printDist(seqDist);
    }
    double start_time = MPI_Wtime();
    parDist = getParallelDijkstras(graph, 0);
    double end_time = MPI_Wtime();
    if (procRank == 0) {
        std::cout << "Parallel: " << end_time - start_time << std::endl;
        printDist(parDist);
        for (int i = 0; i < sizeGraph; i++) {
            ASSERT_EQ(seqDist[i], parDist[i]);
        }
    }
}

TEST(Parallel_Operations_MPI, Test_Size_44) {
    int procRank, procNum;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    int sizeGraph = 44;
    if (sizeGraph < procNum) {
        procNum = sizeGraph;
    }
    std::vector <int> graph(sizeGraph*sizeGraph);
    std::vector <int> seqDist(sizeGraph);
    std::vector <int> parDist(sizeGraph);
    if (procRank == 0) {
        graph = getRandomGraph(sizeGraph);
        double start_time = MPI_Wtime();
        seqDist = getSequentialDijkstras(graph, 0);
        double end_time = MPI_Wtime();
        std::cout << "Seqential: " << end_time - start_time << std::endl;
        printDist(seqDist);
    }
    double start_time = MPI_Wtime();
    parDist = getParallelDijkstras(graph, 0);
    double end_time = MPI_Wtime();
    if (procRank == 0) {
        std::cout << "Parallel: " << end_time - start_time << std::endl;
        printDist(parDist);
        for (int i = 0; i < sizeGraph; i++) {
            ASSERT_EQ(seqDist[i], parDist[i]);
        }
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners &listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
