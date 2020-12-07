// Copyright 2020 Napylov Evgenii
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include "./radix_sort_int_sm.h"

TEST(Radix_sort_MPI, Test_one_value) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = {7};
    std::vector<int> res = RadixSortParallel(vec);
    if (rank == 0) {
        std::qsort(vec.data(), vec.size(), sizeof(int), compare);
        for (int i = 0; i < static_cast<int>(vec.size()); i++) {
            ASSERT_EQ(res[i], vec[i]);
        }
    }
}

TEST(Radix_sort_MPI, Test_identical_values) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = {7, 7, 7, 7, 7};
    std::vector<int> res = RadixSortParallel(vec);
    if (rank == 0) {
        std::qsort(vec.data(), vec.size(), sizeof(int), compare);
        for (int i = 0; i < static_cast<int>(vec.size()); i++) {
            ASSERT_EQ(res[i], vec[i]);
        }
    }
}

TEST(Radix_sort_MPI, Test_sorted) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = {1, 2, 3, 4, 5};
    std::vector<int> res = RadixSortParallel(vec);
    if (rank == 0) {
        std::qsort(vec.data(), vec.size(), sizeof(int), compare);
        for (int i = 0; i < static_cast<int>(vec.size()); i++) {
            ASSERT_EQ(res[i], vec[i]);
        }
    }
}

TEST(Radix_sort_MPI, Test_random_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = RandomVector(100);
    std::vector<int> res = RadixSortParallel(vec);
    if (rank == 0) {
        std::qsort(vec.data(), vec.size(), sizeof(int), compare);
        for (int i = 0; i < static_cast<int>(vec.size()); i++) {
            ASSERT_EQ(res[i], vec[i]);
        }
    }
}

TEST(Radix_sort_MPI, Test_random_1234) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = RandomVector(1234);
    std::vector<int> res = RadixSortParallel(vec);
    if (rank == 0) {
        std::qsort(vec.data(), vec.size(), sizeof(int), compare);
        for (int i = 0; i < static_cast<int>(vec.size()); i++) {
            ASSERT_EQ(res[i], vec[i]);
        }
    }
}

TEST(Radix_sort_MPI, Test_performance_100000) {
    int rank;
    double t0, t1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = RandomVector(100000);
    t0 = MPI_Wtime();
    std::vector<int> res = RadixSortParallel(vec);
    t1 = MPI_Wtime();
    if (rank == 0) {
        std::cout << "par_time: " << t1 - t0 << std::endl;
        std::vector<int> res_seq;
        t0 = MPI_Wtime();
        res_seq = RadixSort(vec);
        t1 = MPI_Wtime();
        std::cout << "seq_time: " << t1 - t0 << std::endl;
        std::qsort(vec.data(), vec.size(), sizeof(int), compare);
        for (int i = 0; i < static_cast<int>(vec.size()); i++) {
            ASSERT_EQ(res[i], vec[i]);
        }
        for (int i = 0; i < static_cast<int>(vec.size()); i++) {
            ASSERT_EQ(res_seq[i], vec[i]);
        }
    }
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
