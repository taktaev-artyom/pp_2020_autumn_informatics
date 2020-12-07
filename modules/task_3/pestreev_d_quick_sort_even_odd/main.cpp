// Copyright 2020 Pestreev Daniil
#include <mpi.h>
#include <gtest/gtest.h>
#include <gtest-mpi-listener.hpp>
#include <math.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>

#include "../../../modules/task_3/pestreev_d_quick_sort_even_odd/quick_sort_even_odd_merge.h"

TEST(TestQuickSort, IntroducedVector) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = {-8098, 1009, -160, -1, 179796, 2,
        -1603, 166, -7396, -199, -18348, -15, 659, 7, 82519};
    double startTime = MPI_Wtime();
    std::vector<int> resP = parallel_sorting(vec);
    double finishTime = MPI_Wtime();
    double paralleltime = finishTime - startTime;
    if (rank == 0) {
        std::vector<int> resS = vec;
        double startTime = MPI_Wtime();
        std::sort(resS.begin(), resS.end());
        double sequentialtime = MPI_Wtime() - startTime;
        std::cout << "Sequential sort:" << sequentialtime << std::endl;
        std::cout << "Parallel sort:" << paralleltime << std::endl;
        ASSERT_EQ(resS, resP);
    }
}

TEST(TestQuickSort, SameNumber1000) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec(1000);
    std::fill(vec.begin(), vec.end(), 8);
    double startTime = MPI_Wtime();
    std::vector<int> resP = parallel_sorting(vec);
    double finishTime = MPI_Wtime();
    double paralleltime = finishTime - startTime;
    if (rank == 0) {
        std::vector<int> resS = vec;
        double startTime = MPI_Wtime();
        std::sort(resS.begin(), resS.end());
        double sequentialtime = MPI_Wtime() - startTime;
        std::cout << "Sequential sort:" << sequentialtime << std::endl;
        std::cout << "Parallel sort:" << paralleltime << std::endl;
        ASSERT_EQ(resS, resP);
    }
}

TEST(TestQuickSort, randomVector3) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vecsize = 3;
    std::vector<int> vec = getRandomVector(vecsize);
    double startTime = MPI_Wtime();
    std::vector<int> resP = parallel_sorting(vec);
    double finishTime = MPI_Wtime();
    double paralleltime = finishTime - startTime;
    if (rank == 0) {
        std::vector<int> resS = vec;
        double startTime = MPI_Wtime();
        std::sort(resS.begin(), resS.end());
        double sequentialtime = MPI_Wtime() - startTime;
        std::cout << "Sequential sort:" << sequentialtime << std::endl;
        std::cout << "Parallel sort:" << paralleltime << std::endl;
        ASSERT_EQ(resS, resP);
    }
}

TEST(TestQuickSort, randomVector131) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vecsize = 131;
    std::vector<int> vec = getRandomVector(vecsize);
    double startTime = MPI_Wtime();
    std::vector<int> resP = parallel_sorting(vec);
    double finishTime = MPI_Wtime();
    double paralleltime = finishTime - startTime;
    if (rank == 0) {
        std::vector<int> resS = vec;
        double startTime = MPI_Wtime();
        std::sort(resS.begin(), resS.end());
        double sequentialtime = MPI_Wtime() - startTime;
        std::cout << "Sequential sort:" << sequentialtime << std::endl;
        std::cout << "Parallel sort:" << paralleltime << std::endl;
        ASSERT_EQ(resS, resP);
    }
}

TEST(TestQuickSort, randomVector100100) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vecsize = 100100;
    std::vector<int> vec = getRandomVector(vecsize);
    double startTime = MPI_Wtime();
    std::vector<int> resP = parallel_sorting(vec);
    double finishTime = MPI_Wtime();
    double paralleltime = finishTime - startTime;
    if (rank == 0) {
        std::vector<int> resS = vec;
        double startTime = MPI_Wtime();
        std::sort(resS.begin(), resS.end());
        double sequentialtime = MPI_Wtime() - startTime;
        std::cout << "Sequential sort:" << sequentialtime << std::endl;
        std::cout << "Parallel sort:" << paralleltime << std::endl;
        ASSERT_EQ(resS, resP);
    }
}

TEST(TestQuickSort, randomVector5234876) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vecsize = 5234876;
    std::vector<int> vec = getRandomVector(vecsize);
    double startTime = MPI_Wtime();
    std::vector<int> resP = parallel_sorting(vec);
    double finishTime = MPI_Wtime();
    double paralleltime = finishTime - startTime;
    if (rank == 0) {
        std::vector<int> resS = vec;
        double startTime = MPI_Wtime();
        std::sort(resS.begin(), resS.end());
        double sequentialtime = MPI_Wtime() - startTime;
        std::cout << "Sequential sort:" << sequentialtime << std::endl;
        std::cout << "Parallel sort:" << paralleltime << std::endl;
        ASSERT_EQ(resS, resP);
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
