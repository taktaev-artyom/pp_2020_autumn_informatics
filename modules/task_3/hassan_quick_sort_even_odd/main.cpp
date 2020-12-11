// Copyright 2020 Hassan EzzAldeen
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./quick_sort_even_odd.h"

TEST(Quick_Sort_Batcher, Test_Quick_Sort_1) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int vec_size = 100;
  std::vector<int> vec_1(vec_size);
  std::vector<int> vec_2(vec_size);

  if (rank == 0) {
    vec_1 = createRandomVector(vec_size);
    vec_2 = vec_1;
    quickSort(&vec_1, 0, 99);
    std::sort(vec_2.begin(), vec_2.end());
    ASSERT_EQ(vec_1, vec_2);
  }
}

TEST(Quick_Sort_Batcher, Test_Quick_Sort_2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int vec_size = 10;
  std::vector<int> vec_1(vec_size);
  std::vector<int> vec_2(vec_size);

  if (rank == 0) {
    vec_1 = { 8, 4, 17, 3, 11, 2, 1, 9, 7, 20 };
    vec_2 = vec_1;
    quickSort(&vec_1, 0, 9);
    std::sort(vec_2.begin(), vec_2.end());
    ASSERT_EQ(vec_1, vec_2);
  }
}

TEST(Quick_Sort_Batcher, Test_Batcher_Sort_1) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int vec_size = 100;
  std::vector<int> vec_1(vec_size);
  std::vector<int> vec_2(vec_size);

  if (rank == 0) {
    vec_1 = createRandomVector(vec_size);
    vec_2 = vec_1;
  }

  quickSortBatcher(&vec_1);
  if (rank == 0) {
    std::sort(vec_2.begin(), vec_2.end());
    ASSERT_EQ(vec_1, vec_2);
  }
}

TEST(Quick_Sort_Batcher, Test_Batcher_Sort_2) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int vec_size = 1000;
  std::vector<int> vec_1(vec_size);
  std::vector<int> vec_2(vec_size);

  if (rank == 0) {
    vec_1 = createRandomVector(vec_size);
    vec_2 = vec_1;
  }

  quickSortBatcher(&vec_1);
  if (rank == 0) {
    std::sort(vec_2.begin(), vec_2.end());
    ASSERT_EQ(vec_1, vec_2);
  }
}

TEST(Quick_Sort_Batcher, TimeTest) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int vecsize = 100000;
  std::vector<int> vec(vecsize);
  std::vector<int> vec2(vecsize);

  if (rank == 0) {
    vec = createRandomVector(vecsize);
    vec2 = vec;
  }
  double startTime = MPI_Wtime();
  quickSortBatcher(&vec);
  double finishTime = MPI_Wtime();

  if (rank == 0) {
    std::cout << "Parallel sort:" <<  finishTime - startTime << std::endl;
    startTime = MPI_Wtime();
    std::sort(vec2.begin(), vec2.end());
    // quickSort(&vec2, 0, 99999);
    double finishTime = MPI_Wtime();
    std::cout << "Quick sort:" << finishTime - startTime << std::endl;
    ASSERT_EQ(vec2, vec);
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
