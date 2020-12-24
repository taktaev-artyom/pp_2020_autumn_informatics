// Copyright 2020 Grigoryan Garry
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>

#include "./sort.h"

TEST(Seq_Sort, sort_vector_positive) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (0 == rank) {
    int size = 9;
    double *vector = new double[size] { 4.0, 9.66057, 1.69, 99.147, 57.3579, 34.23, 44.43, 432.5435, 9.6 };
    bool result = false;
    double start = MPI_Wtime();
    seq_sorting(vector, size);
    double finish = MPI_Wtime();
    result = is_sorted(vector, size);
    ASSERT_EQ(true, result);
    std::cout << finish - start << std::endl;
  }
}

TEST(Seq_Sort, sort_large_vector) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (0 == rank) {
    int size = 1000;
    double *vector = new double[size];
    for (int i = 0; i < size; i++) vector[i] = (size - i) * 0.001 + (i + size / 2);
    bool result = false;
    double start = MPI_Wtime();
    seq_sorting(vector, size);
    result = is_sorted(vector, size);
    double finish = MPI_Wtime();
    ASSERT_EQ(true, result);
    std::cout << finish - start << std::endl;
  }
}

TEST(Seq_Sort, sort_is_correct) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (0 == rank) {
    int size = 5;
    double *result = new double[size] { 1.69, 2.0, 8.66667, 57.3579, 99.147};
    double *vector = new double[size] { 2.0, 8.66667, 1.69, 99.147, 57.3579 };
    double start = MPI_Wtime();
    seq_sorting(vector, size);
    double finish = MPI_Wtime();
    for (int i = 0; i < size; i++) {
      ASSERT_NEAR(result[i], vector[i], 0.01);
    }
    std::cout << finish - start << std::endl;
  }
}

TEST(Par_Sort, sort_vector) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size = 9;
  double *vector = new double[size] { 99.147, 9.66057, 1.69, 4.0, 57.3579, 34.23, 44.43, 432.5435, 96.005 };
  double start = MPI_Wtime();
  par_sorting(&vector, size);
  double finish = MPI_Wtime();
  if (0 == rank) {
    bool result = false;
    result = is_sorted(vector, size);
    ASSERT_EQ(true, result);
    std::cout << finish - start << std::endl;
  }
  delete[] vector;
}

TEST(Par_Sort, sort_large_vector) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int size = 1000;
  double *vector = new double[size];
  for (int i = 0; i < size; i++) vector[i] = (size - i) * 0.001 + (i + size / 2);
  double start = MPI_Wtime();
  par_sorting(&vector, size);
  double finish = MPI_Wtime();
  if (0 == rank) {
    bool result = false;
    result = is_sorted(vector, size);
    ASSERT_EQ(true, result);
    std::cout << finish - start << std::endl;
  }
  delete[] vector;
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
