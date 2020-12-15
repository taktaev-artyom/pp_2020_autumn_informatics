// Copyright 2020 Schekotilova Julia
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./strassen_alg.h"

TEST(TEST_STRASSEN_ALG, Test_addition_of_2_matrix) {
  double eps = 0.001;
  int n = 2;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::vector<double> A = { 1.1, 2.2, 3.3, 4.4};
  std::vector<double> B = { -1.1, -2.2, -3.3, -4.4};
  std::vector<double> C;
  std::vector<double> res = { 0, 0, 0, 0 };
  if (rank == 0) {
    ASSERT_NO_THROW(C = addition_of2m(A, B, n));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        EXPECT_NEAR(C[i * n + j], res[i * n + j], eps);
      }
    }
  }
}

TEST(TEST_STRASSEN_ALG, Test_subtraction_of_2_matrix) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int n = 4;
  std::vector<double> A = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
      16 };
  std::vector<double> B = { 2, 34, 65, 6, 75, 8, 26, 7, 4, 12, 27, 22, 20, 1,
      23, 48 };
  std::vector<double> C = subtraction_of2m(B, A, n);
  std::vector<double> res = { 1, 32, 62, 2, 70, 2, 19, -1, -5, 2, 16, 10, 7,
      -13, 8, 32};
  if (rank == 0) {
    ASSERT_NO_THROW(C = subtraction_of2m(B, A, n););
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        EXPECT_NEAR(C[j + i * n], res[j + i * n], 1e-5);
      }
    }
  }
}

TEST(TEST_STRASSEN_ALG, Test_simple_multiplication) {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size < 2) {
        int n = 4;
        std::vector<double> A = createMatrix(n);
        std::vector<double> B = createMatrix(n);
        std::vector<double> C = createMatrix(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i * n + j] = i + 1;
                B[i * n + j] = i + 1;
                C[i * n + j] = 10 * (i + 1);
            }
        }
        std::vector<double> simple = simple_mult(A, B, n);
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
                ASSERT_TRUE(C[k * n + l] == simple[k * n + l]);
    }
}

TEST(TEST_STRASSEN_ALG, Test_sequential_operations) {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n = 4;
    std::vector<double> A = createMatrix(n);
    std::vector<double> B = createMatrix(n);
    getRandomMatrix(A, n);
    getRandomMatrix(B, n);
    std::vector<double> simple = simple_mult(A, B, n);
    std::vector<double> strassen = getSequentialOperations(A, B, n);
    for (int k = 0; k < n; k++)
      for (int l = 0; l < n; l++)
        ASSERT_TRUE(strassen[k * n + l] == simple[k * n + l]);
}

TEST(TEST_STRASSEN_ALG, Test_parallel_sequential_time) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0) {
    int n = 4;
    std::vector<double> A = createMatrix(n);
    std::vector<double> B = createMatrix(n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i * n + j] = i*15035195.2846+ 1;
        B[i * n + j] = i*1005002 + 1;
      }
    }
    double timeSeque = MPI_Wtime();
    std::vector<double> seq = getSequentialOperations(A, B, n);
    double timeS = MPI_Wtime();
    double timeParal = MPI_Wtime();
    std::vector<double> parallel = getParallelOperations(A, B, n);
    double timeP = MPI_Wtime();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          ASSERT_TRUE(parallel[i * n + j] == seq[i * n + j]);
      }
    }
    printf("The parallel   time: %f seconds\n", timeP-timeParal);
    printf("The sequential time: %f seconds\n", timeS - timeSeque);
  }
}

TEST(TEST_STRASSEN_ALG, Test_parallel_operations) {
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0) {
    int n = 4;
    std::vector<double> A = createMatrix(n);
    std::vector<double> B = createMatrix(n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A[i * n + j] = i * 15035195.2846 + 1;
        B[i * n + j] = i * 1005002 + 1;
      }
    }
    std::vector<double> simple = simple_mult(A, B, n);
    double timeParal = MPI_Wtime();
    std::vector<double> parallel = getParallelOperations(A, B, n);
    double timeP = MPI_Wtime();
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        ASSERT_TRUE(parallel[i * n + j] == simple[i * n + j]);
    printf("The parallel   time: %f seconds\n", timeP - timeParal);
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
