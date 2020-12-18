// Copyright 2020 Zhuravlev Roman
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include "../../../modules/task_3/zhuravlev_r_simpson/Simpson.h"

TEST(Parallel_MPI, Test_Sequantial_Res) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    EXPECT_DOUBLE_EQ(Sequential_Simpson_Double(func, 0, 3, -2, 4, 10, 8), 10.8) << "Wrong res";
  }
}

TEST(Parallel_MPI, Test_Parallel_Res) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  double res = Parallel_Simpson_Double(func, 1, 14, 2, 5, 10, 10);

  if (rank == 0) {
    EXPECT_DOUBLE_EQ(res, 409.5) << "Wrong res";
  }
}

TEST(Parallel_MPI, Test_Exception) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    EXPECT_ANY_THROW(Sequential_Simpson_Double(func, 0, 3, 5, 4, 10, 8));
  }
}

TEST(Parallel_MPI, Test_Parallel_Res_equals_Sequential_Res_small) {
  int rank;
  double t1, t2, t3, t4;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  t1 = MPI_Wtime();
  double parRes = Parallel_Simpson_Double(func, 0, 3, -2, 4, 30, 20);
  t2 = MPI_Wtime();
  t3 = MPI_Wtime();
  double seqRes = Sequential_Simpson_Double(func, 0, 3, -2, 4, 30, 20);
  t4 = MPI_Wtime();
  if (rank == 0) {
    std::cout << "parallel: " << (t2 - t1) << ", seqential: " << (t4 - t3) << std::endl;
    EXPECT_DOUBLE_EQ(seqRes, parRes) << "Wrong res";
  }
}

TEST(Parallel_MPI, Test_Parallel_Res_equals_Sequential_Res_big) {
  int rank;
  double t1, t2, t3, t4;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  t1 = MPI_Wtime();
  double parRes = Parallel_Simpson_Double(func, -10, 31, 1, 15, 100, 80);
  t2 = MPI_Wtime();
  t3 = MPI_Wtime();
  double seqRes = Sequential_Simpson_Double(func, -10, 31, 1, 15, 100, 80);
  t4 = MPI_Wtime();
  if (rank == 0) {
    std::cout << "parallel: " << (t2 - t1) << ", seqential: " << (t4 - t3) << std::endl;
    EXPECT_DOUBLE_EQ(seqRes, parRes) << "Wrong res";
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
