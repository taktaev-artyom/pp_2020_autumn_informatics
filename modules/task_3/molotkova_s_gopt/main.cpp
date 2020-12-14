// Copyright 2020 Molotkova Svetlana
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include "./gopt.h"


double f(double* x) {
  return  3 * *x * *x;
}
double testF(double* _x) {
  double x = *_x;
  return sin(9 * x + 4) * cos(15 * x - 2) + 1.4;
}
double f1(double* x) {
  return exp(-0.5 * *x) * sin(6 * *x - 1.5);
}
double f2(double* _x) {
  double x = *_x;
  return pow(x*x + 5, 1.2);
}

double f3(double* _x) {
  double x = *_x;
  return exp(-0.5 * x);
}

TEST(global_optimization, ts1) {
  const int count_It = 1000;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  StronginMethod ts(5, 10, f1, 1e-5);
  double result = ts.Find_Parallel(count_It);
  if (myrank == 0) {
    double sresult = ts.Find_Sequential(count_It);
    ASSERT_NEAR(result, sresult, 1e-2);
  }
}

TEST(global_optimization, ts2) {
  const int count_It = 1000;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  StronginMethod ts(0.6, 2.2, testF, 1e-5);
  double result = ts.Find_Parallel(count_It);
  if (myrank == 0) {
    double sresult = ts.Find_Sequential(count_It); 
    ASSERT_NEAR(result, sresult, 1e-2);
  }
}

TEST(global_optimization, ts3) {
  const int count_It = 1000;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  StronginMethod ts(-6, 10, f, 1e-5);
  double result = ts.Find_Parallel(count_It);
  if (myrank == 0) {
    double sresult = ts.Find_Sequential(count_It); 
    ASSERT_NEAR(result, sresult, 1e-2);
  }
}

TEST(global_optimization, ts4) {
  const int count_It = 1000;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  StronginMethod ts(5, 10, f2, 1e-5);
  double result = ts.Find_Parallel(count_It);
  if (myrank == 0) {
    double sresult = ts.Find_Sequential(count_It); 
    ASSERT_NEAR(result, sresult, 1e-2);
  }
}

TEST(global_optimization, ts5) {
  const int count_It = 1000;
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  StronginMethod ts(5, 10, f3, 1e-5);
  double result = ts.Find_Parallel(count_It);
  if (myrank == 0) {
    double sresult = ts.Find_Sequential(count_It); 
    ASSERT_NEAR(result, sresult, 1e-2);
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
