// Copyright 2020 Ludina Daria
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include "./ring_topology.h"

TEST(Topology_Ring, Create_Ring_Topology) {
  MPI_Comm ringcomm = createRingcomm(MPI_COMM_WORLD);

  ASSERT_TRUE(ringcomm != MPI_COMM_NULL);
}

TEST(Topology_Ring, Test_Cart) {
  int status;
  MPI_Comm ringcomm = createRingcomm(MPI_COMM_WORLD);
  MPI_Topo_test(ringcomm, &status);

  ASSERT_TRUE(status == MPI_CART);
}

TEST(Topology_Ring, Test_Dims) {
  int ndims;
  MPI_Comm ringcomm = createRingcomm(MPI_COMM_WORLD);
  MPI_Cartdim_get(ringcomm, &ndims);

  ASSERT_EQ(ndims, 1);
}

TEST(Topology_Ring, Test_Periods) {
  int ndims;
  int dims[1], periods[1], coords[1];
  MPI_Comm ringcomm = createRingcomm(MPI_COMM_WORLD);
  MPI_Cartdim_get(ringcomm, &ndims);
  MPI_Cart_get(ringcomm, ndims, dims, periods, coords);

  ASSERT_EQ(periods[0], 1);
}

TEST(Topology_Ring, Test_Send_Array) {
  int rank, size;
  MPI_Comm ringcomm = createRingcomm(MPI_COMM_WORLD);
  MPI_Comm_size(ringcomm, &size);
  MPI_Comm_rank(ringcomm, &rank);

  int size_ar = 4;
  int rank_source = 0;
  int rank_dest = size - 1;

  int* arr = new int[size_ar];
  if (rank == rank_source) {
    for (int i = 0; i < size_ar; i++)
      arr[i] = 2;
  }

  int* arr2 = Send(&arr[0], size_ar, rank_source, rank_dest, ringcomm);
  if (rank == rank_dest) {
    int* res = new int[size_ar] { 2, 2, 2, 2 };
    for (int i = 0; i < size_ar; i++) {
      ASSERT_EQ(res[i], arr2[i]);
    }
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
