// Copyright 2020 Taktaev Artem
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <cmath>
#include "../../../modules/task_3/taktaev_a_strongin_method/strongin_method.h"

TEST(Strongin_Method_MPI, Test_Wrong_Parameters_In_Constructor) {
    ASSERT_ANY_THROW(Strongin X(2, -2, 2, -4, [](double x)->double { return std::sin(x); }));
}

TEST(Strongin_Method_MPI, Test_Correct_Parameters_In_Constructor) {
    ASSERT_NO_THROW(Strongin X(-1, 1, 2, -4, [](double x)->double { return std::sin(x); }));
}

TEST(Strongin_Method_MPI, Test_Wrong_X_In_add_Node) {
    Strongin X(-1, 1, 2, -4, [](double x)->double { return std::sin(x); });
    ASSERT_ANY_THROW(X.addNode(5));
}

TEST(Strongin_Method_MPI, Test_Correct_X_In_add_Node) {
    Strongin X(-1, 1, 2, -4, [](double x)->double { return std::sin(x); });
    ASSERT_NO_THROW(X.addNode(0.5));
}

TEST(Strongin_Method_MPI, Test_Correct_Work_Of_seq_Strongin_Method) {
    Strongin Y(-1, 1, 2, -10, [](double x)->double { return std::sin(x); });
    ASSERT_NEAR(Y.seqStronginSearch(), -1, 1.0e-10);
}

TEST(Strongin_Method_MPI, Test_Correct_Work_Of_par_Strongin_Method) {
    Strongin Z(-1, 1, 2, -10, [](double x)->double { return std::sin(x); });
    ASSERT_NEAR(Z.parStronginSearch(), -1, 1.0e-10);
}

TEST(Strongin_Method_MPI, Test_Eff) {
    Strongin A(-1, 1, 2, -10, [](double x)->double { return std::sin(x); });
    Strongin B(-1, 1, 2, -10, [](double x)->double { return std::sin(x); });
    int proc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    double start, end;
    if (proc_rank == 0) {
        start = MPI_Wtime();
        A.seqStronginSearch();
        end = MPI_Wtime();
        std::cout << "Time sequential = " << end - start << std::endl;
    }
    start = MPI_Wtime();
    B.parStronginSearch();
    end = MPI_Wtime();
    if (proc_rank == 0) {
        std::cout << "Time parallel = " << end - start << std::endl;
    }
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
