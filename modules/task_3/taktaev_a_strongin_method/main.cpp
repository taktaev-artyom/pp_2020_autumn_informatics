// Copyright 2020 Taktaev Artem
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include "../../../modules/task_3/taktaev_a_strongin_method/strongin_method.h"

TEST(Strongin_Method_MPI, Test_Wrong_Parameters_In_Constructor) {
    ASSERT_ANY_THROW(Strongin X(2, -2, 2, -4, [](double x)->double { return sin(x); }));
}

TEST(Strongin_Method_MPI, Test_Correct_Parameters_In_Constructor) {
    ASSERT_NO_THROW(Strongin X(-1, 1, 2, -4, [](double x)->double { return sin(x); }));
}

TEST(Strongin_Method_MPI, Test_Wrong_X_In_add_Node) {
    Strongin X(-1, 1, 2, -4, [](double x)->double { return sin(x); });
    ASSERT_ANY_THROW(X.addNode(5));
}

TEST(Strongin_Method_MPI, Test_Correct_X_In_add_Node) {
    Strongin X(-1, 1, 2, -4, [](double x)->double { return sin(x); });
    ASSERT_NO_THROW(X.addNode(0.5));
}

TEST(Strongin_Method_MPI, Test_Correct_Work_Of_seq_Strongin_Method) {
    Strongin Y(-1, 1, 2, -10, [](double x)->double { return sin(x); });
    ASSERT_NEAR(Y.seqStronginSearch(), -1, 1.0e-10);
}

TEST(Strongin_Method_MPI, Test_Correct_Work_Of_par_Strongin_Method) {
    Strongin Z(-1, 1, 2, -10, [](double x)->double { return sin(x); });
    ASSERT_NEAR(Z.parStronginSearch(), -1, 1.0e-10);
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
