// Copyright 2020 Luckyanchenko Ivan
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <random>
#include <vector>
#include "./integral_simpson.h"

TEST(Integral, Test_func1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_time =  MPI_Wtime();
    double res = get_Paral_Integral(func1, 1, 2, 3, 4, 5, 6, 200);
    double end_time = MPI_Wtime();
    double time_paral = end_time - start_time;
    double start_time1 =  MPI_Wtime();;
    double res1 = get_Integral(func1, 1, 2, 3, 4, 5, 6, 200);
    double end_time1 = MPI_Wtime();
    double time_notParal = end_time1 - start_time1;
    double true_answer = 10.5;
    if (rank == 0) {
        ASSERT_LE(std::abs(res - true_answer), 0.1);
        std::cout <<"parall integral :   " << res << std::endl;
        std::cout <<"not parall integral :   " << res1 << std::endl;
        std::cout << std::fixed << std::showpoint;
        std::cout << std::setprecision(10);
        std::cout << "Time Paral = " << time_paral <<
        std::endl << "Time not Paral  = " << time_notParal << std::endl;
    }
}
TEST(Integral, Test_func1_mega_count) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_time =  MPI_Wtime();
    double res = get_Paral_Integral(func1, 1, 2, 3, 4, 5, 6, 1000);
    double end_time = MPI_Wtime();
    double time_paral = end_time - start_time;
    double start_time1 =  MPI_Wtime();;
    double res1 = get_Integral(func1, 1, 2, 3, 4, 5, 6, 1000);
    double end_time1 = MPI_Wtime();
    double time_notParal = end_time1 - start_time1;
    double true_answer = 10.5;
    if (rank == 0) {
        ASSERT_LE(std::abs(res - true_answer), 0.01);
        std::cout <<"parall integral :   " << res << std::endl;
        std::cout <<"not parall integral :   " << res1 << std::endl;
        std::cout << std::fixed << std::showpoint;
        std::cout << std::setprecision(10);
        std::cout << "Time Paral = " << time_paral <<
        std::endl << "Time not Paral  = " << time_notParal << std::endl;
    }
}
TEST(Integral, Test_func2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_time =  MPI_Wtime();
    double res = get_Paral_Integral(func2, 1, 5, 3, 12, 5, 7, 300);
    double end_time = MPI_Wtime();
    double time_paral = end_time - start_time;
    double start_time1 =  MPI_Wtime();;
    double res1 = get_Integral(func2, 1, 5, 3, 12, 5, 7, 300);
    double end_time1 = MPI_Wtime();
    double time_notParal = end_time1 - start_time1;
    double true_answer = 4.48836967855013;
    if (rank == 0) {
        ASSERT_LE(std::abs(res - true_answer), 0.1);
        std::cout <<"parall integral :   " << res << std::endl;
        std::cout <<"not parall integral :   " << res1 << std::endl;
        std::cout << std::fixed << std::showpoint;
        std::cout << std::setprecision(10);
        std::cout << "Time Paral = " << time_paral <<
        std::endl << "Time not Paral  = " << time_notParal << std::endl;
    }
}
TEST(Integral, Test_func2_mega_count) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_time =  MPI_Wtime();
    double res = get_Paral_Integral(func2, 1, 5, 3, 12, 5, 7, 1000);
    double end_time = MPI_Wtime();
    double time_paral = end_time - start_time;
    double start_time1 =  MPI_Wtime();;
    double res1 = get_Integral(func2, 1, 5, 3, 12, 5, 7, 1000);
    double end_time1 = MPI_Wtime();
    double time_notParal = end_time1 - start_time1;
    double true_answer = 4.48836967855013;
    if (rank == 0) {
        ASSERT_LE(std::abs(res - true_answer), 0.01);
        std::cout <<"parall integral :   " << res << std::endl;
        std::cout <<"not parall integral :   " << res1 << std::endl;
        std::cout << std::fixed << std::showpoint;
        std::cout << std::setprecision(10);
        std::cout << "Time Paral = " << time_paral <<
        std::endl << "Time not Paral  = " << time_notParal << std::endl;
    }
}
TEST(Integral, Test_func3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_time =  MPI_Wtime();
    double res = get_Paral_Integral(func3, 2, 5, 12, 15, 6, 11, 100);
    double end_time = MPI_Wtime();
    double time_paral = end_time - start_time;
    double start_time1 =  MPI_Wtime();;
    double res1 = get_Integral(func3, 2, 5, 12, 15, 6, 11, 100);
    double end_time1 = MPI_Wtime();
    double time_notParal = end_time1 - start_time1;
    double true_answer = -1.28134502288645;
    if (rank == 0) {
        ASSERT_LE(std::abs(res - true_answer), 0.01);
        std::cout <<"parall integral :   " << res << std::endl;
        std::cout <<"not parall integral :   " << res1 << std::endl;
        std::cout << std::fixed << std::showpoint;
        std::cout << std::setprecision(10);
        std::cout << "Time Paral = " << time_paral <<
        std::endl << "Time not Paral  = " << time_notParal << std::endl;
    }
}
TEST(Integral, Test_func3_mega_count) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_time =  MPI_Wtime();
    double res = get_Paral_Integral(func3, 2, 5, 12, 15, 6, 11, 1000);
    double end_time = MPI_Wtime();
    double time_paral = end_time - start_time;
    double start_time1 =  MPI_Wtime();;
    double res1 = get_Integral(func3, 2, 5, 12, 15, 6, 11, 1000);
    double end_time1 = MPI_Wtime();
    double time_notParal = end_time1 - start_time1;
    double true_answer = -1.28134502288645;
    if (rank == 0) {
        ASSERT_LE(std::abs(res - true_answer), 0.01);
        std::cout <<"parall integral :   " << res << std::endl;
        std::cout <<"not parall integral :   " << res1 << std::endl;
        std::cout << std::fixed << std::showpoint;
        std::cout << std::setprecision(10);
        std::cout << "Time Paral = " << time_paral <<
        std::endl << "Time not Paral  = " << time_notParal << std::endl;
    }
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
