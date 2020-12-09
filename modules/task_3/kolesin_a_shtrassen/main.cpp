// Copyright 2020 Kolesin Andrey
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <stdio.h>
#include <random>
#include <vector>
#include <algorithm>
#include <deque>
// #include <chrono>
#include "./shtrassen.h"

TEST(Count_Words, Test_SeqVsSimple) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int n = 8;
        std::vector<int> m1 = getRandomMatrix(n);
        std::vector<int> m2 = getRandomMatrix(n);
        std::vector<int> seq = m1;
        std::vector<int> simp = m1;
        Matrix M1(&m1[0], n);
        Matrix M2(&m2[0], n);
        Matrix Seq(&seq[0], n);
        Matrix Simp(&simp[0], n);
        ShtSeq(M1, M2, Seq);
        SimpleMult(M1, M2, Simp);
        // Seq.print();
        // Simp.print();
        EXPECT_EQ(seq, simp);
    }
}
TEST(Count_Words, Test_ParallVsSeq_16) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 16;
    std::vector<int> m1(n * n);
    std::vector<int> m2(n * n);
    if (rank == 0) {
        m1 = getRandomMatrix(n);
        m2 = getRandomMatrix(n);
    }
    std::vector<int> p(n * n);
    std::vector<int> s(n * n);
    Matrix P(&p[0], n);
    Matrix S(&s[0], n);
    Matrix M1(&m1[0], n);
    Matrix M2(&m2[0], n);
// int ms1 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
    SafeSht(M1, M2, P, MPI_COMM_WORLD);
    if (rank == 0) {
// int ms2 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
        ShtSeq(M1, M2, S);
// int ms3 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
// std::cout<<ms3-ms2<<"   "<<ms2-ms1<<std::endl;
        EXPECT_EQ(p, s);
    }
}
TEST(Count_Words, Test_ParallVsSeq_8) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 8;
    std::vector<int> m1(n * n);
    std::vector<int> m2(n * n);
    if (rank == 0) {
        m1 = getRandomMatrix(n);
        m2 = getRandomMatrix(n);
    }
    std::vector<int> p(n * n);
    std::vector<int> s(n * n);
    Matrix P(&p[0], n);
    Matrix S(&s[0], n);
    Matrix M1(&m1[0], n);
    Matrix M2(&m2[0], n);
// int ms1 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
    SafeSht(M1, M2, P, MPI_COMM_WORLD);
    if (rank == 0) {
// int ms2 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
        ShtSeq(M1, M2, S);
// int ms3 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
// std::cout<<ms3-ms2<<"   "<<ms2-ms1<<std::endl;
        EXPECT_EQ(p, s);
    }
}
TEST(Count_Words, Test_ParallVsSeq_4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 4;
    std::vector<int> m1(n * n);
    std::vector<int> m2(n * n);
    if (rank == 0) {
        m1 = getRandomMatrix(n);
        m2 = getRandomMatrix(n);
    }
    std::vector<int> p(n * n);
    std::vector<int> s(n * n);
    Matrix P(&p[0], n);
    Matrix S(&s[0], n);
    Matrix M1(&m1[0], n);
    Matrix M2(&m2[0], n);
// int ms1 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
    SafeSht(M1, M2, P, MPI_COMM_WORLD);
    if (rank == 0) {
// int ms2 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
        ShtSeq(M1, M2, S);
// int ms3 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
// std::cout<<ms3-ms2<<"   "<<ms2-ms1<<std::endl;
        EXPECT_EQ(p, s);
    }
}
TEST(Count_Words, Test_ParallVsSeq_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 1;
    std::vector<int> m1(n * n);
    std::vector<int> m2(n * n);
    if (rank == 0) {
        m1 = getRandomMatrix(n);
        m2 = getRandomMatrix(n);
    }
    std::vector<int> p(n * n);
    std::vector<int> s(n * n);
    Matrix P(&p[0], n);
    Matrix S(&s[0], n);
    Matrix M1(&m1[0], n);
    Matrix M2(&m2[0], n);
// int ms1 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
    SafeSht(M1, M2, P, MPI_COMM_WORLD);
    if (rank == 0) {
// int ms2 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
        ShtSeq(M1, M2, S);
// int ms3 = std::chrono::duration_cast<std::chrono::milliseconds>
// (std::chrono::system_clock::now().time_since_epoch()).count();
// std::cout<<ms3-ms2<<"   "<<ms2-ms1<<std::endl;
        EXPECT_EQ(p, s);
    }
}
int main(int argc, char **argv) {
    setlocale(LC_ALL, "Russian");
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners &listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
