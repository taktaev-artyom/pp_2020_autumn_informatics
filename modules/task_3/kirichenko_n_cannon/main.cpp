// Copyright 2020 Kirichenko Nikita

#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include "./Cannon.h"

#define EPS 1e-5

TEST(Cannon, DISABLED_Time) {
    int procRank, procCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);

    const size_t blockSize = 200;

    int size = procCount * blockSize;
    double dimProcess = std::sqrt(procCount);
    int intDimProcess = static_cast<int>(dimProcess);
    if (std::fabs(dimProcess - intDimProcess) > EPS)
        return;

    std::vector<double> a, b, c;

    double starttime, endtime;
    if (procRank == 0) {
        a = GetRandomMatrix(size);
        b = GetRandomMatrix(size);

        starttime = MPI_Wtime();
        std::vector<double> ref = MatrixMultiplication(a, b, size);
        endtime = MPI_Wtime();
        std::cout << "Sequential time: " << endtime - starttime << std::endl;

        a = GetRandomMatrix(size);
        b = GetRandomMatrix(size);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    starttime = MPI_Wtime();
    c = Cannon(a, b, size);
    endtime = MPI_Wtime();
    if (procRank == 0)
        std::cout << "Parallel time: " << endtime - starttime << std::endl;
}

TEST(Cannon, IncorrectSize) {
    std::vector<double> a(9, 1);
    std::vector<double> b(16, 1);
    ASSERT_ANY_THROW(Cannon(a, b, 3));
}

TEST(Cannon, SequantalMultuplication) {
    int procRank, procCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);

    if (procRank == 0) {
        const size_t size = 2;
        const size_t countElem = size * size;
        std::vector<double> a = { 0, 1, 2, 3 };
        std::vector<double> b = { 3, 2, 1, 0 };
        std::vector<double> ref = { 1, 0, 9, 4 };
        std::vector<double> c(countElem);
        c = MatrixMultiplication(a, b, size);

        for (size_t i = 0; i < countElem; ++i)
            ASSERT_NEAR(ref[i], c[i], EPS);
    }
}

TEST(Cannon, SmallSize) {
    int procRank, procCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);

    const size_t blockSize = 4;

    int size = procCount * blockSize;
    double dimProcess = std::sqrt(procCount);
    int intDimProcess = static_cast<int>(dimProcess);
    if (std::fabs(dimProcess - intDimProcess) > EPS)
        return;

    std::vector<double> a, b, c;

    if (procRank == 0) {
        a = GetRandomMatrix(size);
        b = GetRandomMatrix(size);
    }

    c = Cannon(a, b, size);

    if (procRank == 0) {
        std::vector<double> ref = MatrixMultiplication(a, b, size);

        const size_t countElem = size * size;
        for (size_t i = 0; i < countElem; ++i)
            ASSERT_NEAR(ref[i], c[i], EPS);
    }
}

TEST(Cannon, MediumSize) {
    int procRank, procCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);

    const size_t blockSize = 11;

    int size = procCount * blockSize;
    double dimProcess = std::sqrt(procCount);
    int intDimProcess = static_cast<int>(dimProcess);
    if (std::fabs(dimProcess - intDimProcess) > EPS)
        return;

    std::vector<double> a, b, c;

    if (procRank == 0) {
        a = GetRandomMatrix(size);
        b = GetRandomMatrix(size);
    }

    c = Cannon(a, b, size);

    if (procRank == 0) {
        std::vector<double> ref = MatrixMultiplication(a, b, size);

        const size_t countElem = size * size;
        for (size_t i = 0; i < countElem; ++i)
            ASSERT_NEAR(ref[i], c[i], EPS);
    }
}

TEST(Cannon, LargeSize) {
    int procRank, procCount;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);

    const size_t blockSize = 23;

    int size = procCount * blockSize;
    double dimProcess = std::sqrt(procCount);
    int intDimProcess = static_cast<int>(dimProcess);
    if (std::fabs(dimProcess - intDimProcess) > EPS)
        return;

    std::vector<double> a, b, c;

    if (procRank == 0) {
        a = GetRandomMatrix(size);
        b = GetRandomMatrix(size);
    }

    c = Cannon(a, b, size);

    if (procRank == 0) {
        std::vector<double> ref = MatrixMultiplication(a, b, size);

        const size_t countElem = size * size;
        for (size_t i = 0; i < countElem; ++i)
            ASSERT_NEAR(ref[i], c[i], EPS);
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
