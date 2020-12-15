// Copyright 2020 Kirillov Konstantin
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include<iostream>
#include "./contrast_enhancement.h"
TEST(Parallel_Operations_MPI, Test_Matrix_7x7) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 7;
    int cols = 7;
    double alpha = 1.5;
    int beta = 1;
    Matrix global_mat(rows, std::vector<int>(cols));
    Matrix local_mat(rows, std::vector<int>(cols));
    global_mat = getRandomMatrix(rows, cols);
    if (procRank == 0) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                local_mat[i][j] = global_mat[i][j];
                local_mat[i][j] = getSequentialContrast(local_mat, i, j, alpha, beta);
            }
        }
    }
    global_mat = getParallelContrast(global_mat, rows, cols, alpha, beta);
    if (procRank == 0) {
        for (int  i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                ASSERT_EQ(local_mat[i][j], global_mat[i][j]);
            }
        }
    }
}
TEST(Parallel_Operations_MPI, Test_Matrix_4x4) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 4;
    int cols = 4;
    double alpha = 1.5;
    int beta = 1;
    Matrix global_mat(rows, std::vector<int>(cols));
    Matrix local_mat(rows, std::vector<int>(cols));
    global_mat = getRandomMatrix(rows, cols);
    if (procRank == 0) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                local_mat[i][j] = global_mat[i][j];
                local_mat[i][j] = getSequentialContrast(local_mat, i, j, alpha, beta);
            }
        }
    }
    global_mat = getParallelContrast(global_mat, rows, cols, alpha, beta);
    if (procRank == 0) {
        for (int  i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                ASSERT_EQ(local_mat[i][j], global_mat[i][j]);
            }
        }
    }
}
TEST(Parallel_Operations_MPI, Test_Matrix_5x5) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 100;
    int cols = 97;
    double alpha = 1.5;
    int beta = 1;
    Matrix global_mat(rows, std::vector<int>(cols));
    Matrix local_mat(rows, std::vector<int>(cols));
    global_mat = getRandomMatrix(rows, cols);
    if (procRank == 0) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                local_mat[i][j] = global_mat[i][j];
                local_mat[i][j] = getSequentialContrast(local_mat, i, j, alpha, beta);
            }
        }
    }
    global_mat = getParallelContrast(global_mat, rows, cols, alpha, beta);
    if (procRank == 0) {
        for (int  i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                ASSERT_EQ(local_mat[i][j], global_mat[i][j]);
            }
        }
    }
}
TEST(Parallel_Operations_MPI, Test_Matrix_9x9) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 9;
    int cols = 9;
    double alpha = 1.5;
    int beta = 1;
    Matrix global_mat(rows, std::vector<int>(cols));
    Matrix local_mat(rows, std::vector<int>(cols));
    global_mat = getRandomMatrix(rows, cols);
    if (procRank == 0) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                local_mat[i][j] = global_mat[i][j];
                local_mat[i][j] = getSequentialContrast(local_mat, i, j, alpha, beta);
            }
        }
    }
    global_mat = getParallelContrast(global_mat, rows, cols, alpha, beta);
    if (procRank == 0) {
        for (int  i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                ASSERT_EQ(local_mat[i][j], global_mat[i][j]);
            }
        }
    }
}
TEST(Parallel_Operations_MPI, Test_Matrix_6x6) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 6;
    int cols = 6;
    double alpha = 1.5;
    int beta = 1;
    Matrix global_mat(rows, std::vector<int>(cols));
    Matrix local_mat(rows, std::vector<int>(cols));
    global_mat = getRandomMatrix(rows, cols);
    if (procRank == 0) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                local_mat[i][j] = global_mat[i][j];
                local_mat[i][j] = getSequentialContrast(local_mat, i, j, alpha, beta);
            }
        }
    }
    global_mat = getParallelContrast(global_mat, rows, cols, alpha, beta);
    if (procRank == 0) {
        for (int  i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                ASSERT_EQ(local_mat[i][j], global_mat[i][j]);
            }
        }
    }
}

int main(int argc, char *argv[]) {
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
