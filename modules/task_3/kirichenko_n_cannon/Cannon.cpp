// Copyright 2020 Kirichenko Nikita
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <cmath>
#include <exception>
#include "../../../modules/task_3/kirichenko_n_cannon/Cannon.h"
#define EPS 1e-5

static unsigned int offset = 0;

std::vector<double> GetRandomMatrix(const size_t size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)) + offset);
    offset += 10;

    const size_t count = size * size;
    std::vector<double> vec(count);
    for (size_t i = 0; i < count; ++i)
        vec[i] = gen() % 20;
    return vec;
}

std::vector<double> MatrixMultiplication(const std::vector<double>& a,
    const std::vector<double>& b, const size_t size) {
    if (a.size() != b.size())
        throw "[ERROR] Matrices A and B have got different size!";

    std::vector<double> c(size * size, 0);
    for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < size; ++j)
            for (size_t k = 0; k < size; ++k)
                c[i * size + j] += a[i * size + k] * b[k * size + j];

    return c;
}

std::vector<double> Cannon(const std::vector<double>& a, const std::vector<double>& b, const size_t size) {
    int procRank, procCount;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    if (procCount == 1)
        return MatrixMultiplication(a, b, size);

    int sizeA = 0;
    int sizeB = 0;
    int globalSize;
    if (procRank == 0) {
        sizeA = static_cast<int>(a.size());
        sizeB = static_cast<int>(b.size());
        globalSize = static_cast<int>(size);
    }
    MPI_Bcast(&sizeA, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sizeB, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&globalSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (sizeA != sizeB) {
        throw "[ERROR] Matrices A and B have got different size!";
    }

    double sqrtProcCount = std::sqrt(procCount);
    if ((sqrtProcCount - std::floor(sqrtProcCount)) >= EPS) {
        throw "[ERROR] Number of processes must be a perfect square!";
    }

    if (globalSize % static_cast<int>(sqrtProcCount) != 0) {
        throw "[ERROR] Incorrect matrices size for Cannon algorithm!";
    }

    std::vector<double> c(globalSize * globalSize);

    int dimProcess = static_cast<int>(sqrtProcCount);
    int dimBlock = globalSize / dimProcess;

    int dims[] = { dimProcess, dimProcess };
    int periods[] = { 1, 1 };
    int reorder = 1;
    MPI_Comm cartComm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartComm);

    int globalSizes[2] = { static_cast<int>(globalSize), static_cast<int>(globalSize) };
    int localSize[2] = { dimBlock, dimBlock };
    int starts[2] = { 0, 0 };
    MPI_Datatype type, subarrtype;
    MPI_Type_create_subarray(2, globalSizes, localSize, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
    MPI_Type_create_resized(type, 0, dimBlock * sizeof(a[0]), &subarrtype);
    MPI_Type_commit(&subarrtype);

    const int blockElemCount = dimBlock * dimBlock;
    std::vector<double> localA(blockElemCount);
    std::vector<double> localB(blockElemCount);
    std::vector<double> localC(blockElemCount, 0);

    int* sendCounts = new int[procCount];
    int* displs = new int[procCount];
    if (procRank == 0) {
        std::fill_n(sendCounts, procCount, 1);

        int disp = 0;
        for (int i = 0; i < dimProcess; i++) {
            for (int j = 0; j < dimProcess; j++) {
                displs[i * dimProcess + j] = disp;
                disp++;
            }
            disp += (dimBlock - 1) * dimProcess;
        }
    }

    MPI_Scatterv(a.data(), sendCounts, displs, subarrtype, &localA[0], blockElemCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b.data(), sendCounts, displs, subarrtype, &localB[0], blockElemCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int coord[2];
    int left, right, up, down;
    MPI_Cart_coords(cartComm, procRank, 2, coord);
    MPI_Cart_shift(cartComm, 1, coord[0], &left, &right);
    MPI_Sendrecv_replace(&localA[0], blockElemCount, MPI_DOUBLE, left, 1, right, 1, cartComm, &status);
    MPI_Cart_shift(cartComm, 0, coord[1], &up, &down);
    MPI_Sendrecv_replace(&localB[0], blockElemCount, MPI_DOUBLE, up, 1, down, 1, cartComm, &status);

    std::vector<double> tmpC(blockElemCount);
    for (int k = 0; k < dimProcess; k++) {
        tmpC = MatrixMultiplication(localA, localB, dimBlock);

        for (int i = 0; i < dimBlock; i++) {
            for (int j = 0; j < dimBlock; j++) {
                localC[i * dimBlock + j] += tmpC[i * dimBlock + j];
            }
        }

        MPI_Cart_shift(cartComm, 1, 1, &left, &right);
        MPI_Cart_shift(cartComm, 0, 1, &up, &down);
        MPI_Sendrecv_replace(&localA[0], blockElemCount, MPI_DOUBLE, left, 1, right, 1, cartComm, &status);
        MPI_Sendrecv_replace(&localB[0], blockElemCount, MPI_DOUBLE, up, 1, down, 1, cartComm, &status);
    }

    MPI_Gatherv(&localC[0], blockElemCount, MPI_DOUBLE, &c[0], sendCounts, displs, subarrtype, 0, MPI_COMM_WORLD);

    return c;
}
