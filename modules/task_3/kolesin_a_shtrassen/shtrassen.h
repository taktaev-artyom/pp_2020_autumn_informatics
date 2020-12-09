// Copyright 2020 Kolesin Andrey
#ifndef MODULES_TASK_3_KOLESIN_A_SHTRASSEN_SHTRASSEN_H_
#define MODULES_TASK_3_KOLESIN_A_SHTRASSEN_SHTRASSEN_H_

#include <mpi.h>
#include <random>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <vector>
#include <deque>
class Matrix {
    int n;
    bool tmp;
    std::deque<int> coords;

 public:
    double *buff;
    int N;
    Matrix(double *_buff, int _n, bool _tmp = false);
    Matrix(Matrix old, int x, int y);
    Matrix(const Matrix &old);
    ~Matrix();
    void print();
    double &getElem(int x, int y);
    void stretch(double *buff);
    void Send(MPI_Comm comm, int sender_rank, int reciver_rank, int marker);
};
std::vector<double> getRandomMatrix(int size = -1);
void ShtSeq(Matrix A, Matrix B, Matrix C);
void SimpleMult(Matrix A, Matrix B, Matrix C);
void Sum(Matrix A, Matrix B, Matrix C);
void Sht(Matrix A, Matrix B, Matrix C, MPI_Comm comm, int level = 0);
void SafeSht(Matrix A, Matrix B, Matrix C, MPI_Comm comm);

#endif  // MODULES_TASK_3_KOLESIN_A_SHTRASSEN_SHTRASSEN_H_
