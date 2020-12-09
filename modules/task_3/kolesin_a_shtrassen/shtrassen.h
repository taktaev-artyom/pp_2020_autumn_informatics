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
    int *buff;
    int N;
    Matrix(int *_buff, int _n, bool _tmp = false);
    Matrix(Matrix old, int x, int y);
    Matrix(const Matrix &old);
    ~Matrix();
    void print();
    int &getElem(int x, int y);
    void stretch(int *buff);
    void Send(MPI_Comm comm, int sender_rank, int reciver_rank, int marker);
};
std::vector<int> getRandomMatrix(int size = -1);
int *Mult(int *M1, int *M2, int n);
void ShtSeq(Matrix A, Matrix B, Matrix C);
void SimpleMult(Matrix A, Matrix B, Matrix C);
void Sum(Matrix A, Matrix B, Matrix C);
void Sht(Matrix A, Matrix B, Matrix C, MPI_Comm comm, int level = 0);
void SafeSht(Matrix A, Matrix B, Matrix C, MPI_Comm comm);

#endif  // MODULES_TASK_3_KOLESIN_A_SHTRASSEN_SHTRASSEN_H_
