// Copyright 2020 Schekotilova Julia
#ifndef MODULES_TASK_3_SCHEKOTILOVA_JU_STRASSEN_ALG_STRASSEN_ALG_H_
#define MODULES_TASK_3_SCHEKOTILOVA_JU_STRASSEN_ALG_STRASSEN_ALG_H_
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <cstdlib>

std::vector<double> createMatrix(const int n);
std::vector<double> getRandomMatrix(std::vector<double> m, int size);

const std::vector<double> simple_mult(const std::vector<double>& m1,
  const std::vector<double>& m2, int size);
const std::vector<double> addition_of2m(const std::vector<double>& matrix1,
  const std::vector<double>& matrix2, int size);
const std::vector<double> addition_of4m(const std::vector<double>& matrix1,
  const std::vector<double>& matrix2, const std::vector<double>& matrix3,
  const std::vector<double>& matrix4, int size);
const std::vector<double> subtraction_of2m(const std::vector<double>& matrix1,
  const std::vector<double>& matrix2, int size);
const std::vector<double> add_of3_sub_of4(const std::vector<double>& matrix1,
  const std::vector<double>& matrix2, const std::vector<double>& matrix3,
  const std::vector<double>& matrix4, int size);
void getFourMatrixBlocks(const std::vector<double>& A, std::vector<double>* A11,
  std::vector<double>* A12, std::vector<double>* A21, std::vector<double>* A22,
  int n_A);
const std::vector<double> getSequentialOperations(const std::vector<double>& matrix1,
  const std::vector<double>& matrix2, int size);
const std::vector<double> getParallelOperations(const std::vector<double>& matr_A,
  const std::vector<double>& matr_B, int n);

#endif  // MODULES_TASK_3_SCHEKOTILOVA_JU_STRASSEN_ALG_STRASSEN_ALG_H_\
