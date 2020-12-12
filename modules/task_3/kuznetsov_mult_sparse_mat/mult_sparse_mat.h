// Copyright 2020 Kuznetsov Nikita
#ifndef MODULES_TASK_3_KUZNETSOV_MULT_SPARSE_MAT_MULT_SPARSE_MAT_H_
#define MODULES_TASK_3_KUZNETSOV_MULT_SPARSE_MAT_MULT_SPARSE_MAT_H_

#include <vector>

std::vector<double> randMat(const int rows, const int cols);

struct sparseMatrix {
  std::vector<double> val;
  std::vector<int> c_rows;
  std::vector<int> c_index;
  int cols, rows, nnz;
  friend const std::vector<double> operator*(const sparseMatrix& A, const sparseMatrix& B);
};

sparseMatrix CCS(const std::vector<double> new_mat, const int _cols, const int _rows);
std::vector<double> matMultiply(sparseMatrix A, sparseMatrix B);

#endif  // MODULES_TASK_3_KUZNETSOV_MULT_SPARSE_MAT_MULT_SPARSE_MAT_H_
