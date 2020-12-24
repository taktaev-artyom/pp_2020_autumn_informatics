// Copyright 2020 Schekotilova Julia
#include "../../modules/task_3/schekotilova_ju_strassen_alg/strassen_alg.h"
#include <vector>

std::vector<double> createMatrix(const int n) {
  std::vector<double> m(n * n);
  return m;
}

std::vector<double> getRandomMatrix(std::vector<double> m, int size) {
  std::mt19937 gen;
  gen.seed(static_cast<unsigned int>(time(0)));
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
        m[i * size + j] = gen() % 100;
    }
  }
  return m;
}

const std::vector<double> simple_mult(const std::vector<double>& m1,
  const std::vector<double>& m2, int size) {
  std::vector<double> res = createMatrix(size);
  int sum;
  for (int i = 0; i < size; i++)
    for (int j = 0; j < size; j++) {
      sum = 0;
      for (int k = 0; k < size; k++)
        sum += m1[i * size + k] * m2[k * size + j];
      res[i * size + j] = sum;
    }
  return res;
}

const std::vector<double> addition_of2m(const std::vector<double>& m1,
    const std::vector<double>& m2, int size) {
  std::vector<double> res = createMatrix(size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
        res[i * size + j] = m1[i * size + j] + m2[i * size + j];
    }
  }
  return res;
}

const std::vector<double> subtraction_of2m(const std::vector<double>& m1,
  const std::vector<double>& m2, int size) {
    std::vector<double> res = createMatrix(size);
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        res[i * size + j] = m1[i * size + j] - m2[i * size + j];
      }
    }
    return res;
}

const std::vector<double> addition_of4m(const std::vector<double>& m1,
  const std::vector<double>& m2, const std::vector<double>& m3,
  const std::vector<double>& m4, int size) {
  std::vector<double> res = createMatrix(size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      res[i * size + j] = m1[i * size + j] + m2[i * size + j] + m3[i * size + j]
          + m4[i * size + j];
      }
  }
  return res;
}


const std::vector<double> add_of3_sub_of4(const std::vector<double>& m1,
  const std::vector<double>& m2, const std::vector<double>& m3,
  const std::vector<double>& m4, int size) {
  std::vector<double> res = createMatrix(size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      res[i * size + j] = m1[i * size + j] + m2[i * size + j] +
          m3[i * size + j] - m4[i * size + j];
    }
  }
  return res;
}

void getFourMatrixBlocks(const std::vector<double>& A, std::vector<double>* A11,
  std::vector<double>* A12, std::vector<double>* A21, std::vector<double>* A22,
    int n_A) {
  int new_sz = n_A / 2;

  (*A11).resize(new_sz * new_sz);
  (*A12).resize(new_sz * new_sz);
  (*A21).resize(new_sz * new_sz);
  (*A22).resize(new_sz * new_sz);

  for (int i = 0; i < new_sz; i++) {
    for (int j = 0; j < new_sz; j++) {
      (*A11)[i * new_sz + j] = A[i * n_A + j];
      (*A12)[i * new_sz + j] = A[i * n_A + new_sz + j];
      (*A21)[i * new_sz + j] = A[(i + new_sz) * n_A + j];
      (*A22)[i * new_sz + j] = A[(i + new_sz) * n_A + new_sz + j];
    }
  }
}

const std::vector<double> getSequentialOperations(const std::vector<double>& m1,
    const std::vector<double>& m2, int n) {
  std::vector<double> res = createMatrix(n);
  if (n <= 4) {
    res = simple_mult(m1, m2, n);
  } else {
    std::vector<double> A(4);
    std::vector<double> B(4);
    std::vector<double> C(4);
    std::vector<double> A11, A12, A21, A22, B11, B12, B21, B22;
    std::vector<double> op1, op2, op3, op4, op5, op6, op7;
    std::vector<double> C11, C12, C21, C22;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++) {
        int new_idx = i * n + j, idx = 2 * i * n + j, sz = n;
        A[new_idx] = A11[idx];
        A[1*n + new_idx] = A12[idx + n];
        A[2*n + new_idx] = A21[idx + sz];
        A[3*n + new_idx] = A22[idx + sz + n];

        B[0 * n + new_idx] = B11[idx];
        B[1 * n + new_idx] = B12[idx + n];
        B[2 * n + new_idx] = B21[idx + sz];
        B[3 * n + new_idx] = B22[idx + sz + n];
      }

    std::vector<double> temp = addition_of2m(A11, A22, n);
    std::vector<double> tmp = addition_of2m(B11, B22, n);
    op1 = getSequentialOperations(temp, tmp, n);

    temp = addition_of2m(A21, A22, n);
    op2 = getSequentialOperations(temp, B11, n);

    temp = subtraction_of2m(B12, B22, n);
    op3 = getSequentialOperations(A11, temp, n);

    temp = subtraction_of2m(B21, B11, n);
    op4 = getSequentialOperations(A22, temp, n);

    temp = addition_of2m(A11, A12, n);
    op5 = getSequentialOperations(temp, B22, n);

    temp = subtraction_of2m(A21, A11, n);
    tmp = addition_of2m(B11, B12, n);
    op6 = getSequentialOperations(temp, tmp, n);

    temp = subtraction_of2m(A12, A22, n);
    tmp = addition_of2m(B21, B22, n);
    op7 = getSequentialOperations(temp, tmp, n);

    C11 = add_of3_sub_of4(op1, op4, op7, op5, n);
    C12 = addition_of2m(op3, op5, n);
    C21 = addition_of2m(op2, op4, n);
    C22 = add_of3_sub_of4(op1, op3, op6, op2, n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          res[i * 2 * n + j] = C11[i * n + j];
          res[i * 2 * n + j + n] = C12[i * n + j];
           res[i * 2 * n + j + 2 * n * n] = C21[i * n + j];
           res[i * 2 * n + j + 2 * n * n + n] = C22[i * n + j];
        }
    }
  }
  return res;
}

const std::vector<double> getParallelOperations(const std::vector<double>& m1,
  const std::vector<double>& m2, int n) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int sz, N;
  std::vector<double> res_m, res1_tmp, res3_tmp, res2_tmp, res4_tmp,
    res5_tmp;
  std::vector<double> A(size);
  std::vector<double> B(size);
  std::vector<double> A11, A12, A21, A22, B11, B12, B21, B22;
  MPI_Status status;
  if (n <= 64) {
    return res_m = getSequentialOperations(m1, m2, n);
  } else {
    if (rank == 0) {
      res_m = createMatrix(n);
      sz = sqrt(size);
      N = n / sz;
      std::vector<double> A11 = createMatrix(N);
      std::vector<double> A12 = createMatrix(N);
      std::vector<double> A21 = createMatrix(N);
      std::vector<double> A22 = createMatrix(N);
      std::vector<double> B11 = createMatrix(N);
      std::vector<double> B12 = createMatrix(N);
      std::vector<double> B21 = createMatrix(N);
      std::vector<double> B22 = createMatrix(N);
      getFourMatrixBlocks(A, &A11, &A12, &A21, &A22, N);
      getFourMatrixBlocks(B, &B11, &B12, &B21, &B22, N);
      MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
      for (int i = 1; i < size; i++) {
          MPI_Send(&A11[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          MPI_Send(&A12[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          MPI_Send(&A21[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          MPI_Send(&A22[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          MPI_Send(&B11[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          MPI_Send(&B12[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          MPI_Send(&B21[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          MPI_Send(&B22[0], N * N, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      }
    }
    if (size != 0) {
      MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
      sz = sqrt(size), N = n / sz;
      std::vector<double> A11 = createMatrix(N);
      std::vector<double> A12 = createMatrix(N);
      std::vector<double> A21 = createMatrix(N);
      std::vector<double> A22 = createMatrix(N);
      std::vector<double> B11 = createMatrix(N);
      std::vector<double> B12 = createMatrix(N);
      std::vector<double> B21 = createMatrix(N);
      std::vector<double> B22 = createMatrix(N);
      for (int i = 0; i < sz; i++) {
        MPI_Recv(&A11[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&A12[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&A21[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&A22[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&B11[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&B12[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&B21[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&B22[0], N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      }
    }
    res2_tmp = getSequentialOperations(A11, B11, N);
    res3_tmp = getSequentialOperations(A12, B12, N);
    res4_tmp = getSequentialOperations(A21, B21, N);
    res5_tmp = getSequentialOperations(A22, B22, N);

    if (size == 4)
      res1_tmp = addition_of2m(res2_tmp, res3_tmp, N);
    if (size == 16)
      res1_tmp = addition_of4m(res2_tmp, res3_tmp, res4_tmp, res5_tmp, N);

    if (rank != 0) {
      MPI_Send(&res1_tmp, N * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
      int coeff = sqrt(size);
      for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
          res_m[coeff * i * N + j] = res1_tmp[i * N + j];
      for (int k = 1; k < size; k++) {
        MPI_Recv(&res1_tmp, N * N, MPI_DOUBLE, k, 0, MPI_COMM_WORLD, &status);
        for (int i = 0; i < N; i++)
          for (int j = 0; j < N; j++)
            res_m[(k / coeff) * N * n + (k % coeff) * N + coeff * i * N + j]
            = res1_tmp[i * N + j];
      }
      return res_m;
    }
  }
  return res_m;
}
