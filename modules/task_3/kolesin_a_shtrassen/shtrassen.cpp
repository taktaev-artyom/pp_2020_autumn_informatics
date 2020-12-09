// Copyright 2020 Kolesin Andrey
#include <deque>
#include <vector>
#include "../../../modules/task_3/kolesin_a_shtrassen/shtrassen.h"

Matrix::Matrix(int *_buff, int _n, bool _tmp) {
    buff = _buff;
    n = _n;
    N = _n;
    tmp = _tmp;
}
Matrix::Matrix(Matrix old, int x, int y) {
    tmp = false;
    buff = old.buff;
    n = old.n;
    N = old.N / 2;
    coords = old.coords;
    coords.push_back(x);
    coords.push_back(y);
}
Matrix::Matrix(const Matrix &old) {
    buff = old.buff;
    n = old.n;
    N = old.N;
    coords = old.coords;
    tmp = false;
}

Matrix::~Matrix() {
    if (tmp) {
        delete[] buff;
    }
}
int &Matrix::getElem(int x, int y) {
    std::deque<int> tmp = coords;
    int start_x = 0, start_y = 0;
    for (int N = n; !tmp.empty(); N /= 2) {
        int X = tmp.front();
        tmp.pop_front();
        int Y = tmp.front();
        tmp.pop_front();
        start_x += X * N / 2;
        start_y += Y * N / 2;
    }
    x += start_x;
    y += start_y;
    return buff[x + n * y];
}
void Matrix::print() {
    std::cout << "++++++++++++++\n";
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            std::cout << getElem(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "-------------" << std::endl;
}
void Matrix::stretch(int *buff) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            buff[i + j * N] = getElem(i, j);
        }
    }
}
void Matrix::Send(MPI_Comm comm, int sender_rank, int reciver_rank, int marker) {
    int rank;
    MPI_Comm_rank(comm, &rank);
    MPI_Request sendrequest, recvrequest;
    MPI_Status status;
    if (reciver_rank == sender_rank) {
        std::cout << "Астанавитес" << std::endl;
        return;
    }
    if (rank == sender_rank) {
        std::vector<int> buff(N * N);
        stretch(&buff[0]);
        MPI_Isend(&buff[0], N * N, MPI_INT, reciver_rank, marker, comm, &sendrequest);
    }
    if (rank == reciver_rank) {
        std::vector<int> buff(N * N);
        MPI_Irecv(&buff[0], N * N, MPI_INT, sender_rank, marker, comm, &recvrequest);
        MPI_Wait(&recvrequest, &status);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                getElem(i, j) = buff[i + j * N];
            }
        }
    }
}

std::vector<int> getRandomMatrix(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> res;
    if (size == -1) {
        size = gen() % 5;
        size = pow(2, size);
    }
    size = size * size;
    for (int i = 0; i < size; i++) {
        res.push_back(gen() % 10);
    }
    return res;
}
void SimpleMult(Matrix A, Matrix B, Matrix C) {
    for (int i = 0; i < C.N; i++) {
        for (int j = 0; j < C.N; j++) {
            C.getElem(i, j) = 0;
            for (int k = 0; k < C.N; k++) {
                C.getElem(i, j) += A.getElem(i, k) * B.getElem(k, j);
            }
        }
    }
}
std::vector<int> getRange(int begin, int end) {
    std::vector<int> res;
    for (int i = begin; i < end; i++) {
        res.push_back(i);
    }
    return res;
}
MPI_Comm getComm(int n, int new_size, MPI_Comm old_comm) {
    MPI_Comm new_comm;
    std::vector<int> range = getRange(new_size * n, new_size * (n + 1));
    MPI_Group new_group, old_group;
    MPI_Comm_group(old_comm, &old_group);

    MPI_Group_incl(old_group, new_size, &range[0], &new_group);
    MPI_Comm_create(old_comm, new_group, &new_comm);
    return new_comm;
}
void calcSubC(int x, int y, Matrix A[2][2], Matrix B[2][2], Matrix C[2][2], MPI_Comm comm, int level) {
    Matrix tmp1(new int[A[0][0].N * A[0][0].N], A[0][0].N, true);
    Matrix tmp2(new int[A[0][0].N * A[0][0].N], A[0][0].N, true);
    Sht(A[x][0], B[0][y], tmp1, comm, level);
    Sht(A[x][1], B[1][y], tmp2, comm, level);
    Sum(tmp1, tmp2, C[x][y]);
}
int pow(int x) {
    int res = 1;
    for (int i = 0; i < x; i++) {
        res *= 4;
    }
    return res;
}

int getPrevPower(int x) {
    int r = -1;
    while (x != 0) {
        x/=4;
        r++;
    }
    return pow(r);
}
void SafeSht(Matrix A, Matrix B, Matrix C, MPI_Comm comm) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (A.N != B.N || B.N != C.N) {
        throw -1;
    }
    MPI_Comm new_comm;
    int reduced_size = size;
    if (A.N*A.N < reduced_size) {
        reduced_size = A.N*A.N;
    }
    reduced_size = getPrevPower(reduced_size);
    int color = rank < reduced_size ? 1 : 0;
    // std::cout<<rank<<" "<<color<<"        "<<reduced_size<<std::endl;
    MPI_Comm_split(comm, color, 0, &new_comm);
    if (color == 1) {
            Sht(A, B, C, new_comm);
    }
    MPI_Comm_free(&new_comm);
}
void Sht(Matrix A, Matrix B, Matrix C, MPI_Comm comm, int level) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    if (size == 1) {
        ShtSeq(A, B, C);
        // C.print();
        return;
    }
    MPI_Comm new_comm;
    int color = rank % 4;
    MPI_Comm_split(comm, color, 0, &new_comm);

    Matrix SubA[2][2] = {{Matrix(A, 0, 0), Matrix(A, 1, 0)},
                         {Matrix(A, 0, 1), Matrix(A, 1, 1)}};
    Matrix SubB[2][2] = {{Matrix(B, 0, 0), Matrix(B, 1, 0)},
                         {Matrix(B, 0, 1), Matrix(B, 1, 1)}};
    Matrix SubC[2][2] = {{Matrix(C, 0, 0), Matrix(C, 1, 0)},
                         {Matrix(C, 0, 1), Matrix(C, 1, 1)}};

    SubA[0][0].Send(comm, 0, 1, 0 + pow(size));
    SubB[0][1].Send(comm, 0, 1, 1 + pow(size));
    SubA[0][1].Send(comm, 0, 1, 2 + pow(size));
    SubB[1][1].Send(comm, 0, 1, 3 + pow(size));

    SubA[1][0].Send(comm, 0, 2, 4 + pow(size));
    SubB[0][0].Send(comm, 0, 2, 5 + pow(size));
    SubA[1][1].Send(comm, 0, 2, 6 + pow(size));
    SubB[1][0].Send(comm, 0, 2, 7 + pow(size));

    SubA[1][0].Send(comm, 0, 3, 8 + pow(size));
    SubB[0][1].Send(comm, 0, 3, 9 + pow(size));
    SubA[1][1].Send(comm, 0, 3, 10 + pow(size));
    SubB[1][1].Send(comm, 0, 3, 11 + pow(size));
    int x, y;
    switch (color) {
    case 0:
        x = 0;
        y = 0;
        break;
    case 1:
        x = 0;
        y = 1;
        break;
    case 2:
        x = 1;
        y = 0;
        break;
    case 3:
        x = 1;
        y = 1;
        break;
    }
    calcSubC(x, y, SubA, SubB, SubC, new_comm, level + 1);
    SubC[0][1].Send(comm, 1, 0, 100 * (1 + pow(size)));
    SubC[1][0].Send(comm, 2, 0, 100 * (2 + pow(size)));
    SubC[1][1].Send(comm, 3, 0, 100 * (3 + pow(size)));
    MPI_Comm_free(&new_comm);
}
void ShtSeq(Matrix A, Matrix B, Matrix C) {
    if (C.N == 1) {
        C.getElem(0, 0) = A.getElem(0, 0) * B.getElem(0, 0);
    } else {
        int halfNSqr = (C.N / 2) * (C.N / 2);
        Matrix tmp1(new int[halfNSqr], C.N / 2, true);
        Matrix tmp2(new int[halfNSqr], C.N / 2, true);
        Matrix A11(A, 0, 0);
        Matrix A12(A, 1, 0);
        Matrix A21(A, 0, 1);
        Matrix A22(A, 1, 1);

        Matrix B11(B, 0, 0);
        Matrix B12(B, 1, 0);
        Matrix B21(B, 0, 1);
        Matrix B22(B, 1, 1);

        Matrix C11(C, 0, 0);
        Matrix C12(C, 1, 0);
        Matrix C21(C, 0, 1);
        Matrix C22(C, 1, 1);

        ShtSeq(A11, B11, tmp1);
        ShtSeq(A12, B21, tmp2);
        Sum(tmp1, tmp2, C11);  // C11 = A11*B11 + A12*B21
        ShtSeq(A11, B12, tmp1);
        ShtSeq(A12, B22, tmp2);
        Sum(tmp1, tmp2, C12);  // C12 = A11*B12 + A12*B22
        ShtSeq(A21, B11, tmp1);
        ShtSeq(A22, B21, tmp2);
        Sum(tmp1, tmp2, C21);  // C21 = A21*B11 + A22*B21
        ShtSeq(A21, B12, tmp1);
        ShtSeq(A22, B22, tmp2);
        Sum(tmp1, tmp2, C22);  // C22 = A21*B12 + A22*B22
    }
}
void Sum(Matrix A, Matrix B, Matrix C) {
    for (int j = 0; j < C.N; j++) {
        for (int i = 0; i < C.N; i++) {
            C.getElem(i, j) = A.getElem(i, j) + B.getElem(i, j);
        }
    }
}
