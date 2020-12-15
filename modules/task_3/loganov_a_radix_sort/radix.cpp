// Copyright 2020 Loganov Andrei
#include "../../../../modules/task_3/loganov_a_radix_sort/radix.h"
#include <stdlib.h>
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <functional>
#include <utility>
#include <string>
std::vector<double> EvenSpliter(std::vector<double> vec1, std::vector<double> vec2, std::vector<double> res) {
    int size1 = static_cast<int>(vec1.size());
    int size2 = static_cast<int>(vec2.size());
    res.resize(size1 + size2);
    int flag = 0;
    if ((size1 == 1) && (size2 == 1)) {
        if (vec1[0] < vec2[0]) {
            res[0] = vec1[0];
            res[1] = vec2[0];
        } else {
            res[0] = vec2[0];
            res[1] = vec1[0];
        }
        flag = 1;
    }
    if (flag == 0) {
        int a = 0;
        int b = 0;
        int i = 0;
        while ((a < size1) && (b < size2)) {
            if (vec1[a] <= vec2[b]) {
                res[i] = vec1[a];
                a += 2;
            } else {
                res[i] = vec2[b];
                b += 2;
            }
            i += 2;
        }
        if (a >= size1) {
            for (int j = b; j < size2; j += 2, i += 2) {
                if (i >= static_cast<int>(res.size())) {
                    res[i - 1] = vec2[j];
                } else {
                    res[i] = vec2[j];
                }
            }
        }
        if (b >= size2) {
            for (int j = a; j < size1; i += 2, j += 2) {
                if (i >= static_cast<int>(res.size())) {
                    res[i - 1] = vec1[j];
                } else {
                    res[i] = vec1[j];
                }
            }
        }
    }
    return res;
}
std::vector<double> OddSpliter(std::vector<double> vec1, std::vector<double> vec2, std::vector<double> res) {
    int size1 = static_cast<int>(vec1.size());
    int size2 = static_cast<int>(vec2.size());
    res.resize(size1 + size2);
    int a = 1;
    int b = 1;
    int i = 1;
    while ((a < size1) && (b < size2)) {
        if (vec1[a] <= vec2[b]) {
            res[i] = vec1[a];
            a += 2;
        } else {
            res[i] = vec2[b];
            b += 2;
        }
        i += 2;
    }
    if (a >= size1) {
        for (int j = b; j < size2; j += 2, i += 2) {
            res[i] = vec2[j];
        }
    }
    if (b >= size2) {
        for (int j = a; j < size1; j += 2, i += 2) {
            res[i] = vec1[j];
        }
    }
    return res;
}
std::vector<double> simpmerg(std::vector<double> res, std::vector<double> chet, std::vector<double> nchet) {
    int  size = static_cast<int>(chet.size());
    res = std::vector<double>(size);
    for (int i = 0; i < static_cast<int>(nchet.size()); i++) {
        if ((i % 2 != 0) && (nchet[i] != 0)) {
            res[i] = nchet[i];
        } else {
            res[i] = chet[i];
        }
    }
    for (int i = 1; i < static_cast<int>(res.size()); i++) {
        if (res[i] < res[i - 1]) {
            std::swap(res[i], res[i - 1]);
        }
    }
    return res;
}
std::vector<double> getRandomVector(int size) {
    std::vector<int> vec(size);
    std::vector<double> vec1(size);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < size; i++) {
        vec[i] = gen() % 10000;
    }
    for (int i = 0; i < size; i++) {
        vec1[i] = vec[i]/100;
    }
    return vec1;
}
bool is_exp_of_2(int n) {
    int k = 1;
    while (k < n) {
        k *= 2;
    }
    if (k == n) {
        return true;
    } else {
        return false;
    }
}
int countBelowPoint(double x) {
    std::string str = std::to_string(x);
    int result = 0;
    if (str.find('.')) {
        int k = static_cast<int>(str.find('.'));
        int size = static_cast<int>(str.size());
        result = size - k - 1;
    }
    return result;
}
int countBeforepoint(int x) {
    int temp = 0;
    while (x > 0) {
        x /= 10;
        temp++;
    }
    return temp;
}
int numbyrank(int dis, double x) {
    if (dis < 0) {
        return uint64_t(x * pow(10, -dis)) % 10;
    }
    double mask = pow(10, dis);
    double tmp = x / mask;
    return static_cast<int>(tmp) % 10;
}
std::vector<double> countingsort2(const std::vector<double>& res, int razrad) {
    typedef  std::vector<std::vector<double>> Matrix;
    std::vector<double> result;
    Matrix t(10);
    for (int i = 0; i < 10; i++) {
        t[i] = std::vector<double>(0);
    }
    for (int i = 0; i < static_cast<int>(res.size()); i++) {
        t[numbyrank(razrad, res[i])].push_back(res[i]);
    }
    for (int i = 0; i < 10; i++) {
        for (auto j : t[i]) {
            result.push_back(j);
        }
    }
    return result;
}
std::vector<double> seqRadixSort(const std::vector<double> tmp) {
    std::vector<double> result(tmp.begin(), tmp.end());
    double maxelem = tmp[0];
    for (int i = 0; i < static_cast<int>(tmp.size()); i++) {
        maxelem = std::max(maxelem, tmp[i]);
    }
    int znbeforePoint = countBeforepoint(maxelem);
    int cpp = 0;
    for (int i = 0; i < static_cast<int>(tmp.size()); i++) {
        cpp = std::max(cpp, countBelowPoint(tmp[i]));
    }
    for (int i = -cpp; i<= znbeforePoint; i++) {
        result = countingsort2(result, i);
    }
    for (int i = 1; i < static_cast<int>(result.size()); i++) {
        if (result[i] < result[i - 1]) {
            std::swap(result[i], result[i - 1]);
        }
    }
    return result;
}
std::vector<double> ParallelSort(std::vector<double> vec) {
    int rank;
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ((size == 1) || (size>static_cast<int>(vec.size()))) {
        vec = seqRadixSort(vec);
        return vec;
    }
    int sizeofvec = vec.size();
    int Delta = sizeofvec / size;
    int ost = sizeofvec % size;
    std::vector<double> localvec(Delta);
    if (rank == 0) {
        localvec.resize(Delta + ost);
        localvec = std::vector < double >(vec.begin(), vec.begin() + ost + Delta);
        for (int i = 1; i < size; i++) {
            MPI_Send(vec.data() + Delta * i + ost, Delta, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Status status;
        MPI_Recv(localvec.data(), Delta, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }
    localvec = seqRadixSort(localvec);
    if (rank % 2 != 0) {
        std::vector<double> result;
        MPI_Send(localvec.data(), Delta, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
        MPI_Status status;
        int count;
        (rank - 1 == 0) ? count = Delta + ost : count = Delta;
        std::vector<double> getvec(count);
        MPI_Recv(getvec.data(), count, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &status);
        result = OddSpliter(getvec, localvec, result);
        MPI_Send(result.data(), static_cast<int>(result.size()), MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
    }
    if ((rank % 2 == 0) && (rank + 1 < size)) {
        std::vector<double> result;
        std::vector<double> getvec(Delta);
        MPI_Status status;
        MPI_Recv(getvec.data(), Delta, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status);
        int count;
        (rank == 0) ? count = Delta + ost : count = Delta;
        MPI_Send(localvec.data(), count, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD);
        result = EvenSpliter(getvec, localvec, result);
        MPI_Status status2;
        std::vector<double> res(static_cast<int>(result.size()));
        MPI_Recv(res.data(), static_cast<int>(result.size()), MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status2);
        localvec = simpmerg(localvec, result, res);
    }
    if ((rank % 2 == 0) && (is_exp_of_2(rank) == false) && (rank != 0)) {
        int off = rank;
        int sdvig = 0;
        while (is_exp_of_2(off) == false) {
            off--;
            sdvig++;
        }
        int sz = static_cast<int>(localvec.size());
        int sz2;
        MPI_Send(&sz, 1, MPI_INT, rank - sdvig, 3, MPI_COMM_WORLD);
        MPI_Status status1;
        MPI_Recv(&sz2, 1, MPI_INT, rank - sdvig, 3, MPI_COMM_WORLD, &status1);
        std::vector<double> getvec(sz2);
        MPI_Send(localvec.data(), sz, MPI_DOUBLE, rank - sdvig, 3, MPI_COMM_WORLD);
        MPI_Status status;
        MPI_Recv(getvec.data(), sz2, MPI_DOUBLE, rank - sdvig, 3, MPI_COMM_WORLD, &status);
        std::vector<double> result(sz + sz2);
        result = OddSpliter(getvec, localvec, result);
        MPI_Send(result.data(), sz + sz2, MPI_DOUBLE, rank - sdvig, 3, MPI_COMM_WORLD);
    }
    if ((is_exp_of_2(rank) == true) && (rank + 2 < size) && (is_exp_of_2(rank + 2) == false) && (rank % 2 == 0)) {
        int off = rank;
        int sdvig = 2;
        while ((is_exp_of_2(off + sdvig) ==false) && ((off + sdvig) < size)) {
            MPI_Status status;
            int sz2;
            MPI_Recv(&sz2, 1, MPI_INT, rank + sdvig, 3, MPI_COMM_WORLD, &status);
            std::vector<double> getvec(sz2);
            int sz = static_cast<int>(localvec.size());
            MPI_Send(&sz, 1, MPI_INT, rank + sdvig, 3, MPI_COMM_WORLD);
            MPI_Status status1;
            MPI_Recv(getvec.data(), sz2, MPI_DOUBLE, rank + sdvig, 3, MPI_COMM_WORLD, &status1);
            std::vector<double> result(sz + sz2);
            result = EvenSpliter(getvec, localvec, result);
            MPI_Send(localvec.data(), sz, MPI_DOUBLE, rank + sdvig, 3, MPI_COMM_WORLD);
            std::vector<double> res(sz + sz2);
            MPI_Status status3;
            MPI_Recv(res.data(), sz + sz2, MPI_DOUBLE, rank + sdvig, 3, MPI_COMM_WORLD, &status3);
            localvec = simpmerg(localvec, result, res);
            sdvig += 2;
        }
    }
    if ((is_exp_of_2(rank) == true) && (rank != 0) && (rank % 2 == 0)) {
        int sz = static_cast<int>(localvec.size());
        MPI_Send(&sz, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(localvec.data(), sz, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        int count = 2;
        while (count < size) {
            MPI_Status status1;
            int sz2;
            MPI_Recv(&sz2, 1, MPI_INT, count, 1, MPI_COMM_WORLD, &status1);
            std::vector<double> getvec(sz2);
            MPI_Status status;
            MPI_Recv(getvec.data(), sz2, MPI_DOUBLE, count, 1, MPI_COMM_WORLD, &status);
            std::vector<double> result(sz2 + static_cast<int>(localvec.size()));
            result = EvenSpliter(getvec, localvec, result);
            result = OddSpliter(getvec, localvec, result);
            int sz = static_cast<int>(localvec.size());
            localvec.resize(sz + sz2);
            localvec = result;
            localvec = simpmerg(localvec, result, result);
            count *= 2;
        }
    }
    return localvec;
}

