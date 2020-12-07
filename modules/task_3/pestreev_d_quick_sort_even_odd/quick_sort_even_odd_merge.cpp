// Copyright 2020 Pestreev Daniil
#include <mpi.h>
#include <iostream>
#include <vector>
#include <random>
#include <limits>
#include <utility>
#include <ctime>
#include "../../../modules/task_3/pestreev_d_quick_sort_even_odd/quick_sort_even_odd_merge.h"

std::vector<std::pair<int, int>> comps;

std::vector<int> getRandomVector(int size) {
    std::mt19937 gen;
    gen.seed(time(0));
    std::vector<int> a(size);

    for (int i = 0; i < size; i++) {
        a[i] = gen();
    }
    return a;
}

void qsort(int* vec, int left, int right) {
    int mid;
    int tmp;
    int l = left;
    int r = right;
    mid = vec[(l + r) / 2];
    do {
        while (vec[l] < mid) l++;
        while (vec[r] > mid) r--;
        if (l <= r) {
            std::swap(vec[l], vec[r]);
            l++;
            r--;
        }
    } while (l < r);
    if (left < r) qsort(vec, left, r);
    if (l < right) qsort(vec, l, right);
}

std::vector<int> quickSortV(const std::vector<int>& vec) {
    int size = vec.size();
    int* arr = new int[size];
    for (int i = 0; i < size; i++) {
        arr[i] = vec[i];
    }
    qsort(arr, 0, size - 1);
    std::vector<int> tmp;
    for (int i = 0; i < size; i++) {
        tmp.push_back(arr[i]);
    }
    return tmp;
}

void recur_merge(const std::vector<int>& left, const std::vector<int>& right) {
    int array_count = left.size() + right.size();
    if (array_count <= 1) {
        return;
    }
    if (array_count == 2) {
        comps.push_back(std::pair<int, int>(left[0], right[0]));
        return;
    }
    std::vector<int> left_odd;
    std::vector<int> left_even;
    std::vector<int> right_odd;
    std::vector<int> right_even;

    for (int i = 0; i < left.size(); i++) {
        if (i % 2) {
            left_even.push_back(left[i]);
        } else {
            left_odd.push_back(left[i]);
        }
    }

    for (int i = 0; i < right.size(); i++) {
        if (i % 2) {
            right_even.push_back(right[i]);
        } else {
            right_odd.push_back(right[i]);
        }
    }

    recur_merge(left_odd, right_odd);
    recur_merge(left_even, right_even);

    std::vector<int> res;

    for (int i = 0; i < left.size(); i++) {
        res.push_back(left[i]);
    }
    for (int i = 0; i < right.size(); i++) {
        res.push_back(right[i]);
    }

    for (int i = 1; i + 1 < array_count; i += 2) {
        comps.push_back(std::pair<int, int>(res[i], res[i + 1]));
    }
}

void networking(const std::vector<int>& arr) {
    int arr_size = arr.size();

    if (arr_size <= 1) {
        return;
    }

    std::vector<int> l(arr.begin(), arr.begin() + arr_size / 2);
    std::vector<int> r(arr.begin() + arr_size / 2, arr.begin() + arr_size);
    networking(l);
    networking(r);
    recur_merge(l, r);
}

void batchers_network(int proc_size) {
    if (proc_size <= 1)
        return;
    std::vector<int> procs;
    for (int i = 0; i < proc_size; i++) {
        procs.push_back(i);
    }

    networking(procs);
}

std::vector<int> parallel_sorting(const std::vector<int>& vec) {
    int proc_rank, proc_size;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    std::vector<int> globalV = vec;
    int vecsizeG = globalV.size();

    if (proc_size < 0) {
        return globalV;
    }

    if (proc_size >= globalV.size() || proc_size == 1) {
        if (proc_rank == 0) {
            globalV = quickSortV(vec);
        }
        return globalV;
    }
    int exitcircle = -1;
    int vecsizechange = 0;
    if (proc_rank == 0) {
        while (exitcircle < 0) {
            exitcircle = -1;
            if (globalV.size() % proc_size) {
                int o = std::numeric_limits<int>::max();
                vecsizechange++;
                globalV.push_back(o);
            } else {
                exitcircle = 1;
            }
        }
    }
    vecsizeG = globalV.size();

    if (proc_rank == 0) {
        for (int i = 1; i < proc_size; i++) {
            MPI_Send(&vecsizeG, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&vecsizeG, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    int vecsizeL = vecsizeG / proc_size;

    std::vector<int> localV(vecsizeL);
    std::vector<int> neighboringV(vecsizeL);
    std::vector<int> tmpV(vecsizeL);
    MPI_Scatter(&globalV[0], vecsizeL, MPI_INT, &localV[0], vecsizeL, MPI_INT, 0, MPI_COMM_WORLD);

    batchers_network(proc_size);
    const std::vector<int> v = localV;
    localV = quickSortV(v);
    for (int i = 0; i < comps.size(); i++) {
        if (proc_rank == comps[i].first) {
            MPI_Send(&localV[0], vecsizeL, MPI_INT, comps[i].second, 0, MPI_COMM_WORLD);
            MPI_Recv(&neighboringV[0], vecsizeL, MPI_INT, comps[i].second, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int localidx = 0;
            int neighboridx = 0;
            for (int tmp_index = 0; tmp_index < vecsizeL; tmp_index++) {
                int local = localV[localidx];
                int neighbor = neighboringV[neighboridx];
                if (local < neighbor) {
                    tmpV[tmp_index] = local;
                    localidx++;
                } else {
                    tmpV[tmp_index] = neighbor;
                    neighboridx++;
                }
            }
            localV = tmpV;
        } else if (proc_rank == comps[i].second) {
            MPI_Recv(&neighboringV[0], vecsizeL, MPI_INT, comps[i].first, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&localV[0], vecsizeL, MPI_INT, comps[i].first, 0, MPI_COMM_WORLD);
            int start = vecsizeL - 1;
            int localidx = start;
            int neighboridx = start;
            for (int tmp_index = start; tmp_index >= 0; tmp_index--) {
                int local = localV[localidx];
                int neighbor = neighboringV[neighboridx];
                if (local > neighbor) {
                    tmpV[tmp_index] = local;
                    localidx--;
                } else {
                    tmpV[tmp_index] = neighbor;
                    neighboridx--;
                }
            }
            localV = tmpV;
        }
    }
    MPI_Gather(&localV[0], vecsizeL, MPI_INT, &globalV[0], vecsizeL, MPI_INT, 0, MPI_COMM_WORLD);

    if (proc_rank == 0 && vecsizechange > 0)
        globalV.erase(globalV.begin() + vecsizeG - vecsizechange, globalV.begin() + vecsizeG);
    return globalV;
}
