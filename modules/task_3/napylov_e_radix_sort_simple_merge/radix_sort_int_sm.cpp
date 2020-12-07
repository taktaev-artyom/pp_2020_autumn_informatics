// Copyright 2020 Napylov Evgenii
#include "../../../modules/task_3/napylov_e_radix_sort_simple_merge/radix_sort_int_sm.h"
#include <mpi.h>
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <random>
#include <ctime>
#include <algorithm>

// #define DEBUG

void print_vec(std::vector<int> vec) {
    for (auto val : vec) {
        std::cout << val << ' ';
    }
    std::cout << std::endl;
}

std::vector<int> RandomVector(int len) {
    static std::mt19937 gen(time(0));
    std::vector<int> result(len);
    std::uniform_int_distribution<int> distr(-10000, 10000-1);
    for (int i = 0; i < len; i++) {
        result[i] = distr(gen);
    }
    return result;
}

std::vector<int> mergeVectors(std::vector<int> vec1, std::vector<int> vec2) {
    // Note: vec1.size() > 0 && vec2.size > 0
    std::vector<int> res(vec1.size() + vec2.size());
    int i = 0;  // index for vec1
    int j = 0;  // index for vec2

    for (int r = 0; r < res.size(); r++) {
        if (i > static_cast<int>(vec1.size()) - 1) {
            res[r] = vec2[j++];
        } else {
            if (j > static_cast<int>(vec2.size()) - 1) {
                res[r] = vec1[i++];
            } else {
                if (vec1[i] < vec2[j]) {
                    res[r] = vec1[i++];
                } else {
                    res[r] = vec2[j++];
                }
            }
        }
    }
    return res;
}

int getMaxDigitCount(std::vector<int> vec) {
    int max = *std::max_element(vec.begin(), vec.end());
    int min = *std::min_element(vec.begin(), vec.end());
    int maxabs;
    max > abs(min) ? maxabs = max : maxabs = min;
    int count = 0;
    while (maxabs != 0) {
        maxabs /= 10;
        count++;
    }
    return count;
}

std::vector<int> RadixSort(std::vector<int> vec) {
    //заводим 10 очередей
    const int max_digit_count = getMaxDigitCount(vec);
    std::vector<std::queue<int>> queue_arr(10);

    for (int i = 0; i < max_digit_count; i++) {
        // раскладываем по очередям
        for (int j = 0; j < static_cast<int>(vec.size()); j++) {
            int digit = abs((vec[j] / static_cast<int>(pow(10, i))) % 10);
            queue_arr[digit].push(vec[j]);
        }

        // достаем из очередей
        int k = 0;
        for (int d = 0; d < 10; d++) {
            while (!queue_arr[d].empty()) {
                vec[k] = queue_arr[d].front();
                queue_arr[d].pop();
                k++;
            }
        }
    }
    std::stack<int> st;
    // последний проход для отрицательных чисел
    for (int j = 0; j < static_cast<int>(vec.size()); j++) {
        if (vec[j] < 0) {
            st.push(vec[j]);
        } else {
            queue_arr[0].push(vec[j]);
        }
    }
    int k = 0;
    while (!st.empty()) {
        vec[k] = st.top();
        st.pop();
        k++;
    }
    while (!queue_arr[0].empty()) {
        vec[k] = queue_arr[0].front();
        queue_arr[0].pop();
        k++;
    }
    return vec;
}

std::vector<int> RadixSortParallel(std::vector<int> vec) {
    if (vec.size() <= 1) return vec;

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int worksize;

    if (size > static_cast<int>(vec.size())) {
        worksize = vec.size();
    } else {
        worksize = size;
    }

    const int count = vec.size() / worksize;
    const int rem = vec.size() % worksize;

    std::vector<int> counts(worksize);

    int local_size;

    if (rank < rem) {
        local_size = count + 1;
    } else {
        local_size = count;
    }

    if (rank == 0) {
        for (int p = 0; p < worksize; p++) {
            counts[p] = count + (p < rem ? 1 : 0);
        }
    }
    MPI_Bcast(counts.data(), worksize, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> local_vec(local_size);

    if (rank == 0) {
        local_vec = std::vector<int>(vec.begin(), vec.begin() + local_size);
        int offset = local_size;
        for (int proc = 1; proc < worksize; proc++) {
            MPI_Send(vec.data() + offset, counts[proc], MPI_INT, proc, 0, MPI_COMM_WORLD);
            offset += counts[proc];
        }
    }
    if (rank > 0 && rank < worksize) {
        MPI_Status st;
        MPI_Recv(local_vec.data(), counts[rank], MPI_INT, 0, 0, MPI_COMM_WORLD, &st);
    }
    /*
        Данные распределены:
            local_vec
            local_vec_size
            counts[]
            worksize
    */

    // сортируем в кажом процессе последовательным алгоритмом
    local_vec = RadixSort(local_vec);

    /*
        Части отсортированы 
        Осталось их слить
    */
    #ifdef DEBUG
    // std::cout << rank << ": "; print_vec(local_vec);
    #endif

    int dist = 1;
    for (int rem_proc = worksize; rem_proc > 1; rem_proc = rem_proc / 2 + rem_proc % 2) {
        if (rank < worksize) {
            int left = rank - dist;
            if (left % (2 * dist) == 0) {
                MPI_Send(&local_size, 1, MPI_INT, left, 0, MPI_COMM_WORLD);
                MPI_Send(local_vec.data(), local_size, MPI_INT, left, 0, MPI_COMM_WORLD);
                #ifdef DEBUG
                std::cout << dist << ". Send from " << rank << " to " << left << " sz " << local_size << std::endl;
                #endif
            }
        }
        if (rank < worksize) {
            int right = rank + dist;
            if ((right < worksize) && (rank % (2 * dist) == 0)) {
                MPI_Status st;
                int recv_count;
                MPI_Recv(&recv_count, 1, MPI_INT, right, 0, MPI_COMM_WORLD, &st);
                std::vector<int> tmp_vec(recv_count);
                MPI_Recv(tmp_vec.data(), recv_count, MPI_INT, right, 0, MPI_COMM_WORLD, &st);
                #ifdef DEBUG
                std::cout << dist << ". Recv from " << right << " with sz " << recv_count << std::endl;
                #endif
                local_vec = mergeVectors(local_vec, tmp_vec);
                local_size = local_size + recv_count;
            }
            dist *= 2;
        }
    }
    return local_vec;
}

int compare(const void *a, const void *b) {
  return (*static_cast<const int*>(a) - *static_cast<const int*>(b));
}
