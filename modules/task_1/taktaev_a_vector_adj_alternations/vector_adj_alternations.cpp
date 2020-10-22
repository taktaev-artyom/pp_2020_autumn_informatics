// Copyright 2020 Taktaev Artem
#include <mpi.h>
#include <random>
#include <vector>
#include "../../../modules/task_1/taktaev_a_vector_adj_alternations/vector_adj_alternations.h"

std::vector<int> createRandomVector(int vec_size) {
    if (vec_size <= 0) throw "Wrong size";
    std::random_device rand_d;
    std::mt19937 gen(rand_d());
    std::vector<int> res_vec(vec_size);
    for (int  i = 0; i < vec_size; i++)
        res_vec[i] = static_cast<int>(gen() % 201) - 100;
    return res_vec;
}

int calculateAdjAlternationsSequential(const std::vector<int> &vec, int inc, int start_index) {
    if (start_index < 1) throw "Wrong start index";
    if (inc < 1) throw "Wrong increment";
    int count = 0;
    for (int i = start_index; i < vec.size(); i = i + inc)
        if (vec[i] * vec[i - 1] < 0) count++;
    return count;
}

int calculateAdjAlternationsParallel(const std::vector<int> &vec, int vec_size) {
    int procNum, procRank;
    if (vec_size <= 0) throw "Wrong size";
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int part_vec_size = vec_size / procNum;
    int tail = vec_size - part_vec_size * procNum;
    if (procRank == 0) {
        for (int i = 1; i < procNum; i++)
            MPI_Send(&vec[0] + part_vec_size * i + tail, part_vec_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    std::vector<int> part_vec(procRank == 0 ? vec_size + tail : vec_size);
    int full_count = 0;
    if (procRank == 0) {
        part_vec = std::vector<int>(vec.begin(), vec.begin() + vec_size + tail);
    full_count += calculateAdjAlternationsSequential(vec, part_vec_size, tail + part_vec_size);
    } else {
        MPI_Status status;
        MPI_Recv(&part_vec[0], part_vec_size, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    int part_count = calculateAdjAlternationsSequential(part_vec, 1, 1);
    MPI_Reduce(&part_count, &full_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return full_count;
}
