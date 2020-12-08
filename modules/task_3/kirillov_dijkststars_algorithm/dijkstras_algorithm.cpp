// Copyright 2020 Kirillov Konstantin
#include <mpi.h>
#include<iomanip>
#include<iostream>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>
#include "../../../modules/task_3/kirillov_dijkststars_algorithm/dijkstrat_algorithm.h"
std::vector<int> getRandomGraph(int sizeGraph) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int>graph(sizeGraph*sizeGraph);
    for (int i = 0; i < sizeGraph; i++) {
        for (int j = 0; j < sizeGraph; j++) {
            if (i == j) {
                graph[i*sizeGraph+j] = 0;
            } else {
                if (j < i) {
                    graph[i*sizeGraph + j] = graph[j*sizeGraph + i];
                } else {
                    graph[i*sizeGraph + j] = gen() % 10;
                }
            }
        }
    }
    return graph;
}

void printGraph(std::vector<int>graph) {
    for (int i = 0; i < sqrt(graph.size()); i++) {
        for (int j = 0; j < sqrt(graph.size()); j++) {
            std::cout <<std::setw(2)<< graph[i* sqrt(graph.size()) + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


std::vector<int> getSequentialDijkstras(std::vector<int>graph, int start) {
    int sizeGraph = sqrt(graph.size());
    std::vector<int>dist(sizeGraph, INT_MAX);
    std::vector<bool>visited(sizeGraph, false);
    int index, u;
    dist[start] = 0;

    for (int i = 0; i < sizeGraph - 1; i++) {
        int min = INT_MAX;
        for (int j = 0; j < sizeGraph; j++) {
            if (!visited[j] && dist[j] <= min) {
                min = dist[j];
                index = j;
            }
        }
        u = index;
        visited[u] = true;
        for (int k = 0; k < sizeGraph; k++) {
            if (!visited[k]&&dist[u] != INT_MAX&&graph[u*sizeGraph + k]
                && dist[u] + graph[u*sizeGraph + k] < dist[k]) {
                dist[k] = dist[u] + graph[u*sizeGraph + k];
            }
        }
    }
    return dist;
}

void printDist(std::vector<int>dist) {
    for (int i = 0; i < dist.size(); i++) {
        std::cout << std::setw(2) << dist[i] << " ";
    }
    std::cout << std::endl;
}


std::vector<int> getParallelDijkstras(std::vector<int>graph, int start) {
    int procNum, procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int sizeGraph = sqrt(graph.size());
    std::vector<bool>visited(sizeGraph, false);
    std::vector<int>dist(sizeGraph, INT_MAX);
    int *global_dist = new int[sizeGraph*procNum];
    int *global_visited = new int[sizeGraph*procNum];

    if (sizeGraph < 2) {
        throw std::runtime_error("Wrong size");
    }
    dist[start] = 0;
    int min = INT_MAX;
    std::vector<int>local_graph(sizeGraph * (sizeGraph / procNum));

    struct {
        int value;
        int index;
    }local_vec, global_vec;

    int delta;
    if (procRank == 0) {
    delta = 0;
    } else {
        delta = sizeGraph % procNum + (sizeGraph / procNum) * procRank;
    }

    MPI_Scatter(&graph[(sizeGraph % procNum) * sizeGraph], (sizeGraph / procNum) * sizeGraph,
    MPI_INT, &local_graph[0], (sizeGraph / procNum) * sizeGraph, MPI_INT, 0, MPI_COMM_WORLD);

    if (procRank == 0) {
        local_graph.insert(
        local_graph.begin(),
        graph.begin(),
        graph.begin() + (sizeGraph % procNum) * sizeGraph);
    }

    int localGraphSize = local_graph.size() / sizeGraph;

    for (int i = 0; i < sizeGraph - 1; i++) {
        local_vec.value = -1;
        local_vec.index = -1;
        for (int j = delta; j < localGraphSize + delta; j++) {
            if (!visited[j] && (local_vec.index == -1 || dist[j] < dist[local_vec.index])) {
                local_vec.index = j;
            local_vec.value = dist[local_vec.index];
            }
        }

        if (local_vec.index == -1) {
            local_vec.value = min;
        }
        MPI_Allreduce(&local_vec, &global_vec, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        if (global_vec.index == -1 || dist[global_vec.index] == min) {
            break;
        }
        visited[global_vec.index] = true;

    for (int k = 0; k < localGraphSize; k++) {
        if (local_graph[global_vec.index + sizeGraph * k] != 0 &&
        dist[global_vec.index] + local_graph[global_vec.index + sizeGraph * k] < dist[k + delta]) {
        dist[k + delta] = dist[global_vec.index] + local_graph[global_vec.index + sizeGraph * k];
        }
    }

    if (procRank == 0) {
        dist.resize(sizeGraph % procNum + sizeGraph / procNum);
        std::vector<int> recv(sizeGraph / procNum);
        for (int i = 1; i < procNum; i++) {
            MPI_Recv(&recv[0], sizeGraph / procNum,
            MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            dist.insert(dist.end(), recv.begin(), recv.end());
        }
    } else {
        for (int i = 1; i < procNum; i++) {
            if (procRank == i) {
                MPI_Send(&dist[sizeGraph % procNum + procRank * (sizeGraph / procNum)], sizeGraph / procNum,
                MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
        }
    }
    MPI_Bcast(&dist[0], dist.size(), MPI_INT, 0, MPI_COMM_WORLD);
    }
    return dist;
}
