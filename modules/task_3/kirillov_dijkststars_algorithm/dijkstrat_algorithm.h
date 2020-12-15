#pragma once
// Copyright 2020 Kirillov Konstantin
#ifndef MODULES_TASK_3_KIRILLOV_DIJKSTSTARS_ALGORITHM_DIJKSTRAT_ALGORITHM_H_
#define MODULES_TASK_3_KIRILLOV_DIJKSTSTARS_ALGORITHM_DIJKSTRAT_ALGORITHM_H_
#include <vector>
std::vector<int> getRandomGraph(int sizeGraph);
void printGraph(std::vector<int>graph);
std::vector<int> getSequentialDijkstras(std::vector<int>graph, int start);
void printDist(std::vector<int>dist);
std::vector<int> getParallelDijkstras(std::vector<int>graph, int start);
#endif  // MODULES_TASK_3_KIRILLOV_DIJKSTSTARS_ALGORITHM_DIJKSTRAT_ALGORITHM_H_
