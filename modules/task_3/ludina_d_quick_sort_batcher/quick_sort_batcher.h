// Copyright 2020 Ludina Daria
#ifndef MODULES_TASK_3_LUDINA_D_QUICK_SORT_BATCHER_QUICK_SORT_BATCHER_H_
#define MODULES_TASK_3_LUDINA_D_QUICK_SORT_BATCHER_QUICK_SORT_BATCHER_H_

#include <mpi.h>
#include <random>
#include <ctime>
#include <vector>
#include <utility>
#include <algorithm>

std::vector<int> createRandomVector(int size_v);
int part(std::vector<int>* vec, int first, int last);
void quickSort(std::vector<int>* vec, int first, int last);
void merge(std::vector<int> vec_up, std::vector<int> vec_down);
void network(std::vector<int> procs);
void quickSortBatcher(std::vector<int>* vec);

#endif  // MODULES_TASK_3_LUDINA_D_QUICK_SORT_BATCHER_QUICK_SORT_BATCHER_H_
