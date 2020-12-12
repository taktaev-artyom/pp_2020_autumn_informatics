// Copyright Hassan EzzAldeen 2020
#ifndef MODULES_TASK_3_HASSAN_QUICK_SORT_EVEN_ODD_QUICK_SORT_EVEN_ODD_H_
#define MODULES_TASK_3_HASSAN_QUICK_SORT_EVEN_ODD_QUICK_SORT_EVEN_ODD_H_

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

#endif  // MODULES_TASK_3_HASSAN_QUICK_SORT_EVEN_ODD_QUICK_SORT_EVEN_ODD_H_
