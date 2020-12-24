// Copyright 2020 Pestreev Daniil
#ifndef MODULES_TASK_3_PESTREEV_D_QUICK_SORT_EVEN_ODD_QUICK_SORT_EVEN_ODD_MERGE_H_
#define MODULES_TASK_3_PESTREEV_D_QUICK_SORT_EVEN_ODD_QUICK_SORT_EVEN_ODD_MERGE_H_

#include <vector>
#include <string>

std::vector<int> getRandomVector(int size);
void Swap(int a, int b);
void qsort(int* vec, int left, int right);
std::vector<int> quickSortV(const std::vector<int>& vec);
void recur_merge(const std::vector<int>& left, const std::vector<int>& right);
void networking(const std::vector<int>& arr);
void batchers_network(int proc_size);
std::vector<int> parallel_sorting(const std::vector<int>& vec);


#endif  // MODULES_TASK_3_PESTREEV_D_QUICK_SORT_EVEN_ODD_QUICK_SORT_EVEN_ODD_MERGE_H_
