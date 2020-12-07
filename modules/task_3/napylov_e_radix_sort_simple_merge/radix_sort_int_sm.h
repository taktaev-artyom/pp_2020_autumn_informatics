// Copyright 2020 Napylov Evgenii
#ifndef MODULES_TASK_3_NAPYLOV_E_RADIX_SORT_SIMPLE_MERGE_RADIX_SORT_INT_SM_H_
#define MODULES_TASK_3_NAPYLOV_E_RADIX_SORT_SIMPLE_MERGE_RADIX_SORT_INT_SM_H_
#include <vector>

std::vector<int> RandomVector(int len);

std::vector<int> mergeVectors(std::vector<int> vec1, std::vector<int> vec2);

int getMaxDigitCount(std::vector<int> vec);

std::vector<int> RadixSort(std::vector<int> vec);

std::vector<int> RadixSortParallel(std::vector<int> vec);

int compare(const void *a, const void *b);

#endif  // MODULES_TASK_3_NAPYLOV_E_RADIX_SORT_SIMPLE_MERGE_RADIX_SORT_INT_SM_H_
