// Copyright 2020 Taktaev Artem
#ifndef MODULES_TASK_1_TAKTAEV_A_VECTOR_ADJ_ALTERNATIONS_VECTOR_ADJ_ALTERNATIONS_H
#define MODULES_TASK_1_TAKTAEV_A_VECTOR_ADJ_ALTERNATIONS_VECTOR_ADJ_ALTERNATIONS_H

#include <vector>

std::vector<int> createRandomVector(int  vec_size);
int calculateAdjAlternationsParallel(const std::vector<int> &vec, int vec_size);
int calculateAdjAlternationsSequential(const std::vector<int> &vec, int inc, int start_index);

#endif  // MODULES_TASK_1_TAKTAEV_A_VECTOR_ADJ_ALTERNATIONS_VECTOR_ADJ_ALTERNATIONS_H