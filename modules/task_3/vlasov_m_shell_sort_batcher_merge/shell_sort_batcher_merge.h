// Copyright 2020 Vlasov Maksim
#ifndef MODULES_TASK_2_VLASOV_M_SHELL_SORT_BATCHER_MERGE_H_
#define MODULES_TASK_2_VLASOV_M_SHELL_SORT_BATCHER_MERGE_H_
#include <vector>

using Vector = std::vector<int>;

Vector createRandomVector(int size);

Vector shellSort(Vector arr);

namespace BatcherMerge {
    Vector parallelSort(Vector arr, std::function<Vector(Vector)> sort_func);
}  // namespace BatcherMerge

#endif  // MODULES_TASK_2_VLASOV_M_SHELL_SORT_BATCHER_MERGE_H_
