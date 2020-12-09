// Copyright 2020 Prokofeva Elizaveta
#ifndef MODULES_TASK_3_PROKOFEVA_E_OPERATOR_SOBEL_OPERATOR_SOBEL_H_
#define MODULES_TASK_3_PROKOFEVA_E_OPERATOR_SOBEL_OPERATOR_SOBEL_H_
using matrix = std::vector<std::vector<int>>;

int clamp(int v, int max, int min);
std::vector<int> getRandomMatrix(int rows, int cols);
std::vector<int> calc_sobel(const std::vector<int> image, int rows, int cols);
int calc_kernel( std::vector<int>& image, int x, int y, int rows, int cols);
std::vector<int> parallel_sobel(std::vector<int>  image, int rows, int cols);

#endif  // MODULES_TASK_3_PROKOFEVA_E_OPERATOR_SOBEL_OPERATOR_SOBEL_H_
