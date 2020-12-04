// Copyright 2020 Prokofeva Elizaveta
#ifndef MODULES_TASK_3_PROKOFEVA_E_OPERATOR_SOBEL_OPERATOR_SOBEL_H_
#define MODULES_TASK_3_PROKOFEVA_E_OPERATOR_SOBEL_OPERATOR_SOBEL_H_

int clamp(int v, int max, int min);
int** calc_sobel(int** image, int rows, int cols);

#endif  // MODULES_TASK_3_PROKOFEVA_E_OPERATOR_SOBEL_OPERATOR_SOBEL_H_
