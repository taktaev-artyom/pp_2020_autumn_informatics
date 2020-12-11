// Copyright 2020 Kirichenko Nikita
#ifndef MODULES_TASK_3_KIRICHENKO_N_CANNON_CANNON_H_
#define MODULES_TASK_3_KIRICHENKO_N_CANNON_CANNON_H_

#include <vector>

std::vector<double> GetRandomMatrix(const size_t size);

std::vector<double> MatrixMultiplication(const std::vector<double>& a, const std::vector<double>& b, const size_t size);
std::vector<double> Cannon(const std::vector<double>& a, const std::vector<double>& b, const size_t size);

#endif  // MODULES_TASK_3_KIRICHENKO_N_CANNON_CANNON_H_#pragma once
