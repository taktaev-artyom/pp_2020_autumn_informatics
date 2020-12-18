// Copyright 2020 Makarov Alexander
#ifndef MODULES_TASK_3_MAKAROV_A_GAUSS_FILTER_COL_GAUSS_FILTER_COL_H_
#define MODULES_TASK_3_MAKAROV_A_GAUSS_FILTER_COL_GAUSS_FILTER_COL_H_

#include <vector>

std::vector<unsigned __int16> generate_image(unsigned int w, unsigned int h);

std::vector<double> createGaussianKernel(double sigma, unsigned int radius);

std::vector<unsigned __int16> gaussFilter(
                               const std::vector<unsigned __int16>& exp_image,
                               unsigned int w, unsigned int h, double sigma,
                               unsigned int radius);

std::vector<unsigned __int16> transpose(
                                   const std::vector<unsigned __int16>& image,
                                   unsigned int w, unsigned int h);

std::vector<unsigned __int16> expand(
                                   const std::vector<unsigned __int16>& image,
                                   unsigned int w, unsigned int h,
                                   unsigned int radius);

std::vector<unsigned __int16> gaussFilterSequential(
                                   const std::vector<unsigned __int16>& image,
                                   unsigned int w, unsigned int h,
                                   double sigma, unsigned int radius);
std::vector<unsigned __int16> gaussFilterParallel(
                                   const std::vector<unsigned __int16>& image,
                                   unsigned int w, unsigned int h,
                                   double sigma, unsigned int radius);

#endif  // MODULES_TASK_3_MAKAROV_A_GAUSS_FILTER_COL_GAUSS_FILTER_COL_H_
