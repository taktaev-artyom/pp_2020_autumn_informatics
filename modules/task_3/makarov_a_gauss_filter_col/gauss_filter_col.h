// Copyright 2020 Makarov Alexander
#ifndef MODULES_TASK_3_MAKAROV_A_GAUSS_FILTER_COL_GAUSS_FILTER_COL_H_
#define MODULES_TASK_3_MAKAROV_A_GAUSS_FILTER_COL_GAUSS_FILTER_COL_H_

#include <vector>

std::vector<unsigned short> generate_image(unsigned int w, unsigned int h);

std::vector<double> createGaussianKernel(double sigma, unsigned int radius);

std::vector<unsigned short> gaussFilter(
                               const std::vector<unsigned short>& exp_image,
                               unsigned int w, unsigned int h, double sigma,
                               unsigned int radius);

std::vector<unsigned short> transpose(
                                   const std::vector<unsigned short>& image,
                                   unsigned int w, unsigned int h);

std::vector<unsigned short> expand(
                                   const std::vector<unsigned short>& image,
                                   unsigned int w, unsigned int h,
                                   unsigned int radius);

std::vector<unsigned short> gaussFilterSequential(
                                   const std::vector<unsigned short>& image,
                                   unsigned int w, unsigned int h,
                                   double sigma, unsigned int radius);
std::vector<unsigned short> gaussFilterParallel(
                                   const std::vector<unsigned short>& image,
                                   unsigned int w, unsigned int h,
                                   double sigma, unsigned int radius);

#endif  // MODULES_TASK_3_MAKAROV_A_GAUSS_FILTER_COL_GAUSS_FILTER_COL_H_
