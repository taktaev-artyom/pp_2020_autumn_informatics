// Copyright 2020 Makarov Alexander
#ifndef MODULES_TASK_2_MAKAROV_A_GAUSS_FILTER_GAUSS_FILTER_H_
#define MODULES_TASK_2_MAKAROV_A_GAUSS_FILTER_GAUSS_FILTER_H_

#include <vector>

std::vector<unsigned short> generate_image(unsigned int w, unsigned int h);

std::vector<unsigned short> gaussFilterSequential(const std::vector<unsigned short>& image, unsigned int w, unsigned int h,
                                                  double sigma, unsigned int radius);
std::vector<unsigned short> gaussFilterParallel(const std::vector<unsigned short>& image, unsigned int w, unsigned int h,
                                                double sigma, unsigned int radius);
												
std::vector<unsigned short> prepare_image(const std::vector<unsigned short>& image, unsigned int w, unsigned int h, 
                                          unsigned int radius);

#endif  // MODULES_TASK_2_MAKAROV_A_GAUSS_FILTER_GAUSS_FILTER_H_
