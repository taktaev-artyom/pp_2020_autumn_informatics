// Copyright 2020 Loganov Andrei
#ifndef  MODULES_TASK_3_LOGANOV_A_RADIX_SORT_RADIX_H_
#define  MODULES_TASK_3_LOGANOV_A_RADIX_SORT_RADIX_H_
#include <vector>
std::vector<double> EvenSpliter(std::vector<double> vec1, std::vector<double> vec2, std::vector<double> res);
std::vector<double> OddSpliter(std::vector<double> vec1, std::vector<double> vec2, std::vector<double> res);
std::vector<double> simpmerg(std::vector<double> res, std::vector<double> chet, std::vector<double> nchet);
std::vector<double> getRandomVector(int size);
bool is_exp_of_2(int n);
int countBelowPoint(double x);
int countBeforepoint(int x);
int numbyrank(int dis, double x);
std::vector<double> countingsort2(const std::vector<double>& res, int razrad);
std::vector<double> seqRadixSort(const std::vector<double> tmp);
std:: vector<double> ParallelSort(std::vector<double> vec);
#endif  //  MODULES_TASK_3_LOGANOV_A_RADIX_SORT_RADIX_H_

