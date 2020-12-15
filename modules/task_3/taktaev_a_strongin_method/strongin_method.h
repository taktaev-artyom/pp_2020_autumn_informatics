// Copyright 2020 Taktaev Artem
#ifndef MODULES_TASK_3_TAKTAEV_A_STRONGIN_METHOD_STRONGIN_METHOD_H_
#define MODULES_TASK_3_TAKTAEV_A_STRONGIN_METHOD_STRONGIN_METHOD_H_
#include <vector>
#include <functional>

class Strongin {
 public:
    std::vector<double> x, z;
    int n;
    double a, b;
    double r;
    int prec;
    std::function<double(double)> f;

    Strongin(double _a, double _b, double _r, int prec, std::function<double(double)> _f);
    void addNode(double xx);
    double calculateM();
    std::vector<double> calculateR(double m);

    double seqStronginSearch();
    double parStronginSearch();
};

#endif  // MODULES_TASK_3_TAKTAEV_A_STRONGIN_METHOD_STRONGIN_METHOD_H_
