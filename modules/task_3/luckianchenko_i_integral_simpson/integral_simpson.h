#ifndef  MODULES_TASK_3_LUCKIANCHENKO_I_INTEGRAL_SIMPSON_INTEGRAL_SIMPSON_H_
#define  MODULES_TASK_3_LUCKIANCHENKO_I_INTEGRAL_SIMPSON_INTEGRAL_SIMPSON_H_
// Copyright 2020 Luckyanchenko Ivan
double func1(double x, double y, double z);
double func2(double x, double y, double z);
double func3(double x, double y, double z);
double get_Integral(double(*f)(double, double, double),
double ax, double bx,
double ay, double by,
double az, double bz, int n);
double formula_simpson(double (*f)(double, double, double), double X, double Y, double a, double b, int n);
double get_Paral_Integral( double (*f)(double, double, double),
                     double ax, double bx,
                     double ay, double by,
                     double az, double bz, int n);

#endif  //   MODULES_TASK_3_LUCKIANCHENKO_I_INTEGRAL_SIMPSON_INTEGRAL_SIMPSON_H_
