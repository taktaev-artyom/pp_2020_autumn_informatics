// Copyright 2020 Luckyanchenko Ivan
#include <mpi.h>
#include <iostream>
#include <numeric>
#include <cmath>
#include "../../../modules/task_3/luckianchenko_i_integral_simpson/integral_simpson.h"

double func1(double x, double y, double z) {
    return  x+y+z;
}
double func2(double x, double y, double z) {
    return x/(x+y+z*y-x*y+x*z);
}
double func3(double x, double y, double z) {
    return cos(x) / (y +z);
}

double get_Integral(double(*f)(double, double, double),
    double ax, double bx,
    double ay, double by,
    double az, double bz,
    int n) {
    double hx = (bx - ax) / n;
    double hy = (by - ay) / n;
    double ans = 0;
    for (int i = 1 ; i < n ; i ++) {
        for (int j = 1 ; j < n ; j ++) {
            double X = ax + i * hx;
            double Y = ay + j * hy;
            if ( ( i % 2 == 0 ) && ( j % 2 == 0 ) )
                ans += 4 * formula_simpson(f, X, Y, az, bz, n);
            if ( ( i % 2 != 0 ) && ( j % 2 != 0 ) )
                ans+=16 * formula_simpson(f, X, Y, az, bz, n);
            if ( ( i % 2 != 0 ) && ( j % 2 == 0 ) )
                ans+= 8 * formula_simpson(f, X, Y, az, bz, n);
            if ( ( i % 2 == 0 ) && ( j % 2 != 0 ) )
                ans += 8 * formula_simpson(f, X, Y, az, bz, n);
        }
    }
    // суммы для первого и последнего икса
    for (int i = 0; i < 2; i++) {
        for (int j = 1; j < n; j++) {
            double Y = ay + j * hy;
            double X = ax + i * (bx - ax);
            if (j % 2 == 0) {
                ans += 2 * formula_simpson(f, X, Y, az, bz, n);
            } else {
                ans += 4 * formula_simpson(f, X, Y, az, bz, n);
            }
        }
    }
    ans += formula_simpson(f, ax, ay, az, bz, n) + formula_simpson(f, bx, by, az, bz, n);
    ans *= (hy*hx)/9;
    return ans;
}

double formula_simpson(double (*f)(double, double, double),
    double X, double Y,
    double a, double b,
    int n) {
    double res = 0;
    double h = (b - a) / n;
    for (int i = 1; i < n; i++) {
        if (i % 2 == 0)
            res += 2 * f(X, Y, a + i * h);
        else
            res += 4 * f(X, Y, a + i * h);
    }
    res += f(X, Y, a) + f(X, Y, b);
    res *= h / 3;
    return res;
}

double get_Paral_Integral( double (*f)(double, double, double),
    double ax, double bx,
    double ay, double by,
    double az, double bz,
    int n) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double ans = 0;
    double integ = 0;
    double X, Y;
    double hx = (bx - ax) / n;
    double hy = (by - ay) / n;
    if (size > 1) {
        if (rank == 0) {
            for (int i = 0; i < 2; i++) {
                for (int j = 1; j < n; j++) {
                    double Y = ay + j * hy;
                    double X = ax + i * (bx - ax);
                    if (j % 2 == 0) {
                        integ += 2 * formula_simpson(f, X, Y, az, bz, n);
                    } else {
                        integ += 4 * formula_simpson(f, X, Y, az, bz, n);
                    }
                }
            }
            integ+= formula_simpson(f, ax, ay, az, bz, n)+ formula_simpson(f, bx, by, az, bz, n);
        }
        if (rank != 0) {
            for (int i = rank; i < n; i += size-1) {
                for (int j = 1; j < n; j++) {
                    Y = ay+j*hy;
                    X = ax+i*hx;
                    if ((i % 2 == 0) && (j % 2 == 0))
                        integ+=4*formula_simpson(f, X, Y, az, bz, n);
                    if ((i % 2 != 0) && (j % 2 != 0)) {
                        integ+=16*formula_simpson(f, X, Y, az, bz, n);
                    }
                    if ((i % 2 != 0) && (j % 2 == 0)) {
                        integ+=8*formula_simpson(f, X, Y, az, bz, n);
                    }
                    if ((i % 2 == 0) && (j % 2 != 0))
                        integ+=8*formula_simpson(f, X, Y, az, bz, n);
                }
            }
        }
        MPI_Reduce(&integ, &ans, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        ans *= (hx*hy)/9;
    } else {
        ans = get_Integral(f, ax, bx, ay, by, az, bz, n);
    }
    return ans;
}
