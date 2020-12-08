// Copyright 2020 Prokofeva Elizaveta
#include <mpi.h>
#include <algorithm>
#include <stdexcept>
#include <ctime>
#include "random"
#include "iostream"
#include "../../../modules/task_3/prokofeva_e_operator_sobel/operator_sobel.h"

int clamp(int v, int max, int min) {
    if (v > max) return max;
    else if (v < min) return min;
    return v;
}
std::vector<int> calc_sobel(const std::vector<int> image, int rows, int cols)
{
    if ((rows < 3) || (cols < 3))
        return image;
    std::vector<int> result_image(rows * cols);
    std::vector<int> vectorx(3 * 3);
    std::vector<int> vectory(3 * 3);
    vectorx = { -1, 0, 1,
                -2, 0, 2 ,
                -1, 0, 1 };

    vectory = { -1, -2, -1,
                 0, 0, 0,
                 1, 2, 1 };
    for (int x = 0; x < rows; x++) {
        for (int y = 0; y < cols; y++) {
            int Gx = 0;
            int Gy = 0;
            int G = 0;
            if (y == 0 || y == cols)
                G = 0;
            else if (x == 0 || x == rows)
                G = 0;
            else {
                for (int i = -1; i <= 1; i++) {
                    for (int j = -1; j <= 1; j++) {
                        int idx = (i + 1) * 3 + j + 1;
                        int x1 = clamp(x + j, rows - 1, 0);
                        int y1 = clamp(y + i, cols - 1, 0);
                        Gx += image[x1 * cols + y1] * vectorx[idx];
                        Gy += image[x1 * cols + y1] * vectory[idx];
                        G = sqrt(pow(Gx, 2) + pow(Gy, 2));
                    }
                }
            }
            G = clamp(G, 255, 0);
            result_image[x * cols + y] = G;
        }
    }
    return result_image;
}

std::vector<int> getRandomMatrix(int rows, int cols) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> matrix(rows*cols);
    for (int i = 0; i < rows * cols; i++)
        matrix[i] = gen() % 255;
    return matrix;
}

