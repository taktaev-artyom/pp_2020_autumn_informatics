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

int** calc_sobel(int** image, int rows, int cols) {
    if ((rows < 3) || (cols < 3))
        return image;
    int G = 0;
    int Gx;
    int Gy;
    int* vectorx;
    int* vectory;
    int** temp_image = new int* [rows];
    for (int i = 0; i < rows; i++) {
        temp_image[i] = new int[cols];
    }
    vectorx = new int[3 * 3];
    vectorx[0] = -1; vectorx[3] = -2; vectorx[6] = -1;
    vectorx[1] = 0; vectorx[4] = 0; vectorx[7] = 0;
    vectorx[2] = 1; vectorx[5] = 2; vectorx[8] = 1;
    vectory = new int[3 * 3];
    vectory[0] = -1; vectory[3] = 0; vectory[6] = 1;
    vectory[1] = -2; vectory[4] = 0; vectory[7] = 2;
    vectory[2] = -1; vectory[5] = 0; vectory[8] = 1;
    for (int y = 0; y < cols; y++) {
        for (int x = 0; x < rows; x++) {
            Gx = 0;
            Gy = 0;
            if (y == 0 || y == cols)
                G = 0;
            else if (x == 0 || x == rows)
                G = 0;
            else {
                for (int i = -1; i <= 1; i++)
                    for (int j = -1; j <= 1; j++) {
                        int idx = (i + 1) * 3 + j + 1;
                        int pix = image[clamp(x + j, rows - 1, 0)][clamp(y + i, cols - 1, 0)];
                        Gx += pix * vectorx[idx];
                        Gy += pix * vectory[idx];
                        G = sqrt(pow(Gx, 2) + pow(Gy, 2));
                        temp_image[x][y] = G;
                    }
            }
            G = clamp(G, 255, 0);
        }
    }
    return temp_image;
}
