// Copyright 2020 Solovev Aleksandr
#ifndef MODULES_TASK_3_SOLOVEV_A_BUILDING_IMAGE
#define MODULES_TASK_3_SOLOVEV_A_BUILDING_IMAGE
#include <mpi.h>
#include <vector>
#include <iostream>
#include <ctime>
#include <random>

struct Point {
    int x, y;
    Point() = default;
    Point(int X, int Y) {
        x = X;
        y = Y;
    }
	 const Point operator=(const Point& point) {
        x = point.x;
        y = point.y;
        return *this;
    }

    bool operator!=(const Point& point) {
        if ((x != point.x)||(y != point.y)){
            return true;
        }
        else{
            return false;
        }
        
    }

    bool operator==(const Point& point) {
        if ((x == point.x)&&(y == point.y)){
            return true;
        }
        else{
            return false;
        }
    }
};

std::vector<Point> interpriate_basic(int ** image, int height, int width);
int searchFirstPoint(const std::vector<Point> &mas);
Point searchFirstPointParallel(const std::vector<Point> &mas, int rank, int procCount);
std::vector<Point> buildConvexHullParallel(const std::vector<Point> &mas);
std::vector<Point> buildConvexHull(const std::vector<Point> &mas);
bool pointCheck(const Point& back, const Point& current, const Point& challenger);
#endif MODULES_TASK_3_SOLOVEV_A_BUILDING_IMAGE