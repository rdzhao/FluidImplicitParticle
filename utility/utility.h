#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <cmath>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <vector>
#include <list>
#include <stack>
#include <map>

#include "Dense"
#include "Sparse"

//using namespace std;

namespace FLIPCore{

const double PI         =   3.1415926535897932384626422832795028841971;
const double TWO_PI     =   6.2831853071795864769252867665590057683943;
const double EPSILON    =   0.000000001;
const double GRAVITY    =   9.8;

typedef Eigen::Vector3d Vector3d;
typedef Eigen::Vector3i Vector3i;
typedef Eigen::VectorXd VectorXd;
typedef Eigen::SparseMatrix<double> SpMatXd;
typedef Eigen::Triplet<double> Triplet;


enum geomtype {SOLID = 2, FLUID = 1, AIR = 0};
enum traversestate {VISITED = 1, UNVISITED = 0};
enum mactype {U = 0, V = 1, W = 2};

}

#endif