#ifndef __KERNEL_H__
#define __KERNEL_H__

#include "utility.h"

namespace FLIPCore{

extern void CheckPoint(bool cls, std::string msg);

extern double LinearKernel(double l);
extern double QuadraticKernel(double l);
extern double CubicKernel(double l);

extern double UniformRandom(double lb, double ub);

}

#endif
