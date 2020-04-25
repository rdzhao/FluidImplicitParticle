#include "kernel.h"

namespace FLIPCore{

void CheckPoint(bool cls, std::string msg){
    if(cls == false){
        std::cout<<"Error: "+msg<<std::endl;
        std::abort();
    }

}

double LinearKernel(double l)
{
    if(abs(l)<1){
        return 1-abs(l);
    }
    else{
        return 0;
    }
}

double QuadraticKernel(double l)
{
    if(abs(l)<0.5){
        return 0.75 - l*l;
    }
    else if(abs(l)>=0.5 && abs(l)<1.5){
        return 0.5*pow((1.5 - abs(l)),2);
    }
    else{
        return 0;
    }
}

double CubicKernel(double l){
    if(abs(l)<1){
        return 0.5*pow(abs(l),3) - pow(abs(l),2) + 2.0/3.0;
    }
    else if(abs(l)>=1 && abs(l)<2){
        return (1.0/6.0)*pow((2 - abs(l)),3);
    }
    else{
        return 0;
    }
}

double UniformRandom(double lb, double ub){
    std::random_device dev;
    std::mt19937_64 eng(dev());
    std::uniform_real_distribution<double> distribution(lb, ub);
    return distribution(eng);
}

}