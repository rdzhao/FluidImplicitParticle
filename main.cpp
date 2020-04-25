#include "simulation/flip.h"

using namespace FLIPCore;

int main(int argc, char* argv[]){
    Vector3i dimension = Vector3i(20,20,20);
    FLIPSim flip(dimension, 12, 0.1);

    flip.init();
    flip.simulate();
    

    return 1;
}