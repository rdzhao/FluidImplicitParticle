#ifndef __FLIP_H__
#define __FLIP_H__

#include "../grid/macgrid.h"
#include "../grid/particlegrid.h"

namespace FLIPCore{

class FLIPSim{
    public:
    FLIPSim();
    FLIPSim(Vector3i dimension, double et, double dt);
    ~FLIPSim();

    void init();
    void simulate();

    protected:
    void initGrid();
    void initParticleGrid();
    void initMACGrid();

    void initTypeGrid();
    void sampleParticles();

    protected:
    void simulateStep(int i);
    void prepareParticleGrid();
    void velocityParticleToGrid();
    void applyExternalForce(); // gravity
    void applyPressureProjection();
    void extrapolateVelocities();
    void advectParticles();
    void adjusteParticles();
    void enforceBoundaryCondition();

    void assembleLaplacianOperator();
    VectorXd computeVelocityDivergence();
    void extrapolateGrid(mactype type);
    void integrationRungeKutta4(Particle& p);
    void subtractMACGrid();
    Vector3d interpolateVelocity(Vector3d p, bool isPIC);
    void adjustPositionToDomain(Vector3d& p);
    
    bool cellIndexValid(Vector3i idx);
    bool macIndexValid(Vector3i idx, mactype type);
    bool isSolid(Vector3i idx);

    // debug functions
    void checkVelocityFieldDivergence();
    void checkParticleVelocity();
    void checkKernels();

    protected:
    double m_delta_t;
    double m_flip_ratio;
    double m_elapsed_time;


    Vector3i m_dimension;
    
    Grid<int> m_grid; // base grid
    Grid<int> m_type_grid; // cell type
    ParticleGrid m_particle_grid;
    MACGrid m_mac_grid;
    MACGrid m_mac_grid_p; // previous macgrid
    MACGrid m_mac_grid_d; // diff macgrid

    SpMatXd m_laplacian_op; // gradient operator
};

}

#endif