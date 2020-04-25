#ifndef __PARTICLEGRID_H__
#define __PARTICLEGRID_H__

#include "grid.h"

namespace FLIPCore{

class Particle{
    public:
    Particle();
    Particle(const Vector3d& p, const Vector3d& u, const double& m, const geomtype& t);
    ~Particle();

    Vector3d& position();
    Vector3d& velocity();

    double& mass();
    int& type(); 

    protected:
    Vector3d m_p; // position
    Vector3d m_u; // velocity
    double m_mass;
    int m_type;
};

class ParticleGrid{
    public:
    ParticleGrid();
    ParticleGrid(Grid<int>* grid, Grid<int>* tg); 
    ~ParticleGrid();

    int particleSize();
    void reserveParticle(int k);
    void insertParticle(Particle p);
    Particle& particleRef(int k);
    std::vector<int>& cell(int k);

    void updateTypeGrid();
    void buildCellParticles();
    void updateParticleGridInfo();

    void writeParticles(std::string sfx="");

    protected:
    Grid<int>* m_base_grid;
    std::vector<Particle> m_particles;

    // better not pointer
    Grid<int>* m_type; 

    // given grid index, map to confined particles
    std::map<int, std::vector<int> > m_cell;
    //              ^^^
    // better using list iterator.
};

}

#endif