#include "particlegrid.h"

namespace FLIPCore{

/*****************************\
 *  Class: Particle
\*****************************/

Particle::Particle()
{
    m_p = Vector3d(0,0,0);
    m_u = Vector3d(0,0,0);
    m_mass = 0;
    m_type = AIR;
}

Particle::Particle(const Vector3d& p, const Vector3d& u, const double& m, const geomtype& t)
{
    m_p = p;
    m_u = u;
    m_mass = m;
    m_type = t;
}

Particle::~Particle()
{

}

Vector3d& Particle::position()
{
    return m_p;
}

Vector3d& Particle::velocity()
{
    return m_u;
}

double& Particle::mass()
{
    return m_mass;
}

int& Particle::type()
{
    return m_type;
}

/*****************************\
 *  Class: ParticleGrid
\*****************************/

ParticleGrid::ParticleGrid()
{
    
}

ParticleGrid::ParticleGrid(Grid<int>* grid, Grid<int>* tg) : 
    m_base_grid(grid), m_type(tg){
    
}

ParticleGrid::~ParticleGrid()
{

}

int ParticleGrid::particleSize()
{
    return m_particles.size();
}

void ParticleGrid::reserveParticle(int k)
{
    m_particles.reserve(k);
}

void ParticleGrid::insertParticle(Particle p)
{
    m_particles.push_back(p);
}

Particle& ParticleGrid::particleRef(int k)
{
    return m_particles[k];
}

std::vector<int>& ParticleGrid::cell(int k)
{
    return m_cell[k];
}

void ParticleGrid::updateTypeGrid()
{
    // assign voxels that are not solid with air
    // then use particle information to assign fluid voxel
    for(int i=0; i<m_base_grid->dimension().x(); ++i){
        for(int j=0; j<m_base_grid->dimension().y(); ++j){
            for(int k=0; k<m_base_grid->dimension().z(); ++k){
                if(m_type->cellRef(i,j,k) != SOLID){
                    m_type->cellRef(i,j,k) = AIR;
                }
            }
        }
    }

    for(int i=0; i<m_particles.size(); ++i){
        Vector3d pos = m_particles[i].position();
        int pi = static_cast<int>(std::floor(pos.x()/m_base_grid->gridSpacing()));
        int pj = static_cast<int>(std::floor(pos.y()/m_base_grid->gridSpacing()));
        int pk = static_cast<int>(std::floor(pos.z()/m_base_grid->gridSpacing()));
    
        //assert(m_base_grid->cellIndexValid(Vector3i(pi,pj,pk)) && "Particle outside domain ...");
        //assert(m_type->cellRef(pi,pj,pk) != SOLID && "Error: Particle in solid ...");
        CheckPoint(m_base_grid->cellIndexValid(Vector3i(pi,pj,pk)), "<ParticleGrid::updateTypeGrid> Particle outside domain ...");
        CheckPoint(m_type->cellRef(pi,pj,pk) != SOLID, "<ParticleGrid::updateTypeGrid> Particle in solid ...");

        m_type->cellRef(pi,pj,pk) = FLUID;
    }
}

void ParticleGrid::buildCellParticles(){
    m_cell.clear();

    for(int i=0; i<m_particles.size(); ++i){
        int x = floor(m_particles[i].position().x()/m_base_grid->gridSpacing());
        int y = floor(m_particles[i].position().y()/m_base_grid->gridSpacing());
        int z = floor(m_particles[i].position().z()/m_base_grid->gridSpacing());
        int idx = x 
                + y*m_base_grid->dimension().x() 
                + z*m_base_grid->dimension().x()*m_base_grid->dimension().y();
        
        if(m_cell.find(idx) == m_cell.end()){
            std::vector<int> vp;
            m_cell[idx] = vp;
            m_cell[idx].push_back(i);
        }
        else{
            m_cell[idx].push_back(i);
        }
    }
}

void ParticleGrid::updateParticleGridInfo(){
    // assign voxels that are not solid with air
    // then use particle information to assign fluid voxel
    //std::cout<<"Step 1 ..."<<std::endl;
    for(int i=0; i<m_base_grid->dimension().x(); ++i){
        for(int j=0; j<m_base_grid->dimension().y(); ++j){
            for(int k=0; k<m_base_grid->dimension().z(); ++k){
                if(m_type->cellRef(i,j,k) != SOLID){
                    m_type->cellRef(i,j,k) = AIR;
                }
            }
        }
    }

    // assign fluid voxels 
    // at the same time, build m_cell map
    //std::cout<<"Step 2 ..."<<std::endl;
    m_cell.clear();
    for(int i=0; i<m_particles.size(); ++i){
        Vector3d& pos = m_particles[i].position();
        int pi = static_cast<int>(std::floor(pos.x()/m_base_grid->gridSpacing()));
        int pj = static_cast<int>(std::floor(pos.y()/m_base_grid->gridSpacing()));
        int pk = static_cast<int>(std::floor(pos.z()/m_base_grid->gridSpacing()));
    
        //assert(m_base_grid->cellIndexValid(Vector3i(pi,pj,pk)) && "Particle outside domain ...");
        //assert(m_type->cellRef(pi,pj,pk) != SOLID && "Error: Particle in solid ...");
        CheckPoint(m_base_grid->cellIndexValid(Vector3i(pi,pj,pk)), "<ParticleGrid::updateParticleGridInfo> Particle outside domain ...");
        CheckPoint(m_type->cellRef(pi,pj,pk) != SOLID, "<ParticleGrid::updateParticleGridInfo> Particle in solid ...");

        m_type->cellRef(pi,pj,pk) = FLUID;

        int idx = pi 
                + pj*m_base_grid->dimension().x() 
                + pk*m_base_grid->dimension().x()*m_base_grid->dimension().y();
        
        if(m_cell.find(idx) == m_cell.end()){
            std::vector<int> vp;
            m_cell[idx] = vp;
            m_cell[idx].push_back(i);
        }
        else{
            m_cell[idx].push_back(i);
        }
    }

}

void ParticleGrid::writeParticles(std::string sfx)
{
    std::string fn = "particles/particles_"+sfx+".vtk";
    std::ofstream out(fn.c_str());

    out<<"# vtk DataFile Version 3.0"<<std::endl;
    out<<"Particles"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET UNSTRUCTURED_GRID"<<std::endl;
    out<<std::endl;

    out<<"POINTS "<<m_particles.size()<<" double"<<std::endl;
    for(int i=0; i<m_particles.size(); ++i){
        out<<m_particles[i].position().x()<<" "
            <<m_particles[i].position().y()<<" "
            <<m_particles[i].position().z()<<std::endl;
    }
    out<<std::endl;

    out<<"CELLS "<<m_particles.size()<<" "<<2*m_particles.size()<<std::endl;
    for(int i=0; i<m_particles.size(); ++i){
        out<<"1 "<<i<<std::endl;
    }

    out<<"CELL_TYPES "<<m_particles.size()<<std::endl;
    for(int i=0; i<m_particles.size(); ++i){
        out<<"1"<<std::endl;
    }

    out.close();
}

} // namespace 