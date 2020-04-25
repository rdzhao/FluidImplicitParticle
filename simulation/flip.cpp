#include "flip.h"

using namespace FLIPCore;

FLIPSim::FLIPSim(){

}

FLIPSim::FLIPSim(Vector3i dimension, double et, double dt) :
    m_dimension(dimension), m_elapsed_time(et), m_delta_t(dt){
        m_flip_ratio = 0.95;
}

FLIPSim::~FLIPSim(){

}

void FLIPSim::init(){
    std::cout<<"Initialize base grid ..."<<std::endl;
    initGrid();
    std::cout<<"Initialize particle grid ..."<<std::endl;
    initParticleGrid();
    std::cout<<"Initialize MAC grid"<<std::endl;
    initMACGrid();

    // debug
    //checkKernels();
}

void FLIPSim::simulate(){
    int steps = m_elapsed_time / m_delta_t;
    std::cout<<m_elapsed_time<<" "<<m_delta_t<<std::endl;
    std::cout<<"Total steps: "<<steps<<std::endl;
    for(int i=0; i<steps; ++i){
        std::cout<<"Simulation step: "<<i<<std::endl;
        simulateStep(i);
    }

    // debug
    //checkParticleVelocity();
}

void FLIPSim::initGrid(){
    m_grid = Grid<int>(m_dimension, 1.0, 0);

    m_grid.writeGrid("test");
}

void FLIPSim::initParticleGrid(){
    m_particle_grid = ParticleGrid(&m_grid, &m_type_grid);
    
    initTypeGrid();
    sampleParticles();

    //m_particle_grid.writeParticles("test");
}

void FLIPSim::initMACGrid(){
    m_mac_grid = MACGrid(&m_grid);
    m_mac_grid_p = MACGrid(&m_grid);
    m_mac_grid_d = MACGrid(&m_grid);

    //std::cout<<"Assign cell based velocity ..."<<std::endl;
    //m_mac_grid.convertMACToVertexBasedVelocity();
    //std::cout<<"Write velocity grid ..."<<std::endl;
    //m_mac_grid.writeVelocityGrid("test");
}

void FLIPSim::initTypeGrid(){
    // hard code solid, fluid and air cells;
    m_type_grid = Grid<int>(m_dimension, 0, 1.0); 
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                // x between 0.3-0.7
                // y all
                // z between 0-0.5
                if(i*1.0/m_dimension.x()>-0.1 
                    && i*1.0/m_dimension.x()<0.4
                    //&& j*1.0/m_dimension.y()>0.2
                    //&& j*1.0/m_dimension.y()<0.7
                    //&& k*1.0/m_dimension.z()>0.5
                    && k*1.0/m_dimension.z()<0.7){
                        m_type_grid.cellRef(i,j,k) = FLUID;
                    }
                else{
                    m_type_grid.cellRef(i,j,k) = AIR;
                }
            }
        }
    }
}

void FLIPSim::sampleParticles(){
    // sample particles in base grid
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                if(m_type_grid.cellRef(i,j,k) == FLUID){
                    Vector3d center = m_grid.getCellCenterCoord(i,j,k);
                    const double gs = m_grid.gridSpacing();
                    Vector3d corner = center - 0.5*Vector3d(gs,gs,gs);

                    for(int l=0; l<2; ++l){
                        for(int m=0; m<2; ++m){
                            for(int n=0; n<2; ++n){
                                Particle p = Particle();
                                double rx = UniformRandom(-0.2, 0.2);
                                double ry = UniformRandom(-0.2, 0.2);
                                double rz = UniformRandom(-0.2, 0.2);
                                //std::cout<<rx<<" "<<ry<<" "<<rz<<std::endl;
                                p.position() = corner 
                                    + Vector3d((0.5*l+0.25)*gs,(0.5*m+0.25)*gs,(0.5*n+0.25)*gs)
                                    + Vector3d(rx,ry,rz);                                
                                p.velocity() = Vector3d(0,0,0);
                                p.mass() = 1.0;

                                m_particle_grid.insertParticle(p);
                            }
                        }
                    }

                }
            }
        }
    }
}

void FLIPSim::simulateStep(int i){
    std::cout<<"Preparing particle grid ..."<<std::endl;
    prepareParticleGrid();
    std::cout<<"Applying external force ..."<<std::endl;
    applyExternalForce();
    //checkParticleVelocity();

    std::cout<<"Splatting particle velocity to grid ..."<<std::endl;
    velocityParticleToGrid();
    
    std::cout<<"Converting ..."<<std::endl;
    m_mac_grid.convertMACToVertexBasedVelocity();
    std::cout<<"Writing ..."<<std::endl;
    //m_mac_grid.writeVelocityGrid("noextra_"+std::to_string(i));

    std::cout<<"Extrapolate MAC velocities ..."<<std::endl;
    extrapolateVelocities();  

    //m_particle_grid.writeParticles("before"+std::to_string(i));
    //m_mac_grid.convertMACToVertexBasedVelocity();
    //m_mac_grid.writeVelocityGrid("test");
    std::cout<<"Converting ..."<<std::endl;
    //m_mac_grid.convertMACToCellBasedVelocity();
    m_mac_grid.convertMACToVertexBasedVelocity();
    //m_mac_grid.checkCellVelocity();
    std::cout<<"Writing ..."<<std::endl;
    //m_mac_grid.writeCellVelocityGrid("before"+std::to_string(i));
    //m_mac_grid.writeVelocityGrid("before"+std::to_string(i));

    std::cout<<"Enforcing boundary conditions ..."<<std::endl;
    enforceBoundaryCondition();

    std::cout<<"Applying pressure projection ..."<<std::endl;
    applyPressureProjection();
    //checkVelocityFieldDivergence();

    std::cout<<"Converting ..."<<std::endl;
    //m_mac_grid.convertMACToCellBasedVelocity();
    m_mac_grid.convertMACToVertexBasedVelocity();
    //m_mac_grid.checkCellVelocity();
    std::cout<<"Writing ..."<<std::endl;
    //m_mac_grid.writeCellVelocityGrid("after"+std::to_string(i));
    m_mac_grid.writeVelocityGrid("after"+std::to_string(i));
    //enforceBoundaryCondition();
    //extrapolateVelocities();
    std::cout<<"Advecting Particles ..."<<std::endl;
    advectParticles();
    //checkParticleVelocity();

    //std::cout<<"Adjusting Particles ..."<<std::endl;
    //adjusteParticles();

    m_particle_grid.writeParticles("after"+std::to_string(i));
}

void FLIPSim::prepareParticleGrid(){
    //m_particle_grid.buildCellParticles();
    //m_particle_grid.updateTypeGrid();
    m_particle_grid.updateParticleGridInfo();
}

void FLIPSim::velocityParticleToGrid(){
    // use quadratic kernel, 2 blocks away
    int nb = 2;

    // MAC grid ux
    std::cout<<"Splatting u to grid ..."<<std::endl;
    for(int i=0; i<m_dimension.x()+1; ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                Vector3d fc = Vector3d(i*m_grid.gridSpacing(), (j+0.5)*m_grid.gridSpacing(), (k+0.5)*m_grid.gridSpacing()); // face center
                std::vector<int> neighbors;

                for(int l=-nb; l<nb; ++l){
                    for(int m=-nb; m<nb+1; ++m){
                        for(int n=-nb; n<nb+1; ++n){
                            int pi = i+l;
                            int pj = j+m;
                            int pk = k+n;
                            int idx = pi + pj*m_dimension.x() + pk*m_dimension.x()*m_dimension.y();

                            if(cellIndexValid(Vector3i(pi, pj, pk))){
                                neighbors.insert(neighbors.end(), m_particle_grid.cell(idx).begin(), m_particle_grid.cell(idx).end());
                            }
                        }
                    }
                }

                //std::cout<<"Neighbor size: "<<neighbors.size()<<std::endl;

                double ux = 0;
                double tw = 0;
                for(int l=0; l<neighbors.size(); ++l){
                    Vector3d d = m_particle_grid.particleRef(neighbors[l]).position()-fc;
                    if(d.norm()<nb*m_grid.gridSpacing()){
                        // within range!
                        double tgs = m_grid.gridSpacing();
                        double w = QuadraticKernel(d.x()/tgs)
                                    *QuadraticKernel(d.y()/tgs)
                                    *QuadraticKernel(d.z()/tgs);
                        //std::cout<<w<<std::endl;
                        //double w = QuadraticKernel(d/(2*m_grid.gridSpacing()));
                        ux += w*m_particle_grid.particleRef(neighbors[l]).velocity().x();
                        tw += w;
                    }
                }
                if(tw == 0){
                    ux = 0;
                    m_mac_grid.gridUM().cellRef(i,j,k) = UNVISITED;
                }
                else{
                    ux /= tw; //average particles per cell.
                    //std::cout<<tw<<std::endl;
                    m_mac_grid.gridUM().cellRef(i,j,k) = VISITED;
                }
                //std::cout<<ux<<std::endl;
                m_mac_grid.gridU().cellRef(i,j,k) = ux;
            }
        }
    }

    // MAC grid uy
    std::cout<<"Splatting v to grid ..."<<std::endl;
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y()+1; ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                Vector3d fc = Vector3d((i+0.5)*m_grid.gridSpacing(), j*m_grid.gridSpacing(), (k+0.5)*m_grid.gridSpacing()); // face center
                std::vector<int> neighbors;

                for(int l=-nb; l<nb+1; ++l){
                    for(int m=-nb; m<nb; ++m){
                        for(int n=-nb; n<nb+1; ++n){
                            int pi = i+l;
                            int pj = j+m;
                            int pk = k+n;
                            int idx = pi + pj*m_dimension.x() + pk*m_dimension.x()*m_dimension.y();

                            if(cellIndexValid(Vector3i(pi,pj,pk))){
                                neighbors.insert(neighbors.end(), m_particle_grid.cell(idx).begin(), m_particle_grid.cell(idx).end());     
                            }

                        }
                    }
                }

                double uy = 0;
                double tw = 0;
                for(int l=0; l<neighbors.size(); ++l){
                    Vector3d d = m_particle_grid.particleRef(neighbors[l]).position()-fc;
                    if(d.norm()<2*m_grid.gridSpacing()){
                        double tgs = m_grid.gridSpacing();
                        double w = QuadraticKernel(d.x()/tgs)
                                    *QuadraticKernel(d.y()/tgs)
                                    *QuadraticKernel(d.z()/tgs);
                        //double w = QuadraticKernel(d/(2*m_grid.gridSpacing()));
                        uy += w*m_particle_grid.particleRef(neighbors[l]).velocity().y();
                        tw += w; 
                    }
                }
                if(tw == 0){
                    uy = 0;
                    m_mac_grid.gridVM().cellRef(i,j,k) = UNVISITED;
                }
                else{
                    uy /= tw; //average particles per cell.
                    m_mac_grid.gridVM().cellRef(i,j,k) = VISITED;
                }
                m_mac_grid.gridV().cellRef(i,j,k) = uy;
            }
        }
    }

    // MAC grid uz
    std::cout<<"Splatting w to grid .. "<<std::endl;
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z()+1; ++k){
                Vector3d fc = Vector3d((i+0.5)*m_grid.gridSpacing(), (j+0.5)*m_grid.gridSpacing(), k*m_grid.gridSpacing()); // face center
                std::vector<int> neighbors;

                for(int l=-nb; l<nb+1; ++l){
                    for(int m=-nb; m<nb+1; ++m){
                        for(int n=-nb; n<nb; ++n){
                            int pi = i+l;
                            int pj = j+m;
                            int pk = k+n;
                            int idx = pi + pj*m_dimension.x() + pk*m_dimension.x()*m_dimension.y();

                            if(cellIndexValid(Vector3i(pi,pj,pk))){
                                neighbors.insert(neighbors.end(), m_particle_grid.cell(idx).begin(), m_particle_grid.cell(idx).end());     
                            }

                        }
                    }
                }

                double uz = 0;
                double tw = 0;
                for(int l=0; l<neighbors.size(); ++l){
                    Vector3d d = m_particle_grid.particleRef(neighbors[l]).position()-fc;
                    if(d.norm()<2*m_grid.gridSpacing()){
                        double tgs = m_grid.gridSpacing();
                        double w = QuadraticKernel(d.x()/tgs)
                                    *QuadraticKernel(d.y()/tgs)
                                    *QuadraticKernel(d.z()/tgs);
                        //double w = QuadraticKernel(d/(2*m_grid.gridSpacing()));
                        uz += w*m_particle_grid.particleRef(neighbors[l]).velocity().z();
                        tw += w; 
                    }
                }

                if(tw == 0){
                    uz = 0;
                    m_mac_grid.gridWM().cellRef(i,j,k) = UNVISITED;
                }
                else{
                    uz /= tw; //average particles per cell.
                    m_mac_grid.gridWM().cellRef(i,j,k) = VISITED;
                }

                m_mac_grid.gridW().cellRef(i,j,k) = uz;
            }
        }
    }

    m_mac_grid_p = m_mac_grid;
}

void FLIPSim::applyExternalForce(){
    // assume there is only gravity as external force.
    // so we only need to update MAC grid UZ component
    //for(int i=0; i<m_dimension.x(); ++i){
    //    for(int j=0; j<m_dimension.y(); ++j){
    //        for(int k=0; k<m_dimension.z()+1; ++k){
    //            m_mac_grid.gridW().cellRef(i,j,k) += m_delta_t*(-GRAVITY);
    //        }
    //    }
    //}

    // update particle directly.
    for(int i=0; i<m_particle_grid.particleSize(); ++i){
        m_particle_grid.particleRef(i).velocity().z() += m_delta_t*(-GRAVITY);
    }
}

void FLIPSim::applyPressureProjection(){
    // vector field decomposition
    std::cout<<"Assembling Laplacian matrix ..."<<std::endl;
    assembleLaplacianOperator();
    std::cout<<"Assembling divergence vector ..."<<std::endl;
    VectorXd vel_div = computeVelocityDivergence();
    //std::cout<<vel_div(0)<<std::endl;

    std::cout<<"Poisson equation ..."<<std::endl;
    Eigen::ConjugateGradient<SpMatXd> solver;
    //Eigen::SimplicialLDLT<SpMatXd> solver;
    std::cout<<"Prefactoring Laplacian matrix ..."<<std::endl;
    solver.compute(m_delta_t*m_laplacian_op);
    //assert(solver.info() == Eigen::Success && "Compute Fail ...");
    CheckPoint(solver.info() == Eigen::Success, "<FLIPSim::applyPressureProjection> Linear compute fail ...");
    std::cout<<"Solving Poisson equation ..."<<std::endl;
    VectorXd solution = solver.solve(vel_div);
    //assert(solver.info() == Eigen::Success && "Solve Fail ...");
    CheckPoint(solver.info() == Eigen::Success, "<FLIPSim::applyPressureProjection> Linear solve fail ...");

    std::cout<<"Eliminate gradient field to enforce incompressibility ..."<<std::endl;
    // modify u
    for(int i=0; i<m_dimension.x()+1; ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                if(i ==0 || i == m_dimension.x()){
                    m_mac_grid.gridU().cellRef(i,j,k) = 0;
                }
                else{
                    int s_idx = (i-1) + j*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                    int l_idx = i + j*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                    m_mac_grid.gridU().cellRef(i,j,k) -= m_delta_t*(solution(l_idx) - solution(s_idx))/m_grid.gridSpacing(); 
                }
            }
        }
    }    

    // modify v
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y()+1; ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                if(j == 0 || j == m_dimension.y()){
                    m_mac_grid.gridV().cellRef(i,j,k) = 0;
                }
                else{
                    int s_idx = i + (j-1)*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                    int l_idx = i + (j)*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                    m_mac_grid.gridV().cellRef(i,j,k) -= m_delta_t*(solution(l_idx) - solution(s_idx))/m_grid.gridSpacing(); 
                }
            }
        }
    }

    // modify w
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z()+1; ++k){
                if(k ==0 || k == m_dimension.z()){
                    m_mac_grid.gridW().cellRef(i,j,k) = 0;
                }
                else{
                    int s_idx = i + j*m_dimension.x() + (k-1)*m_dimension.x()*m_dimension.y();
                    int l_idx = i + j*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                    m_mac_grid.gridW().cellRef(i,j,k) -= m_delta_t*(solution(l_idx) - solution(s_idx))/m_grid.gridSpacing(); 
                }
            }
        }
    }

}

void FLIPSim::advectParticles(){
    subtractMACGrid();

    for(int i=0; i<m_particle_grid.particleSize(); ++i){
        integrationRungeKutta4(m_particle_grid.particleRef(i));

        m_particle_grid.particleRef(i).velocity() = 
            m_flip_ratio*(m_particle_grid.particleRef(i).velocity()+interpolateVelocity(m_particle_grid.particleRef(i).position(), false))
            + (1-m_flip_ratio)*interpolateVelocity(m_particle_grid.particleRef(i).position(), true);
    }
}

void FLIPSim::adjusteParticles(){
    /*
    // currently only report error 
    for(int i=0; i<m_particle_grid.particleSize(); ++i){
        Vector3d pos = m_particle_grid.particleRef(i).position();
        int pi = static_cast<int>(std::floor(pos.x()/m_grid.gridSpacing()));
        int pj = static_cast<int>(std::floor(pos.y()/m_grid.gridSpacing()));
        int pk = static_cast<int>(std::floor(pos.z()/m_grid.gridSpacing()));
        //std::cout<<pi<<" "<<pj<<" "<<pk<<std::endl;
        //std::cout<<m_type_grid.cellRef(Vector3i(pi,pj,pk))<<std::endl;
        //CheckPoint(m_grid.cellIndexValid(Vector3i(pi,pj,pk)), "<adjustParticles> Particle outside domain ...");
        //CheckPoint(m_type_grid.cellRef(Vector3i(pi,pj,pk)) != SOLID, "<adjustParticles> Particle in solid ...");
    
        if(!m_grid.cellIndexValid(Vector3i(pi,pj,pk))){
            adjustPositionToDomain(pos, Vector3i(pi,pj,pk));

            //if(pk == -1){
            //    std::cout<<"Check: "<<pos<<std::endl;
            //}
            m_particle_grid.particleRef(i).position() = pos;
        }
    }
    */
}

void FLIPSim::enforceBoundaryCondition(){
    // x face
    for(int i=0; i<m_dimension.x()+1; ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                if(i == 0 || i == m_dimension.x()){
                    m_mac_grid.gridU().cellRef(i,j,k) = 0;
                }
                else{
                    if(isSolid(Vector3i(i-1,j,k))*isSolid(Vector3i(i,j,k))){
                        m_mac_grid.gridU().cellRef(i,j,k) = 0;
                    }
                }
            }
        }
    }

    // y face
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y()+1; ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                if(j == 0 || j == m_dimension.y()){
                    m_mac_grid.gridV().cellRef(i,j,k) = 0;
                }
                else{
                    if(isSolid(Vector3i(i,j-1,k))*isSolid(Vector3i(i,j,k))){
                        m_mac_grid.gridV().cellRef(i,j,k) = 0;
                    }
                }
            }
        }
    }

    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z()+1; ++k){
                if(k == 0 || k == m_dimension.z()){
                    m_mac_grid.gridW().cellRef(i,j,k) = 0;
                }
                else{
                    if(isSolid(Vector3i(i,j,k-1))*isSolid(Vector3i(i,j,k))){
                        m_mac_grid.gridW().cellRef(i,j,k) = 0;
                    }
                }
            }
        }
    }
}

void FLIPSim::extrapolateVelocities(){
    extrapolateGrid(U);
    extrapolateGrid(V);
    extrapolateGrid(W);
}

void FLIPSim::assembleLaplacianOperator(){
    //int numF = 0;
    //int numA = 0;
    int dim = m_dimension.x()*m_dimension.y()*m_dimension.z(); 
    std::vector<Triplet> triplets;
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                if(m_type_grid.cellRef(i,j,k) == FLUID){
                    //++numF;              
                    int idx = i + j*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                    std::vector<Vector3i> neighbors;
                    neighbors.push_back(Vector3i(i-1,j,k));
                    neighbors.push_back(Vector3i(i+1,j,k));
                    neighbors.push_back(Vector3i(i,j-1,k));
                    neighbors.push_back(Vector3i(i,j+1,k));
                    neighbors.push_back(Vector3i(i,j,k-1));
                    neighbors.push_back(Vector3i(i,j,k+1));

                    int adjvalid = 0;
                    for(int l=0;l<neighbors.size(); ++l){
                        Vector3i& nv = neighbors[l];
                        int nidx = nv.x() + nv.y()*m_dimension.x() + nv.z()*m_dimension.x()*m_dimension.y();
                    
                        if(cellIndexValid(nv)){
                            if(m_type_grid.cellRef(nv) == FLUID){
                                triplets.push_back(Triplet(idx, nidx, 1.0/std::pow(m_grid.gridSpacing(),2)));
                                ++adjvalid;
                            }
                            else if(m_type_grid.cellRef(nv) == AIR){
                                ++adjvalid;
                            }
                        }
                    }
                    triplets.push_back(Triplet(idx, idx, -1.0*adjvalid/std::pow(m_grid.gridSpacing(),2)));
                }
                else if(m_type_grid.cellRef(i,j,k) == AIR){
                    //++numA;
                    int idx = i + j*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                    triplets.push_back(Triplet(idx, idx, 1));
                }
                else if(m_type_grid.cellRef(i,j,k) == SOLID){
                    // TODO
                }
            }
        }
    }    
    
    //std::cout<<"NumFA: "<<numF<<" "<<numA<<std::endl;
    m_laplacian_op.resize(dim, dim);
    m_laplacian_op.setFromTriplets(triplets.begin(), triplets.end());
}

VectorXd FLIPSim::computeVelocityDivergence(){
    int dim = m_dimension.x()*m_dimension.y()*m_dimension.z(); 
    VectorXd div_vel;
    div_vel.resize(dim);
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                int idx = i + j*m_dimension.x() + k*m_dimension.x()*m_dimension.y();
                if(m_type_grid.cellRef(i,j,k) == FLUID){
                    div_vel(idx) = m_mac_grid.gridU().cellRef(i+1,j,k)
                                - m_mac_grid.gridU().cellRef(i,j,k)
                                + m_mac_grid.gridV().cellRef(i,j+1,k)
                                - m_mac_grid.gridV().cellRef(i,j,k)
                                + m_mac_grid.gridW().cellRef(i,j,k+1)
                                - m_mac_grid.gridW().cellRef(i,j,k);
                    div_vel(idx) /= m_grid.gridSpacing();
                }
                else if(m_type_grid.cellRef(i,j,k) == AIR){
                    div_vel(idx) = 0;
                }
                else if(m_type_grid.cellRef(i,j,k) == SOLID){
                    // TODO
                }
            }
        }
    }
    return div_vel;
}

Vector3d FLIPSim::interpolateVelocity(Vector3d p, bool isPIC){
    int pi = static_cast<int>(std::floor(p.x()/m_grid.gridSpacing()));
    int pj = static_cast<int>(std::floor(p.y()/m_grid.gridSpacing()));
    int pk = static_cast<int>(std::floor(p.z()/m_grid.gridSpacing()));
    //std::cout<<p.x()<<" "<<p.y()<<" "<<p.z()<<std::endl;
    //std::cout<<pi<<" "<<pj<<" "<<pk<<std::endl;

    double ru = (p.x()-pi*m_grid.gridSpacing())/m_grid.gridSpacing();
    double rv = (p.y()-pj*m_grid.gridSpacing())/m_grid.gridSpacing();
    double rw = (p.z()-pk*m_grid.gridSpacing())/m_grid.gridSpacing();

    double pu, pv, pw;
    if(!isPIC){
        pu = (1-ru)*m_mac_grid_d.gridU().cellRef(pi,pj,pk) + ru*m_mac_grid_d.gridU().cellRef(pi+1,pj,pk);
        pv = (1-rv)*m_mac_grid_d.gridV().cellRef(pi,pj,pk) + rv*m_mac_grid_d.gridV().cellRef(pi,pj+1,pk);
        pw = (1-rw)*m_mac_grid_d.gridW().cellRef(pi,pj,pk) + rw*m_mac_grid_d.gridW().cellRef(pi,pj,pk+1);
    }
    else{
        pu = (1-ru)*m_mac_grid.gridU().cellRef(pi,pj,pk) + ru*m_mac_grid.gridU().cellRef(pi+1,pj,pk);
        pv = (1-rv)*m_mac_grid.gridV().cellRef(pi,pj,pk) + rv*m_mac_grid.gridV().cellRef(pi,pj+1,pk);
        pw = (1-rw)*m_mac_grid.gridW().cellRef(pi,pj,pk) + rw*m_mac_grid.gridW().cellRef(pi,pj,pk+1);
    }
    return Vector3d(pu, pv, pw);

    /*
    int nb = 2;

    int pi = static_cast<int>(std::floor(p.x()/m_grid.gridSpacing()));
    int pj = static_cast<int>(std::floor(p.y()/m_grid.gridSpacing()));
    int pk = static_cast<int>(std::floor(p.z()/m_grid.gridSpacing()));

    double u = 0;
    for(int i=-nb+1; i<nb+1; ++i){
        for(int j=-nb; j<nb+1; ++j){
            for(int k=-nb; k<nb+1; ++k){
                Vector3i idx = Vector3i(pi+i, pj+j, pk+k);
                if(m_mac_grid.velocityGridIndexValid(idx, U)){
                    Vector3d fc = Vector3d(idx.x()*m_grid.gridSpacing(), (idx.y()+0.5)*m_grid.gridSpacing(), (idx.z()+0.5)*m_grid.gridSpacing());
                    Vector3d d = fc - p;
                    double wt = QuadraticKernel(d.x()/m_grid.gridSpacing())
                                *QuadraticKernel(d.y()/m_grid.gridSpacing())
                                *QuadraticKernel(d.z()/m_grid.gridSpacing());
                    u += wt*m_mac_grid.gridU().cellRef(idx);
                }
            }
        }
    }

    double v = 0;
    for(int i=-nb; i<nb+1; ++i){
        for(int j=-nb+1; j<nb+1; ++j){
            for(int k=-nb; k<nb+1; ++k){
                Vector3i idx = Vector3i(pi+i, pj+j, pk+k);
                if(m_mac_grid.velocityGridIndexValid(idx, V)){
                    Vector3d fc = Vector3d((idx.x()+0.5)*m_grid.gridSpacing(), idx.y()*m_grid.gridSpacing(), (idx.z()+0.5)*m_grid.gridSpacing());
                    Vector3d d = fc - p;
                    double wt = QuadraticKernel(d.x()/m_grid.gridSpacing())
                                *QuadraticKernel(d.y()/m_grid.gridSpacing())
                                *QuadraticKernel(d.z()/m_grid.gridSpacing());
                    v += wt*m_mac_grid.gridV().cellRef(idx);
                }
            }
        }
    }

    double w = 0;
    for(int i=-nb; i<nb+1; ++i){
        for(int j=-nb; j<nb+1; ++j){
            for(int k=-nb+1; k<nb+1; ++k){
                Vector3i idx = Vector3i(pi+i, pj+j, pk+k);
                if(m_mac_grid.velocityGridIndexValid(idx, W)){
                    Vector3d fc = Vector3d((idx.x()+0.5)*m_grid.gridSpacing(), (idx.y()+0.5)*m_grid.gridSpacing(), idx.z()*m_grid.gridSpacing());
                    Vector3d d = fc - p;
                    double wt = QuadraticKernel(d.x()/m_grid.gridSpacing())
                                *QuadraticKernel(d.y()/m_grid.gridSpacing())
                                *QuadraticKernel(d.z()/m_grid.gridSpacing());
                    w += wt*m_mac_grid.gridW().cellRef(idx);
                }
            }
        }
    }

    return Vector3d(u, v, w);
    */
}

void FLIPSim::adjustPositionToDomain(Vector3d& p)
{
    double margin = 0.001*m_grid.gridSpacing();
    Vector3d mincoord(margin, margin, margin);
    Vector3d maxcoord(m_dimension.x()*m_grid.gridSpacing()-margin, 
                    m_dimension.y()*m_grid.gridSpacing()-margin, 
                    m_dimension.z()*m_grid.gridSpacing()-margin);

    p = Vector3d(std::max(mincoord.x(), p.x()), std::max(mincoord.y(), p.y()), std::max(mincoord.z(), p.z()));
    p = Vector3d(std::min(maxcoord.x(), p.x()), std::min(maxcoord.y(), p.y()), std::min(maxcoord.z(), p.z()));
}

void FLIPSim::extrapolateGrid(mactype type){
    std::cout<<"Initialize wavefront .."<<std::endl;
    int dim;
    Grid<int> distgrid;
    int ub, vb, wb;
    Grid<double>* velg;
    Grid<int>* velmg;
    if(type == U){
        dim = (m_dimension.x()+1)*m_dimension.y()*m_dimension.z();
        distgrid = Grid<int>(Vector3i(m_dimension.x()+1, m_dimension.y(), m_dimension.z()), m_grid.gridSpacing(), dim);
        ub=m_dimension.x()+1;
        vb=m_dimension.y();
        wb=m_dimension.z();
        velg = &m_mac_grid.gridU();
        velmg = &m_mac_grid.gridUM();
    }
    else if(type == V){
        dim = m_dimension.x()*(m_dimension.y()+1)*m_dimension.z();
        distgrid = Grid<int>(Vector3i(m_dimension.x(), m_dimension.y()+1, m_dimension.z()), m_grid.gridSpacing(), dim);
        ub=m_dimension.x();
        vb=m_dimension.y()+1;
        wb=m_dimension.z();
        velg = &m_mac_grid.gridV();
        velmg = &m_mac_grid.gridVM();
    }
    else if(type == W){
        dim = m_dimension.x()*m_dimension.y()*(m_dimension.z()+1);
        distgrid = Grid<int>(Vector3i(m_dimension.x(), m_dimension.y(), m_dimension.z()+1), m_grid.gridSpacing(), dim);
        ub=m_dimension.x();
        vb=m_dimension.y();
        wb=m_dimension.z()+1;
        velg = &m_mac_grid.gridW();
        velmg = &m_mac_grid.gridWM();
    }
    else{
        //assert(false && "Error: mactype not UVW!");
        CheckPoint(false, "<FLIPSim::extrapolateGrid> MAC type not UVW ...");
    }

    std::list<Vector3i> wavefront;
    for(int i=0; i<ub; ++i){
        for(int j=0; j<vb; ++j){
            for(int k=0; k<wb; ++k){
                if(velmg->cellRef(i,j,k) == VISITED){
                    distgrid.cellRef(i,j,k) = 0;
                }
                else{
                    distgrid.cellRef(i,j,k) = dim;

                    std::vector<Vector3i> neighbors;
                    neighbors.push_back(Vector3i(i-1, j, k));
                    neighbors.push_back(Vector3i(i+1, j, k));
                    neighbors.push_back(Vector3i(i, j-1, k));
                    neighbors.push_back(Vector3i(i, j+1, k));
                    neighbors.push_back(Vector3i(i, j, k-1));
                    neighbors.push_back(Vector3i(i, j, k+1));

                    bool gotcha = false;
                    for(int l=0; l<neighbors.size(); ++l){
                        if(macIndexValid(neighbors[l], type)){
                            if(velmg->cellRef(neighbors[l]) == VISITED){
                                gotcha = true;
                            }
                        }
                    }
                    if(gotcha){
                        wavefront.push_back(Vector3i(i,j,k));
                        distgrid.cellRef(i,j,k) = 1;
                    }
                }
            }
        }
    }

    std::cout<<"Propagate to the whole grid ..."<<std::endl;
    while(!wavefront.empty()){
        Vector3i pos = wavefront.front();
        wavefront.pop_front();
        double vel = 0;
        int count = 0;

        std::vector<Vector3i> neighbors;
        neighbors.push_back(Vector3i(pos.x()-1, pos.y(), pos.z()));
        neighbors.push_back(Vector3i(pos.x()+1, pos.y(), pos.z()));
        neighbors.push_back(Vector3i(pos.x(), pos.y()-1, pos.z()));
        neighbors.push_back(Vector3i(pos.x(), pos.y()+1, pos.z()));
        neighbors.push_back(Vector3i(pos.x(), pos.y(), pos.z()-1));
        neighbors.push_back(Vector3i(pos.x(), pos.y(), pos.z()+1));

        for(int i=0; i<neighbors.size(); ++i){
            if(macIndexValid(neighbors[i], type)){
                if(distgrid.cellRef(neighbors[i])<distgrid.cellRef(pos)){
                    vel += velg->cellRef(neighbors[i]);
                    ++count;
                }
                else{
                    if(distgrid.cellRef(neighbors[i]) == dim){
                        distgrid.cellRef(neighbors[i]) = distgrid.cellRef(pos)+1;
                        wavefront.push_back(neighbors[i]);
                    }
                }
            }
        }
        //assert(count != 0 && "Error: Wavefront doest not have visited neighbors ...");
        CheckPoint(count != 0, "<FLIPSim::extrapolateGrid> Wavefront doest not have visited neighbors ...");
        vel/=count;

        velg->cellRef(pos) = vel;
    }
}

void FLIPSim::subtractMACGrid()
{
    for(int i=0; i<m_dimension.x()+1; ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                m_mac_grid_d.gridU().cellRef(i,j,k) = 
                    m_mac_grid.gridU().cellRef(i,j,k)
                    - m_mac_grid_p.gridU().cellRef(i,j,k);
            }
        }
    }

    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y()+1; ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                m_mac_grid_d.gridV().cellRef(i,j,k) =
                    m_mac_grid.gridV().cellRef(i,j,k)
                    - m_mac_grid_p.gridV().cellRef(i,j,k);
            }
        }
    }

    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z()+1; ++k){
                m_mac_grid_d.gridW().cellRef(i,j,k) = 
                    m_mac_grid.gridW().cellRef(i,j,k)
                    - m_mac_grid_p.gridW().cellRef(i,j,k);
            }
        }
    }
}

void FLIPSim::integrationRungeKutta4(Particle& p){
    //std::cout<<"In RK4 ..."<<std::endl;
    Vector3d v1 = m_flip_ratio*(p.velocity()+interpolateVelocity(p.position(), false)) + (1-m_flip_ratio)*interpolateVelocity(p.position(), true);
    Vector3d k1 = m_delta_t*v1;

    Vector3d p2 = p.position()+0.5*k1;
    adjustPositionToDomain(p2);
    Vector3d v2 = m_flip_ratio*(p.velocity()+interpolateVelocity(p2, false)) + (1-m_flip_ratio)*interpolateVelocity(p2, true);
    Vector3d k2 = m_delta_t*v2;
    
    Vector3d p3 = p.position()+0.5*k2;
    adjustPositionToDomain(p3);
    Vector3d v3 = m_flip_ratio*(p.velocity()+interpolateVelocity(p3, false)) + (1-m_flip_ratio)*interpolateVelocity(p3, true);
    Vector3d k3 = m_delta_t*v3;
    
    Vector3d p4 = p.position()+k3;
    adjustPositionToDomain(p4);
    Vector3d v4 = m_flip_ratio*(p.velocity()+interpolateVelocity(p4, false)) + (1-m_flip_ratio)*interpolateVelocity(p4, true);
    Vector3d k4 = m_delta_t*v4; 

    p.position() += (k1 + 2*k2 + 2*k3 + k4) / 6;
    adjustPositionToDomain(p.position());
}

bool FLIPSim::cellIndexValid(Vector3i idx){
    if(idx.x()>=0 && idx.x()<m_dimension.x()
        && idx.y()>=0 && idx.y()<m_dimension.y()
        && idx.z()>=0 && idx.z()<m_dimension.z()){
            return true;
        }
    else{
        return false;
    }
}

bool FLIPSim::macIndexValid(Vector3i idx, mactype type){
    bool valid = false;

    if(type == U){
        if(idx.x()>=0 && idx.x()<m_dimension.x()+1
            && idx.y()>=0 && idx.y()<m_dimension.y()
            && idx.z()>=0 && idx.z()<m_dimension.z()){
                valid = true;
        }
    }
    else if(type == V){
        if(idx.x()>=0 && idx.x()<m_dimension.x()
            && idx.y()>=0 && idx.y()<m_dimension.y()+1
            && idx.z()>=0 && idx.z()<m_dimension.z()){
                valid = true;
        }
    }
    else if(type == W){
        if(idx.x()>=0 && idx.x()<m_dimension.x()
            && idx.y()>=0 && idx.y()<m_dimension.y()
            && idx.z()>=0 && idx.z()<m_dimension.z()+1){
                valid = true;
        }
    }
    else{
        //assert(0 && "Error: Type not supported!");
        CheckPoint(false, "<FLIPSim::macIndexValid> MAC type not UVW ...");
    }

    return valid;
}

bool FLIPSim::isSolid(Vector3i idx){   
    return (m_type_grid.cellRef(idx) == SOLID);
}

/******************************\
 * Debug Functions
\******************************/

void FLIPSim::checkVelocityFieldDivergence(){
    for(int i=0; i<m_dimension.x(); ++i){
        for(int j=0; j<m_dimension.y(); ++j){
            for(int k=0; k<m_dimension.z(); ++k){
                double div = m_mac_grid.gridU().cellRef(i+1,j,k)
                            - m_mac_grid.gridU().cellRef(i,j,k)
                            + m_mac_grid.gridV().cellRef(i,j+1,k)
                            - m_mac_grid.gridV().cellRef(i,j,k)
                            + m_mac_grid.gridW().cellRef(i,j,k+1)
                            - m_mac_grid.gridW().cellRef(i,j,k);
                
                if(div>EPSILON){
                    std::cout<<i<<" "<<j<<" "<<k<<" divergence: "<<div<<std::endl;
                }
            }
        }
    }
}

void FLIPSim::checkParticleVelocity()
{
    std::cout<<"Displaying particle velocities ..."<<std::endl;
    for(int i=0; i<m_particle_grid.particleSize(); ++i){
        if(i<10){
            Vector3d vel = m_particle_grid.particleRef(i).velocity();
            std::cout<<vel.x()<<" "<<vel.y()<<" "<<vel.z()<<std::endl;
        }
    }
}

void FLIPSim::checkKernels()
{
    std::cout<<QuadraticKernel(0.1)<<std::endl;
    std::cout<<QuadraticKernel(0.3)<<std::endl;
    std::cout<<QuadraticKernel(0.5)<<std::endl;
    std::cout<<QuadraticKernel(0.7)<<std::endl;
    std::cout<<QuadraticKernel(0.9)<<std::endl;
}