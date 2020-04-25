#include "macgrid.h"

namespace FLIPCore{

MACGrid::MACGrid()
{

}

MACGrid::MACGrid(Grid<int>* grid)
{
    m_base_grid = grid;

    Vector3i dim_u = m_base_grid->dimension()+Vector3i(1,0,0);
    Vector3i dim_v = m_base_grid->dimension()+Vector3i(0,1,0);
    Vector3i dim_w = m_base_grid->dimension()+Vector3i(0,0,1);
    Vector3i dim_uvw = m_base_grid->dimension()+Vector3i(1,1,1);

    m_u = Grid<double>(dim_u, m_base_grid->gridSpacing(), 0);
    m_v = Grid<double>(dim_v, m_base_grid->gridSpacing(), 0);
    m_w = Grid<double>(dim_w, m_base_grid->gridSpacing(), 0);
    m_u_m = Grid<int>(dim_u, m_base_grid->gridSpacing(), 0);
    m_v_m = Grid<int>(dim_v, m_base_grid->gridSpacing(), 0);
    m_w_m = Grid<int>(dim_w, m_base_grid->gridSpacing(), 0);
    
    m_velocity = Grid<Vector3d>(dim_uvw, m_base_grid->gridSpacing(), Vector3d());
    m_cell_velocity = Grid<Vector3d>(m_base_grid->dimension(), m_base_grid->gridSpacing(), Vector3d());
}

MACGrid::~MACGrid()
{
    
}

Grid<int>* MACGrid::baseGrid()
{
    return m_base_grid;
}

Grid<double>& MACGrid::gridU()
{
    return m_u;
}

Grid<double>& MACGrid::gridV()
{
    return m_v;
}

Grid<double>& MACGrid::gridW()
{
    return m_w;
}

Grid<int>& MACGrid::gridUM()
{
    return m_u_m;
}

Grid<int>& MACGrid::gridVM()
{
    return m_v_m;
}

Grid<int>& MACGrid::gridWM()
{
    return m_w_m;
}

Grid<double>& MACGrid::gridD()
{   
    return m_d;
}

Grid<double>& MACGrid::gridP()
{
    return m_p;
}

void MACGrid::convertMACToVertexBasedVelocity()
{
    // test
    // assign m_velocity directly
    for(int i=0; i<m_base_grid->dimension().x()+1; ++i){
        for(int j=0; j<m_base_grid->dimension().y()+1; ++j){
            for(int k=0; k<m_base_grid->dimension().z()+1; ++k){
                // average u
                double au=0;
                int nu=0;
                for(int m=-1; m<1; ++m){
                    for(int n=-1; n<1; ++n){
                        int pi = i;
                        int pj = j+m;
                        int pk = k+n;
                    
                        if(velocityGridIndexValid(Vector3i(pi,pj,pk), 0)){
                            au += m_u.cellRef(pi,pj,pk);
                            ++nu;
                        }
                    }
                }
                au /= nu;

                // averrage v
                double av=0;
                int nv=0;
                for(int m=-1; m<1; ++m){
                    for(int n=-1; n<1; ++n){
                        int pi = i+m;
                        int pj = j;
                        int pk = k+n;
                    
                        if(velocityGridIndexValid(Vector3i(pi,pj,pk), 1)){
                            av += m_v.cellRef(pi,pj,pk);
                            ++nv;
                        }
                    }
                }
                av /= nv;

                // average w
                double aw=0;
                int nw=0;
                for(int m=-1; m<1; ++m){
                    for(int n=-1; n<1; ++n){
                        int pi = i+m;
                        int pj = j+n;
                        int pk = k;
                    
                        if(velocityGridIndexValid(Vector3i(pi,pj,pk), 2)){
                            aw += m_w.cellRef(pi,pj,pk);
                            ++nw;
                        }
                    }
                }
                aw /= nw;

                m_velocity.cellRef(i,j,k) = Vector3d(au, av, aw);
            }
        }
    }
}

void MACGrid::convertMACToCellBasedVelocity()
{
    for(int i=0; i<m_base_grid->dimension().x(); ++i){
        for(int j=0; j<m_base_grid->dimension().y(); ++j){
            for(int k=0; k<m_base_grid->dimension().z(); ++k){
                double au = 0;
                int nu = 0;
                for(int l=0; l<2; ++l){
                    int pi = i+l;
                    int pj = j;
                    int pk = k;
                
                    if(velocityGridIndexValid(Vector3i(pi,pj,pk), 0)){
                        au+=m_u.cellRef(pi,pj,pk);
                        ++nu;
                    }
                }
                au/=nu;

                double av = 0;
                int nv = 0;
                for(int l=0; l<2; ++l){
                    int pi = i;
                    int pj = j+l;
                    int pk = k;
                
                    if(velocityGridIndexValid(Vector3i(pi,pj,pk), 1)){
                        av+=m_v.cellRef(pi,pj,pk);
                        ++nv;
                    }
                }
                av/=nv;

                double aw = 0;
                int nw = 0;
                for(int l=0; l<2; ++l){
                    int pi = i;
                    int pj = j;
                    int pk = k+l;
                
                    if(velocityGridIndexValid(Vector3i(pi,pj,pk), 2)){
                        aw+=m_w.cellRef(pi,pj,pk);
                        ++nw;
                    }
                }
                aw/=nw;

                m_cell_velocity.cellRef(i,j,k) = Vector3d(au, av, aw);
            }
        }
    }
}

bool MACGrid::velocityGridIndexValid(Vector3i idx, int type)
{
    bool valid = false;

    // type - 0:u 1:v 2:w
    if(type == 0){ // u grid
        if(idx.x()>=0 && idx.x()<m_base_grid->dimension().x()+1
            && idx.y()>=0 && idx.y()<m_base_grid->dimension().y()
            && idx.z()>=0 && idx.z()<m_base_grid->dimension().z()){
                valid = true;
            }
    }
    if(type == 1){
        if(idx.x()>=0 && idx.x()<m_base_grid->dimension().x()
            && idx.y()>=0 && idx.y()<m_base_grid->dimension().y()+1
            && idx.z()>=0 && idx.z()<m_base_grid->dimension().z()){
                valid = true;
            }
    }
    if(type == 2){
        if(idx.x()>=0 && idx.x()<m_base_grid->dimension().x()
            && idx.y()>=0 && idx.y()<m_base_grid->dimension().y()
            && idx.z()>=0 && idx.z()<m_base_grid->dimension().z()+1){
                valid = true;
            }
    }

    return valid;
}

void MACGrid::writeVelocityGrid(std::string sfx)
{
    std::string fn = "velocitygrid/velocitygrid_"+sfx+".vtk";
    std::ofstream out(fn.c_str());

    out<<"# vtk DataFile Version 3.1"<<std::endl;
    out<<"Grid with Velocity"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET RECTILINEAR_GRID"<<std::endl;
    out<<"DIMENSIONS "
        <<m_base_grid->dimension().x()+1<<" "
        <<m_base_grid->dimension().y()+1<<" "
        <<m_base_grid->dimension().z()+1<<std::endl;
    
    out<<"X_COORDINATES "<<m_base_grid->dimension().x()+1<<" double"<<std::endl;
    for(int i=0; i<m_base_grid->dimension().x()+1; ++i){
        out<<i*m_base_grid->gridSpacing()<<" ";
    }
    out<<std::endl;
    out<<"Y_COORDINATES "<<m_base_grid->dimension().y()+1<<" double"<<std::endl;
    for(int i=0; i<m_base_grid->dimension().y()+1; ++i){
        out<<i*m_base_grid->gridSpacing()<<" ";
    }
    out<<std::endl;
    out<<"Z_COORDINATES "<<m_base_grid->dimension().z()+1<<" double"<<std::endl;
    for(int i=0; i<m_base_grid->dimension().z()+1; ++i){
        out<<i*m_base_grid->gridSpacing()<<" ";
    }
    out<<std::endl;

    int cellNum = (m_base_grid->dimension().x()+1)
                *(m_base_grid->dimension().y()+1)
                *(m_base_grid->dimension().z()+1);
    out<<"POINT_DATA "<<cellNum<<std::endl;
    out<<"SCALARS L double 1"<<std::endl;
    out<<"LOOKUP_TABLE table"<<std::endl;
    for(int k=0; k<m_base_grid->dimension().z()+1; ++k){
        for(int j=0; j<m_base_grid->dimension().y()+1; ++j){
            for(int i=0 ;i<m_base_grid->dimension().x()+1; ++i){
                out<<m_velocity.cellRef(i,j,k).norm()<<std::endl;
            }
        }
    }
    out<<"VECTORS U double"<<std::endl;
    for(int k=0; k<m_base_grid->dimension().z()+1; ++k){
        for(int j=0; j<m_base_grid->dimension().y()+1; ++j){
            for(int i=0 ;i<m_base_grid->dimension().x()+1; ++i){
                out<<m_velocity.cellRef(i,j,k).x()<<" "
                    <<m_velocity.cellRef(i,j,k).y()<<" "
                    <<m_velocity.cellRef(i,j,k).z()<<std::endl;
            }
        }
    }

    out.close();
}

void MACGrid::writeCellVelocityGrid(std::string sfx)
{
    std::string fn = "velocitygrid/cellvelocitygrid_"+sfx+".vtk";
    std::ofstream out(fn.c_str());

    out<<"# vtk DataFile Version 3.1"<<std::endl;
    out<<"Grid with Velocity"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET RECTILINEAR_GRID"<<std::endl;
    out<<"DIMENSIONS "
        <<m_base_grid->dimension().x()+1<<" "
        <<m_base_grid->dimension().y()+1<<" "
        <<m_base_grid->dimension().z()+1<<std::endl;
    
    out<<"X_COORDINATES "<<m_base_grid->dimension().x()+1<<" double"<<std::endl;
    for(int i=0; i<m_base_grid->dimension().x()+1; ++i){
        out<<i*m_base_grid->gridSpacing()<<" ";
    }
    out<<std::endl;
    out<<"Y_COORDINATES "<<m_base_grid->dimension().y()+1<<" double"<<std::endl;
    for(int i=0; i<m_base_grid->dimension().y()+1; ++i){
        out<<i*m_base_grid->gridSpacing()<<" ";
    }
    out<<std::endl;
    out<<"Z_COORDINATES "<<m_base_grid->dimension().z()+1<<" double"<<std::endl;
    for(int i=0; i<m_base_grid->dimension().z()+1; ++i){
        out<<i*m_base_grid->gridSpacing()<<" ";
    }
    out<<std::endl;

    int cellNum = (m_base_grid->dimension().x())
                *(m_base_grid->dimension().y())
                *(m_base_grid->dimension().z());
    out<<"CELL_DATA "<<cellNum<<std::endl;
    out<<"SCALARS L double 1"<<std::endl;
    out<<"LOOKUP_TABLE table"<<std::endl;
    for(int k=0; k<m_base_grid->dimension().z(); ++k){
        for(int j=0; j<m_base_grid->dimension().y(); ++j){
            for(int i=0 ;i<m_base_grid->dimension().x(); ++i){
                out<<static_cast<float>(m_cell_velocity.cellRef(i,j,k).norm())<<std::endl;
            }
        }
    }
    out<<"VECTORS U double"<<std::endl;
    for(int k=0; k<m_base_grid->dimension().z(); ++k){
        for(int j=0; j<m_base_grid->dimension().y(); ++j){
            for(int i=0 ;i<m_base_grid->dimension().x(); ++i){
                out<<static_cast<float>(m_cell_velocity.cellRef(i,j,k).x())<<" "
                    <<static_cast<float>(m_cell_velocity.cellRef(i,j,k).y())<<" "
                    <<static_cast<float>(m_cell_velocity.cellRef(i,j,k).z())<<std::endl;
            }
        }
    }

    out.close();
}

void MACGrid::checkCellVelocity()
{
    double maxspeed = 0;
    for(int k=0; k<m_base_grid->dimension().z(); ++k){
        for(int j=0; j<m_base_grid->dimension().y(); ++j){
            for(int i=0 ;i<m_base_grid->dimension().x(); ++i){
                if(m_cell_velocity.cellRef(i,j,k).norm()>maxspeed){
                    maxspeed = m_cell_velocity.cellRef(i,j,k).norm();
                    //std::cout<<m_cell_velocity.cellRef(i,j,k).norm()<<std::endl;
                }
            }
        }
    }

    std::cout<<"Max grid speed: "<<maxspeed<<std::endl;
}

}