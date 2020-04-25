#ifndef __MACGRID_H__
#define __MACGRID_H__

#include "grid.h"

namespace FLIPCore{
 
class MACGrid{
    public:
    MACGrid();
    MACGrid(Grid<int>* grid);
    ~MACGrid();

    Grid<int>* baseGrid();

    Grid<double>& gridU();
    Grid<double>& gridV(); 
    Grid<double>& gridW();

    Grid<int>& gridUM();
    Grid<int>& gridVM(); 
    Grid<int>& gridWM();

    Grid<double>& gridD();
    Grid<double>& gridP();

    void convertMACToVertexBasedVelocity();
    void convertMACToCellBasedVelocity();

    bool velocityGridIndexValid(Vector3i idx, int type);

    // I/O functions
    void writeVelocityGrid(std::string sfx="");
    void writeCellVelocityGrid(std::string sfx="");

    // debug
    void checkCellVelocity();

    protected:
    Grid<int>* m_base_grid;

    // velocities in 3 directions
    Grid<double> m_u; 
    Grid<double> m_v;
    Grid<double> m_w;
    Grid<Vector3d> m_velocity; // approximated vertex based velocity
    Grid<Vector3d> m_cell_velocity; // approximated cell based velocity

    // visited or not (for extrapolation)
    Grid<int> m_u_m;
    Grid<int> m_v_m;
    Grid<int> m_w_m;

    Grid<double> m_d; // divergence grid
    Grid<double> m_p; // pressure grid
};

}

#endif