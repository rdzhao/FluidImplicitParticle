#ifndef __GRID_H__
#define __GRID_H__

#include "../utility/utility.h"
#include "../utility/kernel.h"

namespace FLIPCore{

template<typename T> 
class Grid{
    public:
    Grid(); 
    Grid(const Vector3i& dimensions, const double& grid_spacing, const T& background);
    ~Grid();

    void setDimension(const Vector3i& dimensions);
    void setGridSpacing(const T& grid_spacing);

    T& cellRef(const Vector3i& index);
    T& cellRef(const int& x, const int& y, const int& z);
    Vector3i& dimension();
    double& gridSpacing();
    T& background();

    Vector3d getCellCenterCoord(const Vector3i& index);
    Vector3d getCellCenterCoord(const int& x, const int& y, const int& z); 

    void clear();

    bool cellIndexValid(Vector3i idx);

    // I/O functions
    void writeGrid(std::string sfx = "");

    protected:
    std::vector<std::vector<std::vector<T> > > m_raw_grid;
    Vector3i m_dimensions;
    double m_grid_spacing;
    T m_background;

    // we assume that the grid always starts with origin of coordinate system.
};

}

#include "grid.inl"

#endif