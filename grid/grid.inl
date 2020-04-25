namespace FLIPCore{

template<typename T>
Grid<T>::Grid()
{
    m_dimensions = Vector3i(0,0,0);
    m_background = T();
    m_grid_spacing = 0;
}

template<typename T>
Grid<T>::Grid(const Vector3i& dimensions, const double& grid_spacing, const T& background){
    m_dimensions = dimensions;
    m_background = background;
    m_grid_spacing = grid_spacing;

    m_raw_grid.reserve(m_dimensions.x());
    for(int i=0; i<m_dimensions.x(); ++i){
        std::vector<std::vector<T> > block_2d;
        block_2d.reserve(m_dimensions.y());
        
        for(int j=0; j<m_dimensions.y(); ++j){
            std::vector<T> block_1d;
            block_1d.reserve(m_dimensions.z());
            
            for(int k=0; k<m_dimensions.z(); ++k){
                block_1d.push_back(background);
            }

            block_2d.push_back(block_1d);
        }

        m_raw_grid.push_back(block_2d);
    }
}

template<typename T>
Grid<T>::~Grid()
{

}

template<typename T>
void Grid<T>::setDimension(const Vector3i& dimensions){
    m_dimensions = dimensions;

    m_raw_grid.reserve(m_dimensions.x());
    for(int i=0; i<m_dimensions.x(); ++i){
        std::vector<std::vector<T> > block_2d;
        block_2d.reserve(m_dimensions.y());
        
        for(int j=0; j<m_dimensions.y(); ++j){
            std::vector<T> block_1d;
            block_1d.reserve(m_dimensions.z());
            
            for(int k=0; k<m_dimensions.z(); ++k){
                block_1d.push_back(m_background);
            }

            block_2d.push_back(block_1d);
        }

        m_raw_grid.push_back(block_2d);
    }
}

template<typename T>
void Grid<T>::setGridSpacing(const T& grid_spacing){
    m_grid_spacing = grid_spacing;
}
template<typename T>
T& Grid<T>::cellRef(const Vector3i& index)
{
    return m_raw_grid[index.x()][index.y()][index.z()];
}

template<typename T>
T& Grid<T>::cellRef(const int& x, const int& y, const int& z)
{
    return m_raw_grid[x][y][z];
}

template<typename T>
Vector3d Grid<T>::getCellCenterCoord(const Vector3i& index)
{
    return Vector3d((index.x()+0.5)*m_grid_spacing, (index.y()+0.5)*m_grid_spacing, (index.z()+0.5)*m_grid_spacing);
}

template<typename T>
Vector3d Grid<T>::getCellCenterCoord(const int& x, const int& y, const int& z)
{
    return Vector3d((x+0.5)*m_grid_spacing, (y+0.5)*m_grid_spacing, (z+0.5)*m_grid_spacing);
}

template<typename T>
Vector3i& Grid<T>::dimension()
{
    return m_dimensions;
}

template<typename T>
double& Grid<T>::gridSpacing()
{
    return m_grid_spacing;
}

template<typename T>
T& Grid<T>::background()
{
    return m_background;
}

template<typename T>
void Grid<T>::clear()
{
    m_raw_grid.clear();
    m_dimensions = Vector3i(0,0,0);
    m_grid_spacing = 0;
    m_background = T();
}
template<typename T>
bool Grid<T>::cellIndexValid(Vector3i idx)
{
    if(idx.x()>=0 && idx.x()<m_dimensions.x()
        && idx.y()>=0 && idx.y()<m_dimensions.y()
        && idx.z()>=0 && idx.z()<m_dimensions.z()){
            return true;
        }
    else{
        return false;
    }
}

template<typename T>
void Grid<T>::writeGrid(std::string sfx)
{
    std::string fn = "basegrid_"+sfx+".vtk";
    std::ofstream out(fn.c_str());

    out<<"# vtk DataFile Version 3.1"<<std::endl;
    out<<"FLIP Simluation Base Grid"<<std::endl;
    out<<"ASCII"<<std::endl;
    out<<"DATASET RECTILINEAR_GRID"<<std::endl;
    out<<"DIMENSIONS "
        <<m_dimensions.x()+1<<" "
        <<m_dimensions.y()+1<<" "
        <<m_dimensions.z()+1<<std::endl;
    
    out<<"X_COORDINATES "<<m_dimensions.x()+1<<" double"<<std::endl;
    for(int i=0; i<m_dimensions.x()+1; ++i){
        out<<i*m_grid_spacing<<" ";
    }
    out<<std::endl;
    out<<"Y_COORDINATES "<<m_dimensions.y()+1<<" double"<<std::endl;
    for(int i=0; i<m_dimensions.y()+1; ++i){
        out<<i*m_grid_spacing<<" ";
    }
    out<<std::endl;
    out<<"Z_COORDINATES "<<m_dimensions.z()+1<<" double"<<std::endl;
    for(int i=0; i<m_dimensions.z()+1; ++i){
        out<<i*m_grid_spacing<<" ";
    }
    out<<std::endl;

    out.close();
}


}