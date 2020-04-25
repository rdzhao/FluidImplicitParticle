cmake_minimum_required(VERSION 3.1)

project(FLIPSim)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

include_directories(Eigen)

set(Headers
    grid/grid.h
    grid/macgrid.h
    grid/particlegrid.h
    simulation/flip.h
    utility/utility.h
    utility/kernel.h
)

set(Sources
    main.cpp
    grid/grid.inl
    grid/macgrid.cpp
    grid/particlegrid.cpp
    simulation/flip.cpp
    utility/kernel.cpp
)

add_executable(FLIPSim ${Headers} ${Sources})

