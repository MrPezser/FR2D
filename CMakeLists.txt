cmake_minimum_required(VERSION 3.25)
project(FR2D)

set(CMAKE_CXX_STANDARD 23)

add_executable(FR2D main.cpp basis.cpp basis.h EulerFlux.cpp EulerFlux.h indexing.h SailorMath.cpp SailorMath.h SpatialDiscretization.cpp SpatialDiscretization.h Setup2D.cpp Setup2D.h output.cpp output.h)
