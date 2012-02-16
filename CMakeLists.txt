project(voxel CXX C)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -ffast-math -march=native")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -ffast-math -march=native")

find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# We need to boost libraries.
find_package(Boost REQUIRED COMPONENTS filesystem iostreams system thread)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})

find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIR})
add_definitions(${PNG_DEFINITIONS})

find_package(PkgConfig REQUIRED)

pkg_check_modules(GLIB REQUIRED glib-2.0)
include_directories(${GLIB_INCLUDE_DIRS})
link_directories(${GLIB_LIBRARY_DIRS})

pkg_check_modules(GTHREAD REQUIRED gthread-2.0)
include_directories(${GTHREAD_INCLUDE_DIRS})
link_directories(${GTHREAD_LIBRARY_DIRS})

pkg_check_modules(EIGEN eigen3 REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})
link_directories(${EIGEN_LIBRARY_DIRS})

# Third-party libraries
add_subdirectory(third-party)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/third-party")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/third-party/c10t/src")

# Generic utility functions
add_library(util
	util/common.c
	util/sampling.c
	util/vector.c
)

add_subdirectory(mcvoxel)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(map2octree tools/map2octree.cpp)
target_link_libraries(map2octree mcvoxel)

add_executable(simple_path_trace tools/simple_path_trace.cpp)
target_link_libraries(simple_path_trace mcvoxel)

add_executable(sample_sky tools/sample_sky.cpp)
target_link_libraries(sample_sky mcvoxel)

add_executable(voxel tools/voxel.cpp)
target_link_libraries(voxel mcvoxel)
