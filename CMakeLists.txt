project(voxel CXX C)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

# We need to boost libraries.
find_package(Boost REQUIRED COMPONENTS filesystem system)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})

# Generic third-party library includes
include_directories(third-party)

# Add an Octree implementation.
include_directories(third-party/octree)

# Build a demo of using the octree code.
add_executable(octree-demo
	third-party/octree/demo.cpp
)

# Code from c10t to load level files
include_directories(third-party/c10t/src)
add_subdirectory(third-party/c10t/src/mc)
add_subdirectory(third-party/c10t/src/nbt)

# HACK: Create an 'image' library with just enough of c10t in it to link.
add_library(c10t-image third-party/c10t/src/image/color.cpp)

add_executable(voxel
	main.cpp
)

# Fast Ray/Axis-Aligned Bounding Box Overlap Tests using Ray Slopes
# http://jgt.akpeters.com/papers/EisemannEtAl07/
add_library(rayslope
	third-party/rayslope/aabox.cpp
	third-party/rayslope/ray.cpp
	third-party/rayslope/slope.cpp
	third-party/rayslope/slopeint_div.cpp
	third-party/rayslope/slopeint_mul.cpp
)

target_link_libraries(voxel
	c10t-mc c10t-image
	rayslope
	${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}
	${ZLIB_LIBRARY})
