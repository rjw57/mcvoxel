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
find_package(Boost REQUIRED COMPONENTS filesystem iostreams system)
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

# Generic third-party library includes
include_directories(third-party)

# Code from c10t to load level files
include_directories(third-party/c10t/src)
add_subdirectory(third-party/c10t/src/mc)
add_subdirectory(third-party/c10t/src/nbt)

# HACK: Create an 'image' library with just enough of c10t in it to link.
add_library(c10t-image third-party/c10t/src/image/color.cpp)

# HDR image support
add_library(hdrloader third-party/hdrloader/hdrloader.cpp)

# Generic utility functions
add_library(util
	util/common.c
	util/sampling.c
	util/vector.c
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_executable(voxel
	io.cpp
	main.cpp
	octree.cpp
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
	c10t-mc c10t-image rayslope hdrloader util
	${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_IOSTREAMS_LIBRARY}
	${ZLIB_LIBRARY} ${GLIB_LIBRARIES} ${GTHREAD_LIBRARIES} ${PNG_LIBRARIES}
)
