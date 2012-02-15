# Generic third-party library includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Code from c10t to load level files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/c10t/src")
add_subdirectory(c10t/src/mc)
add_subdirectory(c10t/src/nbt)

# HACK: Create an 'image' library with just enough of c10t in it to link.
add_library(c10t-image c10t/src/image/color.cpp)

# HDR image support
add_library(hdrloader hdrloader/hdrloader.cpp)

# Fast Ray/Axis-Aligned Bounding Box Overlap Tests using Ray Slopes
# http://jgt.akpeters.com/papers/EisemannEtAl07/
add_library(rayslope
	rayslope/aabox.cpp
	rayslope/ray.cpp
	rayslope/slope.cpp
	rayslope/slopeint_div.cpp
	rayslope/slopeint_mul.cpp
)

