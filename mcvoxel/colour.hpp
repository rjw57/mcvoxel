#ifndef MC_VOXEL_COLOUR_HPP__
#define MC_VOXEL_COLOUR_HPP__

#include <Eigen/Dense>

namespace mcvoxel
{

/// @brief Convert a RGB triplet into a luminance value.
///
/// @param v
inline float rgb2y(const Eigen::Vector3f& v)
{
	return v.dot(Eigen::Vector3f(0.299f, 0.114f, 0.587f));
}

}

#endif // MC_VOXEL_COLOUR_HPP__
