#ifndef MC_VOXEL_SAMPLING_HPP__
#define MC_VOXEL_SAMPLING_HPP__

#include <Eigen/Dense>

namespace mcvoxel
{

/// @brief Draw a real uniformly from the interval [\p first, \p last).
///
/// @param first
/// @param last
float uniform_real(float first = 0.f, float last = 1.f);

/// @brief Draw a direction uniformly from the unit sphere.
Eigen::Vector3f uniform_direction();

}

#endif // MC_VOXEL_SAMPLING_HPP__