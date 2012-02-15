#ifndef MC_VOXEL_SKY_HPP__
#define MC_VOXEL_SKY_HPP__

#include <boost/utility.hpp>
#include <Eigen/Dense>
#include <hdrloader/hdrloader.h>

namespace mcvoxel
{

/// @brief A sky represented by an HDR light probe image.
class sky : boost::noncopyable
{
public:
	sky() : max_lum_(0.f), integral_(0.f) { sky_probe_.cols = NULL; }

	sky(const char* filename) { load_hdr_probe(filename); }

	/// @brief Load a HDR sky probe image from a filename.
	///
	/// @param filename
	void load_hdr_probe(const char* filename);

	/// @brief Sample a direction <em>to</em> the sky with an associated chrominance.
	///
	/// This will sample a direction to the light as if the light-probe luminance were a spherical PDF. The
	/// chrominance (pixel value divided by luminance) of the direction sampled is returned in \p chrominance.
	///
	/// @param to_dir
	/// @param chrominance
	void sample_direction(Eigen::Vector3f& to_dir, Eigen::Vector3f& chrominance) const;

	/// @brief Return the pixel value in the specified direction.
	///
	/// @param direction
	Eigen::Vector3f value_in_direction(const Eigen::Vector3f& direction) const;

protected:
	HDRLoaderResult sky_probe_;

	float max_lum_;

	float integral_;

	/// @brief Return the solid angle of the pixel (\p x, \p y)
	///
	/// @param x
	/// @param y
	float pixel_solid_angle(int x, int y) const;
};

}

#endif // MC_VOXEL_SKY_HPP__
