#include <boost/assert.hpp>
#include <cmath>
#include <hdrloader/hdrloader.h>
#include <limits>
#include <stdexcept>

#include "colour.hpp"
#include "sampling.hpp"
#include "sky.hpp"

namespace mcvoxel
{

/// @brief Map from normalised texture co-ordinates to the spherical co-ordinate for that point.
///
/// @param x
/// @param y
/// @param r_theta
/// @param r_phi
inline void normalised_texture_to_spherical(float x, float y, float& r_theta, float& r_phi)
{
	r_theta = y * 3.14159f;
	r_phi = ((x * 2.f) - 1.f) * 3.14159f;
}

/// @brief Map from a spherical co-ordinate to the normalised pixel co-ordinate for the corresponding pixel.
///
/// @param theta
/// @param phi
/// @param r_x
/// @param r_y
inline void spherical_to_normalised_texture(float theta, float phi, float& r_x, float& r_y)
{
	r_x = std::max(0.f, std::min(1.f, 0.5f * (1.f + (phi / 3.14159f))));
	r_y = std::max(0.f, std::min(1.f, theta / 3.14159f));
}

void sky::load_hdr_probe(const char* filename)
{
	if(!HDRLoader::load(filename, sky_probe_))
		throw std::runtime_error("failed to load light probe");

	// find sky spherical integral
	double integral = 0.f;
	max_lum_ = 0.f;
	for(int sy=0; sy<sky_probe_.height; ++sy)
	{
		for(int sx=0; sx<sky_probe_.width; ++sx)
		{
			float *p_sky_pix = &(sky_probe_.cols[3*(sx+(sy*sky_probe_.width))]);
			Eigen::Vector3f pixel(p_sky_pix[0], p_sky_pix[1], p_sky_pix[2]);
			float sky_lum = rgb2y(pixel);

			max_lum_ = std::max(max_lum_, sky_lum);
			integral += pixel_solid_angle(sx, sy) * sky_lum;
		}
	}
	integral_ = integral;
}

void sky::sample_direction(Eigen::Vector3f& to_dir, Eigen::Vector3f& chrominance) const
{
	// sanity check: have we a light probe loaded?
	BOOST_ASSERT(sky_probe_.cols != NULL);

	// we sample via rejection sampling: consider the uniform spherical PDF g(t,p) = 1/(4*pi).
	// since we know that f(t,p) <= max_lum/sky_norm, then f(t,p) <= 4*pi*max_lum/sky_norm * g(t,p).
	// and hence f(t,p) <= M g(x) for some constant M.
	//
	// sampling from g() is easy we may therefore use rejection sampling to sample from f
	for(int try_idx=0; try_idx<16; ++try_idx)
	{
		// sample from g()
		to_dir = uniform_direction();

		// calculate f * sky_norm
		Eigen::Vector3f sky_pixel = value_in_direction(to_dir);
		float f = rgb2y(sky_pixel);
		chrominance = sky_pixel / std::max(std::numeric_limits<float>::epsilon(), f);

		// accept? (u < f / (M g(x))) === (u < sky_norm * f / max_lum)
		if(uniform_real() < (f / max_lum_))
		{
			// accepted
			return;
		}
	}

	// if we failed, at least return a zero chrominance
	chrominance = Eigen::Vector3f(0.f,0.f,0.f);
}

float sky::pixel_solid_angle(int x, int y) const
{
	BOOST_ASSERT(sky_probe_.cols != NULL);

	float pix_dx = 1.f / static_cast<float>(sky_probe_.width);
	float pix_dy = 1.f / static_cast<float>(sky_probe_.height);

	float min_x = static_cast<float>(x) * pix_dx;
	float min_y = static_cast<float>(y) * pix_dy;

	float max_x = static_cast<float>(x+1) * pix_dx;
	float max_y = static_cast<float>(y+1) * pix_dy;

	float min_theta, min_phi;
	normalised_texture_to_spherical(min_x, min_y, min_theta, min_phi);

	float max_theta, max_phi;
	normalised_texture_to_spherical(max_x, max_y, max_theta, max_phi);

	return fabs(cos(min_theta) - cos(max_theta)) * fabs(max_phi - min_phi);
}

Eigen::Vector3f sky::value_in_direction(const Eigen::Vector3f& direction) const
{
	BOOST_ASSERT(sky_probe_.cols != NULL);

	float m = direction.norm();
	float theta = acos(direction[1]/m);
	float phi = atan2(direction[2], direction[0]);

	float sky_x, sky_y;
	spherical_to_normalised_texture(theta, phi, sky_x, sky_y);

	int pix_x = static_cast<int>(sky_x * sky_probe_.width) % sky_probe_.width;
	int pix_y = static_cast<int>(sky_y * sky_probe_.height) % sky_probe_.height;

	float *p_sky_pix = &(sky_probe_.cols[3*(pix_x+(pix_y*sky_probe_.width))]);
	return Eigen::Vector3f(p_sky_pix[0], p_sky_pix[1], p_sky_pix[2]);
}

}
