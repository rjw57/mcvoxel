#include <cmath>
#include <iostream>
#include <stdexcept>

#include "sky.hpp"

#include <hdrloader/hdrloader.h>

namespace sky
{

// utility functions
inline void normalised_texture_to_spherical(float x, float y, float& r_theta, float& r_phi);
inline void spherical_to_normalised_texture(float theta, float phi, float& r_x, float& r_y);

sky::sky(const pixel_f32& sky_constant_)
{
	set_constant(sky_constant_);
}

sky::~sky()
{ }

sky::sky(const char* hdr_filename)
{
	load_hdr(hdr_filename);
}

void sky::load_hdr(const char* hdr_filename)
{
	HDRLoaderResult sky_light_probe;

	if(!HDRLoader::load(hdr_filename, sky_light_probe))
	{
		throw std::runtime_error("failed to load light probe");
	}

	sky_image_.width = sky_light_probe.width;
	sky_image_.height = sky_light_probe.height;
	sky_image_.pixels.clear();
	sky_image_.pixels.reserve(sky_image_.width * sky_image_.height);

	for(off_t idx=0; idx<3*sky_image_.width*sky_image_.height; idx+=3)
	{
		float *p_pix = sky_light_probe.cols + idx;
		sky_image_.pixels.push_back(pixel_f32(p_pix[0], p_pix[1], p_pix[2]));
	}

	std::cout << "Loaded " << sky_image_.width << "x" << sky_image_.height
		<< " sky light probe image from " << hdr_filename << ".\n";

	update_integral();
}

void sky::set_constant(const data::pixel<float>& sky_constant)
{
	sky_image_.width = sky_image_.height = 0;
	sky_image_.pixels.clear();
	sky_constant_ = sky_constant;
	update_integral();
}

const data::pixel<float>& sky::pixel(int x, int y) const
{
	if((width() == 0) || (height() == 0))
		return sky_constant_;

	x %= width();
	if(x < 0)
		x += width();

	y %= height();
	if(y < 0)
		y += height();

	return sky_image_.at(x, y);
}

void sky::update_integral()
{
	if((width() == 0) || (height() == 0))
	{
		const float const_lum = data::rgb2y(sky_constant_);
		sky_lum_integral_ = const_lum * 4.f * 3.14159;
		sky_max_lum_ = const_lum;
		return;
	}

	// find sky spherical integral - since we're adding lots of small values together, do the addition
	// using double precision.
	double integral = 0.;
	sky_max_lum_ = 0.f;
	for(int sy=0; sy<sky_image_.height; ++sy)
	{
		for(int sx=0; sx<sky_image_.width; ++sx)
		{
			const pixel_f32& sky_pix(pixel(sx, sy));
			float sky_lum = data::rgb2y(sky_pix);
			sky_max_lum_ = std::max(sky_max_lum_, sky_lum);
			integral += pixel_solid_angle(sx, sy) * sky_lum;
		}
	}
	sky_lum_integral_ = integral;
}

float sky::pixel_solid_angle(int x, int y) const
{
	if((sky_image_.width == 0) || (sky_image_.height == 0))
		return 4.f * 3.14159f;

	float pix_dx = 1.f / static_cast<float>(sky_image_.width);
	float pix_dy = 1.f / static_cast<float>(sky_image_.height);

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

const data::pixel<float>& sky::sample(float i, float j, float k) const
{
	if((width() == 0) || (height() == 0))
		return sky_constant_;

	float m = sqrt(i*i+j*j+k*k);
	float theta = acos(j/m);
	float phi = atan2(k, i);

	float sky_x, sky_y;
	spherical_to_normalised_texture(theta, phi, sky_x, sky_y);

	int pix_x = static_cast<int>(sky_x * width()) % width();
	int pix_y = static_cast<int>(sky_y * height()) % height();

	if(pix_x < 0)
		pix_x += width();
	if(pix_y < 0)
		pix_y += height();

	return sky_image_.at(pix_x, pix_y);
}

// utility functions

inline void normalised_texture_to_spherical(float x, float y, float& r_theta, float& r_phi)
{
	r_theta = y * 3.14159f;
	r_phi = ((x * 2.f) - 1.f) * 3.14159f;
}

inline void spherical_to_normalised_texture(float theta, float phi, float& r_x, float& r_y)
{
	r_x = std::max(0.f, std::min(1.f, 0.5f * (1.f + (phi / 3.14159f))));
	r_y = std::max(0.f, std::min(1.f, theta / 3.14159f));
}

}
