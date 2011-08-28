#ifndef SKY_HPP
#define SKY_HPP

#include <boost/utility.hpp>

#include <rayslope/ray.h>

#include "datamodel.hpp"
#include "util/vector.h"

namespace sky
{

class sky : public boost::noncopyable
{
	public:
		// by default have uniformly luminous sky
		sky(const data::pixel<float>& sky_constant = data::pixel<float>(1,1,1));
		~sky();

		sky(const char* hdr_filename);

		// load sky from HDR image
		void load_hdr(const char* hdr_filename);

		// set sky to be constant colout
		void set_constant(const data::pixel<float>& sky_constant);

		// solid angle of a pixel in the image
		float pixel_solid_angle(int x, int y) const;

		// sample sky
		const data::pixel<float>& sample(float i, float j, float k) const;
		const data::pixel<float>& sample(const ray& r) const { return sample(r.i, r.j, r.k); }
		const data::pixel<float>& sample(const Vector& v) const
		{
			return sample(pt_vector_get_x(v), pt_vector_get_y(v), pt_vector_get_z(v));
		}

		// convenience accessors
		int width() const { return sky_image_.width; }
		int height() const { return sky_image_.height; }
		const data::image<data::pixel<float> >& image() const { return sky_image_; }

		float max_lum() const { return sky_max_lum_; }
		float lum_integral() const { return sky_lum_integral_; }

	protected:
		typedef data::pixel<float> pixel_f32;
		typedef data::image<pixel_f32> image_f32;

		image_f32 sky_image_;
		pixel_f32 sky_constant_; // if image not loaded, use this constant illumination for sky.

		float sky_lum_integral_; // integral of sky luminance
		float sky_max_lum_; // maximum sky luminance

		const data::pixel<float>& pixel(int x, int y) const;
		void update_integral();
};

}

#endif // SKY_HPP


