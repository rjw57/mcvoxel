#include <cmath>
#include <iterator>
#include <deque>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <rayslope/ray.h>

#include "scene.hpp"

namespace scene
{

// draw a sample uniformly distributed on [0, 1).
static float uniform_real();

scene::scene()
	: current_x_(0.f), current_y_(0.f)
	, r1_(0.1f), r2_(1.f), log_r2_over_r1_(log(r2_/r1_))
{
	initialise(1, 1);
	set_camera(0.f, 0.f, 0.f, 0.f, 0.f);
}

scene::~scene()
{ }

void scene::initialise(int w, int h)
{
	samples_.resize(w, h);

	current_x_ = 0.5f * samples_.width;
	current_y_ = 0.5f * samples_.height;

	r2_ = sqrt(0.05f * samples_.width * samples_.height); // 5% of image area
	log_r2_over_r1_ = log(r2_/r1_);
}

void scene::draw()
{
	peturb_image_loc(current_x_, current_y_, current_x_, current_y_);

	// make an eye ray
	ray eye_ray;
	make_eye_ray(current_x_, current_y_, eye_ray);

	std::deque<world_position> eye_path;
	trace_ray(eye_ray, 7, std::back_inserter(eye_path));

	pixel_f32 sample(0,0,0);
	if(eye_path.size() > 0)
	{
		sample = pixel_f32(1,1,1) * eye_path.size() / 7.f;
	}

	samples_.at(floor(current_x_), floor(current_y_)).record(sample);
}

void scene::peturb_image_loc(float x, float y, float& new_x, float& new_y) const
{
	do
	{
		// choose angle and distance to peturb location
		float theta = 2.f * 3.14159f * uniform_real();
		float r = r2_ * exp(-log_r2_over_r1_ * uniform_real());
		new_x = x + r * cos(theta);
		new_y = y + r * sin(theta);
	}
	while((new_x < 0) || (new_y < 0) || (new_x >= samples_.width) || (new_y >= samples_.height));
}

void scene::set_camera(float x, float y, float z, float yaw, float pitch)
{
	// convert yaw and pitch to radians
	cam_pitch_ = pitch * (2.f * 3.14159f / 360.f);
	cam_yaw_ = yaw * (2.f * 3.14159f / 360.f);

	// pre-calculate trig. values
	cos_pitch_ = cos(cam_pitch_); sin_pitch_ = sin(cam_pitch_);
	cos_yaw_ = cos(cam_yaw_); sin_yaw_ = sin(cam_yaw_);

	// save camera origin
	cam_x_ = x; cam_y_ = y; cam_z_ = z;
}

void scene::make_eye_ray(float fx, float fy, ray& out_ray) const
{
	fx -= (samples_.width)>>1; fy -= (samples_.height)>>1;

	float i(-fx), j(fy), k(samples_.height);

	float new_j = cos_pitch_*j - sin_pitch_*k, new_k = sin_pitch_*j + cos_pitch_*k;
	j = new_j; k = new_k;

	float new_i = cos_yaw_*i - sin_yaw_*k, new_k2 = sin_yaw_*i + cos_yaw_*k;
	i = new_i; k = new_k2;

	float mag = sqrt(i*i + j*j + k*k);
	make_ray(cam_x_, cam_y_, cam_z_, i/mag, j/mag, k/mag, &out_ray);
}

//// UTILITY FUNCTIONS ////

// use a fixed seed so results are replicable
static boost::mt19937 rng(42);
static boost::uniform_real<> uni_dist(0,1);
static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni_gen(rng, uni_dist);

static float uniform_real()
{
	return uni_gen();
}

}
