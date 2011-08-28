#include <cmath>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "scene.hpp"

namespace scene
{

// draw a sample uniformly distributed on [0, 1).
static float uniform_real();

scene::scene()
	: current_x_(0.f), current_y_(0.f)
	, r1_(0.1f), r2_(1.f), log_r2_over_r1_(log(r2_/r1_))
{ }

scene::~scene()
{ }

void scene::initialise(image& im)
{
	current_x_ = 0.5f * im.width;
	current_y_ = 0.5f * im.height;

	r2_ = sqrt(0.05f * im.width * im.height); // 5% of image area
	log_r2_over_r1_ = log(r2_/r1_);
}

void scene::draw(scene::image& im)
{
	peturb_image_loc(current_x_, current_y_, im, current_x_, current_y_);
	im.at(floor(current_x_), floor(current_y_)).record(pixel_f32(1,1,1));
}

void scene::peturb_image_loc(float x, float y, const scene::image& im, float& new_x, float& new_y)
{
	do
	{
		// choose angle and distance to peturb location
		float theta = 2.f * 3.14159f * uniform_real();
		float r = r2_ * exp(-log_r2_over_r1_ * uniform_real());
		new_x = x + r * cos(theta);
		new_y = y + r * sin(theta);
	}
	while((new_x < 0) || (new_y < 0) || (new_x >= im.width) || (new_y >= im.height));
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
