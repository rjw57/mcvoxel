#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>

#include "sampling.hpp"

namespace mcvoxel
{
static boost::mt19937 rng_(std::time(0));

float uniform_real(float first, float last)
{
	boost::uniform_real<> uni_dist(first, last);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(rng_, uni_dist);
	return uni();
}

Eigen::Vector3f uniform_direction()
{
	float u(uniform_real()), v(uniform_real());
	float theta = 2.f * 3.14159f * u;
	float phi = acos(2.f * v - 1.f);
	float st = sin(theta), ct = cos(theta);
	float sp = sin(phi), cp = cos(phi);
	return Eigen::Vector3f(ct*sp, st*sp, cp);
}

}
