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

Eigen::Vector3f cosine_weighted_hemisphere_direction(const Eigen::Vector3f& normal)
{
	// method: raise a uniform sampled disk into the hemisphere.
	float u1 = uniform_real();
	float r = sqrt(u1);
	float theta = 2.f * M_PI * uniform_real();
	float x = r * cos(theta), y = r * sin(theta);

	// a dirty trick for rotating the sample
	Eigen::Vector3f norm_normal(normal.normalized()), h(norm_normal);
	if((fabs(h[0]) <= fabs(h[1])) && (fabs(h[0]) <= fabs(h[2])))
	{
		h[0] = 1.f;
	}
       	else if((fabs(h[1]) <= fabs(h[0])) && (fabs(h[1]) <= fabs(h[2])))
	{
		h[1] = 1.f;
	}
	else
	{
		h[2] = 1.f;
	}

	Eigen::Vector3f xv = h.cross(normal).normalized();
	Eigen::Vector3f yv = xv.cross(normal).normalized();

	return x*xv + y*yv + sqrt(std::max(0.f, 1.f-u1)) * norm_normal;
}

}
