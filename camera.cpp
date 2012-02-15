#include <boost/assert.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "camera.hpp"

namespace mcvoxel
{

void camera::get_frame(Eigen::Vector3f& centre,
	 	       Eigen::Vector3f& look_at,
 		       Eigen::Vector3f& up,
 		       Eigen::Vector3f& right) const
{
	// Centre is the origin of the camera-space
	centre = transform() * Eigen::Vector3f(0.f, 0.f, 0.f);

	// Up is the +ve y-axis.
	up = transform() * Eigen::Vector3f(0.f, 1.f, 0.f) - centre;

	// Right is the +ve x-axis.
	right = transform() * Eigen::Vector3f(1.f, 0.f, 0.f) - centre;

	// Look at is the is the -ve z-axis.
	look_at = transform() * Eigen::Vector3f(0.f, 0.f, -1.f) - centre;
}

void camera::set_centre(Eigen::Vector3f centre)
{
	// Get translation from centre to old centre
	Eigen::Vector3f delta(centre - this->centre());

	// Translate the co-ordinate system appropriately.
	translate(delta);
}

bool camera::pixel_coord(const Eigen::Vector3f& direction, float& x, float& y) const
{
	Eigen::Vector3f centre, look_at, up, right;
	get_frame(centre, look_at, up, right);

	// check direction actually intersects plane
	if(direction.dot(look_at) <= 0.f)
		return false;

	// extract components of direction in camera co-ordinate system
	float la(look_at.dot(direction)), u(up.dot(direction)), r(right.dot(direction));

	// eye ray direction is f * look_at + x * right + y * up. If the look_at component is la, we need to scale each
	// component by f / la to get the final result.
	float f(focal_length());
	x = (f/la) * r;
	y = (f/la) * u;

	return true;
}

}
