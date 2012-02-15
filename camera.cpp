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

}
