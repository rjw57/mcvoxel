#ifndef MC_VOXEL_RAY_HPP__
#define MC_VOXEL_RAY_HPP__

#include <Eigen/Dense>

namespace mcvoxel
{

class ray
{
public:
	const Eigen::Vector3f origin() const { return origin_; }

	const Eigen::Vector3f direction() const { return direction_; }

	void set_origin(const Eigen::Vector3f& o) { origin_ = o; }

	void set_direction(const Eigen::Vector3f& d) { direction_ = d.normalized(); }

	/// @brief Create a ray pointing from the origin along the -ve z-axis.
	ray()
		: origin_(0.f,0.f,0.f)
		, direction_(0.f,0.f,-1.f)
	{ }

	ray(const Eigen::Vector3f origin,
	    const Eigen::Vector3f direction)
		: origin_(origin)
		, direction_(direction.normalized())
	{ }

	ray(const ray& r)
		: origin_(r.origin_)
		, direction_(r.direction_)
	{ }

	const ray& operator = (const ray& r)
	{
		origin_ = r.origin_;
		direction_ = r.direction_;
		return *this;
	}

protected:
	Eigen::Vector3f origin_;
	Eigen::Vector3f direction_;
};

}

#endif // MC_VOXEL_RAY_HPP__
