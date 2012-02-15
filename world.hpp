#ifndef MC_VOXEL_WORLD_HPP___
#define MC_VOXEL_WORLD_HPP___

#include <boost/utility.hpp>
#include <deque>

#include "octree.hpp"
#include "ray.hpp"

namespace mcvoxel
{

struct surface_location
{
	Eigen::Vector3f position;
};

/// @brief A world composed of one or more crystalised octrees.
class world : boost::noncopyable
{
public:
	world() { }

	/// @brief Load the world from a serialised and compressed octree representation.
	///
	/// @param filename
	world(const char* filename) { load_from(filename); }

	/// @brief Load the world from a serialised and compressed octree representation.
	///
	/// @param filename
	void load_from(const char* filename);

	/// @brief Attempt to intersect the ray \p r with the world returning the world-space location for the
	/// intersection, the surface normal and the surface material in \p intersection. It there is no intersection
	/// return \c false otherwise return \c true.
	///
	/// @param r
	/// @param intersection
	bool intersect(const ray& r, surface_location& intersection) const;

protected:
	std::deque<octree::crystalised_octree> trees_;
};

}

#endif // MC_VOXEL_WORLD_HPP___
