#include <boost/foreach.hpp>
#include <io.hpp>
#include <iterator>

#include "world.hpp"

namespace mcvoxel
{

void world::load_from(const char* filename)
{
	trees_.clear();
	io::load_crystal_octrees(filename, std::back_inserter(trees_));
}

bool world::intersect(const ray& r, surface_location& intersection) const
{
	// construct a rayslope ray corresponding to r
	::ray rs_ray;
	make_ray(r.origin()[0], r.origin()[1], r.origin()[2],
		 r.direction()[0], r.direction()[1], r.direction()[2],
		 &rs_ray);

	// intersect with each tree
	bool found_intersection = false;
	BOOST_FOREACH(const octree::crystalised_octree& tree, trees_)
	{
		octree::sub_location sub_loc;
		if(!tree.ray_intersect<data::block>(rs_ray, sub_loc))
			continue;

		// Found a location. Extract the position
		Eigen::Vector3f pos(sub_loc.coords[0], sub_loc.coords[1], sub_loc.coords[2]);

		// Calculate normal
		octree::location node_centre_l(sub_loc.node_extent.midpoint());
		Eigen::Vector3f node_centre(node_centre_l.x, node_centre_l.y, node_centre_l.z);
		Eigen::Vector3f proto_norm(pos-node_centre);

		if((fabs(proto_norm[0]) >= fabs(proto_norm[1])) && (fabs(proto_norm[0]) >= fabs(proto_norm[2])))
		{
			proto_norm[1] = proto_norm[2] = 0.f;
		}
		else if((fabs(proto_norm[1]) >= fabs(proto_norm[0])) && (fabs(proto_norm[1]) >= fabs(proto_norm[2])))
		{
			proto_norm[0] = proto_norm[2] = 0.f;
		}
		else
		{
			proto_norm[0] = proto_norm[1] = 0.f;
		}
		Eigen::Vector3f norm(proto_norm.normalized());

		// Have we already found a location?
		if(found_intersection)
		{
			// Work out if this is closer
			float new_pos_dist = (pos - r.origin()).norm();
			float old_pos_dist = (intersection.position - r.origin()).norm();

			// If not, skip
			if(old_pos_dist < new_pos_dist)
				continue;
		}

		// Update the intersection record
		intersection.position = pos;
		intersection.normal = norm;
		found_intersection = true;
	}

	return found_intersection;
}

}
