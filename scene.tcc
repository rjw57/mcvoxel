#ifndef INSIDE_SCENE_HPP
#  error "This file must only be included by scene.hpp"
#endif

#include "util/sampling.h"
#include "util/vector.h"

namespace scene
{

template<typename OutputIterator>
void scene::trace_ray(ray& in_r, int max_bounces, OutputIterator out)
{
	world_position bounce_pos;
	octree::sub_location sub_loc;
	data::block block;

	ray r;
	make_ray(in_r.x, in_r.y, in_r.z, in_r.i, in_r.j, in_r.k, &r);

	// don't bounce more than necessary
	while(max_bounces > 0)
	{
		// if we don't hit, terminate here
		if(!world::cast_ray(world, r, sub_loc, block,
					bounce_pos.norm_x, bounce_pos.norm_y, bounce_pos.norm_z))
			return;

		// create out new position structure and add it to the iterator
		bounce_pos.pos_x = sub_loc.coords[0];
		bounce_pos.pos_y = sub_loc.coords[1];
		bounce_pos.pos_z = sub_loc.coords[2];

		bounce_pos.node_ext = sub_loc.node_extent;

		bounce_pos.node_block = block;

		*out++ = bounce_pos;
		--max_bounces;

		// if we no more bounces, we're done
		if(max_bounces <= 0)
			continue;

		// otherwise, choose a new ray direction and bounce it
		Vector normal = pt_vector_make(bounce_pos.norm_x, bounce_pos.norm_x, bounce_pos.norm_z, 0.f);
		Vector new_dir = pt_sampling_hemisphere(normal);

		make_ray(bounce_pos.pos_x, bounce_pos.pos_y, bounce_pos.pos_z,
				pt_vector_get_x(new_dir),
				pt_vector_get_y(new_dir),
				pt_vector_get_z(new_dir),
				&r);
	}
}

}
