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
		bounce new_bounce;

		// if we don't hit, terminate here
		if(!bounce_ray(r, new_bounce))
			return;

		*out++ = new_bounce;
		--max_bounces;

		make_ray(
				new_bounce.where.pos.x,
				new_bounce.where.pos.y,
				new_bounce.where.pos.z,
				new_bounce.bounce_dir.x,
				new_bounce.bounce_dir.y,
				new_bounce.bounce_dir.z,
				&r);
	}
}

}
