#ifndef MC_VOXEL_TRACE_HPP__
#define MC_VOXEL_TRACE_HPP__

#include "ray.hpp"
#include "world.hpp"

namespace mcvoxel
{

/// @brief Trace a ray into the world propagating it through multiple bounces to form a light path.
///
/// The underlying PDF sampled from is proportional to the luminance transfer along the path.
///
/// @tparam OutputIterator An implementation of the \c OutputIterator concept for surface_location instances.
/// @param r The initial ray
/// @param w The world to cast rays into
/// @param max_bounces The maximum number of bounces to write to \p output_bounces
/// @param output_final_ray The final ray direction sampled
/// @param output_bounces Where to write the output bounces to
template<typename OutputIterator>
void trace_path(const ray& r,
		const world& w,
     		size_t max_bounces,
		ray& output_final_ray,
		OutputIterator output_bounces)
{
	ray next_ray(r);

	for(size_t n_bounces = 0; n_bounces < max_bounces; ++n_bounces, ++output_bounces)
	{
		surface_location intersection;
		if(!w.intersect(next_ray, intersection))
			break;

		*output_bounces = intersection;

		next_ray = ray(intersection.position,
				  cosine_weighted_hemisphere_direction(intersection.normal));
	}

	output_final_ray = next_ray;
}

}

#endif // MC_VOXEL_TRACE_HPP__
