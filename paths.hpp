#ifndef MC_VOXEL_PATHS_HPP___
#define MC_VOXEL_PATHS_HPP___

#include <boost/tuple/tuple.hpp>
#include "octree.hpp"

namespace mcvoxel
{

typedef boost::tuple<float, float, float> vector3;

typedef octree::sub_location surface_location;

/// @brief A vertex on an eye or light path consists of where the bounce point
/// is and the likelihood of that bounce having occurred.
struct path_vertex : public boost::tuple<surface_location, float>
{
	path_vertex(const surface_location& loc, float lik)
		: boost::tuple<surface_location, float>(loc, lik)
	{ }

	path_vertex(const path_vertex& v)
		: boost::tuple<surface_location, float>(v)
	{ }

	const path_vertex operator = (const path_vertex& v)
	{ boost::tuple<surface_location, float>::operator =(v); return *this; }

	const surface_location& where() const { return get<0>(); }

	float likelihood() const { return get<1>(); }
};

template<typename InputIterator, typename OutputIterator>
void trace_path(InputIterator first_tree,
		InputIterator last_tree,
		vector3 start_loc, vector3 start_dir,
		OutputIterator output_path_vertices);

}

#define WITHIN_MC_VOXEL_PATHS_HPP__
#include "paths.tcc"
#undef WITHIN_MC_VOXEL_PATHS_HPP__
#endif // MC_VOXEL_PATHS_HPP___
