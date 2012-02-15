#ifndef WITHIN_MC_VOXEL_PATHS_HPP__
#error "this file must only be included from within paths.hpp"
#endif

namespace mcvoxel
{

template<typename InputIterator, typename OutputIterator>
void trace_path(InputIterator first_tree,
		InputIterator last_tree,
		vector3 start_loc, vector3 start_dir,
		OutputIterator output_path_vertices);

}
