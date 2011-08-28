#include <algorithm>
#include <functional>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <nbt/nbt.hpp>
#include <mc/blocks.hpp>
#include <mc/level.hpp>
#include <mc/utils.hpp>
#include <mc/world.hpp>

#include "octree.hpp"
#include "world.hpp"

namespace world
{

// internal struct to load a level
struct load_level : public std::unary_function<const mc::utils::level_coord&, void>
{
	load_level(boost::shared_ptr<mc::region> region, octree::octree<data::block>& octree)
		: region(region), octree(octree)
	{ }

	void operator() (const mc::utils::level_coord& coord)
	{
		boost::shared_ptr<mc::level_info> level_info(new mc::level_info(region, coord));
		mc::level level(level_info);
		mc::dynamic_buffer region_buffer(mc::region::CHUNK_MAX);
		level.read(region_buffer);

		boost::shared_ptr<nbt::ByteArray> blocks(level.get_blocks());

		int32_t x(0), y(0), z(0);

		int32_t start_x(coord.get_x() * 16), start_z(coord.get_z() * 16);
		octree::location offset(octree.first_loc());

		for(int32_t idx(0); idx < blocks->length; ++idx)
		{
			// if appropriate, insert this block into the octree
			Byte block_id = blocks->values[idx];
			if(block_id != mc::Air)
			{
				octree::location loc(x + start_x, y, z + start_z);
				octree.set(loc + offset, data::block(block_id));
			}

			// go to next co-ord
			++y;
			if(y >= 128)
			{
				y = 0; ++z;
				if(z >= 16)
				{
					z = 0; ++x;
				}
			}
		}
	}

	boost::shared_ptr<mc::region>   region;
	octree::octree<data::block>&        octree;
};

void load_world(const char* filename, world& whence)
{
	mc::world world(filename);
	whence.clear();

	mc::region_iterator region_iterator(world.get_iterator());
	std::deque<mc::utils::level_coord> level_coords;
	while(region_iterator.has_next())
	{
		boost::shared_ptr<mc::region> region(region_iterator.next());
		std::cout << "loaded region from: " << region->get_path() << std::endl;

		level_coords.clear();
		region->read_header();
		region->read_coords(level_coords);
		std::cout << " - " << level_coords.size() << " level coords." << std::endl;

		mc::utils::level_coord rc = mc::utils::path_to_region_coord(region->get_path());
		octree::location start_coord = octree::location(rc.get_x() * 16, 0, rc.get_z() * 16);

		octree::octree<data::block> loaded_tree(9, start_coord, mc::Air);
		std::for_each(level_coords.begin(), level_coords.end(), load_level(region, loaded_tree));
		std::cout << " - raw node count " << loaded_tree.node_count() << " nodes." << std::endl;
		std::cout << "   ... compacting" << std::endl;
		loaded_tree.compact();
		std::cout << " - new node count " << loaded_tree.node_count() << " nodes." << std::endl;

		// crystalise into a read-only form
		whence.push_back(octree::crystalised_octree(loaded_tree));
	}
}

void save_cached_world(std::ostream& os, const world& w)
{
	namespace bio = boost::iostreams;
	bio::filtering_ostream output;
	output.push(bio::zlib_compressor());
	output.push(os);

	uint32_t n_trees = w.size();
	io::nbo::write(output, n_trees);
	BOOST_FOREACH(const octree::crystalised_octree& tree, w)
	{
		tree.serialise(output);
	}
}

void load_cached_world(std::istream& is, world& w)
{
	namespace bio = boost::iostreams;
	bio::filtering_istream input;
	input.push(bio::zlib_decompressor());
	input.push(is);

	w.clear();

	uint32_t n_trees;
	io::nbo::read(input, n_trees);
	std::cout << "tree count: " << n_trees << std::endl;
	for(size_t i=0; i<n_trees; ++i)
	{
		w.push_back(octree::crystalised_octree(0, octree::location(0,0,0)));
		w.back().deserialise(input);
		std::cout << " - tree " << i+1 << " loaded." << std::endl;
	}
}

}