#include <algorithm>
#include <deque>
#include <functional>
#include <iostream>
#include <stdexcept>

#include <boost/foreach.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <nbt/nbt.hpp>
#include <mc/blocks.hpp>
#include <mc/level.hpp>
#include <mc/utils.hpp>
#include <mc/world.hpp>

#include <mcvoxel/datamodel.hpp>
#include <mcvoxel/io.hpp>
#include <mcvoxel/octree.hpp>

const int CHUNK_SIZE = 16;

/// @brief A functor which loads levels into an octree reference.
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

/// @brief The main program code itself.
struct main_program
{
	std::vector<octree::crystalised_octree>    crystal_octrees;

	int operator() (int argc, char** argv)
	{
		if(argc != 3)
		{
			std::cerr << "usage: " << argv[0] << " <path to world> <output.oct>" << std::endl;
			return 1;
		}

		mc::world world(argv[1]);

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
			crystal_octrees.push_back(octree::crystalised_octree(loaded_tree));
		}

		std::cout << "writing to " << argv[2] << std::endl;
		namespace bio = boost::iostreams;
		bio::filtering_ostream output;
		output.push(bio::zlib_compressor());
		output.push(bio::file_sink(argv[2]));

		uint32_t n_trees = crystal_octrees.size();
		io::nbo::write(output, n_trees);
		BOOST_FOREACH(const octree::crystalised_octree& tree, crystal_octrees)
		{
			tree.serialise(output);
		}

		return 0;
	}
};

int main(int argc, char** argv)
{
	main_program prog;
	return prog(argc, argv);
}
