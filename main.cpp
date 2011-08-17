#include <algorithm>
#include <functional>
#include <iostream>
#include <deque>

#include <octree.h>

#include <nbt/nbt.hpp>
#include <mc/blocks.hpp>
#include <mc/level.hpp>
#include <mc/utils.hpp>
#include <mc/world.hpp>

template<typename T, int AS>
class RayCastingOctree : public Octree<T, AS>
{
	public:

	RayCastingOctree(int size, const T &emptyVal) : Octree<T, AS>(size, emptyVal)
	{ }

	RayCastingOctree(const Octree<T, AS>& o) : Octree<T, AS>(o)
	{ }

	~RayCastingOctree()
	{ }
};

typedef RayCastingOctree<uint8_t, 16> block_octree;

struct load_level : public std::unary_function<const mc::utils::level_coord&, void>
{
	load_level(boost::shared_ptr<mc::region> region, block_octree& octree)
		: region(region), octree(octree)
	{ }

	void operator() (const mc::utils::level_coord& coord)
	{
		mc::level level(boost::shared_ptr<mc::level_info>(new mc::level_info(region, coord)));
		mc::dynamic_buffer region_buffer(mc::region::CHUNK_MAX);
		level.read(region_buffer);

		boost::shared_ptr<nbt::ByteArray> blocks(level.get_blocks());

		int32_t x(coord.get_x() << 4), y(0), z(coord.get_z() << 4);
		for(size_t idx(0); idx < blocks->length; ++idx)
		{
			// if appropriate, insert this block into the octree
			Byte block_id = blocks->values[idx];
			if(block_id != mc::Air)
			{
				if((x < 0) || (y < 0) || (z < 0))
					break;
				if((x >= octree.size()) || (y >= octree.size()) || (z >= octree.size()))
					break;
				octree.set(x,y,z,block_id);
			}

			// go to next co-ord
			++y;
			if(y > 128)
			{
				y = 0; ++z;
				if(z > 16)
				{
					z = 0; ++x;
				}
			}
		}
	}

	boost::shared_ptr<mc::region>   region;
	block_octree&                   octree;
};

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		std::cerr << "usage: " << argv[0] << " <path to world>" << std::endl;
		return 1;
	}

	block_octree octree(4096, mc::Air);

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

		std::for_each(level_coords.begin(), level_coords.end(), load_level(region, octree));
	}

	std::cout << "Octree has a total of " << octree.nodes() << " nodes." << std::endl;

	return 0;
}
