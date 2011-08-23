#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <boost/foreach.hpp>

#include <nbt/nbt.hpp>
#include <mc/blocks.hpp>
#include <mc/level.hpp>
#include <mc/utils.hpp>
#include <mc/world.hpp>

#include <rayslope/aabox.h>
#include <rayslope/ray.h>
#include <rayslope/slope.h>
#include <rayslope/slopeint_mul.h>

#include "datamodel.hpp"
#include "io.hpp"
#include "octree.hpp"

typedef octree::tree<uint8_t> block_tree;

struct load_level : public std::unary_function<const mc::utils::level_coord&, void>
{
	load_level(boost::shared_ptr<mc::region> region, block_tree& tree)
		: region(region), tree(tree)
	{ }

	void operator() (const mc::utils::level_coord& coord)
	{
		boost::shared_ptr<mc::level_info> level_info(new mc::level_info(region, coord));
		mc::level level(level_info);
		mc::dynamic_buffer region_buffer(mc::region::CHUNK_MAX);
		level.read(region_buffer);

		boost::shared_ptr<nbt::ByteArray> blocks(level.get_blocks());

		int32_t x(0), y(0), z(0);

		int32_t start_x(level_info->get_x() * 16), start_z(level_info->get_z() * 16);

		for(int32_t idx(0); idx < blocks->length; ++idx)
		{
			// if appropriate, insert this block into the octree
			Byte block_id = blocks->values[idx];
			if(block_id != mc::Air)
			{
				octree::location loc(x + start_x, y, z + start_z);
				if(!tree.root().contains(loc))
					continue;
				tree.set(loc, block_id);
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
	block_tree&                     tree;
};

struct intersection_record
{
	octree::node<uint8_t>* node;
	float                  t;
	bool operator < (const intersection_record& rhs) const { return t < rhs.t; }
};

octree::node<uint8_t>* cast_ray(ray* r, const std::vector<block_tree>& trees, float *t)
{
	float min_t = 0.f;
	octree::node<uint8_t>* closest_node = NULL;

	BOOST_FOREACH(const block_tree& tree, trees)
	{
		float temp_t;
		octree::node<uint8_t>* node = tree.ray_intersect(r, &temp_t);
		if(node == NULL)
			continue;

		if((closest_node == NULL) || (temp_t < min_t))
		{
			closest_node = node;
			min_t = temp_t;
		}
	}

	*t = min_t;
	return closest_node;
}

inline bool is_surrounded(const block_tree& tree, const octree::location& loc)
{
	if(tree.get(loc.x-1, loc.y, loc.z) == mc::Air) return false;
	if(tree.get(loc.x+1, loc.y, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y-1, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y+1, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y, loc.z-1) == mc::Air) return false;
	if(tree.get(loc.x, loc.y, loc.z+1) == mc::Air) return false;
	return true;
}

void sparsify(block_tree& tree)
{
	octree::location first_loc = tree.root().min_loc;
	octree::location last_loc(first_loc);

	last_loc.x += tree.root().size();
	last_loc.y += tree.root().size();
	last_loc.z += tree.root().size();

	for(int x=first_loc.x+1; x<last_loc.x-1; ++x)
		for(int y=1; y<128-1; ++y)
			for(int z=first_loc.z+1; z<last_loc.z-1; ++z)
			{
				if(is_surrounded(tree, octree::location(x, y, z)))
					tree.set(x,y,z,0xff);
			}
}

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		std::cerr << "usage: " << argv[0] << " <path to world>" << std::endl;
		return 1;
	}

	std::vector<block_tree> trees;
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
		trees.push_back(block_tree(9, start_coord, mc::Air));

		std::for_each(level_coords.begin(), level_coords.end(), load_level(region, trees.back()));
		std::cout << " - raw node count " << trees.back().nodes() << " nodes." << std::endl;
		std::cout << "   ... sparsifying" << std::endl;
		sparsify(trees.back());
		std::cout << "   ... compacting" << std::endl;
		trees.back().compact();
		std::cout << " - compacted and sparsified node count " << trees.back().nodes() << " nodes." << std::endl;
	}

	std::ofstream output("foo.dat");
	BOOST_FOREACH(const block_tree& tree, trees)
	{
		output << tree;
	}

	const int w=640, h=480;
	boost::shared_ptr<data::pixel> pixels(new data::pixel[w*h]);

	int32_t x(0), y(h);
	for(uint32_t idx=0; idx<w*h; ++idx)
	{
		data::pixel* out = pixels.get() + idx;
		float fx(x), fy(y), t(0.f);
		//block_pos pos;
		ray r;

		if((idx & 0xff) == 0)
		{
			std::cout << 100*idx/(w*h) << "%\r" << std::flush;
		}

		fx -= w>>1; fy -= h>>1;

		float i(fx), j(fy), k(0.5f*h);

		float mag = sqrt(i*i + j*j + k*k);

		make_ray(50, 90, -10, i/mag, j/mag, k/mag, &r);
		octree::node<uint8_t>* int_node = cast_ray(&r, trees, &t);
		if(int_node != NULL)
		{
			uint8_t block_id = int_node->value;
			switch(block_id)
			{
				case mc::Stone:
					out->r = out->g = out->b = 0x33;
					break;
				case mc::Grass:
					out->r = 0x00; out->g = 0x7f; out->b = 0x00;
					break;
				case mc::Wood:
				case mc::Log:
					out->r = 0x30; out->g = 0x30; out->b = 0x00;
					break;
				case mc::Dirt:
					out->r = 0x60; out->g = 0x60; out->b = 0x00;
					break;
				case mc::Sand:
					out->r = 0xe0; out->g = 0xe0; out->b = 0x00;
					break;
				case mc::Water:
				case mc::StationaryWater:
					out->r = 0x00; out->g = 0x00; out->b = 0x80;
					break;
				case mc::Leaves:
					out->r = 0x00; out->g = 0xff; out->b = 0x00;
					break;
				default:
					//std::cout << "unknown block: 0x" << std::hex << (int) block_id << std::endl;
					out->r = out->g = out->b = 0x7f;
					break;
			}
		}
		else
		{
			out->r = out->g = out->b = 0x00;
		}

		++x;
		if(x >= w)
		{
			x = 0; --y;
		}
	}

	std::ofstream output_fstream("output.ppm");
	io::write_ppm(output_fstream, pixels.get(), w, h);

	return 0;
}
