#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <iostream>
#include <stdexcept>

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

		int32_t start_x((level_info->get_x() * 16) + 1024), start_z((level_info->get_z() * 16) + 1024);

		for(int32_t idx(0); idx < blocks->length; ++idx)
		{
			// if appropriate, insert this block into the octree
			Byte block_id = blocks->values[idx];
			if(block_id != mc::Air)
			{
				if((x < 0) || (y < 0) || (z < 0))
					break;
				if((x >= tree.size()) || (y >= tree.size()) || (z >= tree.size()))
					break;
				tree.set(start_x + x, y, start_z + z, block_id);
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

struct out_pixel
{
	uint8_t r,g,b;
};

void write_ppm(std::ostream& os, out_pixel* src, uint32_t w, uint32_t h)
{
    // write the header and record that we're writing in big endian order
    os << "P6" << std::endl << w << " " << h << std::endl << 255 << std::endl;
    os.write(reinterpret_cast<char*>(src), 3*w*h);
}

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		std::cerr << "usage: " << argv[0] << " <path to world>" << std::endl;
		return 1;
	}

	block_tree tree(12, mc::Air);

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

		std::for_each(level_coords.begin(), level_coords.end(), load_level(region, tree));
		std::cout << " - raw node count " << tree.nodes() << " nodes." << std::endl;
		tree.compact();
		std::cout << " - compacted node count " << tree.nodes() << " nodes." << std::endl;
	}

	std::cout << "Octree has a total of " << tree.nodes() << " nodes." << std::endl;

	const int w=640, h=480;
	boost::shared_ptr<out_pixel> pixels(new out_pixel[w*h]);

	int32_t x(0), y(h);
	for(uint32_t idx=0; idx<w*h; ++idx)
	{
		out_pixel* out = pixels.get() + idx;
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

		make_ray(920, 80, 900, i/mag, j/mag, k/mag, &r);
		octree::node<uint8_t>* int_node = tree.ray_intersect(&r, &t);
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
	write_ppm(output_fstream, pixels.get(), w, h);

	return 0;
}
