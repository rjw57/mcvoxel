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

const int CHUNK_SIZE = 16;

struct load_level : public std::unary_function<const mc::utils::level_coord&, void>
{
	load_level(boost::shared_ptr<mc::region> region, octree::octree<uint8_t>& octree)
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
				octree.set(loc + offset, block_id);
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
	octree::octree<uint8_t>&        octree;
};

template<typename T>
const T* cast_ray(const ray& r, const std::vector< octree::octree<T> >& trees, octree::sub_location& out_sub_loc)
{
	const T* closest_node_p = NULL;
	float min_dist_sq = 0.f;

	BOOST_FOREACH(const octree::octree<T>& tree, trees)
	{
		octree::sub_location temp_sub_loc;

		if(!tree.ray_intersect(r, temp_sub_loc))
			continue;

		float dx = (temp_sub_loc.coords[0] - r.x);
		float dy = (temp_sub_loc.coords[1] - r.y);
		float dz = (temp_sub_loc.coords[2] - r.z);
		float dist_sq = dx*dx + dy*dy + dz*dz;

		if((closest_node_p == NULL) || (dist_sq < min_dist_sq))
		{
			closest_node_p = &(tree.get(temp_sub_loc.node_extent.loc));
			min_dist_sq = dist_sq;
			out_sub_loc = temp_sub_loc;
		}
	}

	return closest_node_p;
}


template<typename T>
bool is_surrounded(octree::octree<T>& tree, const octree::location& loc)
{
	if(tree.get(loc.x-1, loc.y, loc.z) == mc::Air) return false;
	if(tree.get(loc.x+1, loc.y, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y-1, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y+1, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y, loc.z-1) == mc::Air) return false;
	if(tree.get(loc.x, loc.y, loc.z+1) == mc::Air) return false;
	return true;
}

template<typename T>
void sparsify(octree::octree<T>& tree)
{
	octree::location first_loc = tree.first_loc();
	octree::location last_loc(first_loc);

	last_loc.x += tree.size();
	last_loc.y += tree.size();
	last_loc.z += tree.size();

	for(int x=first_loc.x+1; x<last_loc.x-1; ++x)
		for(int y=1; y<128-1; ++y)
			for(int z=first_loc.z+1; z<last_loc.z-1; ++z)
			{
				if(is_surrounded(tree, octree::location(x, y, z)))
					tree.set(x,y,z,0xff);
			}
}

struct main_program
{
	std::vector< octree::octree<uint8_t> > octrees;

	float light_x, light_y, light_z;

	int operator() (int argc, char** argv)
	{
		if(argc != 2)
		{
			std::cerr << "usage: " << argv[0] << " <path to world>" << std::endl;
			return 1;
		}

#if 1
		mc::world world(argv[1]);

		std::cout << "node size: " << sizeof(octree::octree<uint8_t>::branch_or_leaf_node_t) << " bytes" << std::endl;

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
			octree::location start_coord = octree::location(rc.get_x() * 8, 0, rc.get_z() * 8);

			octrees.push_back(octree::octree<uint8_t>(9, start_coord, mc::Air));

			std::for_each(level_coords.begin(), level_coords.end(), load_level(region, octrees.back()));
			std::cout << " - raw node count " << octrees.back().node_count() << " nodes." << std::endl;
#if 0
			std::cout << "   ... sparsifying" << std::endl;
			sparsify(octrees.back());
#endif
			std::cout << "   ... compacting" << std::endl;
			octrees.back().compact();
			std::cout << " - new node count " << octrees.back().node_count() << " nodes." << std::endl;
		}

		{
			std::ofstream output("foo.dat");
			output << octrees.size() << std::endl;
			BOOST_FOREACH(const octree::octree<uint8_t>& tree, octrees)
			{
				output << tree;
			}
			output.close();
		}
#else
		octrees.clear();
		{
			std::ifstream input("foo.dat");
			size_t n_trees;
			input >> n_trees;
			std::cout << "tree count: " << n_trees << std::endl;
			for(size_t i=0; i<n_trees; ++i)
			{
				octrees.push_back(octree::octree<uint8_t>(9));
				input >> octrees.back();
				std::cout << " - tree " << i+1 << " loaded with "
					<< octrees.back().node_count() << " nodes." << std::endl;
			}
			input.close();
		}
#endif

		const int w=640, h=480;
		boost::shared_ptr< data::pixel<uint8_t> > pixels(new data::pixel<uint8_t>[w*h]);

		light_x = 1.f;
		light_y = 4.f;
		light_z = -1.f;
		float light_mag = sqrt(light_x*light_x + light_y*light_y + light_z*light_z);

		light_x /= light_mag;
		light_y /= light_mag;
		light_z /= light_mag;

		int32_t x(0), y(h);
		for(uint32_t idx=0; idx<w*h; ++idx)
		{
			data::pixel<uint8_t>* out = pixels.get() + idx;

			data::pixel<float> fout;

#if 1
			float fx(x), fy(y);
			fout = sample(fx, fy, w, h) * 2;
			fout = fout + sample(fx - 0.25f, fy, w, h);
			fout = fout + sample(fx + 0.25f, fy, w, h);
			fout = fout + sample(fx, fy - 0.25f, w, h);
			fout = fout + sample(fx, fy + 0.25f, w, h);
			fout = fout / 6;
#else
			fout = sample(x, y, w, h);
#endif

			out->r = std::max(0, std::min(0xff, static_cast<int>(fout.r)));
			out->g = std::max(0, std::min(0xff, static_cast<int>(fout.g)));
			out->b = std::max(0, std::min(0xff, static_cast<int>(fout.b)));

			if((idx & 0xff) == 0)
			{
				std::cout << 100*idx/(w*h) << "%\r" << std::flush;
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

	data::pixel<float> sample(float fx, float fy, int w, int h)
	{
		data::pixel<float> output;

		ray r;

		fx -= w>>1; fy -= h>>1;

		//float i(fx), k(fy), j(-0.5f*h);
		float i(fx), j(fy), k(0.5f*h);

		float mag = sqrt(i*i + j*j + k*k);

		make_ray(40, 87.8f, 10, i/mag, j/mag, k/mag, &r);
		//make_ray(0, 1200, 0, i/mag, j/mag, k/mag, &r);

		octree::sub_location node_sub_loc;
		const uint8_t* block_id_p = cast_ray(r, octrees, node_sub_loc);

		if(block_id_p != NULL)
		{
			const octree::extent& node_ext = node_sub_loc.node_extent;

			float mid_x = static_cast<float>(node_ext.loc.x) + 0.5f * static_cast<float>(node_ext.size);
			float mid_y = static_cast<float>(node_ext.loc.y) + 0.5f * static_cast<float>(node_ext.size);
			float mid_z = static_cast<float>(node_ext.loc.z) + 0.5f * static_cast<float>(node_ext.size);

			float hit_x = node_sub_loc.coords[0];
			float hit_y = node_sub_loc.coords[1];
			float hit_z = node_sub_loc.coords[2];

			float normal_x = hit_x - mid_x;
			float normal_y = hit_y - mid_y;
			float normal_z = hit_z - mid_z;

			float mag_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
			mag_normal = 0.5f * static_cast<float>(node_ext.size);
			normal_x /= mag_normal;
			normal_y /= mag_normal;
			normal_z /= mag_normal;

			if((abs(normal_y) >= abs(normal_x)) && (abs(normal_y) >= abs(normal_z)))
			{
				normal_y = normal_y > 0.f ? 1.f : -1.f;
				normal_x = normal_z = 0.f;
			}
			else if((abs(normal_x) >= abs(normal_y)) && (abs(normal_x) >= abs(normal_z)))
			{
				normal_x = normal_x > 0.f ? 1.f : -1.f;
				normal_y = normal_z = 0.f;
			}
			else if((abs(normal_z) >= abs(normal_x)) && (abs(normal_z) >= abs(normal_y)))
			{
				normal_z = normal_z > 0.f ? 1.f : -1.f;
				normal_x = normal_y = 0.f;
			}

			float luminance = 0.2f + 0.8f * std::max(0.f, light_x*normal_x + light_y*normal_y + light_z*normal_z);

			uint8_t block_id = *block_id_p;
			switch(block_id)
			{
				case mc::Stone:
					output.r = output.g = output.b = 0x33;
					break;
				case mc::Grass:
					output.r = 0x00; output.g = 0x7f; output.b = 0x00;
					break;
				case mc::Wood:
				case mc::Log:
					output.r = 0x30; output.g = 0x30; output.b = 0x00;
					break;
				case mc::Dirt:
					output.r = 0x60; output.g = 0x60; output.b = 0x00;
					break;
				case mc::Sand:
					output.r = 0xe0; output.g = 0xe0; output.b = 0x00;
					break;
				case mc::Water:
				case mc::StationaryWater:
					output.r = 0x00; output.g = 0x00; output.b = 0x80;
					break;
				case mc::Leaves:
					output.r = 0x00; output.g = 0xff; output.b = 0x00;
					break;
				default:
					//std::cout << "unknown block: 0x" << std::hex << (int) block_id << std::endl;
					output.r = output.g = output.b = 0x7f;
					break;
			}

			//output.r = 127 + (127*normal_x);
			//output.g = 127 + (127*normal_y);
			//output.b = 127 + (127*normal_z);

			output.r *= luminance;
			output.g *= luminance;
			output.b *= luminance;
		}
		else
		{
			output.r = 0x1e;
			output.g = 0x90;
			output.b = 0xff;
		}

		return output;
	}
};

int main(int argc, char** argv)
{
	main_program prog;
	return prog(argc, argv);
}
