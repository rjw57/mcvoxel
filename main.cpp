#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <iostream>
#include <stdexcept>

#include <octree.h>

#include <nbt/nbt.hpp>
#include <mc/blocks.hpp>
#include <mc/level.hpp>
#include <mc/utils.hpp>
#include <mc/world.hpp>

#include <rayslope/aabox.h>
#include <rayslope/ray.h>
#include <rayslope/slope.h>
#include <rayslope/slopeint_mul.h>

struct block_pos
{
	int x, y, z;
};

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

	bool ray_int(ray* r, float* t, block_pos* pos)
	{
		return ray_node_int(r, this->root(), 0, 0, 0, this->size(), t, pos);
	}

	protected:

	typedef typename Octree<T, AS>::Node OctreeNode;
	typedef typename Octree<T, AS>::Branch OctreeBranch;
	typedef typename Octree<T, AS>::Aggregate OctreeAggregate;
	typedef typename Octree<T, AS>::Leaf OctreeLeaf;

	bool ray_node_int(ray* r, OctreeNode* n, int x0, int y0, int z0, int size, float* t, block_pos* pos)
	{
		switch(n->type())
		{
			case Octree<T, AS>::BranchNode:
				return ray_branch_node_int(r, reinterpret_cast<OctreeBranch*>(n), x0, y0, z0, size, t, pos);
			case Octree<T, AS>::AggregateNode:
				return ray_agg_node_int(r, reinterpret_cast<OctreeAggregate*>(n), x0, y0, z0, size, t, pos);
			case Octree<T, AS>::LeafNode:
				return ray_leaf_node_int(r, reinterpret_cast<OctreeLeaf*>(n), x0, y0, z0, size, t, pos);
			default:
				throw std::runtime_error("unknown node type.");
		}
	}

	bool ray_branch_node_int(ray* r, OctreeBranch* n, int x0, int y0, int z0, int size, float* t, block_pos* pos)
	{
		aabox box;
		make_aabox(x0, y0, z0, x0+size, y0+size, z0+size, &box);
		if(!slope(r, &box))
			return false;

		//std::cout << "in branch @ (" << x0 << "," << y0 << "," << z0 << "), size: " << size << std::endl;

		size >>= 1;

		bool found_intersection = false;
		float temp_t;
		block_pos temp_pos;

		for(int dx=0; dx<=1; ++dx)
			for(int dy=0; dy<=1; ++dy)
				for(int dz=0; dz<=1; ++dz)
				{
					OctreeNode* child_node(n->child(dx, dy, dz));
					if(child_node == NULL)
						continue;

					if(!ray_node_int(r, child_node, x0 + dx*size, y0 + dy*size, z0 + dz*size,
								size, &temp_t, &temp_pos))
						continue;

					if(found_intersection)
					{
						if(temp_t < *t)
						{
							*t = temp_t;
							*pos = temp_pos;
						}
					}
					else
					{
						found_intersection = true;
						*t = temp_t;
						*pos = temp_pos;
					}
				}

		return found_intersection;
	}

	bool ray_agg_node_int(ray* r, OctreeAggregate* n, int x0, int y0, int z0, int size, float* t, block_pos* pos)
	{
		aabox box;

		make_aabox(x0, y0, z0, x0+size, y0+size, z0+size, &box);
		if(!slope(r, &box))
			return false;

		//std::cout << "in aggregate @ (" << x0 << "," << y0 << "," << z0 << "), size: " << size << std::endl;

		bool found_intersection = false;
		float temp_t;

		for(int dx=0; dx<size; ++dx)
			for(int dy=0; dy<size; ++dy)
				for(int dz=0; dz<size; ++dz)
				{
					if(n->value(dx, dy, dz) == this->emptyValue())
						continue;

					make_aabox(x0+dx, y0+dy, z0+dz, x0+dx+1, y0+dy+1, z0+dz+1, &box);
					if(!slopeint_mul(r, &box, &temp_t))
						continue;

					if(found_intersection)
					{
						if(temp_t < *t)
						{
							pos->x = x0 + dx;
							pos->y = y0 + dy;
							pos->z = z0 + dz;
							*t = temp_t;
						}
					}
					else
					{
						found_intersection = true;
						pos->x = x0 + dx;
						pos->y = y0 + dy;
						pos->z = z0 + dz;
						*t = temp_t;
					}
				}

		//if(found_intersection)
		//	std::cout << " - intersection @ (" << best_x << "," << best_y << "," << best_z << ")" << std::endl;
		return found_intersection;
	}

	bool ray_leaf_node_int(ray* r, OctreeLeaf* n, int x0, int y0, int z0, int size, float* t, block_pos* pos)
	{
		if(n->value() == this->emptyValue())
			return false;

		//std::cout << "look at leaf @ (" << x0 << "," << y0 << "," << z0 << "), size: " << size << std::endl;

		aabox box;
		make_aabox(x0, y0, z0, x0 + size, y0 + size, z0 + size, &box);
		if(slopeint_mul(r, &box, t))
		{
			pos->x = x0;
			pos->y = y0;
			pos->z = z0;
			return true;
		}

		return false;
	}
};

typedef RayCastingOctree<uint8_t, 8> block_octree;

struct load_level : public std::unary_function<const mc::utils::level_coord&, void>
{
	load_level(boost::shared_ptr<mc::region> region, block_octree& octree)
		: region(region), octree(octree)
	{ }

	bool is_surrounded(nbt::Byte* bp, int x, int y, int z)
	{
		bool surrounded = true;
		for(int dx=-1; dx <= 1; dx += 2)
			for(int dy=-1; dy <= 1; dy += 2)
				for(int dz=-1; dz <= 1; dz += 2)
				{
					if(bp[y + (128*z) + (128*16*x)] == mc::Air)
						surrounded = false;
				}
		return surrounded;
	}

	void operator() (const mc::utils::level_coord& coord)
	{
		boost::shared_ptr<mc::level_info> level_info(new mc::level_info(region, coord));
		mc::level level(level_info);
		mc::dynamic_buffer region_buffer(mc::region::CHUNK_MAX);
		level.read(region_buffer);

		boost::shared_ptr<nbt::ByteArray> blocks(level.get_blocks());

		int32_t x(0), y(0), z(0);

		int32_t start_x(level_info->get_x() + 1024), start_z(level_info->get_z() + 1024);
		//std::cout << " - (" << start_x << "," << start_z << ")" << std::endl;

		// convert surrounded blocks into air
		for(y=1; y<127; ++y)
			for(z=1; z<15; ++z)
				for(x=1; x<15; ++x)
				{
					if(is_surrounded(blocks->values, x, y, z))
					{
						blocks->values[y + (128*z) + (128*16*x)] = mc::Air;
					}
				}

		x = y = z = 0;

		for(int32_t idx(0); idx < blocks->length; ++idx)
		{
			// if appropriate, insert this block into the octree
			Byte block_id = blocks->values[idx];
			if(block_id != mc::Air)
			{
				if((x < 0) || (y < 0) || (z < 0))
					break;
				if((x >= octree.size()) || (y >= octree.size()) || (z >= octree.size()))
					break;
				octree.set(start_x + x, y, start_z + z, block_id);
				//std::cout << "(" << x << "," << y << "," << z << ") ";
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
	block_octree&                   octree;
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

	const int w=640, h=480;
	boost::shared_ptr<out_pixel> pixels(new out_pixel[w*h]);

	int32_t x(0), y(h);
	for(uint32_t idx=0; idx<w*h; ++idx)
	{
		out_pixel* out = pixels.get() + idx;
		float fx(x), fy(y), t(0.f);
		block_pos pos;
		ray r;

		fx -= w>>1; fy -= h>>1;

		float i(fx), j(fy), k(0.33f*w);

		float mag = sqrt(i*i + j*j + k*k);

		make_ray(1044, 100, 1000, i/mag, j/mag, k/mag, &r);
		if(octree.ray_int(&r, &t, &pos))
		{
			out->r = (octree.at(pos.x, pos.y, pos.x) << 3) && 0xff;
			out->g = pos.y;
			out->b = t;
			//std::cout << t << std::endl;
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
