#ifndef WORLD_HPP
#define WORLD_HPP

#include <iostream>

#include "octree.hpp"

namespace world
{

// a world is just a set of octrees
typedef std::vector<octree::crystalised_octree> world;

void load_world(const char* filename, world& w);
void save_cached_world(std::ostream& os, const world& w);
void load_cached_world(std::istream& is, world& w);

// cast a ray into the world. If we hit a non-air block, return true and set out_sub_loc and out_block as
// appropriate
bool cast_ray(const world& w, const ray& r, octree::sub_location& out_sub_loc, data::block& out_block,
		float& nx, float& ny, float &nz);

// as above, but don't return position
bool cast_ray(const world& w, const ray& r);

}

#endif // WORLD_HPP

