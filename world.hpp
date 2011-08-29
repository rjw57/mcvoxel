#ifndef WORLD_HPP
#define WORLD_HPP

#include <iostream>

#include "octree.hpp"

namespace world
{

// a world is just a set of octrees
typedef std::vector<octree::crystalised_octree> world;

void load(const char* filename, world& w);
void save_cached(std::ostream& os, const world& w);
void load_cached(std::istream& is, world& w);

// work out the minimum and maximum locations in the world
void extent(const world& w, octree::location& first_loc, octree::location& last_loc);

// cast a ray into the world. If we hit a non-air block, return true and set out_sub_loc and out_block as
// appropriate
bool cast_ray(const world& w, const ray& r, octree::sub_location& out_sub_loc, data::block& out_block, data::vec3_f32& normal);

// as above, but don't return position
bool cast_ray(const world& w, const ray& r);

}

#endif // WORLD_HPP

