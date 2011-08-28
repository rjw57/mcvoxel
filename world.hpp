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

}

#endif // WORLD_HPP

