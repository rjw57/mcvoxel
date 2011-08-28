#ifndef SCENE_HPP
#define SCENE_HPP

#include <boost/utility.hpp>

#include "world.hpp"
#include "sky.hpp"

namespace scene
{

struct scene : public boost::noncopyable
{
	// the world we loaded
	world::world world;

	// our sky dome
	sky::sky sky;
};

}

#endif // SCENE_HPP
