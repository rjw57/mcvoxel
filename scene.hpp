#ifndef SCENE_HPP
#define SCENE_HPP

#include <boost/utility.hpp>

#include "datamodel.hpp"
#include "sky.hpp"
#include "world.hpp"

namespace scene
{

struct scene : public boost::noncopyable
{
	typedef data::pixel<float>               pixel_f32;
	typedef data::sample_recorder<pixel_f32> sample_recorder;
	typedef data::image<sample_recorder>     image;

	// the world we loaded
	world::world world;

	// our sky dome
	sky::sky sky;

	// draw a sample
	void draw(image& im);
};

}

#endif // SCENE_HPP
