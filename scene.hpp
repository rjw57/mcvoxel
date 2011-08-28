#ifndef SCENE_HPP
#define SCENE_HPP

#include <boost/utility.hpp>

#include "datamodel.hpp"
#include "sky.hpp"
#include "world.hpp"

namespace scene
{

class scene : public boost::noncopyable
{
	public:
		typedef data::pixel<float>               pixel_f32;
		typedef data::sample_recorder<pixel_f32> sample_recorder;
		typedef data::image<sample_recorder>     image;

		// the world we loaded
		world::world world;

		// our sky dome
		sky::sky sky;

		scene();
		~scene();

		// initialise chains for sampling into image
		void initialise(image& im);

		// draw a sample
		void draw(image& im);

	protected:

		// current image plane co-ordinate
		float current_x_, current_y_;
};

}

#endif // SCENE_HPP
