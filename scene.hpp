#ifndef SCENE_HPP
#define SCENE_HPP

#include <boost/utility.hpp>
#include <rayslope/ray.h>

#include "datamodel.hpp"
#include "octree.hpp"
#include "sky.hpp"
#include "world.hpp"
#include "util/vector.h"

namespace scene
{

struct world_position
{
	Vector           position;
	octree::location node_loc;
	octree::extent   node_ext;
};

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

		// initialise camera position
		void set_camera(float x, float y, float z, float yaw, float pitch);

		// initialise chains for sampling into image
		void initialise(int w, int h);

		// draw a sample
		void draw();

		// accessor
		const image& samples() const { return samples_; }

	protected:
		// samples currently drawn
		image samples_;

		// current image plane co-ordinate
		float current_x_, current_y_;

		// camera position, pose and cached parameters
		float cam_x_, cam_y_, cam_z_;
		float cam_pitch_, cos_pitch_, sin_pitch_;
		float cam_yaw_, cos_yaw_, sin_yaw_;

		// image plane peturbation parameters
		float r1_, r2_, log_r2_over_r1_;

		// choose a new image location, given the existing one
		void peturb_image_loc(float x, float y, float& new_x, float& new_y) const;

		// convert an image plane co-ordinate into an eye ray
		void make_eye_ray(float x, float y, ray& out_ray) const;
};

}

#endif // SCENE_HPP
