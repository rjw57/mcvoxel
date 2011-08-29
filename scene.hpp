#ifndef SCENE_HPP
#define SCENE_HPP

#include <deque>
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
	float            pos_x, pos_y, pos_z;
	float            norm_x, norm_y, norm_z;

	octree::extent   node_ext;
	data::block      node_block;
};

// a path from a light source to the eye, optionally starting at the sky
struct path
{
	typedef std::deque<world_position> vertex_collection;
	typedef std::deque<bool>           flag_collection;

	// the vertices of the path starting from the light, ending at the bounce before the eye
	vertex_collection          vertices;

	// optimisation: there must be one element in this collection for each element in the vertex
	// collection plus one. If known_visible[i] is true, then we know that vertices[i] can see the vertex
	// (or sky) immediately before it in the path. If the last element is true, we know that the final
	// vertex in the path can see the eye.
	flag_collection            known_visible;

	// sky support, valid if from_sky is true
	bool                       from_sky;
	float                      sky_x, sky_y, sky_z; // direction _from_ sky

	// eye position is assumed to be (cam_x_, cam_y_, cam_z_)
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

		// choose a random sky ray - it is not guaranteed to intersect the world
		void make_sky_ray(ray& out_ray) const;

		// given a path through the world, compute the luminance transfer for that path
		// N.B. this will modify the known_visible flags in path as appropriate
		float luminance_transfer(path& p) const;

		// return true iff a and b are visible to each other
		bool visible(const Vector& a, const Vector& b) const;

		// return true iff the sky is visible from p in the direction dir
		bool sky_visible(const Vector& p, const Vector& dir) const;

		// trace a ray through the world for a maximum of max_bounces
		// add a world_position to out for every bounce
		template<typename OutputIterator>
		void trace_ray(ray& r, int max_bounces, OutputIterator out);
};

}

#define INSIDE_SCENE_HPP
#include "scene.tcc"
#undef INSIDE_SCENE_HPP

#endif // SCENE_HPP
