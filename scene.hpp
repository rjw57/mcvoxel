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
	// position on world and normal at this point
	data::vec3_f32  pos, norm;

	// extent and content of octree node at this point
	octree::extent   node_ext;
	data::block      node_block;
};

// a record of one bounce
struct bounce
{
	world_position where;        // where the bounce ocurred
	data::vec3_f32 incident_dir; // incident direction (pointing into surface)
	data::vec3_f32 bounce_dir;   // bounce direction (pointing out of surface)
	float          brdf;         // the value of the BRDF for this bounce
};

// a path from a light source to the eye. The path either starts at the sky (in which case from_sky is true)
// or it starts from the first element of bounces (and hence the incident_dir field in that element is
// ignored).
struct path
{
	typedef std::deque<bounce> bounce_collection;

	// the bounces of the path starting from the light, ending at the bounce before the eye
	bounce_collection          bounces;

#if 0
	// optimisation: there must be one element in this collection for each element in the vertex
	// collection plus one. If known_visible[i] is true, then we know that vertices[i] can see the vertex
	// (or sky) immediately before it in the path. If the last element is true, we know that the final
	// vertex in the path can see the eye.
	flag_collection            known_visible;
#endif

	// if true, the incident_direction of the first bounce is from the sky. Otherwise, ignore the sky.
	bool                       from_sky;

	// eye position: included for completeness although for the moment this is just the scene camera
	// location.
	data::vec3_f32             eye_pos;
};

// a structure giving the current state of a MH chain
struct mh_state
{
	float image_x, image_y; // location in image plane
	path  light_path;       // light transport path through world
	float light_path_lik;   // the likelihood of the path == luminance_transfer(path)
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
		void set_camera(const data::vec3_f32& origin, float yaw, float pitch);

		// initialise chains for sampling into image
		void initialise(int w, int h);

		// draw a sample
		void draw();

		// accessor
		const image& samples() const { return samples_; }

	protected:
		// samples currently drawn
		image samples_;

		// the current state of the MH chain
		mh_state current_state_;

		// camera position, pose and cached parameters
		data::vec3_f32 camera_origin_;
		float cam_pitch_, cos_pitch_, sin_pitch_;
		float cam_yaw_, cos_yaw_, sin_yaw_;

		// image plane perturbation parameters
		float r1_, r2_, log_r2_over_r1_;

		// choose a new image location, given the existing one
		void perturb_image_loc(float x, float y, float& new_x, float& new_y) const;

		// convert an image plane co-ordinate into an eye ray
		void make_eye_ray(float x, float y, ray& out_ray) const;

		// choose a random sky ray - it is not guaranteed to intersect the world
		void make_sky_ray(ray& out_ray) const;

		// return true iff a and b are visible to each other
		bool visible(const Vector& a, const Vector& b) const;

		// return true iff the sky is visible from p in the direction dir
		bool sky_visible(const Vector& p, const Vector& dir) const;

		// given a light ray, see if we intersect the world (return false if we do not). If we
		// intersect, sample a bounce ray and return it's location, direction and probability
		// in out_bounce
		bool bounce_ray(ray& incident, bounce& bounce) const;

		// trace a ray through the world for a maximum of max_bounces
		// add a bounce to out for every bounce
		template<typename OutputIterator>
		void trace_ray(ray& r, int max_bounces, OutputIterator out);

		// MH-related methods

		// generate an initial light -> eye path using bidirectional path tracing
		void generate_initial_path(path& out_p) const;

		// given a path through the world, compute the luminance transfer for that path
		// i.e., the value of f(X) for that path.
		float luminance_transfer(const path& p) const;
};

}

#define INSIDE_SCENE_HPP
#include "scene.tcc"
#undef INSIDE_SCENE_HPP

#endif // SCENE_HPP
