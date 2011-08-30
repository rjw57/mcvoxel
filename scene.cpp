#include <cfloat>
#include <cmath>
#include <deque>
#include <iterator>
#include <utility>
#include <boost/foreach.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <rayslope/ray.h>

#include "scene.hpp"

namespace scene
{

// draw a sample uniformly distributed on [0, 1).
static float uniform_real();

scene::scene()
	: r1_(0.1f), r2_(1.f), log_r2_over_r1_(log(r2_/r1_))
{
	initialise(1, 1);
	set_camera(data::vec3_f32(0.f, 0.f, 0.f), 0.f, 0.f);
}

scene::~scene()
{ }

void scene::initialise(int w, int h)
{
	image_.resize(w, h);
	n_samples_ = 0;

	// set the initial MH state to sensible defaults
	current_state_.image_x = w >> 1;
	current_state_.image_y = h >> 1;

	current_state_.light_path.bounces.clear();
	current_state_.light_path.eye_pos = camera_origin_;
	current_state_.light_path.from_sky = false;
	current_state_.light_path_lik = 0.f;

	// pre-cache the image plane perturbation parameters
	r2_ = sqrt(0.05f * image_.width * image_.height); // 5% of image area
	log_r2_over_r1_ = log(r2_/r1_);
}

bool scene::bounce_ray(ray& incident, bounce& bounce) const
{
	world_position bounce_pos;
	octree::sub_location sub_loc;
	data::block block;
	data::vec3_f32 norm;

	// if we don't hit, terminate here
	if(!world::cast_ray(world, incident, sub_loc, block, norm))
		return false;

	// start filling in the bounce structure
	bounce.where.pos = data::vec3_f32(sub_loc.coords[0], sub_loc.coords[1], sub_loc.coords[2]);
	bounce.where.norm = norm;
	bounce.where.node_ext = sub_loc.node_extent;
	bounce.where.node_block = block;

	// choose a new ray direction
	//bounce.bounce_dir = pt_sampling_cosine(norm);
	bounce.bounce_dir = pt_sampling_hemisphere(norm);

	return true;
}

bool scene::generate_initial_path(path& p, float ix, float iy) const
{
	typedef std::deque<bounce> subpath_t;
	const int max_eye_path_len = 1;
	const int max_sky_path_len = 0;

	// clear the input path
	p.bounces.clear();
	p.from_sky = false;
	p.eye_pos = camera_origin_;

	// construct a path from the eye
	ray eye_ray; subpath_t eye_path;
	make_eye_ray(ix, iy, eye_ray);
	trace_ray(eye_ray, max_eye_path_len, std::back_inserter(eye_path));

	// do we only see sky?
	if(eye_path.size() == 0)
	{
		p.from_sky = true; p.sky_dir = pt_vector_neg(pt_vector_make(eye_ray.i, eye_ray.j, eye_ray.k, 0.f));
		return true;
	}

	// we need to reverse the eye ray so that it starts at the last bounce and heads towards the eye
	subpath_t rev_eye_path;

	// firstly, just copy the bounces
	std::reverse_copy(eye_path.begin(), eye_path.end(), std::back_inserter(rev_eye_path));
	
	// now fix up the bounce directions
	subpath_t::iterator now(rev_eye_path.begin());
	subpath_t::iterator next(rev_eye_path.begin()); ++next;
	for(; next != rev_eye_path.end(); ++next, ++now)
	{
		Vector now_pos = now->where.pos;
		Vector next_pos = next->where.pos;
		now->bounce_dir = pt_vector_normalise3(pt_vector_sub(next_pos, now_pos));
	}

	// last bounce is to eye
	rev_eye_path.back().bounce_dir = pt_vector_normalise3(pt_vector_sub(p.eye_pos, rev_eye_path.back().where.pos));
	
	// construct a path from the sky
	ray sky_ray; subpath_t sky_path;
	make_sky_ray(sky_ray);
	trace_ray(sky_ray, max_sky_path_len, std::back_inserter(sky_path));
	
	// record sky direction in path
	p.from_sky = true;
	p.sky_dir = data::vec3_f32(sky_ray.i, sky_ray.j, sky_ray.k); // vector pointing _from_ sky _to_ world/eye

	if(sky_path.size() > 0)
	{
		// now try all pairs of start/end bounces
		for(size_t n_sky = 1; n_sky <= sky_path.size(); ++n_sky)
		{
			for(size_t n_eye = 1; n_eye <= rev_eye_path.size(); ++n_eye)
			{
				// go from begining of sky path n_sky bounces along
				subpath_t::iterator sky_it(sky_path.begin()); std::advance(sky_it, n_sky-1);

				// go from n_eye bounces into eye path to end
				subpath_t::iterator eye_it(rev_eye_path.begin()); std::advance(eye_it, n_eye-1);

				if(visible(eye_it->where.pos, sky_it->where.pos))
				{
					// mutually visible! woo

					// need to fix the last sky bounce we added
					Vector first_eye = eye_it->where.pos;
					Vector last_sky = sky_it->where.pos;

					p.bounces.clear();
					std::copy(sky_path.begin(), ++sky_it, std::back_inserter(p.bounces));
					std::copy(eye_it, rev_eye_path.end(), std::back_inserter(p.bounces));

					p.bounces[n_sky-1].bounce_dir = pt_vector_normalise3(pt_vector_sub(first_eye, last_sky));

					return true;
				}
			}
		}
	}
	
	// can the last bounce on the eye path see the sky?
	if(!sky_visible(eye_path.back().where.pos, pt_vector_neg(p.sky_dir)))
		return false;

	// yes, copy reverse eye path -> path and return success.
	p.bounces = rev_eye_path;
	return true;
}

void scene::draw()
{
	if(n_samples_ == 0)
	{
		// initialise a path for the MH transport
		while(!generate_initial_path(current_state_.light_path, current_state_.image_x, current_state_.image_y)) { /* nop */ }
		current_state_.light_path_lik = luminance_transfer(current_state_.light_path);
	}

	typedef std::deque<bounce> subpath_t;

	perturb_image_loc(current_state_.image_x, current_state_.image_y, current_state_.image_x, current_state_.image_y);

	int px(floor(current_state_.image_x)), py(floor(current_state_.image_y));

	while(!generate_initial_path(current_state_.light_path, current_state_.image_x, current_state_.image_y)) { /* nop */ }
	current_state_.light_path_lik = luminance_transfer(current_state_.light_path);

	pixel_f32 sample_value(0,0,0);

	//bidirectional_mutate(current_state_.light_path);

	const path& current_path(current_state_.light_path);
	if(current_path.from_sky)
	{
		pixel_f32 sky_sample = sky.sample(-current_path.sky_dir.x, -current_path.sky_dir.y, -current_path.sky_dir.z);
		sample_value = sky_sample * current_state_.light_path_lik;
	}

	++n_samples_;

	// scale by number of pixels
	sample_value *= image_.width * image_.height;
	image_.at(px, py) += sample_value;
}

void scene::set_camera(const data::vec3_f32& origin, float yaw, float pitch)
{
	// convert yaw and pitch to radians
	cam_pitch_ = pitch * (2.f * 3.14159f / 360.f);
	cam_yaw_ = (180.f + yaw) * (2.f * 3.14159f / 360.f);

	// pre-calculate trig. values
	cos_pitch_ = cos(cam_pitch_); sin_pitch_ = sin(cam_pitch_);
	cos_yaw_ = cos(cam_yaw_); sin_yaw_ = sin(cam_yaw_);

	// save camera origin
	camera_origin_ = origin;
}

void scene::make_eye_ray(float fx, float fy, ray& out_ray) const
{
	fx -= (image_.width)>>1; fy -= (image_.height)>>1;

	float i(-fx), j(-fy), k(image_.height);

	float new_j = cos_pitch_*j - sin_pitch_*k, new_k = sin_pitch_*j + cos_pitch_*k;
	j = new_j; k = new_k;

	float new_i = cos_yaw_*i - sin_yaw_*k, new_k2 = sin_yaw_*i + cos_yaw_*k;
	i = new_i; k = new_k2;

	float mag = sqrt(i*i + j*j + k*k);
	make_ray(camera_origin_.x, camera_origin_.y, camera_origin_.z, i/mag, j/mag, k/mag, &out_ray);
}

void scene::make_sky_ray(ray& out_ray) const
{
	octree::location world_first(0,0,0), world_last(0,0,0);
	world::extent(world, world_first, world_last);

	float dx = world_last.x - world_first.x;
	float dy = world_last.y - world_first.y;
	float dz = world_last.z - world_first.z;
	float world_radius = sqrt(dx*dx + dy*dy + dz*dz);

	Vector sky_dir = pt_sampling_sphere();

	float aim_x = world_first.x + dx * uniform_real();
	float aim_y = world_first.y + dy * uniform_real();
	float aim_z = world_first.z + dz * uniform_real();

	float sky_x = aim_x + 32.f * world_radius * - pt_vector_get_x(sky_dir);
	float sky_y = aim_y + 32.f * world_radius * - pt_vector_get_y(sky_dir);
	float sky_z = aim_z + 32.f * world_radius * - pt_vector_get_z(sky_dir);

	make_ray(sky_x, sky_y, sky_z,
			pt_vector_get_x(sky_dir),
			pt_vector_get_y(sky_dir),
			pt_vector_get_z(sky_dir),
			&out_ray);
}

#define ALWAYS_VISIBLE 1

float scene::luminance_transfer(const path& p) const
{
	// assert(p.vertices.size() + 1 == p.known_visible.size());
	if(p.bounces.size() == 0)
	{
		// FIXME: check visibilirt
		return 1.f;
	}

	// the full transfer is built up multiplicatively
	float transfer = 1.f;

	// there needs to be at least one bounce
	assert(p.bounces.size() > 0);

	// handle the last bounce -> eye
	{
		const bounce &last_bounce(p.bounces.back());

		// extract location and normal of bounce
		Vector n = last_bounce.where.norm;

#if !ALWAYS_VISIBLE
		Vector eye_pos = p.eye_pos;
		Vector p = last_bounce.where.pos;
		// check visibility
		if(!visible(eye_pos, p))
			return  0.f;
#endif

		// fold in transfer
		transfer *= pt_vector_get_w(pt_vector_dot3(n, last_bounce.bounce_dir));
	}

	// now we need to handle the sky->first bounce assuming luminance came from the sky
	if(p.from_sky)
	{
		// find first bounce
		const bounce& first_bounce(p.bounces.front());

		// extract position and normal of bounce
		Vector norm = first_bounce.where.norm;

		// direction _to_ sky from bounce
		Vector sky_dir = pt_vector_neg(p.sky_dir);

#if !ALWAYS_VISIBLE
		Vector pos = first_bounce.where.pos;
		// if no sky visibility, no path
		if(!sky_visible(pos, sky_dir))
			return 0.f;
#endif

		transfer *= pt_vector_get_w(pt_vector_dot3(norm, sky_dir));
	}

	// handle the bounce-to-bounce differential beam throughput, starting from the first bounce
	if(p.bounces.size() > 1)
	{
		path::bounce_collection::const_iterator b_it(p.bounces.begin());

		// work out the position and normal of the first vertex.
		Vector start_point = b_it->where.pos, start_normal = b_it->where.norm;
		Vector start_bounce = b_it->bounce_dir;

		// for the remaining vertices...
		BOOST_FOREACH(const bounce& bounce, std::make_pair(++b_it, p.bounces.end()))
		{
			Vector end_point = bounce.where.pos, end_normal = bounce.where.norm;

#if !ALWAYS_VISIBLE
			// if we ever fail a visibility test, the entire path is dark
			if(!visible(start_point, end_point))
				return 0.f;
#endif

			// compute 1/||end - start||^2
			Vector delta = pt_vector_sub(end_point, start_point);
			float recip_delta_sq = pt_vector_get_w(pt_vector_w_reciprocal(pt_vector_dot3(delta, delta)));

			// what is the normalised delta from start -> end, i.e. this bounce's incident dir
			Vector end_bounce = bounce.bounce_dir;

			// work out cosine delta to normal for start and end
			transfer *= pt_vector_get_w(pt_vector_abs(pt_vector_mul(
							pt_vector_dot3(start_normal, start_bounce),
							pt_vector_dot3(end_normal, pt_vector_neg(start_bounce)))));

			// multiply in 1/||end - start||^2
			transfer *= recip_delta_sq;

			// prepare for the next iteration
			start_point = end_point;
			start_normal = end_normal;
			start_bounce = end_bounce;
		}
	}

	return transfer;
}

void scene::perturb_image_loc(float x, float y, float& new_x, float& new_y) const
{
	do
	{
		// choose angle and distance to perturb location
		float theta = 2.f * 3.14159f * uniform_real();
		float r = r2_ * exp(-log_r2_over_r1_ * uniform_real());
		new_x = x + r * cos(theta);
		new_y = y + r * sin(theta);
	}
	while((new_x < 0) || (new_y < 0) || (new_x >= image_.width) || (new_y >= image_.height));
}

inline float p_d(int len)
{
	assert(len > 0);

	if(len == 1)
		return 0.25f;
	if(len == 2)
		return 0.5f;

	return 1.f / static_cast<float>(1 << len);
}

inline std::pair<int, float> sample_new_subpath_len(int old_len)
{
	assert(old_len > 0);

	// hack: consider other lengths
	float v = uniform_real() * 0.8f;
	if(v <= 0.15f)
		return std::make_pair(old_len - 1, 0.15/0.8f);
	if(v <= 0.65f)
		return std::make_pair(old_len, 0.5f/0.8f);
	return std::make_pair(old_len + 1, 0.15f/0.8f);
}

inline std::pair<int, float> sample_subpath_len(int max_len)
{
	assert(max_len > 0);

	float norm = 0.f;
	for(int len=1; len<=max_len; ++len)
	{
		norm += p_d(len);
	}

	float v = uniform_real() * norm;
	float sum = 0.f;

	for(int len=1; len<=max_len; ++len)
	{
		float p = p_d(len);
		sum += p;
		if(v <= sum)
			return std::make_pair(len, p);
	}

	assert(false && "unreachable");
}

bool scene::bidirectional_mutate(path& p) const
{
	// how many vertices are there in this path? we don't count the eye but we do count the sky
	int n_vertices = p.bounces.size() + p.from_sky ? 1 : 0;

	// cannot mutate zero-length path (eye-only)
	if(n_vertices < 1)
		return false;

	// vertices are numbered from 0 (light) to k (immediately before eye). n_vertices = k+1
	// choose which subpath to delete
	std::pair<int, float> sp_len = sample_subpath_len(n_vertices);

	if(sp_len.first == n_vertices)
	{
		// throw away entire path and start again
		// ...
	}

	return false;
}

bool scene::visible(const Vector& a, const Vector& b) const
{
	// compute delta vector from a to b
	Vector delta = pt_vector_sub(b, a);
	Vector norm_delta = pt_vector_normalise3(delta);
	float mag_delta = pt_vector_get_w(pt_vector_w_sqrt(pt_vector_dot3(delta, delta)));

	// create a ray from a to b and b to a
	ray a_to_b, b_to_a;
	make_ray(pt_vector_get_x(a), pt_vector_get_y(a), pt_vector_get_z(a),
			pt_vector_get_x(norm_delta), pt_vector_get_y(norm_delta), pt_vector_get_z(norm_delta),
			&a_to_b);
	make_ray(pt_vector_get_x(b), pt_vector_get_y(b), pt_vector_get_z(b),
			-pt_vector_get_x(norm_delta), -pt_vector_get_y(norm_delta), -pt_vector_get_z(norm_delta),
			&b_to_a);

	// by default assume we don't intersect - this might be the case when a
	// or b is not on the surface and so one of the rays may fly striaght
	// off into space.
	float mag_a2b_hit_delta = FLT_MAX, mag_b2a_hit_delta = FLT_MAX;

	// cast it...
	octree::sub_location a2bl, b2al; data::block bl; data::vec3_f32 normal;

	if(world::cast_ray(world, a_to_b, a2bl, bl, normal))
	{
		// where did we hit and how far away was that?
		Vector a2b_hit = pt_vector_make(a2bl.coords[0], a2bl.coords[1], a2bl.coords[2], 0.f);
		Vector a2b_hit_delta = pt_vector_sub(a2b_hit, a);
		mag_a2b_hit_delta = pt_vector_get_w(pt_vector_w_sqrt(pt_vector_dot3(a2b_hit_delta, a2b_hit_delta)));
	}
	
	if(world::cast_ray(world, b_to_a, b2al, bl, normal))
	{
		// where did we hit and how far away was that?
		Vector b2a_hit = pt_vector_make(b2al.coords[0], b2al.coords[1], b2al.coords[2], 0.f);
		Vector b2a_hit_delta = pt_vector_sub(b2a_hit, b);
		mag_b2a_hit_delta = pt_vector_get_w(pt_vector_w_sqrt(pt_vector_dot3(b2a_hit_delta, b2a_hit_delta)));
	}

	// now the logic here is that if both of these deltas is >0.5 x the delta from a to b, there can't be anything
	// in the way
	return (mag_a2b_hit_delta > 0.5f*mag_delta) && (mag_b2a_hit_delta > 0.5f*mag_delta);
}

bool scene::sky_visible(const Vector& p, const Vector& dir) const
{
	// create a ray
	ray test_ray;
	make_ray(pt_vector_get_x(p), pt_vector_get_y(p), pt_vector_get_z(p),
			pt_vector_get_x(dir), pt_vector_get_y(dir), pt_vector_get_z(dir),
			&test_ray);

	// try to intersect world
	return !world::cast_ray(world, test_ray);
}

//// UTILITY FUNCTIONS ////

// use a fixed seed so results are replicable
static boost::mt19937 rng(42);
static boost::uniform_real<> uni_dist(0,1);
static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni_gen(rng, uni_dist);

static float uniform_real()
{
	return uni_gen();
}

}
