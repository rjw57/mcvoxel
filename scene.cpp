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
	samples_.resize(w, h);

	// set the initial MH state to sensible defaults
	current_state_.image_x = w >> 1;
	current_state_.image_y = h >> 1;

	current_state_.light_path.bounces.clear();
	current_state_.light_path.eye_pos = camera_origin_;
	current_state_.light_path.from_sky = false;
	current_state_.light_path_lik = 0.f;

	// pre-cache the image plane perturbation parameters
	r2_ = sqrt(0.05f * samples_.width * samples_.height); // 5% of image area
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

	bounce.incident_dir = data::vec3_f32(incident.i, incident.j, incident.k);

	// Lambertian BRDF (cosine term is wrapped in luminance transfer function)
	bounce.brdf = 1.f;

	// choose a new ray direction
	bounce.bounce_dir = pt_sampling_cosine(norm);

	return true;
}

void scene::generate_initial_path(path& out_p) const
{
}

void scene::draw()
{
	typedef std::deque<bounce> subpath_t;

	perturb_image_loc(current_state_.image_x, current_state_.image_y, current_state_.image_x, current_state_.image_y);
	sample_recorder& record(samples_.at(floor(current_state_.image_x), floor(current_state_.image_y)));

	const int max_eye_path_len = 1;
	const int max_sky_path_len = 0;

	// make an eye ray path
	ray eye_ray;
	subpath_t eye_path;

	make_eye_ray(current_state_.image_x, current_state_.image_y, eye_ray);
	trace_ray(eye_ray, max_eye_path_len, std::back_inserter(eye_path));

	if(eye_path.size() == 0)
	{
		// we hit sky
		record(sky.sample(eye_ray));
		return;
	}

	// make a sky ray path
	ray sky_ray;
	subpath_t sky_path;

	make_sky_ray(sky_ray);
	trace_ray(sky_ray, max_sky_path_len, std::back_inserter(sky_path));

	pixel_f32 sky_sample = sky.sample(-sky_ray.i, -sky_ray.j, -sky_ray.k);

	// now for each possible path made by joining these together, work out the luminance - we always use
	// the first vertex of the eye ray and zero or more vertices of the sky ray
	for(size_t n_sky=0; n_sky <= sky_path.size(); ++n_sky)
	{
		// must have at least one item from the eye path since it makes sure the eye ray corresponds
		// to this pixel.
		for(size_t n_eye=1; n_eye <= eye_path.size(); ++n_eye)
		{
			path p;

			// we know ahead of time how many items will be in these collections. Should we move
			// to using vectors, we could make use of that knowledge here.
			// p.vertices.reserve(n_sky + n_eye);
			// p.known_visible.reserve(1 + n_sky + n_eye);

			// fill in sky part of the path
			p.from_sky = true;
			data::vec3_f32 sky_dir = data::vec3_f32(sky_ray.i, sky_ray.j, sky_ray.k);

			// fill in the eye part of the path
			p.eye_pos = camera_origin_;

			// copy the appropriate number of sky_ray items
			subpath_t::const_iterator sky_first(sky_path.begin());
			subpath_t::const_iterator sky_last(sky_path.begin()); std::advance(sky_last, n_sky);
			std::copy(sky_first, sky_last, std::back_inserter(p.bounces));

			// we want to copy the first n_eye items from the eye ray but insert them backwards
			subpath_t::const_iterator eye_first(eye_path.begin());
			subpath_t::const_iterator eye_last(eye_path.begin()); std::advance(eye_last, n_eye);
			std::reverse_copy(eye_first, eye_last, std::back_inserter(p.bounces));

			// we need to fix up the bounce directions at the join
			off_t first_eye_idx = n_sky;
			bounce& first_eye(p.bounces[first_eye_idx]);
			if(first_eye_idx > 0)
			{
				bounce& last_sky(p.bounces[first_eye_idx-1]);
				Vector sky_to_eye = pt_vector_normalise3(
						pt_vector_sub(first_eye.where.pos, last_sky.where.pos));
				first_eye.incident_dir = sky_to_eye;
				last_sky.bounce_dir = sky_to_eye;

				// FIXME: Update BRDF when we support them
			}
			else
			{
				first_eye.incident_dir = sky_dir;

				// FIXME: Update BRDF when we support them
			}

			// record the sky sample times luminance transfer
			record(sky_sample * luminance_transfer(p));
		}
	}
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
	while((new_x < 0) || (new_y < 0) || (new_x >= samples_.width) || (new_y >= samples_.height));
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
	fx -= (samples_.width)>>1; fy -= (samples_.height)>>1;

	float i(-fx), j(-fy), k(samples_.height);

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

float scene::luminance_transfer(const path& p) const
{
	// assert(p.vertices.size() + 1 == p.known_visible.size());

	// the full transfer is built up multiplicatively
	float transfer = 1.f;

	// Create a vector for the eye position
	Vector eye_pos = p.eye_pos;

	// there needs to be at least one bounce
	assert(p.bounces.size() > 0);

	// handle the last bounce -> eye
	{
		const bounce &last_bounce(p.bounces.back());

		// extract location and normal of bounce
		Vector p = last_bounce.where.pos, n = last_bounce.where.norm;

		// check visibility
		if(!visible(eye_pos, p))
			return  0.f;

		// fold in transfer
		transfer *= pt_vector_get_w(pt_vector_dot3(n, last_bounce.bounce_dir));
	}

	// now we need to handle the sky->first bounce assuming luminance came from the sky
	if(p.from_sky)
	{
		// find first bounce
		const bounce& first_bounce(p.bounces.front());

		// extract position and normal of bounce
		Vector p = first_bounce.where.pos, n = first_bounce.where.norm;

		// direction _to_ sky from bounce
		Vector sky_dir = pt_vector_neg(first_bounce.incident_dir);

		// if no sky visibility, no path
		if(!sky_visible(p, sky_dir))
			return 0.f;

		transfer *= pt_vector_get_w(pt_vector_dot3(n, sky_dir));
	}

	// handle the bounce-to-bounce differential beam throughput, starting from the first bounce
	if(p.bounces.size() > 1)
	{
		path::bounce_collection::const_iterator b_it(p.bounces.begin());

		// work out the position and normal of the first vertex.
		Vector start_point = b_it->where.pos, start_normal = b_it->where.norm;

		// for the remaining vertices...
		BOOST_FOREACH(const bounce& bounce, std::make_pair(++b_it, p.bounces.end()))
		{
			Vector end_point = bounce.where.pos, end_normal = bounce.where.norm;

			// if we ever fail a visibility test, the entire path is dark
			if(!visible(start_point, end_point))
				return 0.f;

			// compute 1/||end - start||^2
			Vector delta = pt_vector_sub(end_point, start_point);
			float recip_delta_sq = pt_vector_get_w(pt_vector_w_reciprocal(pt_vector_dot3(delta, delta)));

			// what is the normalised delta from start -> end, i.e. this bounce's incident dir
			Vector end_incident_dir = bounce.incident_dir;

			// work out cosine delta to normal for start and end
			transfer *= pt_vector_get_w(pt_vector_abs(pt_vector_mul(
							pt_vector_dot3(start_normal, end_incident_dir),
							pt_vector_dot3(end_normal, pt_vector_neg(end_incident_dir)))));

			// multiply in 1/||end - start||^2
			transfer *= recip_delta_sq;

			// prepare for the next iteration
			start_point = end_point;
			start_normal = end_normal;
		}
	}

	return transfer;
}

bool scene::visible(const Vector& a, const Vector& b) const
{
	// compute delta vector from a to b
	Vector delta = pt_vector_sub(b, a);
	Vector norm_delta = pt_vector_normalise3(delta);
	float mag_delta = pt_vector_get_w(pt_vector_w_sqrt(pt_vector_dot3(delta, delta)));

	// create a ray
	ray test_ray;
	make_ray(pt_vector_get_x(a), pt_vector_get_y(a), pt_vector_get_z(a),
			pt_vector_get_x(norm_delta), pt_vector_get_y(norm_delta), pt_vector_get_z(norm_delta),
			&test_ray);

	// cast it...
	octree::sub_location sl; data::block bl;
	data::vec3_f32 normal;
	if(!world::cast_ray(world, test_ray, sl, bl, normal))
		return false; // (shouldn't really happen I think)

	// where did we hit?
	Vector hit = pt_vector_make(sl.coords[0], sl.coords[1], sl.coords[2], 0.f);

	// how far away was that?
	Vector hit_delta = pt_vector_sub(hit, a);
	float mag_hit_delta = pt_vector_get_w(pt_vector_w_sqrt(pt_vector_dot3(hit_delta, hit_delta)));

	// bit of a HACK here...
	return mag_hit_delta >= mag_delta;
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
