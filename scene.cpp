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
	: current_x_(0.f), current_y_(0.f)
	, r1_(0.1f), r2_(1.f), log_r2_over_r1_(log(r2_/r1_))
{
	initialise(1, 1);
	set_camera(data::vec3_f32(0.f, 0.f, 0.f), 0.f, 0.f);
}

scene::~scene()
{ }

void scene::initialise(int w, int h)
{
	samples_.resize(w, h);

	current_x_ = 0.5f * samples_.width;
	current_y_ = 0.5f * samples_.height;

	r2_ = sqrt(0.05f * samples_.width * samples_.height); // 5% of image area
	log_r2_over_r1_ = log(r2_/r1_);
}

void scene::draw()
{
	typedef std::deque<world_position> subpath_t;

	peturb_image_loc(current_x_, current_y_, current_x_, current_y_);
	sample_recorder& record(samples_.at(floor(current_x_), floor(current_y_)));

	const int max_eye_path_len = 1;
	const int max_sky_path_len = 0;

	// make an eye ray path
	ray eye_ray;
	subpath_t eye_path;

	make_eye_ray(current_x_, current_y_, eye_ray);
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
			p.sky_dir = data::vec3_f32(sky_ray.i, sky_ray.j, sky_ray.k);

			// fill in the eye part of the path
			p.eye_pos = camera_origin_;

			// copy the appropriate number of sky_ray items
			subpath_t::const_iterator sky_first(sky_path.begin());
			subpath_t::const_iterator sky_last(sky_path.begin()); std::advance(sky_last, n_sky);
			std::copy(sky_first, sky_last, std::back_inserter(p.vertices));

			// we know that all the sky path vertices can see each other
			for(size_t i=0; i < n_sky; ++i)
				p.known_visible.push_back(true);

			// we don't know if the join between the sky path and eye path is visible
			p.known_visible.push_back(false);

			// we know that all the eye path vertices can see each other
			for(size_t i=0; i < n_eye; ++i)
				p.known_visible.push_back(true);

			// we want to copy the first n_eye items from the eye ray but insert them backwards
			subpath_t::const_iterator eye_first(eye_path.begin());
			subpath_t::const_iterator eye_last(eye_path.begin()); std::advance(eye_last, n_eye);
			std::reverse_copy(eye_first, eye_last, std::back_inserter(p.vertices));

			// record the sky sample times luminance transfer
			record(sky_sample * luminance_transfer(p));
		}
	}
}

void scene::peturb_image_loc(float x, float y, float& new_x, float& new_y) const
{
	do
	{
		// choose angle and distance to peturb location
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
	assert(p.vertices.size() + 1 == p.known_visible.size());

	float contribution = 1.f;

	// Create a vector for the eye position
	Vector eye_pos = p.eye_pos;

	// if no intermediate vertices, just do an eye sky test
	if(p.vertices.size() == 0)
	{
		if(!p.from_sky)
			return 0.;

		// work out direction _to_ sky from eye.
		Vector sky_dir = pt_vector_neg(p.sky_dir);

		if(p.known_visible.front())
			return 1.f;

		if(sky_visible(eye_pos, sky_dir))
			return 1.f;

		return 0.f;
	}

	// handle the last vertex -> eye
	const world_position& last_v(p.vertices.back());

	// work out the position and normal of the first vertex
	Vector pos = last_v.pos;
	Vector normal = last_v.norm;

	// if eye not visible, the whole path is invalid
	if(!p.known_visible.back() && !visible(eye_pos, pos))
		return  0.f;

	// fold in last vertex -> eye contribution
	Vector to_eye = pt_vector_normalise3(pt_vector_sub(eye_pos, pos));
	contribution *= pt_vector_get_w(pt_vector_dot3(normal, to_eye));

	// now we need to handle the sky->first vertex assuming luminance came from the sky
	if(p.from_sky)
	{
		// find first vertex
		const world_position& v(p.vertices.front());

		// work out the position and normal of the first vertex
		Vector pos = v.pos;
		Vector normal = v.norm;

		// direction _to_ sky from vertex
		Vector sky_dir = pt_vector_neg(p.sky_dir);

		// if no sky visibility, no path
		if(!p.known_visible.front() && !sky_visible(pos, sky_dir))
			return 0.f;

		contribution *= pt_vector_get_w(pt_vector_dot3(normal, sky_dir));
	}

	// handle the vertex-to-vertex differential beam throughput, starting from the first vertex
	if(p.vertices.size() > 1)
	{
		path::vertex_collection::const_iterator v_it(p.vertices.begin());
		path::flag_collection::const_iterator known_visible_it(p.known_visible.begin());

		// work out the position and normal of the first vertex.
		Vector start_point = v_it->pos;
		Vector start_normal = v_it->norm;

		// for the remaining vertices...
		BOOST_FOREACH(const world_position& wp, std::make_pair(++v_it, p.vertices.end()))
		{
			++known_visible_it;

			Vector end_point = wp.pos;
			Vector end_normal = wp.norm;

			// if we ever fail a visibility test, the entire path is dark
			if(!(*known_visible_it) && !visible(start_point, end_point))
				return 0.f;

			float vertex_contribution = 1.f;

			// compute 1/||end - start||^2
			Vector delta = pt_vector_sub(end_point, start_point);
			float recip_delta_sq = pt_vector_get_w(pt_vector_w_reciprocal(pt_vector_dot3(delta, delta)));

			// what is the normalised delta from start -> end
			Vector norm_delta = pt_vector_normalise3(delta);

			// work out cosine delta to normal for start and end
			vertex_contribution *= pt_vector_get_w(pt_vector_abs(pt_vector_mul(
							pt_vector_dot3(start_normal, norm_delta),
							pt_vector_dot3(end_normal, pt_vector_neg(norm_delta)))));

			// multiply in 1/||end - start||^2
			vertex_contribution *= recip_delta_sq;

			// fold this vertex's contribution into the whole
			contribution *= vertex_contribution;

			// prepare for the next iteration
			start_point = end_point;
			start_normal = end_normal;
		}
	}

	return contribution;
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
