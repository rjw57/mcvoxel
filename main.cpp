#include <algorithm>
#include <cmath>
#include <deque>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <glib.h>

#include <boost/foreach.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <nbt/nbt.hpp>
#include <mc/blocks.hpp>
#include <mc/level.hpp>
#include <mc/utils.hpp>
#include <mc/world.hpp>

#include <rayslope/aabox.h>
#include <rayslope/ray.h>
#include <rayslope/slope.h>
#include <rayslope/slopeint_mul.h>

#include <hdrloader/hdrloader.h>

#include <util/sampling.h>

#include "datamodel.hpp"
#include "io.hpp"
#include "octree.hpp"

const int CHUNK_SIZE = 16;

static float uniform_real()
{
	static boost::mt19937 generator(42);
	static boost::uniform_real<> uni_dist(0,1);
	static boost::variate_generator<boost::mt19937&, boost::uniform_real<> > uni(generator, uni_dist);
	return uni();
}

static void sample_spherical_ray(float x, float y, float z, ray* r)
{
	float u(uniform_real()), v(uniform_real());
	float theta = 2.f * 3.14159f * u;
	float phi = acos(2.f * v - 1.f);
	float st = sin(theta), ct = cos(theta);
	float sp = sin(phi), cp = cos(phi);
	make_ray(x, y, z, ct*sp, st*sp, cp, r);
}

template<typename T>
T rgb2y(const data::pixel<T>& p)
{
	return 0.299f*p.r + 0.114f*p.g + 0.587f*p.b;
}

struct load_level : public std::unary_function<const mc::utils::level_coord&, void>
{
	load_level(boost::shared_ptr<mc::region> region, octree::octree<data::block>& octree)
		: region(region), octree(octree)
	{ }

	void operator() (const mc::utils::level_coord& coord)
	{
		boost::shared_ptr<mc::level_info> level_info(new mc::level_info(region, coord));
		mc::level level(level_info);
		mc::dynamic_buffer region_buffer(mc::region::CHUNK_MAX);
		level.read(region_buffer);

		boost::shared_ptr<nbt::ByteArray> blocks(level.get_blocks());

		int32_t x(0), y(0), z(0);

		int32_t start_x(coord.get_x() * 16), start_z(coord.get_z() * 16);
		octree::location offset(octree.first_loc());

		for(int32_t idx(0); idx < blocks->length; ++idx)
		{
			// if appropriate, insert this block into the octree
			Byte block_id = blocks->values[idx];
			if(block_id != mc::Air)
			{
				octree::location loc(x + start_x, y, z + start_z);
				octree.set(loc + offset, data::block(block_id));
			}

			// go to next co-ord
			++y;
			if(y >= 128)
			{
				y = 0; ++z;
				if(z >= 16)
				{
					z = 0; ++x;
				}
			}
		}
	}

	boost::shared_ptr<mc::region>   region;
	octree::octree<data::block>&        octree;
};


template<typename T>
bool is_surrounded(octree::octree<T>& tree, const octree::location& loc)
{
	if(tree.get(loc.x-1, loc.y, loc.z) == mc::Air) return false;
	if(tree.get(loc.x+1, loc.y, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y-1, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y+1, loc.z) == mc::Air) return false;
	if(tree.get(loc.x, loc.y, loc.z-1) == mc::Air) return false;
	if(tree.get(loc.x, loc.y, loc.z+1) == mc::Air) return false;
	return true;
}

template<typename T>
void sparsify(octree::octree<T>& tree)
{
	octree::location first_loc = tree.first_loc();
	octree::location last_loc(first_loc);

	last_loc.x += tree.size();
	last_loc.y += tree.size();
	last_loc.z += tree.size();

	for(int x=first_loc.x+1; x<last_loc.x-1; ++x)
		for(int y=1; y<128-1; ++y)
			for(int z=first_loc.z+1; z<last_loc.z-1; ++z)
			{
				if(is_surrounded(tree, octree::location(x, y, z)))
					tree.set(x,y,z,0xff);
			}
}

static float gamma_est(float sample_mean, float sample_log_mean)
{
	if(sample_mean <= 0.f)
		return sample_mean;

	float s = log(sample_mean) - sample_log_mean;

	float alpha = (s-3.f)*(s-3.f)+24.f*s;
	if(alpha < 0.f)
		return sample_mean;

	float k = (3.f - s + sqrt(alpha)) / (12.f*s);
	if(k < 1.f)
		return sample_mean;

	float theta = sample_mean / k;

	return (k-1.f) * theta;
}

struct main_program
{
	std::vector<octree::crystalised_octree>    crystal_octrees;
	HDRLoaderResult                            sky_light_probe;
	float                                      sky_max_lum; // for importance sampling
	float                                      sky_integral; // for importance sampling

	data::image<data::pixel<float> >           terrain_colour;
	data::image<float>                         terrain_alpha;

	typedef data::pixel<uint8_t> pixel_u8;
	typedef data::pixel<float> pixel_f32;
	typedef data::sample_recorder<pixel_f32> sample_recorder;

	bool cast_ray(const ray& r, pixel_f32& transmission) const
	{
		octree::sub_location temp_sub_loc;
		data::block temp_hit_block;
		return cast_ray(r, temp_sub_loc, temp_hit_block, transmission);
	}

	bool cast_ray(const ray& r, octree::sub_location& out_sub_loc, data::block& out_block, pixel_f32& transmission) const
	{
		ray cont_ray(r);
		transmission = pixel_f32(1,1,1);

		while(cast_ray_raw(cont_ray, out_sub_loc, out_block))
		{
			// we hit something, is it transparent?
			if((out_block.id != mc::Torch) && (out_block.id != mc::Water) && (out_block.id != mc::StationaryWater) && (out_block.id != mc::Glass))
				return true;

			// handle transparent blocks

			float hit_x = out_sub_loc.coords[0];
			float hit_y = out_sub_loc.coords[1];
			float hit_z = out_sub_loc.coords[2];
			const octree::extent& node_ext = out_sub_loc.node_extent;

			// work out where the ray emerges
			ray emerge_ray;
			make_ray(hit_x+cont_ray.i*node_ext.size*4.f, hit_y+cont_ray.j*node_ext.size*4.f, hit_z+cont_ray.k*node_ext.size*4.f,
					-cont_ray.i, -cont_ray.j, -cont_ray.k, &emerge_ray);
			aabox node_box(node_ext.make_aabox());

			float emerge_distance = 0.f;
			const float epsilon = 1e-2f;

			if(!slopeint_mul(&emerge_ray, &node_box, &emerge_distance))
			{
				// we didn't hit, this is usually due to hitting just on a corner/edge, just nudge us slightly...
				emerge_distance = epsilon;
				std::cout << '.' << std::flush;

				make_ray(hit_x+epsilon*cont_ray.i, hit_y+epsilon*cont_ray.j, hit_z+epsilon*cont_ray.k,
						cont_ray.i, cont_ray.j, cont_ray.k, &cont_ray);
			}
			else
			{
				float emerge_x = emerge_ray.x + emerge_ray.i*emerge_distance;
				float emerge_y = emerge_ray.y + emerge_ray.j*emerge_distance;
				float emerge_z = emerge_ray.z + emerge_ray.k*emerge_distance;

				make_ray(emerge_x-epsilon*cont_ray.i, emerge_y-epsilon*cont_ray.j, emerge_z-epsilon*cont_ray.k,
						cont_ray.i, cont_ray.j, cont_ray.k, &emerge_ray);

				float dx = emerge_ray.x - hit_x;
				float dy = emerge_ray.y - hit_y;
				float dz = emerge_ray.z - hit_z;
				emerge_distance = sqrt(dx*dx + dy*dy + dz*dz);

				cont_ray = emerge_ray;
			}
			
			make_ray(hit_x+epsilon*cont_ray.i, hit_y+epsilon*cont_ray.j, hit_z+epsilon*cont_ray.k,
					cont_ray.i, cont_ray.j, cont_ray.k, &cont_ray);
			//make_ray(hit_x, hit_y, hit_z, cont_ray.i, cont_ray.j, cont_ray.k, &cont_ray);

			if((out_block.id == mc::Water) || (out_block.id == mc::StationaryWater))
			{
				pixel_f32 block_t(powf(0.8, emerge_distance), powf(0.9, emerge_distance), powf(0.95, emerge_distance));
				transmission = transmission * block_t;
			}
		}

		return false;
	}

	bool cast_ray_raw(const ray& r, octree::sub_location& out_sub_loc, data::block& out_block) const
	{
		float min_dist_sq = -1.f;

		BOOST_FOREACH(const octree::crystalised_octree& tree, crystal_octrees)
		{
			octree::sub_location temp_sub_loc;

			if(!tree.ray_intersect<data::block>(r, temp_sub_loc))
				continue;

			float dx = (temp_sub_loc.coords[0] - r.x);
			float dy = (temp_sub_loc.coords[1] - r.y);
			float dz = (temp_sub_loc.coords[2] - r.z);
			float dist_sq = dx*dx + dy*dy + dz*dz;

			if((min_dist_sq < 0.f) || (dist_sq < min_dist_sq))
			{
				min_dist_sq = dist_sq;
				out_sub_loc = temp_sub_loc;
				out_block = data::block(tree.get(temp_sub_loc.node_extent.loc));
			}
		}

		return (min_dist_sq >= 0.f);
	}

	// draw a sample from f(t,p) = sky_lum(t,p)/[\int sky_lum(t,p) dt dp] == sky_lum(t,p)/sky_norm
	// return sky(t,p) / sky_lum;
	data::pixel<float> sample_sky_pdf(float x, float y, float z, ray& out_ray) const
	{
		// we sample via rejection sampling: consider the uniform spherical PDF g(t,p) = 1/(4*pi).
		// since we know that f(t,p) <= sky_max_lum/sky_norm, then f(t,p) <= 4*pi*sky_max_lum/sky_norm * g(t,p).
		// and hence f(t,p) <= M g(x) for some constant M.
		//
		// sampling from g() is easy we may therefore use rejection sampling to sample from f
		for(int try_idx=0; try_idx<16; ++try_idx)
		{
			// sample from g()
			sample_spherical_ray(x, y, z, &out_ray);

			// calculate f * sky_norm
			data::pixel<float> sky = sample_sky(out_ray);
			float f = rgb2y(sky);

			// accept? (u < f / (M g(x))) === (u < sky_norm * f / sky_max_lum)
			if(uniform_real() < (f / sky_max_lum))
			{
				// accepted
				return sky / f;
			}
		}

		return data::pixel<float>(0,0,0);
	}

	data::pixel<float> sample_sky(const Vector& v, float* sa = NULL) const
	{
		ray temp_ray;
		make_ray(0, 0, 0, pt_vector_get_x(v), pt_vector_get_y(v), pt_vector_get_z(v), &temp_ray);
		return sample_sky(temp_ray, sa);
	}

	static float sky_pixel_solid_angle(int x, int y, int w, int h)
	{
		float pix_dx = 1.f / static_cast<float>(w);
		float pix_dy = 1.f / static_cast<float>(h);

		float min_x = static_cast<float>(x) * pix_dx;
		float min_y = static_cast<float>(y) * pix_dy;

		float max_x = static_cast<float>(x+1) * pix_dx;
		float max_y = static_cast<float>(y+1) * pix_dy;

		float min_theta, min_phi;
		normalised_texture_to_spherical(min_x, min_y, min_theta, min_phi);

		float max_theta, max_phi;
		normalised_texture_to_spherical(max_x, max_y, max_theta, max_phi);

		return fabs(cos(min_theta) - cos(max_theta)) * fabs(max_phi - min_phi);
	}

	static void normalised_texture_to_spherical(float x, float y, float& r_theta, float& r_phi)
	{
		r_theta = y * 3.14159f;
		r_phi = ((x * 2.f) - 1.f) * 3.14159f;
	}

	static void spherical_to_normalised_texture(float theta, float phi, float& r_x, float& r_y)
	{
		r_x = std::max(0.f, std::min(1.f, 0.5f * (1.f + (phi / 3.14159f))));
		r_y = std::max(0.f, std::min(1.f, theta / 3.14159f));
	}

	data::pixel<float> sample_sky(const ray& r, float* sa = NULL) const
	{
#if 0
		float c1 = fabs(r.i), c2 = fabs(r.k), c3 = r.j;
		if(c1 < 0.f)
			c1 = 0.f;
		if(c2 < 0.f)
			c2 = 0.f;
		if(c3 < 0.f)
			c3 = 0.f;

		//c1 = 0.f;
		//c2 = 0.f;

		float strength = 30.f;
		return (data::pixel<float>(1,0,0) * powf(c1,strength)
			+ data::pixel<float>(0,1,0) * powf(c2,strength)
			+ data::pixel<float>(0,0,1) * powf(c3,strength));
#endif

		float m = sqrt(r.i*r.i+r.j*r.j+r.k*r.k);
		float theta = acos(r.j/m);
		float phi = atan2(r.k, r.i);

		float sky_x, sky_y;
		spherical_to_normalised_texture(theta, phi, sky_x, sky_y);

		if(sa != NULL)
		{
			*sa = sky_pixel_solid_angle(sky_x, sky_y, sky_light_probe.width, sky_light_probe.height);
		}

		data::pixel<float> output;
		//output.r = sky_x; output.g = sky_y; output.b = 0.f;

		int pix_x = static_cast<int>(sky_x * sky_light_probe.width) % sky_light_probe.width;
		int pix_y = static_cast<int>(sky_y * sky_light_probe.height) % sky_light_probe.height;

		float *p_sky_pix = &(sky_light_probe.cols[3*(pix_x+(pix_y*sky_light_probe.width))]);
		output.r = p_sky_pix[0];
		output.g = p_sky_pix[1];
		output.b = p_sky_pix[2];

		return output;
	}

	void sample_pixel(sample_recorder& record, float fx, float fy, int w, int h) const
	{
		ray r;

		fx -= w>>1; fy -= h>>1;

		//float i(fx), k(fy), j(-0.5f*h);
		float i(-fx), j(fy), k(h);

		float pitch = 10.f * (2.f*3.14159f/360.f);
		float yaw = 0.f * (2.f*3.14159f/360.f);

		float cp = cos(pitch), sp = sin(pitch);
		float new_j = cp*j - sp*k, new_k = sp*j + cp*k;
		j = new_j; k = new_k;

		float cy = cos(yaw), sy = sin(yaw);
		float new_i = cy*i - sy*k, new_k2 = sy*i + cy*k;
		i = new_i; k = new_k2;

		float mag = sqrt(i*i + j*j + k*k);
		//make_ray(107.2, 77.8f, 87.1, i/mag, j/mag, k/mag, &r);
		//make_ray(107.2, 75.8f, 102.1, i/mag, j/mag, k/mag, &r);
		//make_ray(-10, 400, 102, i/mag, j/mag, k/mag, &r);
		//make_ray(107.2, 81.8, 102.1, i/mag, j/mag, k/mag, &r);
		//make_ray(200.f, 107.f, -100.f, i/mag, j/mag, k/mag, &r);

		//make_ray(-20, 73, 122, i/mag, j/mag, k/mag, &r);
		//make_ray(-100, 130, -200, i/mag, j/mag, k/mag, &r);
		make_ray(0, 128, 0, i/mag, j/mag, k/mag, &r);
		//make_ray(-100, 500, -100, i/mag, j/mag, k/mag, &r);

		// have 2 bounces of indirect illumination
		sample_ray(record, r, 1);
	}

	enum block_side { TOP, BOTTOM, SIDE };

	data::pixel<float> sample_terrain(int bx, int by, float tx, float ty) const
	{
		assert((bx >= 0) && (bx < 16));
		assert((by >= 0) && (by < 16));

		int px = (bx << 4) + static_cast<int>(floor(tx * 16.f));
		int py = (by << 4) + 15 - static_cast<int>(floor(ty * 16.f));

		px = std::max(0, std::min(terrain_colour.width-1, px));
		py = std::max(0, std::min(terrain_colour.height-1, py));

		return terrain_colour.at(px, py);
	}

	data::pixel<float> block_surface_colour(
			const data::block& block, const octree::sub_location hit_loc,
			float nx, float ny, float nz) const
	{
		data::pixel<float> output;

		// return data::pixel<float>(1,1,1);

		// extract the fractional part of the intersection co-ord (this will be the texture)
		float tmp, tx, ty;
		float fx(fabs(modff(hit_loc.coords[0], &tmp)));
		float fy(fabs(modff(hit_loc.coords[1], &tmp)));
		float fz(fabs(modff(hit_loc.coords[2], &tmp)));

		if(hit_loc.coords[0] < 0.f)
			fx = 1.f - fx;
		if(hit_loc.coords[1] < 0.f)
			fy = 1.f - fy;
		if(hit_loc.coords[2] < 0.f)
			fz = 1.f - fz;

		// work out which side of the block we are and the texture co-ord.
		block_side side = SIDE;
		if(ny > 0.95f)
		{
			side = TOP;
			tx = fx; ty = fz;
		}
		else if(ny < -0.95f)
		{
			side = BOTTOM;
			tx = fx; ty = fz;
		}
		else if(fabs(nx) > 0.95f)
		{
			side = SIDE;
			tx = fz; ty = fy;
		}
		else if(fabs(nz) > 0.95f)
		{
			side = SIDE;
			tx = fx; ty = fy;
		}
		else
		{
			// ????
			tx = ty = 0.f;
		}

		switch(block.id)
		{
			case mc::Stone:
				output = sample_terrain(1, 0, tx, ty);
				break;
			case mc::Dirt:
				output = sample_terrain(2, 0, tx, ty);
				break;
			case mc::Grass:
				output = sample_terrain(3, 0, tx, ty);
				if(side == TOP)
					output = sample_terrain(0, 0, tx, ty) * data::pixel<float>(0.f, 1.f, 0.f);
				break;
			case mc::Wood:
				output = sample_terrain(4, 0, tx, ty);
				break;
			case mc::Cobblestone:
				output = sample_terrain(0, 1, tx, ty);
				break;
			case mc::Log:
				output = sample_terrain(4, 1, tx, ty);
				break;
			case mc::Sand:
				output = sample_terrain(2, 1, tx, ty);
				break;
			case mc::CoalOre:
				output = sample_terrain(2, 2, tx, ty);
				break;
			case mc::Water:
			case mc::StationaryWater:
				output = data::pixel<float>(0.f, 0.f, 0.5f);
				break;
			case mc::Leaves:
				output = sample_terrain(5, 3, tx, ty) * data::pixel<float>(0.f, 1.f, 0.f);
				break;
			default:
				//std::cout << "unknown block: 0x" << std::hex << (int) block_id << std::endl;
				output = data::pixel<float>(0.5f, 0.5f, 0.5f);
				break;
		}

		return output;
	}

	void sample_ray(sample_recorder& record, ray& r, int recurse_depth, const pixel_f32& sample_filter_ = pixel_f32(1,1,1)) const
	{
		octree::sub_location node_sub_loc;
		data::block hit_block;
		pixel_f32 ray_transmission;

		if(!cast_ray(r, node_sub_loc, hit_block, ray_transmission))
		{
			// we didn't hit the world, sample the sky
			record(sample_filter_ * ray_transmission * sample_sky(r));
			return;
		}

		// update the sample filter with the tranmission of this ray
		pixel_f32 sample_filter = sample_filter_ * ray_transmission;

		const octree::extent& node_ext = node_sub_loc.node_extent;

		float mid_x = static_cast<float>(node_ext.loc.x) + 0.5f * static_cast<float>(node_ext.size);
		float mid_y = static_cast<float>(node_ext.loc.y) + 0.5f * static_cast<float>(node_ext.size);
		float mid_z = static_cast<float>(node_ext.loc.z) + 0.5f * static_cast<float>(node_ext.size);

		float hit_x = node_sub_loc.coords[0];
		float hit_y = node_sub_loc.coords[1];
		float hit_z = node_sub_loc.coords[2];

		float to_obs_x = - hit_x + r.x;
		float to_obs_y = - hit_y + r.y;
		float to_obs_z = - hit_z + r.z;
		float mag_to_obs = sqrt(to_obs_x*to_obs_x + to_obs_y*to_obs_y + to_obs_z*to_obs_z);
		to_obs_x /= mag_to_obs;
		to_obs_y /= mag_to_obs;
		to_obs_z /= mag_to_obs;

		float normal_x = hit_x - mid_x;
		float normal_y = hit_y - mid_y;
		float normal_z = hit_z - mid_z;

		float mag_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
		normal_x /= mag_normal;
		normal_y /= mag_normal;
		normal_z /= mag_normal;

		// convert the spherical normal into a cubical one...
		float abs_x = fabs(normal_x), abs_y = fabs(normal_y), abs_z = fabs(normal_z);
		float almost_one = 1.f;

		if((almost_one*abs_x > abs_y) && (almost_one*abs_x > abs_z))
		{
			normal_x = normal_x > 0.f ? 1.f : -1.f;
			normal_y = normal_z = 0.f;
		}
		else if((almost_one*abs_y > abs_x) && (almost_one*abs_y > abs_z))
		{
			normal_y = normal_y > 0.f ? 1.f : -1.f;
			normal_x = normal_z = 0.f;
		}
		else if((almost_one*abs_z > abs_x) && (almost_one*abs_z > abs_y))
		{
			normal_z = normal_z > 0.f ? 1.f : -1.f;
			normal_x = normal_y = 0.f;
		}

		// handle luminous blocks
		const float glow_scale = 1.f * std::max(0.f,
				to_obs_x*normal_x + to_obs_y*normal_y + to_obs_z*normal_z);
		switch(hit_block.id)
		{
			case mc::Torch:
				record(sample_filter * data::pixel<float>(1,1,1) * glow_scale);
				return; break;

			case mc::Lava:
			case mc::StationaryLava:
				record(sample_filter * data::pixel<float>(0.25,0.1,0.01) * glow_scale);
				return; break;
		}

		data::pixel<float> surface_colour = block_surface_colour(hit_block, node_sub_loc,
				normal_x, normal_y, normal_z);

		// for the sky (or, rather, it's luminosity) we want to evaluate the integral
		// 
		// I = \int sky(t,p) * max(0, cos((t,p) . normal)) dt dp [t == theta, p == phi]
		// 
		// let g1(t,p) be the spherical PDF sky(t,p)/[\int sky(t,p) dt dp] == sky(t,p)/sky_norm. then
		//
		// I = sky_norm \int g1(t,p) * max(0,cos((t,p) . normal)) dt dp
		// 
		// let g2(t,p) = max(0, cos((t,p) . normal)) / \int max(0, cos((t,p) /. normal)) dt dp = max(0, cos((t,p) . normal)) / pi
		//
		// then
		//
		// I = pi * sky_norm \int g1(t,p) * g2(t,p) dt dp where g1() and g2() are pdfs
		//
		// We can therefore evaluate the integral via
		// importance sampling: viewing it as the expectation
		// of g1 w.r.t. g2 or the expectation of g2 w.r.t. g1.
		//
		// We can draw samples from g2 analytically and can draw samples from g1 via rejection sampling.

		// FIXME: is there some normalisation for the lambertian lighting BRDF?
		const float lambertian_norm = 1.f; // 3.14159f;

		// calculate expectation of g2 w.r.t. g1

		// dont waste the opportunity just because we end up sampling a sky ray behind us
		float g2 = 0.f;
		for(int tries = 0; (g2 <= 0.f) && (tries < 4); ++tries)
		{
			pixel_f32 sky_sample(0,0,0);

			// sample g1
			ray g1_ray;
			data::pixel<float> sky_colour = sample_sky_pdf(hit_x, hit_y, hit_z, g1_ray);

			// calculate g2()
			float g2 = g1_ray.i * normal_x + g1_ray.j * normal_y + g1_ray.k * normal_z;
			if(g2 > 0.f)
			{
				// can we see the sky?
				pixel_f32 to_sky_transmission;
				if(!cast_ray(g1_ray, to_sky_transmission))
				{
					sky_sample = sky_colour * g2 * sky_integral * to_sky_transmission;
				}
			}

			// record the sky sample (including black samples)
			record(sample_filter * sky_sample * surface_colour * lambertian_norm);
		}
		
		// calculate expectation of sky w.r.t. g2 (avoids having to divide and multiply by sky_norm)
	
		// sample g2
		Vector g2_ray_dir = pt_sampling_cosine(
				pt_vector_make(normal_x, normal_y, normal_z, 0.f));

		ray g2_ray;
		make_ray(hit_x, hit_y, hit_z,
				pt_vector_get_x(g2_ray_dir),
				pt_vector_get_y(g2_ray_dir),
				pt_vector_get_z(g2_ray_dir),
				&g2_ray);

		// we need to decide whether to recursively sample the world as well
		if(recurse_depth > 0)
		{
			// FIXME: Think about this... is it biasing the result?
			pixel_f32 new_filter = sample_filter * surface_colour * lambertian_norm;
			sample_ray(record, g2_ray, recurse_depth - 1, new_filter);
		}
		else
		{
			// by default we see blackness
			data::pixel<float> sample(0,0,0);

			// can we see the sky?
			pixel_f32 to_sky_transmission;
			if(!cast_ray(g2_ray, to_sky_transmission))
			{
				// calculate g1() * sky_norm
				sample = sample_sky(g2_ray);
			}

			record(sample_filter * sample * surface_colour * lambertian_norm);
		}
	}

	int operator() (int argc, char** argv)
	{
		if(argc != 4)
		{
			std::cerr << "usage: " << argv[0] << " <path to world> <sky HDR image> <output PPM>" << std::endl;
			return 1;
		}

		io::read_png("terrain.png", terrain_colour, terrain_alpha);
		std::cout << "Read terrain at " << terrain_colour.width << " x " << terrain_colour.height << std::endl;

		if(!HDRLoader::load(argv[2], sky_light_probe))
		{
			std::cerr << "error loading sky light probe." << std::endl;
			return 1;
		}

		std::cout << "Loaded skybox size: " << sky_light_probe.width << "x" << sky_light_probe.height << std::endl;
		// find sky spherical integral
		double integral = 0.f;
		sky_max_lum = 0;
		for(int sy=0; sy<sky_light_probe.height; ++sy)
		{
			for(int sx=0; sx<sky_light_probe.width; ++sx)
			{
				float *p_sky_pix = &(sky_light_probe.cols[3*(sx+(sy*sky_light_probe.width))]);
				data::pixel<float> sky_pix(p_sky_pix[0], p_sky_pix[1], p_sky_pix[2]);
				float sky_lum = rgb2y(sky_pix);

				sky_max_lum = std::max(sky_max_lum, sky_lum);
				integral += sky_pixel_solid_angle(sx, sy, sky_light_probe.width, sky_light_probe.height) * sky_lum;
			}
		}
		sky_integral = integral;
		std::cout << "sky light probe integral: " << integral << std::endl;
		std::cout << "sky maximum luminosity: " << sky_max_lum << std::endl;

#if 0
		mc::world world(argv[1]);

		mc::region_iterator region_iterator(world.get_iterator());
		std::deque<mc::utils::level_coord> level_coords;
		while(region_iterator.has_next())
		{
			boost::shared_ptr<mc::region> region(region_iterator.next());
			std::cout << "loaded region from: " << region->get_path() << std::endl;

			level_coords.clear();
			region->read_header();
			region->read_coords(level_coords);
			std::cout << " - " << level_coords.size() << " level coords." << std::endl;

			mc::utils::level_coord rc = mc::utils::path_to_region_coord(region->get_path());
			octree::location start_coord = octree::location(rc.get_x() * 16, 0, rc.get_z() * 16);

			octree::octree<data::block> loaded_tree(9, start_coord, mc::Air);
			std::for_each(level_coords.begin(), level_coords.end(), load_level(region, loaded_tree));
			std::cout << " - raw node count " << loaded_tree.node_count() << " nodes." << std::endl;
#if 0
			std::cout << "   ... sparsifying" << std::endl;
			sparsify(loaded_tree);
#endif
			std::cout << "   ... compacting" << std::endl;
			loaded_tree.compact();
			std::cout << " - new node count " << loaded_tree.node_count() << " nodes." << std::endl;

			// crystalise into a read-only form
			crystal_octrees.push_back(octree::crystalised_octree(loaded_tree));
		}

		{
			namespace bio = boost::iostreams;
			bio::filtering_ostream output;
			output.push(bio::zlib_compressor());
			output.push(bio::file_sink("foo.dat"));

			uint32_t n_trees = crystal_octrees.size();
			io::nbo::write(output, n_trees);
			BOOST_FOREACH(const octree::crystalised_octree& tree, crystal_octrees)
			{
				tree.serialise(output);
			}
		}
#endif

		crystal_octrees.clear();
		{
			namespace bio = boost::iostreams;
			bio::filtering_istream input;
			input.push(bio::zlib_decompressor());
			input.push(bio::file_source("foo.dat"));

			uint32_t n_trees;
			io::nbo::read(input, n_trees);
			std::cout << "tree count: " << n_trees << std::endl;
			for(size_t i=0; i<n_trees; ++i)
			{
				crystal_octrees.push_back(octree::crystalised_octree(0, octree::location(0,0,0)));
				crystal_octrees.back().deserialise(input);
				std::cout << " - tree " << i+1 << " loaded." << std::endl;
			}
		}

		const int w=850, h=480;

		std::vector<pixel_u8>        pixels(w*h);
		std::vector<sample_recorder> recorders(w*h);

		const int32_t n_samples = 2048;
		for(int32_t sample_idx=0; sample_idx<n_samples; ++sample_idx)
		{
			std::cout << "pass " << sample_idx+1 << "/" << n_samples << std::endl;

			for(int32_t idx=0; idx<w*h; ++idx)
			{
				int x = idx % w;
				int y = h - 1 - idx / w;

				float fx(x), fy(y);
				sample_pixel(recorders[idx], fx+uniform_real()-0.5f, fy+uniform_real()-0.5f, w, h);
			}

			// rationale: a uniformly bright sky (with integral 4*pi) will result in a lambertian surface
			// having brightness pi => want to scale brightness by 1/pi == 4/integral
			float lum_scale = 4.f / sky_integral;
			for(int32_t idx=0; idx<w*h; ++idx)
			{
				const pixel_f32& mean = recorders[idx].sample_mean;
				pixel_u8& out = pixels[idx];

				pixel_f32 fout = mean * lum_scale;

				fout.r = 255.f * sqrt(fout.r);
				fout.g = 255.f * sqrt(fout.g);
				fout.b = 255.f * sqrt(fout.b);

				out.r = std::max(0, std::min(0xff, static_cast<int>(fout.r)));
				out.g = std::max(0, std::min(0xff, static_cast<int>(fout.g)));
				out.b = std::max(0, std::min(0xff, static_cast<int>(fout.b)));
			}

			std::ofstream output_fstream(argv[3]);
			io::write_ppm(output_fstream, &(pixels[0]), w, h);
		}

		return 0;
	}
};

int main(int argc, char** argv)
{
	g_thread_init(NULL);
	main_program prog;
	return prog(argc, argv);
}
