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

	data::image<data::pixel<float> >           terrain_colour;
	data::image<float>                         terrain_alpha;

	bool cast_ray(const ray& r) const
	{
		octree::sub_location temp_sub_loc;
		data::block temp_hit_block;
		return cast_ray(r, temp_sub_loc, temp_hit_block);
	}

	bool cast_ray(const ray& r, octree::sub_location& out_sub_loc, data::block& out_block) const
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

	data::pixel<float> sample_sky(const Vector& v) const
	{
		ray temp_ray;
		make_ray(0, 0, 0, pt_vector_get_x(v), pt_vector_get_y(v), pt_vector_get_z(v), &temp_ray);
		return sample_sky(temp_ray);
	}

	data::pixel<float> sample_sky(const ray& r) const
	{
		float m = sqrt(r.i*r.i+r.j*r.j+r.k*r.k);
		float theta = acos(r.j/m);
		float phi = atan2(r.k, r.i);

		float sky_x = std::max(0.f, std::min(1.f, 0.5f * (1.f + (phi / 3.14159f))));
		float sky_y = std::max(0.f, std::min(1.f, theta / 3.14159f));

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

	data::pixel<float> sample_pixel(float fx, float fy, int w, int h) const
	{
		ray r;

		fx -= w>>1; fy -= h>>1;

		//float i(fx), k(fy), j(-0.5f*h);
		float i(fx), j(fy), k(h);

		float pitch = 50.f * (2.f*3.14159f/360.f);
		float yaw = -20.f * (2.f*3.14159f/360.f);

		float cp = cos(pitch), sp = sin(pitch);
		float new_j = cp*j - sp*k, new_k = sp*j + cp*k;
		j = new_j; k = new_k;

		float cy = cos(yaw), sy = sin(yaw);
		float new_i = cy*i - sy*k, new_k2 = sy*i + cy*k;
		i = new_i; k = new_k2;

		float mag = sqrt(i*i + j*j + k*k);
		//make_ray(107.2, 77.8f, 87.1, i/mag, j/mag, k/mag, &r);
		//make_ray(107.2, 75.8f, 102.1, i/mag, j/mag, k/mag, &r);
		//make_ray(0, 800, 0, i/mag, j/mag, k/mag, &r);
		//make_ray(107.2, 81.8, 102.1, i/mag, j/mag, k/mag, &r);
		make_ray(0.f, 260.f, 0.f, i/mag, j/mag, k/mag, &r);

		// have 2 bounces of indirect illumination
		return sample_ray(r, 2);
	}

	enum block_side { TOP, BOTTOM, SIDE };

	data::pixel<float> sample_terrain(int bx, int by, float tx, float ty) const
	{
		assert((bx >= 0) && (bx < 16));
		assert((by >= 0) && (by < 16));

		int px = (bx << 4) + static_cast<int>(floor(tx * 16.f));
		int py = (by << 4) + 15 - static_cast<int>(floor(ty * 16.f));

		return terrain_colour.at(px, py);
	}

	data::pixel<float> block_surface_colour(
			const data::block& block, const octree::sub_location hit_loc,
			float nx, float ny, float nz) const
	{
		data::pixel<float> output;

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
			case mc::Log:
				output = sample_terrain(4, 1, tx, ty);
				break;
			case mc::Sand:
				output = sample_terrain(2, 1, tx, ty);
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

	data::pixel<float> sample_ray(ray& r, int recurse_depth) const
	{
		data::pixel<float> output;

		octree::sub_location node_sub_loc;
		data::block hit_block;

		if(cast_ray(r, node_sub_loc, hit_block))
		{
			const octree::extent& node_ext = node_sub_loc.node_extent;

			float mid_x = static_cast<float>(node_ext.loc.x) + 0.5f * static_cast<float>(node_ext.size);
			float mid_y = static_cast<float>(node_ext.loc.y) + 0.5f * static_cast<float>(node_ext.size);
			float mid_z = static_cast<float>(node_ext.loc.z) + 0.5f * static_cast<float>(node_ext.size);

			float hit_x = node_sub_loc.coords[0];
			float hit_y = node_sub_loc.coords[1];
			float hit_z = node_sub_loc.coords[2];

			float to_obs_x = hit_x - r.x;
			float to_obs_y = hit_y - r.y;
			float to_obs_z = hit_z - r.z;
			float mag_to_obs = sqrt(to_obs_x*to_obs_x + to_obs_y*to_obs_y + to_obs_z*to_obs_z);
			to_obs_x /= mag_to_obs;
			to_obs_y /= mag_to_obs;
			to_obs_z /= mag_to_obs;

			float normal_x = hit_x - mid_x;
			float normal_y = hit_y - mid_y;
			float normal_z = hit_z - mid_z;

			// convert the spherical normal into a cubical one...
			float abs_x = fabs(normal_x), abs_y = fabs(normal_y), abs_z = fabs(normal_z);
			float almost_one = 0.9999f;

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
			else
			{
				float mag_normal = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z);
				normal_x /= mag_normal;
				normal_y /= mag_normal;
				normal_z /= mag_normal;
			}

			data::pixel<float> surface_colour = block_surface_colour(hit_block, node_sub_loc,
					normal_x, normal_y, normal_z);

			output.r = output.g = output.b = 0.f;

			// bounce samples
			const int n_iterations = 1;
			int samples_drawn = 0;
			for(int iteration_idx = 0; iteration_idx < n_iterations; ++iteration_idx)
			{
#if 1
				// firstly, we'll importance sample the sky. We choose a direction by rejection sampling
				// cosine weighted samples. We'll give up if we get to a maximum # of samples with no
				// result.
				Vector sky_ray_dir;
				bool sky_accepted = false;
				data::pixel<float> sky_sample;

				const int max_sky_samples = 256;
				for(int sky_sample_idx = 0; !sky_accepted && (sky_sample_idx < max_sky_samples); ++sky_sample_idx)
				{
					sky_ray_dir = pt_sampling_cosine(
						pt_vector_make(normal_x, normal_y, normal_z, 0.f));
					sky_sample = sample_sky(sky_ray_dir);

					float sky_lum = rgb2y(sky_sample);
					float prob_accept = sky_lum / sky_max_lum;

					if(uniform_real() < prob_accept)
					{
						// samples are already normalise w.r.t luminosity
						sky_sample = sky_sample / prob_accept;
						sky_accepted = true;
					}
				}

				if(sky_accepted)
				{
					// make a ray pointing from our intersection point to our sampled sky
					// direction
					ray sky_ray;
					make_ray(hit_x, hit_y, hit_z,
							pt_vector_get_x(sky_ray_dir),
							pt_vector_get_y(sky_ray_dir),
							pt_vector_get_z(sky_ray_dir),
							&sky_ray);

					// we see black if we intersect the world
					if(!cast_ray(sky_ray))
					{
						output = output + (surface_colour * sky_sample * 0.5f);
					}
				}
				++samples_drawn;
#endif

#if 1
				// we need to decide whether to recursively sample the world as well
				if(recurse_depth > 0)
				{
					// sample a ray in a weighted cosine centred on the normal
					Vector bounce_ray_dir = pt_sampling_cosine(
							pt_vector_make(normal_x, normal_y, normal_z, 0.f));

					ray bounce_ray;
					make_ray(hit_x, hit_y, hit_z,
							pt_vector_get_x(bounce_ray_dir),
							pt_vector_get_y(bounce_ray_dir),
							pt_vector_get_z(bounce_ray_dir),
							&bounce_ray);

					data::pixel<float> sample = sample_ray(bounce_ray, recurse_depth - 1);
					output = output + (surface_colour * sample * 0.5f);
					++samples_drawn;
				}
#endif
			}
			output = output / samples_drawn;
		}
		else
		{
			// we didn't hit the world, sample the sky
			output = sample_sky(r);
		}

		return output;
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
		// find sky max luminosity
		sky_max_lum = 0;
		for(int i=0; i<sky_light_probe.width*sky_light_probe.height; ++i)
		{
			float *sky_pix = sky_light_probe.cols + 3 * i;
			data::pixel<float> p;
			p.r = sky_pix[0];
			p.g = sky_pix[1];
			p.b = sky_pix[2];
			sky_max_lum = std::max(sky_max_lum, rgb2y(p));
		}
		std::cout << "Maximum luminosity: " << sky_max_lum << std::endl;

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
			octree::location start_coord = octree::location(rc.get_x() * 8, 0, rc.get_z() * 8);

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

		//const int w=850, h=480;
		const int w=1280, h=720;
		boost::shared_ptr< data::pixel<uint8_t> > pixels(new data::pixel<uint8_t>[w*h]);
		boost::shared_ptr< data::pixel<float> > float_pixels(new data::pixel<float>[w*h]);
		boost::shared_ptr< data::pixel<float> > float_pixels_sq(new data::pixel<float>[w*h]);
		boost::shared_ptr< data::pixel<float> > float_pixels_log(new data::pixel<float>[w*h]);
		boost::shared_ptr< float > luminance_pixels(new float[w*h]);
		boost::shared_ptr< float > luminance_sq_pixels(new float[w*h]);
		boost::shared_ptr< int > n_samples_pixels(new int[w*h]);

		for(int32_t idx=0; idx<w*h; ++idx)
		{
			data::pixel<float>* out = float_pixels.get() + idx;
			out->r = out->g = out->b = 0.f;

			data::pixel<float>* out_sq = float_pixels_sq.get() + idx;
			out_sq->r = out_sq->g = out_sq->b = 0.f;

			data::pixel<float>* out_log = float_pixels_log.get() + idx;
			out_log->r = out_log->g = out_log->b = 0.f;

			(luminance_pixels.get())[idx] = (luminance_sq_pixels.get())[idx] = 0.f;
			(n_samples_pixels.get())[idx] = 0;
		}

		const int32_t n_samples = 512;
		for(int32_t sample_idx=0; sample_idx<n_samples; ++sample_idx)
		{
			std::cout << "pass " << sample_idx+1 << "/" << n_samples << std::endl;

			float max_sigma = 0.f;
			for(int32_t idx=0; idx<w*h; ++idx)
			{
				int *p_n_samples = n_samples_pixels.get() + idx;
				if(*p_n_samples == 0)
					continue;

				float *p_luminance = luminance_pixels.get() + idx;
				float *p_luminance_sq = luminance_sq_pixels.get() + idx;

				float mu = *p_luminance / *p_n_samples;
				float local_variance = (*p_luminance_sq / *p_n_samples) - mu*mu;
				float local_sigma = 0.f;

				if(local_variance > 0.f)
					local_sigma = sqrt(local_variance);

				// divide by local mean luminance
				if(mu > 0.f)
					local_sigma /= mu;

				max_sigma = std::max(max_sigma, local_sigma);
			}

#			pragma omp parallel for schedule(dynamic, 1024)
			for(int32_t idx=0; idx<w*h; ++idx)
			{
				int x = idx % w;
				int y = h - 1 - idx / w;

				data::pixel<float>* out = float_pixels.get() + idx;
				data::pixel<float>* out_sq = float_pixels_sq.get() + idx;
				data::pixel<float>* out_log = float_pixels_log.get() + idx;

				int *p_n_samples = n_samples_pixels.get() + idx;
				float *p_luminance = luminance_pixels.get() + idx;
				float *p_luminance_sq = luminance_sq_pixels.get() + idx;

				float local_sigma = 100.f;
				if(*p_n_samples > 0)
				{
					float mu = *p_luminance / *p_n_samples;
					float local_variance = (*p_luminance_sq / *p_n_samples) - mu*mu;
					local_sigma = sqrt(local_variance);

					if(local_variance > 0.f)
						local_sigma = sqrt(local_variance);

					// divide by local mean luminance
					if(mu > 0.f)
						local_sigma /= mu;

					if(max_sigma > 0.f)
						local_sigma /= max_sigma;
				}

				if((sample_idx < 8) || (uniform_real() <= 0.25f + 0.75f * local_sigma))
				{
					float fx(x), fy(y);
					data::pixel<float> pixel_value =
						sample_pixel(fx+uniform_real()-0.5f, fy+uniform_real()-0.5f, w, h);

					*out = *out + pixel_value;
					*out_sq = *out_sq + pixel_value * pixel_value;

					if(out->r > 0.f)
						out_log->r += log(out->r);
					if(out->g > 0.f)
						out_log->g += log(out->g);
					if(out->b > 0.f)
						out_log->b += log(out->b);

					*p_luminance += rgb2y(pixel_value);
					*p_luminance_sq += rgb2y(pixel_value) * rgb2y(pixel_value);
					*p_n_samples += 1;
				}
			}

			float max_lum = 0.f;
			for(int32_t idx=0; idx<w*h; ++idx)
			{
				data::pixel<float> fout = *(float_pixels.get() + idx);
				int *p_n_samples = n_samples_pixels.get() + idx;

				if(*p_n_samples > 0)
				{
					// a better estimator based on assuming the distribution is gamma
					data::pixel<float> mean = fout / (*p_n_samples);
					max_lum = std::max(max_lum, rgb2y(mean));
				}
			}

			max_lum = sqrt(max_lum);

			for(int32_t idx=0, x=0, y=0; idx<w*h; ++idx)
			{
				data::pixel<float> fout = *(float_pixels.get() + idx);
				data::pixel<float> fout_sq = *(float_pixels_sq.get() + idx);
				data::pixel<float> fout_log = *(float_pixels_log.get() + idx);

				data::pixel<uint8_t>* out = pixels.get() + idx;
				int *p_n_samples = n_samples_pixels.get() + idx;

				if(*p_n_samples > 0)
				{
					// a better estimator based on assuming the distribution is gamma
					data::pixel<float> mean = fout / (*p_n_samples);
					data::pixel<float> mean_log = fout_log / (*p_n_samples);

					fout.r = gamma_est(mean.r, mean_log.r);
					fout.g = gamma_est(mean.g, mean_log.g);
					fout.b = gamma_est(mean.b, mean_log.b);

					/*
					data::pixel<float> mean_sq = fout_sq / (*p_n_samples);
					data::pixel<float> variance = mean_sq - mean*mean;
					if(rgb2y(variance) > 0.f)
					{
						data::pixel<float> theta = variance / mean;
						data::pixel<float> k = mean / theta;

						if((rgb2y(k) >= 1.f))
							fout = (k-1) * theta;
					}
					*/

					fout = mean;
				}

				fout.r = 255.f * sqrt(fout.r) / max_lum;
				fout.g = 255.f * sqrt(fout.g) / max_lum;
				fout.b = 255.f * sqrt(fout.b) / max_lum;

				out->r = std::max(0, std::min(0xff, static_cast<int>(fout.r)));
				out->g = std::max(0, std::min(0xff, static_cast<int>(fout.g)));
				out->b = std::max(0, std::min(0xff, static_cast<int>(fout.b)));

				//out->r = out->g = out->b = (250*(*p_n_samples))/(sample_idx+1);
#if 0
				float *p_luminance = luminance_pixels.get() + idx;
				float *p_luminance_sq = luminance_sq_pixels.get() + idx;

				float mu = *p_luminance / *p_n_samples
				float local_variance = (*p_luminance_sq / *p_n_samples) - mu*mu;
				float local_sigma = sqrt(local_variance);
				out->r = out->g = out->b = (250*local_sigma/max_sigma);
#endif

				++x;
				if(x >= w)
				{
					x = 0; --y;
				}
			}

			std::ofstream output_fstream(argv[3]);
			io::write_ppm(output_fstream, pixels.get(), w, h);
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
