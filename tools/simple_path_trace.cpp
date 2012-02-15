#include <boost/foreach.hpp>
#include <mcvoxel/camera.hpp>
#include <cmath>
#include <deque>
#include <Eigen/Dense>
#include <iterator>
#include <iostream>
#include <mcvoxel/octree.hpp>
#include <mcvoxel/sky.hpp>
#include <mcvoxel/sampling.hpp>
#include <stdlib.h>
#include <vector>
#include <mcvoxel/world.hpp>

namespace mcvoxel
{

template<typename OutputIterator>
void trace(const ray& r, const world& w, OutputIterator bounces, int max_bounces = 5)
{
	ray current_ray(r);
	for(int n_bounces = 0; n_bounces < max_bounces; ++n_bounces, ++bounces)
	{
		surface_location intersection;
		if(!w.intersect(current_ray, intersection))
			return;
		*bounces = intersection;

		// choose new ray
		current_ray = ray(intersection.position,
				  cosine_weighted_hemisphere_direction(intersection.normal));
	}
}

}

int main(int argc, char** argv)
{
	if(argc != 4)
	{
		std::cerr << "syntax: " << argv[0] << " input.oct sky.hdr output.ppm" << std::endl;
		return EXIT_FAILURE;
	}

	// Load the world
	std::cout << "Loading word from: " << argv[1] << std::endl;
	mcvoxel::world world(argv[1]);

	// Load the sky
	std::cout << "Loading sky from: " << argv[2] << std::endl;
	mcvoxel::sky sky(argv[2]);

	// Size of output
	const int w(850), h(480);

	// Create a camera
	mcvoxel::camera camera;
	camera.set_centre(Eigen::Vector3f(0.f, 128.f, 0.f));
	camera.set_focal_length(h);
	camera.yaw_left(190.f * (2.f * M_PI / 360.f));
	camera.pitch_up(-20.f * (2.f * M_PI / 360.f));

	// What is the ray corresponding to (0,0)
	mcvoxel::ray origin_ray(camera.eye_ray(0.f,0.f));
	std::cout << "Camera origin: " << origin_ray.origin().transpose() << '\n';
	std::cout << "   Looking at: " << origin_ray.direction().transpose() << '\n';

	// Create a collection of pixel samples which is 3x(w*h)
	Eigen::ArrayXXf samples(3, w*h);

	for(int j=0; j<32; ++j)
	{
		std::cout << "j: " << j << std::endl;
		for(int i=0; i<w*h; ++i)
		{
			// choose some pixel
			float x(mcvoxel::uniform_real(-0.5f*w, 0.5f*w));
			float y(mcvoxel::uniform_real(-0.5f*h, 0.5f*h));
			int px = static_cast<int>(floor(x + 0.5f*w));
			int py = static_cast<int>(floor(-y + 0.5f*h));

			// in case of rounding error
			if((px < 0) || (px >= w) || (py < 0) || (py >= h))
				continue;

			// sample an eye ray
			mcvoxel::ray starting_ray(camera.eye_ray(x,y));
			Eigen::Vector3f sample(0.f,0.f,0.f);

			// trace a path
			std::deque<mcvoxel::surface_location> bounces;
			mcvoxel::trace(starting_ray, world, std::back_inserter(bounces));

			int n_sky_samples = 6;

			if(!bounces.empty())
			{
				const mcvoxel::surface_location& last_bounce(bounces.back());
				for(int i=0; i<n_sky_samples; ++i)
				{
					// choose a sky direction and chrominance
					Eigen::Vector3f sky_dir, sky_chrominance;
					sky.sample_direction(sky_dir, sky_chrominance);

					// don't bother if the sky direction faces away from us
					if(sky_dir.dot(last_bounce.normal) <= 0.f)
						continue;

					// form a ray from last bounce to sky
					mcvoxel::ray final_ray(last_bounce.position, sky_dir);

					// if we don't hit anything...
					mcvoxel::surface_location tmp_loc;
					if(!world.intersect(final_ray, tmp_loc))
					{
						// increment sample
						sample += sky_chrominance * last_bounce.normal.dot(sky_dir);
					}
				}
			}
			else
			{
				sample += n_sky_samples * sky.value_in_direction(starting_ray.direction());
			}

			int pidx = py * w + px;
			samples.matrix().col(pidx) += sample;
		}

		// Poor-man's tone-mapping
		Eigen::ArrayXXf tone_mapped_samples((samples / std::max(1e-3f, samples.maxCoeff())).cwiseSqrt());

		std::vector<data::pixel<uint8_t> > pixels(w*h, data::pixel<uint8_t>(0,0,0));
		for(int i=0; i<w*h; ++i)
		{
			pixels[i].r = static_cast<uint8_t>(255.f * tone_mapped_samples(0,i));
			pixels[i].g = static_cast<uint8_t>(255.f * tone_mapped_samples(1,i));
			pixels[i].b = static_cast<uint8_t>(255.f * tone_mapped_samples(2,i));
		}
		std::ofstream output(argv[3]);
		io::write_ppm(output, &(pixels[0]), w, h);
	}

	return EXIT_SUCCESS;
}
