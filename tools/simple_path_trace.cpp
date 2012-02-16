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
#include <mcvoxel/trace.hpp>
#include <mcvoxel/world.hpp>
#include <stdlib.h>
#include <vector>

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
	camera.set_centre(Eigen::Vector3f(64.f, 100.f, 512.f));
	camera.set_focal_length(h);
	camera.yaw_left(10.f * (2.f * M_PI / 360.f));
	camera.pitch_up(-10.f * (2.f * M_PI / 360.f));

	// What is the ray corresponding to (0,0)
	mcvoxel::ray origin_ray(camera.eye_ray(0.f,0.f));
	std::cout << "Camera origin: " << origin_ray.origin().transpose() << '\n';
	std::cout << "   Looking at: " << origin_ray.direction().transpose() << '\n';

	// Create a collection of pixel samples which is 4x(w*h). The bottom row giving the sample count for that pixel.
	Eigen::ArrayXXf samples(4, w*h);

	for(int j=0; j<512; ++j)
	{
		float spp = samples.matrix().row(3).sum() / static_cast<float>(w*h);
		std::cout << "mean samples per-pixel: " << static_cast<size_t>(spp) << std::endl;

#		pragma omp parallel for schedule(dynamic, w)
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

			// sample a path from the eye
			mcvoxel::ray starting_ray(camera.eye_ray(x,y));
			const size_t max_bounces = 4;
			mcvoxel::ray finishing_ray;
			std::deque<mcvoxel::surface_location> bounces;
			mcvoxel::trace_path(starting_ray, world, max_bounces,
					    finishing_ray, std::back_inserter(bounces));

			// This is the sample we build up
			Eigen::Vector3f sample(0.f,0.f,0.f);
			int n_samples = 0;

			// Did we escape the world before the maximum number of bounces?
			if(bounces.size() < max_bounces)
			{
				sample += sky.value_in_direction(finishing_ray.direction());
				++n_samples;
			}

			if(!bounces.empty())
			{
				// Draw a light direction and try to form an explicit light path
				Eigen::Vector3f sky_dir, sky_sample;
				sky.sample_direction(sky_dir, sky_sample);

				BOOST_FOREACH(const mcvoxel::surface_location& loc, bounces)
				{
					// skip if the light faces away from the bounce surface
					float cos_theta(loc.normal.dot(sky_dir));
					if(cos_theta < 0.f)
						continue;

					// do we hit something on the way?
					if(world.intersect(mcvoxel::ray(loc.position, sky_dir)))
						continue;

					// no! add the sample
					sample += sky_sample * (cos_theta / M_PI);
				}

				n_samples += bounces.size();
			}

			int pidx = py * w + px;
#			pragma omp critical
			{
				samples.matrix().col(pidx).head<3>() += sample;
				samples(3, pidx) += static_cast<float>(n_samples);
			}
		}

		// normalise samples
		Eigen::ArrayXXf normalised_samples(samples.topRows<3>());
		for(int i=0; i<3; ++i)
			normalised_samples.matrix().row(i).cwiseQuotient(samples.matrix().row(3));

		// Poor-man's tone-mapping
		Eigen::ArrayXXf tone_mapped_samples((normalised_samples /
						     std::max(1e-3f, normalised_samples.maxCoeff())).cwiseSqrt());

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
