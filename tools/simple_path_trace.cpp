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

	for(int j=0; j<64; ++j)
	{
		std::cout << "iteration: " << j << std::endl;
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
			mcvoxel::ray finishing_ray;
			Eigen::Vector3f sample(0.f,0.f,0.f);

			// trace a path
			std::deque<mcvoxel::surface_location> bounces;
			const size_t max_bounces = 5;
			Eigen::Vector3f normalisation;
			mcvoxel::trace_path(starting_ray, world, max_bounces,
					    normalisation, finishing_ray, std::back_inserter(bounces));

			// check final ray doesn't intersect if necessary
			mcvoxel::surface_location intersection;
			if((bounces.size() < max_bounces) || !world.intersect(finishing_ray, intersection))
				sample += sky.value_in_direction(finishing_ray.direction()).cwiseProduct(normalisation);

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
