#include <boost/foreach.hpp>
#include <camera.hpp>
#include <cmath>
#include <deque>
#include <Eigen/Dense>
#include <iterator>
#include <iostream>
#include <octree.hpp>
#include <stdlib.h>
#include <vector>
#include <world.hpp>

typedef Eigen::Vector3f pixel_sample;
typedef Eigen::aligned_allocator<Eigen::Vector3f> pixel_sample_allocator;
typedef std::deque<pixel_sample, pixel_sample_allocator> pixel_sample_deque;

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		std::cerr << "syntax: " << argv[0] << " input.oct output.ppm" << std::endl;
		return EXIT_FAILURE;
	}

	// Load the world
	std::cout << "Loading word from: " << argv[1] << std::endl;
	mcvoxel::world world(argv[1]);

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

	for(int py=0; py<h; ++py)
	{
		float y(-static_cast<float>(py) + 0.5f * static_cast<float>(h));
		for(int px=0; px<w; ++px)
		{
			float x(static_cast<float>(px) - 0.5f * static_cast<float>(w));
			int idx = py*w + px;

			mcvoxel::ray eye_ray(camera.eye_ray(x,y));
			mcvoxel::surface_location intersection;
			if(!world.intersect(eye_ray, intersection))
				continue;

			float delta = std::max(0.f, (intersection.position - eye_ray.origin()).norm() - 128.f);

			samples(0, idx) = delta;
			samples(1, idx) = delta;
			samples(2, idx) = delta;
		}
	}

	std::cout << "Rendering done." << std::endl;

	// Poor-man's tone-mapping
	Eigen::ArrayXXf tone_mapped_samples((samples / std::max(1e-3f, samples.maxCoeff())).cwiseSqrt());

	std::vector<data::pixel<uint8_t> > pixels(w*h, data::pixel<uint8_t>(0,0,0));
	for(int i=0; i<w*h; ++i)
	{
		pixels[i].r = static_cast<uint8_t>(255.f * tone_mapped_samples(0,i));
		pixels[i].g = static_cast<uint8_t>(255.f * tone_mapped_samples(1,i));
		pixels[i].b = static_cast<uint8_t>(255.f * tone_mapped_samples(2,i));
	}
	std::ofstream output(argv[2]);
	io::write_ppm(output, &(pixels[0]), w, h);

	return EXIT_SUCCESS;
}
