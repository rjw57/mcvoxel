#include <mcvoxel/camera.hpp>
#include <mcvoxel/datamodel.hpp>
#include <Eigen/Dense>
#include <mcvoxel/io.hpp>
#include <iostream>
#include <mcvoxel/sky.hpp>
#include <stdlib.h>
#include <vector>

int main(int argc, char** argv)
{
	if(argc != 3)
	{
		std::cerr << "syntax: " << argv[0] << " input_sky.hdr output.ppm" << std::endl;
		return EXIT_FAILURE;
	}

	// Load the sky
	mcvoxel::sky sky(argv[1]);

	// Size of output
	const int w(850), h(480);

	// Create a camera
	mcvoxel::camera camera;
	camera.set_centre(Eigen::Vector3f(0.f, 128.f, 0.f));
	camera.set_focal_length(h);
	camera.yaw_left(190.f * (2.f * M_PI / 360.f));
	//camera.pitch_up(-20.f * (2.f * M_PI / 360.f));
	camera.pitch_up(20.f * (2.f * M_PI / 360.f));

	// What is the ray corresponding to (0,0)
	mcvoxel::ray origin_ray(camera.eye_ray(0.f,0.f));
	std::cout << "Camera origin: " << origin_ray.origin().transpose() << '\n';
	std::cout << "   Looking at: " << origin_ray.direction().transpose() << '\n';

	// Create a collection of pixel samples which is 3x(w*h)
	Eigen::ArrayXXf samples(3, w*h);

	// Get camera vectors
	Eigen::Vector3f centre, look_at, up, right;
	camera.get_frame(centre, look_at, up, right);

	// Draw samples from the sky.
	for(int j=0; j<15; ++j)
	{
		std::cout << "j: " << j << std::endl;
		for(int i=0; i<w*h; ++i)
		{
			Eigen::Vector3f direction, chrominance;
			sky.sample_direction(direction, chrominance);

			// Get the image-plane co-ordinates for that ray
			float x, y;
			if(!camera.pixel_coord(direction, x, y))
				continue;

			// convert to image pixel co-ordinate
			x = floor(x + 0.5f * w); y = floor(-y + 0.5f * h);

			if((x < 0.f) || (x >= w) || (y < 0.f) || (y >= h))
				continue;

			int idx = static_cast<int>(y) * w + static_cast<int>(x);
			samples.matrix().col(idx) += chrominance;
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
