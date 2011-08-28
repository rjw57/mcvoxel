// Voxel renderer main program.

#include <error.h>
#include <stdint.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <boost/program_options.hpp>

#include "scene.hpp"
#include "sky.hpp"
#include "world.hpp"

namespace voxel
{

namespace po = boost::program_options;

struct main
{
	typedef data::pixel<uint8_t> pixel_u8;
	typedef data::pixel<float> pixel_f32;
	typedef data::sample_recorder<pixel_f32> sample_recorder;

	// our copy of the parameters passed to main.
	int argc; char** argv;

	// parsed command line parameters
	po::variables_map options;

	// description of the command line options
	po::options_description cmdline_opts;

	// the scene we want to render
	scene::scene scene;

	main(int argc, char** argv)
		: argc(argc), argv(argv)
	{ }

	~main()
	{ }

	// load the specified world and save a cache if requested
	void world_io()
	{
		if(options.count("world"))
		{
			world::load_world(options["world"].as<std::string>().c_str(), scene.world);

			// save a cache if we were asked to
			if(options.count("cached-world"))
			{
				std::ofstream output(options["cached-world"].as<std::string>().c_str());
				world::save_cached_world(output, scene.world);
			}
		}
		else if(options.count("cached-world"))
		{
			std::ifstream input(options["cached-world"].as<std::string>().c_str());
			world::load_cached_world(input, scene.world);
		}
		else
		{
			throw std::runtime_error("must specify at least one of world or cached world");
		}
	}

	void write_output() const
	{
		if(!options.count("output"))
			throw std::runtime_error("no output stem specified.");

		std::string out_stem(options["output"].as<std::string>());

		const scene::scene::image& samples(scene.samples());

		int w(samples.width), h(samples.height);

		// rationale: a uniformly bright sky (with integral 4*pi) will result in a lambertian surface
		// having brightness pi => want to scale brightness by 1/pi == 4/integral
		float lum_scale = 4.f / scene.sky.lum_integral();
		data::image<pixel_u8> tone_mapped(samples.width, samples.height);
		for(int32_t idx=0; idx<w*h; ++idx)
		{
			const pixel_f32& mean = samples.pixels[idx].sample_mean;
			pixel_u8& out = tone_mapped.pixels[idx];

			pixel_f32 fout = mean * lum_scale;

			fout.r = 255.f * sqrt(fout.r);
			fout.g = 255.f * sqrt(fout.g);
			fout.b = 255.f * sqrt(fout.b);

			out.r = std::max(0, std::min(0xff, static_cast<int>(fout.r)));
			out.g = std::max(0, std::min(0xff, static_cast<int>(fout.g)));
			out.b = std::max(0, std::min(0xff, static_cast<int>(fout.b)));
		}

		std::string tone_mapped_fname(out_stem + ".ppm");
		{
			std::ofstream output_fstream(tone_mapped_fname.c_str());
			io::write_ppm(output_fstream, &(tone_mapped.pixels[0]), w, h);
			output_fstream.close();
		}
		std::cout << "wrote tone-mapped output to: " << tone_mapped_fname << "\n";
	}

	int operator () ()
	{
		try {
			parse_command_line();
		} catch (const po::unknown_option& e) {
			error(1, 0, "error: failed to parse command line: %s", e.what());
		}

		// print usage
		if(options.count("help"))
		{
			std::cout << cmdline_opts << "\n";
			return 0;
		}

		if(!options.count("output"))
		{
			error(1, 0, "must specify output stem via --output");
		}

		// try to load the sky
		try {
			if(options.count("sky"))
			{
				scene.sky.load_hdr(options["sky"].as<std::string>().c_str());
			}
		} catch (const std::exception& e) {
			error(1, 0, "error: failed to load sky dome: %s", e.what());
		}

		// FIXME: some debug output
		std::cout << "sky light probe integral: " << scene.sky.lum_integral() << std::endl;
		std::cout << "sky maximum luminosity: " << scene.sky.max_lum() << std::endl;

		// try to load the world
		try {
			world_io();
		} catch (const std::exception& e) {
			error(1, 0, "error: failed to load world: %s", e.what());
		}

		// set output size and initialise the scene
		scene.initialise(848, 480);
		scene.set_camera(0, 67, 0, 0, 0);

		try
		{
			const off_t n_passes = options["num-passes"].as<int>();
			const off_t n_samples_per_pass = scene.samples().width * scene.samples().height;

			for(off_t pass_idx = 0; pass_idx < n_passes; ++pass_idx)
			{
				std::cout << "pass " << pass_idx + 1 << "/" << n_passes << std::endl;

				for(off_t sample_idx=0; sample_idx < n_samples_per_pass; ++sample_idx)
				{
					scene.draw();
				}

				write_output();
			}
		}
		catch (const std::exception& e)
		{
			error(1, 0, "error: %s", e.what());
		}

		return 0;
	}

	void parse_command_line()
	{
		// set up command line options
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h", "print a brief usage summary")
			("output,o", po::value<std::string>(), "write output to files whose name begin with arg")
			("num-passes,n", po::value<int>()->default_value(64), "number of passes (default 64)")
		;
		cmdline_opts.add(generic);

		po::options_description world("Minecraft world options");
		world.add_options()
			("world,w", po::value<std::string>(), "load Minecraft world file")
			("cached-world,c", po::value<std::string>(), "save/load cached world")
		;
		cmdline_opts.add(world);

		po::options_description scene("Scene description options");
		scene.add_options()
			("sky,s", po::value<std::string>(), "load HDR sky probe image")
		;
		cmdline_opts.add(scene);

		// parse command line
		po::store(po::command_line_parser(argc, argv).
				options(cmdline_opts).run(), options);
		po::notify(options);
	}
};

}

int main(int argc, char** argv)
{
	// make an instance of the main program class, run it and return.
	return voxel::main(argc, argv)();
}
