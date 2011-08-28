// Voxel renderer main program.

#include <error.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <boost/program_options.hpp>

#include "sky.hpp"
#include "world.hpp"

namespace voxel
{

namespace po = boost::program_options;

struct main
{
	// our copy of the parameters passed to main.
	int argc; char** argv;

	// parsed command line parameters
	po::variables_map options;

	// description of the command line options
	po::options_description cmdline_opts;

	// the world we loaded
	world::world world;

	// our sky dome
	sky::sky sky;

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
			world::load_world(options["world"].as<std::string>().c_str(), world);

			// save a cache if we were asked to
			if(options.count("cached-world"))
			{
				std::ofstream output(options["cached-world"].as<std::string>().c_str());
				world::save_cached_world(output, world);
			}
		}
		else if(options.count("cached-world"))
		{
			std::ifstream input(options["cached-world"].as<std::string>().c_str());
			world::load_cached_world(input, world);
		}
		else
		{
			throw std::runtime_error("must specify at least one of world or cached world");
		}
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

		// try to load the sky
		try {
			if(options.count("sky"))
			{
				sky.load_hdr(options["sky"].as<std::string>().c_str());
			}
		} catch (const std::exception& e) {
			error(1, 0, "error: failed to load sky dome: %s", e.what());
		}

		// FIXME: some debug output
		std::cout << "sky light probe integral: " << sky.lum_integral() << std::endl;
		std::cout << "sky maximum luminosity: " << sky.max_lum() << std::endl;

		// try to load the world
		try {
			world_io();
		} catch (const std::exception& e) {
			error(1, 0, "error: failed to load world: %s", e.what());
		}

		return 0;
	}

	void parse_command_line()
	{
		// set up command line options
		po::options_description generic("Generic options");
		generic.add_options()
			("help,h", "print a brief usage summary")
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
