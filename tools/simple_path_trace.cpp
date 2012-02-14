#include <boost/foreach.hpp>
#include <deque>
#include <Eigen/Dense>
#include <iterator>
#include <io.hpp>
#include <octree.hpp>
#include <stdlib.h>
#include <vector>

typedef Eigen::Vector3f pixel_sample;
typedef Eigen::aligned_allocator<Eigen::Vector3f> pixel_sample_allocator;
typedef std::deque<pixel_sample, pixel_sample_allocator> pixel_sample_deque;

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		std::cerr << "syntax: " << argv[0] << " input.oct" << std::endl;
		return EXIT_FAILURE;
	}

	// Load the octrees
	std::deque<octree::crystalised_octree> trees;
	io::load_crystal_octrees(argv[1], std::back_inserter(trees));
	BOOST_FOREACH(const octree::crystalised_octree& tree, trees)
	{
		std::cout << "tree first loc: " << tree.extent_().loc << " + " << tree.extent_().size << '\n';
	}

	// Size of output
	const int w(1280), h(720);

	// Create a collection of pixel samples which is 3x(w*h)
	Eigen::ArrayXXf samples(3, w*h);

	float max_sample = samples.maxCoeff();

	std::vector<data::pixel<uint8_t> > pixels(w*h, data::pixel<uint8_t>(0,0,0));

	return EXIT_SUCCESS;
}
