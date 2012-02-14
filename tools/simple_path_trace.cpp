#include <boost/foreach.hpp>
#include <deque>
#include <iterator>
#include <io.hpp>
#include <octree.hpp>
#include <stdlib.h>

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		std::cerr << "syntax: " << argv[0] << " input.oct" << std::endl;
		return EXIT_FAILURE;
	}

	std::deque<octree::crystalised_octree> trees;
	io::load_crystal_octrees(argv[1], std::back_inserter(trees));

	BOOST_FOREACH(const octree::crystalised_octree& tree, trees)
	{
		std::cout << "tree first loc: " << tree.extent_().loc << " + " << tree.extent_().size << '\n';
	}


	return EXIT_SUCCESS;
}
