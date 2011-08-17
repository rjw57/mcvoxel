#include <iostream>

#include <mc/world.hpp>

int main(int argc, char** argv)
{
	if(argc != 2)
	{
		std::cerr << "usage: " << argv[0] << " <path to world>" << std::endl;
		return 1;
	}

	mc::world world(argv[1]);

	mc::region_iterator region_iterator(world.get_iterator());
	while(region_iterator.has_next())
	{
		std::cout << "hello" << std::endl;
		region_iterator.next();
	}

	return 0;
}
