#ifndef IO_HPP
#define IO_HPP

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <iostream>

#include "datamodel.hpp"
#include "octree.hpp"

namespace io
{

void write_ppm(std::ostream& os, data::pixel<uint8_t>* src, uint32_t w, uint32_t h);

void read_png(const char* filename, data::image<data::pixel<float> >& im, data::image<float>& alpha);

// generic writing and reading in network byte order
namespace nbo
{
	std::ostream& write(std::ostream& os, uint32_t val);
	std::ostream& write(std::ostream& os, uint32_t val);
	std::ostream& write(std::ostream& os, int32_t val);
	std::ostream& write(std::ostream& os, int16_t val);

	std::istream& read(std::istream& is, uint32_t& r_val);
	std::istream& read(std::istream& is, uint16_t& r_val);
	std::istream& read(std::istream& is, int32_t& r_val);
	std::istream& read(std::istream& is, int16_t& r_val);
}

/// @brief Load a set of octree::crystalised_octree instances in from a filename.
///
/// @tparam OutputIterator
/// @param filename
/// @param output
template<typename OutputIterator>
void load_crystal_octrees(
		const char* filename,
		OutputIterator output)
{
	namespace bio = boost::iostreams;
	bio::filtering_istream input;
	input.push(bio::zlib_decompressor());
	input.push(bio::file_source(filename));

	uint32_t n_trees;
	io::nbo::read(input, n_trees);
	std::cout << "tree count: " << n_trees << std::endl;
	for(size_t i=0; i<n_trees; ++i, ++output)
	{
		octree::crystalised_octree tree(0, octree::location(0,0,0));
		tree.deserialise(input);
		*output = tree;
		std::cout << " - tree " << i+1 << " loaded." << std::endl;
	}
}

}

namespace octree { struct location; }

//template< typename T >
//std::ostream& operator << (std::ostream& os, const octree::tree<T>& tree);

std::ostream& operator << (std::ostream& os, const octree::location& loc);

#define INSIDE_IO_HPP
#include "io.tcc"
#undef INSIDE_IO_HPP

#endif // IO_HPP
