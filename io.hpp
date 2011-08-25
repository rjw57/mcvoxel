#ifndef IO_HPP
#define IO_HPP

#include <iostream>

#include "datamodel.hpp"

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

}

namespace octree { struct location; }

//template< typename T >
//std::ostream& operator << (std::ostream& os, const octree::tree<T>& tree);

std::ostream& operator << (std::ostream& os, const octree::location& loc);

#define INSIDE_IO_HPP
#include "io.tcc"
#undef INSIDE_IO_HPP

#endif // IO_HPP
