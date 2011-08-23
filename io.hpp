#ifndef IO_HPP
#define IO_HPP

#include <iostream>

#include "datamodel.hpp"


namespace io
{

void write_ppm(std::ostream& os, data::pixel<uint8_t>* src, uint32_t w, uint32_t h);

}

namespace octree { struct location; }

//template< typename T >
//std::ostream& operator << (std::ostream& os, const octree::tree<T>& tree);

std::ostream& operator << (std::ostream& os, const octree::location& loc);

#define INSIDE_IO_HPP
#include "io.tcc"
#undef INSIDE_IO_HPP

#endif // IO_HPP
