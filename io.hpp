#ifndef IO_HPP
#define IO_HPP

#include <iostream>

#include "datamodel.hpp"
#include "octree.hpp"

namespace io
{

void write_ppm(std::ostream& os, data::pixel* src, uint32_t w, uint32_t h);

}

template< typename T >
std::ostream& operator << (std::ostream& os, const octree::tree<T>& tree);

#define INSIDE_IO_HPP
#include "io.tcc"
#undef INSIDE_IO_HPP

#endif // IO_HPP
