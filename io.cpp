#include "io.hpp"
#include "octree.hpp"

namespace io
{

void write_ppm(std::ostream& os, data::pixel<uint8_t>* src, uint32_t w, uint32_t h)
{
    // write the header and record that we're writing in big endian order
    os << "P6" << std::endl << w << " " << h << std::endl << 255 << std::endl;
    os.write(reinterpret_cast<char*>(src), 3*w*h);
}

}

std::ostream& operator << (std::ostream& os, const octree::location& loc)
{
	os << "(" << loc.x << "," << loc.y << "," << loc.z << ")";
	return os;
}
