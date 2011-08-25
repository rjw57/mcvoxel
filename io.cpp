#include <arpa/inet.h>

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

namespace nbo
{

std::ostream& write(std::ostream& os, uint32_t val)
{
	uint32_t nval = htonl(val);
	os.write(reinterpret_cast<char*>(&nval), sizeof(uint32_t));
	return os;
}

std::ostream& write(std::ostream& os, uint16_t val)
{
	uint16_t nval = htons(val);
	os.write(reinterpret_cast<char*>(&nval), sizeof(uint16_t));
	return os;
}

std::istream& read(std::istream& is, uint32_t& r_val)
{
	uint32_t nval;
	is.read(reinterpret_cast<char*>(&nval), sizeof(uint32_t));
	r_val = ntohl(nval);
	return is;
}

std::istream& read(std::istream& is, uint16_t& r_val)
{
	uint16_t nval;
	is.read(reinterpret_cast<char*>(&nval), sizeof(uint16_t));
	r_val = ntohs(nval);
	return is;
}

}

}

std::ostream& operator << (std::ostream& os, const octree::location& loc)
{
	os << "(" << loc.x << "," << loc.y << "," << loc.z << ")";
	return os;
}
