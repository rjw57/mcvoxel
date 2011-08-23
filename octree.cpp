#include "octree.hpp"

namespace octree
{

location::location(long x, long y, long z)
	: x(x), y(y), z(z)
{ }

location::location(const location& loc)
	: x(loc.x), y(loc.y), z(loc.z)
{ }

const location& location::operator = (const location& loc)
{
	x = loc.x; y = loc.y; z = loc.z;
	return *this;
}

location location::operator + (const location& rhs) const
{
	return location(x+rhs.x, y+rhs.y, z+rhs.z);
}

location location::operator - (const location& rhs) const
{
	return location(x-rhs.x, y-rhs.y, z-rhs.z);
}

}
