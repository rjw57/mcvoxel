#ifndef DATAMODEL_HPP
#define DATAMODEL_HPP

#include <stdint.h>

namespace data
{

template<typename T>
struct pixel
{
	T r,g,b;

	pixel<T> operator + (const pixel<T>& rhs)
	{
		pixel<T> rv;
		rv.r = r + rhs.r;
		rv.g = g + rhs.g;
		rv.b = b + rhs.b;
		return rv;
	}

	pixel<T> operator - (const pixel<T>& p)
	{
		pixel<T> rv;
		rv.r = r - p.r;
		rv.g = g - p.g;
		rv.b = b - p.b;
		return rv;
	}

	pixel<T> operator / (const pixel<T>& p)
	{
		pixel<T> rv;
		rv.r = r / p.r;
		rv.g = g / p.g;
		rv.b = b / p.b;
		return rv;
	}

	pixel<T> operator * (const pixel<T>& p)
	{
		pixel<T> rv;
		rv.r = r * p.r;
		rv.g = g * p.g;
		rv.b = b * p.b;
		return rv;
	}

	pixel<T> operator + (const T& v)
	{
		pixel<T> rv;
		rv.r = r + v;
		rv.g = g + v;
		rv.b = b + v;
		return rv;
	}

	pixel<T> operator - (const T& v)
	{
		pixel<T> rv;
		rv.r = r - v;
		rv.g = g - v;
		rv.b = b - v;
		return rv;
	}

	pixel<T> operator / (const T& v)
	{
		pixel<T> rv;
		rv.r = r / v;
		rv.g = g / v;
		rv.b = b / v;
		return rv;
	}

	pixel<T> operator * (const T& v)
	{
		pixel<T> rv;
		rv.r = r * v;
		rv.g = g * v;
		rv.b = b * v;
		return rv;
	}
};

}

#endif // DATAMODEL_HPP
