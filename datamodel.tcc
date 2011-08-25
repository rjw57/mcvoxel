#ifndef INSIDE_DATAMODEL_HPP
#  error "This file must only be included by datamodel.hpp"
#endif

namespace data
{

template<typename T>
pixel<T> pixel<T>::operator + (const pixel<T>& p)
{
	pixel<T> rv;
	rv.r = r + p.r;
	rv.g = g + p.g;
	rv.b = b + p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator - (const pixel<T>& p)
{
	pixel<T> rv;
	rv.r = r - p.r;
	rv.g = g - p.g;
	rv.b = b - p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator / (const pixel<T>& p)
{
	pixel<T> rv;
	rv.r = r / p.r;
	rv.g = g / p.g;
	rv.b = b / p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator * (const pixel<T>& p)
{
	pixel<T> rv;
	rv.r = r * p.r;
	rv.g = g * p.g;
	rv.b = b * p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator + (const T& v)
{
	pixel<T> rv;
	rv.r = r + v;
	rv.g = g + v;
	rv.b = b + v;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator - (const T& v)
{
	pixel<T> rv;
	rv.r = r - v;
	rv.g = g - v;
	rv.b = b - v;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator / (const T& v)
{
	pixel<T> rv;
	rv.r = r / v;
	rv.g = g / v;
	rv.b = b / v;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator * (const T& v)
{
	pixel<T> rv;
	rv.r = r * v;
	rv.g = g * v;
	rv.b = b * v;
	return rv;
}

template<typename T>
image<T>::image()
	: width(0), height(0), pixels()
{ }

template<typename T>
image<T>::image(const image& im)
	: width(im.width), height(im.height), pixels(im.pixels)
{ }

template<typename T>
const image<T>& image<T>::operator = (const image<T>& im)
{
	width = im.width;
	height = im.height;
	pixels = im.pixels;
}

template<typename T>
const T& image<T>::at(int32_t x, int32_t y) const
{
	return pixels.at(x + y*width);
}

template<typename T>
T& image<T>::at(int32_t x, int32_t y)
{
	return pixels.at(x + y*width);
}

}
