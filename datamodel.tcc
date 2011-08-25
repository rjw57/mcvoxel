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

}
