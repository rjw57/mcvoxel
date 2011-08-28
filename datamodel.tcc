#ifndef INSIDE_DATAMODEL_HPP
#  error "This file must only be included by datamodel.hpp"
#endif

namespace data
{

template<typename T>
pixel<T> pixel<T>::operator + (const pixel<T>& p) const
{
	pixel<T> rv;
	rv.r = r + p.r;
	rv.g = g + p.g;
	rv.b = b + p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator - (const pixel<T>& p) const
{
	pixel<T> rv;
	rv.r = r - p.r;
	rv.g = g - p.g;
	rv.b = b - p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator / (const pixel<T>& p) const
{
	pixel<T> rv;
	rv.r = r / p.r;
	rv.g = g / p.g;
	rv.b = b / p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator * (const pixel<T>& p) const
{
	pixel<T> rv;
	rv.r = r * p.r;
	rv.g = g * p.g;
	rv.b = b * p.b;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator + (const T& v) const
{
	pixel<T> rv;
	rv.r = r + v;
	rv.g = g + v;
	rv.b = b + v;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator - (const T& v) const
{
	pixel<T> rv;
	rv.r = r - v;
	rv.g = g - v;
	rv.b = b - v;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator / (const T& v) const
{
	pixel<T> rv;
	rv.r = r / v;
	rv.g = g / v;
	rv.b = b / v;
	return rv;
}

template<typename T>
pixel<T> pixel<T>::operator * (const T& v) const
{
	pixel<T> rv;
	rv.r = r * v;
	rv.g = g * v;
	rv.b = b * v;
	return rv;
}

template<typename T>
T rgb2y(const data::pixel<T>& p)
{
	return 0.299f*p.r + 0.114f*p.g + 0.587f*p.b;
}

// IMAGE

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

// SAMPLE RECORDER

template<typename T>
sample_recorder<T>::sample_recorder(const T& init_mean, const T& init_sq_mean)
	: sample_mean(init_mean), sample_sq_mean(init_sq_mean), n_samples(0)
{ }

template<typename T>
sample_recorder<T>::~sample_recorder()
{ }

template<typename T>
sample_recorder<T>::sample_recorder(const sample_recorder<T>& sr)
	: sample_mean(sr.sample_mean), sample_sq_mean(sr.sample_sq_mean), n_samples(sr.n_samples)
{ }

template<typename T>
const sample_recorder<T>& sample_recorder<T>::operator = (const sample_recorder<T>& sr)
{
	sample_mean = sr.sample_mean;
	sample_sq_mean = sr.sample_sq_mean;
	n_samples = sr.n_samples;
}

template<typename T>
void sample_recorder<T>::record(const T& sample)
{
	float f_n_samples = static_cast<float>(n_samples);
	const T sample_sq = sample * sample;

	sample_mean = sample_mean * (f_n_samples / (1.f + f_n_samples));
	sample_mean = sample_mean + (sample / (1.f + f_n_samples));

	sample_sq_mean = sample_sq_mean * (f_n_samples / (1.f + f_n_samples));
	sample_sq_mean = sample_sq_mean + (sample_sq / (1.f + f_n_samples));

	++n_samples;
}

}
