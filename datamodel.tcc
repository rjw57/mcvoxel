#ifndef INSIDE_DATAMODEL_HPP
#  error "This file must only be included by datamodel.hpp"
#endif

namespace data
{

template<typename T>
vec3<T>::vec3()
	: x(), y(), z()
{ }

template<typename T>
vec3<T>::vec3(const T& x, const T& y, const T& z)
	: x(x), y(y), z(z)
{ }

template<typename T>
vec3<T>::vec3(const Vector& v)
	: x(pt_vector_get_x(v)), y(pt_vector_get_y(v)), z(pt_vector_get_z(v))
{ }

template<typename T>
vec3<T> vec3<T>::operator+=(const vec3<T>& p)
{
	x += p.x; y += p.y; z += p.z;
	return *this;
}

template<typename T>
vec3<T> vec3<T>::operator-=(const vec3<T>& p)
{
	x -= p.x; y -= p.y; z -= p.z;
	return *this;
}

template<typename T>
vec3<T> vec3<T>::operator*=(const T& v)
{
	x *= v; y *= v; z *= v;
	return *this;
}

template<typename T>
vec3<T> vec3<T>::operator/=(const T& v)
{
	x /= v; y /= v; z /= v;
	return *this;
}

template<typename T>
vec3<T>::operator Vector () const
{
	return pt_vector_make(x, y, z, 0.f);
}

// PIXEL

template<typename T>
pixel<T>::pixel()
	: r(), g(), b()
{ }

template<typename T>
pixel<T>::pixel(const T& r, const T& g, const T& b)
	: r(r), g(g), b(b)
{ }

template<typename T>
pixel<T>::pixel(const Vector& v)
	: r(pt_vector_get_x(v)), g(pt_vector_get_y(v)), b(pt_vector_get_z(v))
{ }

template<typename T>
pixel<T> pixel<T>::operator+=(const pixel<T>& p)
{
	r += p.r; g += p.g; b += p.b;
	return *this;
}

template<typename T>
pixel<T> pixel<T>::operator-=(const pixel<T>& p)
{
	r -= p.r; g -= p.g; b -= p.b;
	return *this;
}

template<typename T>
pixel<T> pixel<T>::operator*=(const T& v)
{
	r *= v; g *= v; b *= v;
	return *this;
}

template<typename T>
pixel<T> pixel<T>::operator/=(const T& v)
{
	r /= v; g /= v; b /= v;
	return *this;
}

template<typename T>
pixel<T>::operator Vector () const
{
	return pt_vector_make(r, g, b, 0.f);
}

template<typename T>
T rgb2y(const data::pixel<T>& p)
{
	return 0.299f*p.r + 0.114f*p.g + 0.587f*p.b;
}

// IMAGE

template<typename T>
image<T>::image()
{
	reset();
}

template<typename T>
image<T>::image(int w, int h)
{
	resize(w, h);
}

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
void image<T>::resize(int w, int h)
{
	assert((w >= 0) && (h >= 0));
	width = w; height = h;
	pixels.resize(w*h);
}

template<typename T>
void image<T>::reset()
{
	width = height = 0;
	pixels.clear();
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
	return *this;
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
