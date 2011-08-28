#ifndef DATAMODEL_HPP
#define DATAMODEL_HPP

#include <stdint.h>
#include <vector>

#include <mc/blocks.hpp>

namespace data
{

struct block
{
	uint8_t id, data;

	block() : id(mc::Air), data(0) { }
	~block() { }

	// copy and assignment constructors
	block(uint8_t id, uint8_t data = 0) : id(id), data(data) { }
	const block& operator = (const block& b) { id = b.id; data = b.data; return *this; }

	// packing and unpacking into an int32. FIXME in C++11, make this conversion operator explicit.
	operator int32_t () const { return static_cast<int32_t>(id) | (static_cast<int32_t>(data)<<8); }
	explicit block(int32_t i) : id(i & 0xff), data((i >> 8) & 0xff) { }

	// comparison
	bool operator == (const block& b) const { return (b.id == id) && (b.data == data); }
	bool operator != (const block& b) const { return ! operator == (b); }

	bool is_transparent() const { return id == mc::Air; }
};

template<typename T>
struct pixel
{
	T r,g,b;

	pixel() : r(), g(), b() { }
	pixel(const T& r, const T& g, const T& b) : r(r), g(g), b(b) { }

	pixel<T> operator + (const pixel<T>& p) const;
	pixel<T> operator - (const pixel<T>& p) const;
	pixel<T> operator / (const pixel<T>& p) const;
	pixel<T> operator * (const pixel<T>& p) const;

	pixel<T> operator + (const T& v) const;
	pixel<T> operator - (const T& v) const;
	pixel<T> operator / (const T& v) const;
	pixel<T> operator * (const T& v) const;
};

template<typename T>
T rgb2y(const data::pixel<T>& p);

template<typename T>
struct image
{
	int32_t        width, height;
	std::vector<T> pixels;

	image();
	image(int w, int h);
	image(const image& im);
	const image<T>& operator = (const image<T>& im);

	void resize(int w, int h);
	void reset();

	const T& at(int32_t x, int32_t y) const;
	T& at(int32_t x, int32_t y);
};

template<typename T>
struct sample_recorder
{
	T sample_mean;
	T sample_sq_mean;
	long n_samples;

	sample_recorder(const T& init_mean = T(), const T& init_sq_mean = T());
	~sample_recorder();

	// copy and asignment
	sample_recorder(const sample_recorder<T>& sr);
	const sample_recorder<T>& operator = (const sample_recorder<T>& sr);

	// record sample
	void record(const T& sample);
	void operator() (const T& sample) { record(sample); }
};

}

#define INSIDE_DATAMODEL_HPP
#include "datamodel.tcc"
#undef INSIDE_DATAMODEL_HPP

#endif // DATAMODEL_HPP
