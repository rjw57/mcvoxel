#ifndef DATAMODEL_HPP
#define DATAMODEL_HPP

#include <stdint.h>

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

	pixel<T> operator + (const pixel<T>& p);
	pixel<T> operator - (const pixel<T>& p);
	pixel<T> operator / (const pixel<T>& p);
	pixel<T> operator * (const pixel<T>& p);

	pixel<T> operator + (const T& v);
	pixel<T> operator - (const T& v);
	pixel<T> operator / (const T& v);
	pixel<T> operator * (const T& v);
};

template<typename T>
struct image
{
	int32_t width, height;
};

}

#define INSIDE_DATAMODEL_HPP
#include "datamodel.tcc"
#undef INSIDE_DATAMODEL_HPP

#endif // DATAMODEL_HPP
