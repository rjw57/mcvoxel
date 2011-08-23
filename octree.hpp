#ifndef _OCTREE_HPP
#define _OCTREE_HPP

#include <iostream>

#include <boost/variant.hpp>

#include <rayslope/ray.h>

namespace octree
{

struct location
{
	long x, y, z;

	location(long x, long y, long z);
	location(const location& loc);
	const location& operator = (const location& loc);

	location operator + (const location& rhs) const;
	location operator - (const location& rhs) const;
};

template<typename T>
struct branch_node
{
	typedef branch_node<T> branch_node_t;
	typedef boost::variant< T, boost::recursive_wrapper< branch_node_t > > branch_or_leaf_node_t;
	branch_or_leaf_node_t children[8];

	// branch nodes must be given a default value.
	branch_node(const T& default_value = T());

	~branch_node();
};

struct extent
{
	location loc;
	long size;

	extent()
		: loc(0,0,0), size(0)
	{ }

	extent(const location& loc, const long& size)
		: loc(loc), size(size)
	{ }
};

struct sub_location
{
	float  coords[3];
	extent node_extent;
};

template<typename T>
class octree
{
	public:

	typedef branch_node<T>                                branch_node_t;
	typedef T                                             leaf_node_t;
	typedef typename branch_node_t::branch_or_leaf_node_t branch_or_leaf_node_t;

	octree(int log_2_size, const location& first_loc = location(0,0,0), const T& default_value = T());
	octree(const octree<T>& tree);
	octree(std::istream& os);
	~octree();
	const octree<T>& operator = (const octree<T>& tree);

	const T& get(const location& loc) const;
	const T& get(long x, long y, long z) const { return get(location(x,y,z)); }

	void set(const location& loc, const T& val);
	void set(long x, long y, long z, const T& val ) { set(location(x,y,z), val); }

	const location& first_loc() const { return first_loc_; }
	long size() const { return 1l << log_2_size_; }
	long log_2_size() const { return log_2_size_; }
	const branch_or_leaf_node_t& root() const { return root_node_; }

	long node_count() const;
	void compact();

	bool ray_intersect(const ray& r, sub_location& out_sub_loc) const;

	std::ostream& serialise(std::ostream& os) const;
	std::istream& deserialise(std::istream& os);

	protected:

	static void compact_node(branch_or_leaf_node_t& n);

	const branch_or_leaf_node_t* get_leaf(const location& loc, location& leaf_loc, long& leaf_size) const;
	branch_or_leaf_node_t* get_leaf(const location& loc, location& leaf_loc, long& leaf_size);

	int                   log_2_size_;
	location              first_loc_;
	branch_or_leaf_node_t root_node_;
};

}

template<typename T>
std::istream& operator >> (std::istream& is, typename octree::octree<T>& tree);

template<typename T>
std::ostream& operator << (std::istream& os, const typename octree::octree<T>& tree);

#define INSIDE_OCTREE_HPP
#include "octree.tcc"
#undef INSIDE_OCTREE_HPP

#endif // _OCTREE_HPP
