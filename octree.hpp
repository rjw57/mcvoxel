#ifndef _OCTREE_HPP
#define _OCTREE_HPP

#include <algorithm>
#include <cassert>
#include <deque>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/variant.hpp>

#include <rayslope/aabox.h>
#include <rayslope/ray.h>
#include <rayslope/slope.h>
#include <rayslope/slopeint_mul.h>

#include <mc/blocks.hpp>

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
	branch_node(const T& default_value);

	~branch_node();
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
	~octree();
	const octree<T>& operator = (const octree<T>& tree);

	const T& get(const location& loc) const;
	const T& get(long x, long y, long z) const { return get(location(x,y,z)); }

	void set(const location& loc, const T& val);
	void set(long x, long y, long z, const T& val ) { set(location(x,y,z), val); }

	const location& first_loc() const { return first_loc_; }
	long size() const { return 1l << log_2_size_; }
	long node_count() const;

	void compact();

	bool ray_intersect(const ray& r, float& out_distance, location& out_loc, long& out_size) const;

	protected:

	static void compact_node(branch_or_leaf_node_t& n);

	const branch_or_leaf_node_t* get_leaf(const location& loc, location& leaf_loc, long& leaf_size) const;
	branch_or_leaf_node_t* get_leaf(const location& loc, location& leaf_loc, long& leaf_size);

	int                   log_2_size_;
	location              first_loc_;
	branch_or_leaf_node_t root_node_;
};

}

#define INSIDE_OCTREE_HPP
#include "octree.tcc"
#undef INSIDE_OCTREE_HPP

#endif // _OCTREE_HPP
