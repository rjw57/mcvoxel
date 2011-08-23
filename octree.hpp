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

	location(long x, long y, long z)
		: x(x), y(y), z(z)
	{ }

	location(const location& loc)
		: x(loc.x), y(loc.y), z(loc.z)
	{ }

	const location& operator = (const location& loc)
	{
		x = loc.x; y = loc.y; z = loc.z;
		return *this;
	}
};

inline location operator + (const location& lhs, const location& rhs)
{
	return location(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
}

inline location operator - (const location& lhs, const location& rhs)
{
	return location(lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z);
}

template< typename T >
struct node
{
	location              min_loc;
	unsigned int          log_2_size;
	T                     value;
	std::vector<node<T>*> children; // empty => leaf node

	node(const location& min_loc, unsigned int log_2_size, const T& value)
		: min_loc(min_loc), log_2_size(log_2_size), value(value)
	{ }

	node(const node<T>& n)
		: min_loc(0,0,0), log_2_size(0)
	{
		assign_from(n);
	}

	~node()
	{
		join();
	}

	const node<T>& operator = (const node<T>& n)
	{
		assign_from(n);
		return *this;
	}

	unsigned long size() const
	{
		return 1ul << log_2_size;
	}

	bool contains(const location& point) const
	{
		long sz = 1l << log_2_size;
		if((point.x < min_loc.x) || (point.x >= min_loc.x + sz))
			return false;
		if((point.y < min_loc.y) || (point.y >= min_loc.y + sz))
			return false;
		if((point.z < min_loc.z) || (point.z >= min_loc.z + sz))
			return false;
		return true;
	}

	bool is_leaf() const
	{
		return children.size() == 0;
	}

	node<T>* child_containing(const location& point) const
	{
		size_t idx = index_of_child_containing(point);
		node<T>* child_ptr = children[idx];
		return child_ptr;
	}

	node<T>& ensure_child_containing(const location& point)
	{
		size_t idx = index_of_child_containing(point);
		node<T>* child_ptr = children[idx];

		if(child_ptr != NULL)
			return *child_ptr;

		children[idx] = new node<T>(child_min_loc(idx), log_2_size-1, value);

		return *(children[idx]);
	}

	void split()
	{
		if(!is_leaf())
			throw std::runtime_error("cannot split a non-leaf node");

		if(log_2_size == 0)
			throw std::runtime_error("cannot split a leaf node one unit in size");

		children.clear();
		children.reserve(8);

		for(size_t i=0; i<8; ++i)
		{
			children.push_back(NULL);
		}
	}

	void compact()
	{
		// leaf nodes need no compacting
		if(is_leaf())
			return;

		// compact all children
		for(size_t i=0; i<8; ++i)
		{
			if(children[i] != NULL)
				children[i]->compact();
		}

		int n_default_nodes(0), n_leaf_nodes(0), n_branch_nodes(0);
		for(size_t i=0; i<8; ++i)
		{
			if(children[i] == NULL)
			{
				++n_default_nodes;
			}
			else if(children[i]->is_leaf())
			{
				++n_leaf_nodes;
			}
			else
			{
				++n_branch_nodes;
			}
		}

		// if we're all default value nodes, merge
		if(n_default_nodes == 8)
		{
			join();

			return;
		}

		// if we're all leaf nodes, merge if we have the same values for each
		if(n_leaf_nodes == 8)
		{
			T match_val(children[0]->value);
			bool all_match = true;
			for(size_t i=1; (i<8) && all_match; ++i)
			{
				if(children[i]->value != match_val)
					all_match = false;
			}

			if(all_match)
			{
				value = match_val;
				join();
			}

			return;
		}

		// if we're all leaf nodes or default, merge if we have the
		// same values for each
		if((n_default_nodes + n_leaf_nodes) == 8)
		{
			T match_val(value);
			bool all_match = true;
			for(size_t i=0; (i<8) && all_match; ++i)
			{
				if(children[i] == NULL)
					continue;
				if(children[i]->value != match_val)
					all_match = false;
			}

			if(all_match)
				join();

			return;
		}
	}

	const T& child_value(size_t idx)
	{
		if(is_leaf())
			return value;
		if(children[idx] == NULL)
			return value;
		return children[idx]->value;
	}

	unsigned long nodes_count() const
	{
		if(is_leaf())
			return 1;

		unsigned long count = 1;

		for(size_t i=0; i<8; ++i)
		{
			if(children[i] != NULL)
				count += children[i]->nodes_count();
		}

		return count;
	}

	unsigned long nodes_count_at_size(unsigned int log_2_node_size) const
	{
		if(is_leaf())
			return (log_2_size == log_2_node_size) ? 1 : 0;

		unsigned long count = 0;

		for(size_t i=0; i<8; ++i)
		{
			if(children[i] != NULL)
				count += children[i]->nodes_count_at_size(log_2_node_size);
		}

		return count;
	}

	struct child_intersection
	{
		node<T>* node_ptr;
		aabox    box;
		float    t;

		bool operator < (const child_intersection& rhs) const { return t < rhs.t; }
	};

	node<T>* ray_intersect(ray* r, float* t)
	{
		aabox node_box;
		make_aabox(&node_box);

		// if we don't intersect this node at all, return NULL
		if(!slopeint_mul(r, &node_box, t))
			return NULL;

		// if this node is a leaf node, return it as the intersection
		if(is_leaf())
		{
			if(value == mc::Air)
				return NULL;
			return this;
		}

		// try each child in turn
		std::vector<child_intersection> intersections;
		intersections.reserve(8);

		for(size_t i=0; i<8; ++i)
		{
			child_intersection record;
			make_child_aabox(i, &(record.box));

			if(slopeint_mul(r, &(record.box), &(record.t)))
			{
				record.node_ptr = children[i];
				intersections.push_back(record);
			}
		}

		if(intersections.size() == 0)
			return NULL;

		// sort intersections by distance
		std::sort(intersections.begin(), intersections.end());

		// try intersections, nearest first
		BOOST_FOREACH(const child_intersection& intersection, std::make_pair(intersections.begin(), intersections.end()))
		{
			if(intersection.node_ptr == NULL)
			{
				if(value == mc::Air)
					continue;
				*t = intersection.t;
				return this;
			}

			node<T>* child_intersection = intersection.node_ptr->ray_intersect(r, t);
			if(child_intersection != NULL)
			{
				return child_intersection;
			}
		}

		// return no intersection
		return NULL;
	}

	protected:

	void assign_from(const node<T>& n)
	{
		min_loc = n.min_loc;
		log_2_size = n.log_2_size;
		value = n.value;
		children = n.children;

		BOOST_FOREACH(node<T>*& child_entry, std::make_pair(children.begin(), children.end()))
		{
			if(child_entry == NULL)
				continue;
			child_entry = new node<T>(*child_entry);
		}
	}

	void make_aabox(aabox* box) const
	{
		::make_aabox(min_loc.x, min_loc.y, min_loc.z,
				min_loc.x+size(), min_loc.y+size(), min_loc.z+size(),
				box);
	}

	void make_child_aabox(size_t index, aabox* box) const
	{
		assert(log_2_size > 0);

		location child_loc(child_min_loc(index));
		unsigned long child_size = 1ul << (log_2_size-1);

		::make_aabox(child_loc.x, child_loc.y, child_loc.z,
				child_loc.x+child_size, child_loc.y+child_size, child_loc.z+child_size,
				box);
	}

	void join()
	{
		if(is_leaf())
			return;

		BOOST_FOREACH(node<T>* node_ptr, children) {
			if(node_ptr != NULL)
				delete node_ptr;
		}
		children.clear();
	}

	location child_min_loc(size_t index) const
	{
		location node_loc(min_loc);
		long half_size = size() >> 1;

		if(!!(index & 0x1))
			node_loc.x += half_size;
		if(!!(index & 0x2))
			node_loc.y += half_size;
		if(!!(index & 0x4))
			node_loc.z += half_size;

		return node_loc;
	}

	size_t index_of_child_containing(const location& point) const
	{
		if(!contains(point))
			throw std::invalid_argument("point not contained within node.");

		if(is_leaf())
			throw std::runtime_error("leaf nodes have no children.");

		location sub_location(point - min_loc);
		long half_size = size() >> 1;
		assert(half_size > 0);

		int dx = sub_location.x < half_size ? 0 : 1;
		int dy = sub_location.y < half_size ? 0 : 1;
		int dz = sub_location.z < half_size ? 0 : 1;

		return dx + (dy << 1) + (dz << 2);
	}
};

template<typename T>
class tree
{
	public:

	tree(unsigned int log_2_size, const T& default_value = T())
		: root_(new node<T>(location(0,0,0), log_2_size, default_value))
	{ }

	tree(unsigned int log_2_size, const location& min_loc, const T& default_value = T())
		: root_(new node<T>(min_loc, log_2_size, default_value))
	{ }

	tree(const tree& t)
		: root_(new node<T>(*t.root_))
	{ }

	const tree& operator = (const tree& t)
	{
		root_ = std::auto_ptr< node<T> >(new node<T>(*t.root_));
		return *this;
	}

	const node<T>& root() const { return *root_; }

	unsigned long nodes() const { return root_->nodes_count(); }

	unsigned long nodes_at_size(unsigned int log_2_size) const
	{
		return root_->nodes_count_at_size(log_2_size);
	}

	long size() const { return root_->size(); }

	node<T>* ray_intersect(ray* r, float* t) const {
		ray trans_ray(*r);
		location min_loc(root_->min_loc);
		make_ray(r->x-min_loc.x, r->y-min_loc.y, r->z-min_loc.z, r->i, r->j, r->k, &trans_ray);
		return root_->ray_intersect(&trans_ray, t);
	}

	void compact() { root_->compact(); }

	const T& get(long x, long y, long z) const { return get(location(x,y,z)); }

	const T& get(const location& loc) const {
		// find the leaf node which exists for this location
		node<T> *next_node = root_.get(), *leaf_node = NULL;
		do
		{
			leaf_node = next_node;

			if(leaf_node->is_leaf())
			{
				next_node = NULL;
			}
			else
			{
				next_node = leaf_node->child_containing(loc);
			}
		}
		while(next_node != NULL);

		return leaf_node->value;
	}

	void set(long x, long y, long z, const T& val)
	{
		set(location(x,y,z), val);
	}

	void set(const location& loc, const T& val)
	{
		// find the leaf node which exists for this location
		node<T> *next_node = root_.get(), *leaf_node = NULL;
		do
		{
			leaf_node = next_node;

			if(leaf_node->is_leaf())
			{
				next_node = NULL;
			}
			else
			{
				next_node = leaf_node->child_containing(loc);
			}
		}
		while(next_node != NULL);

		// if the leaf node already has the right value, we need do nothing
		if(leaf_node->value == val)
			return;

		// we need to recurse down to the lowest level of the tree and set the node
		while((leaf_node->log_2_size != 0) && (leaf_node->value != val))
		{
			if(leaf_node->is_leaf())
				leaf_node->split();
			leaf_node = &(leaf_node->ensure_child_containing(loc));
		}

		assert((leaf_node->log_2_size == 0) || (leaf_node->value == val));
		leaf_node->value = val;
	}

	protected:

	std::auto_ptr< node<T> > root_;
};

}

#endif // _OCTREE_HPP
