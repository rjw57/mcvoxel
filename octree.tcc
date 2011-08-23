#ifndef INSIDE_OCTREE_HPP
#  error "This file must only be included by octree.hpp"
#endif

#include <algorithm>
#include <deque>
#include <utility>

#include <boost/foreach.hpp>

#include <rayslope/aabox.h>
#include <rayslope/ray.h>
#include <rayslope/slope.h>
#include <rayslope/slopeint_mul.h>

#include "io.hpp"

namespace octree
{

// BRANCH NODE MEMBERS

template<typename T>
branch_node<T>::branch_node(const T& default_value)
{
	for(size_t i=0; i<8; ++i)
	{
		children[i] = default_value;
	}
}

template<typename T>
branch_node<T>::~branch_node()
{ }

// INTERNAL UTILITY FUNCTIONS

inline bool node_contains(const location& loc, const location& branch_loc, long branch_size)
{
	if((loc.x < branch_loc.x) || (loc.x >= branch_loc.x + branch_size)) return false;
	if((loc.y < branch_loc.y) || (loc.y >= branch_loc.y + branch_size)) return false;
	if((loc.z < branch_loc.z) || (loc.z >= branch_loc.z + branch_size)) return false;
	return true;
}

inline size_t index_of_child_containing(const location& loc, const location& branch_loc, long branch_size)
{
	assert(branch_size > 1);
	assert(node_contains(loc, branch_loc, branch_size));

	long child_size = branch_size >> 1;
	int dx = (loc.x >= branch_loc.x + child_size) ? 1 : 0;
	int dy = (loc.y >= branch_loc.y + child_size) ? 1 : 0;
	int dz = (loc.z >= branch_loc.z + child_size) ? 1 : 0;

	return dx + (dy<<1) + (dz<<2);
}

inline location location_of_child(size_t child_idx, const location& branch_loc, long branch_size)
{
	assert(branch_size > 1);

	location child_loc(branch_loc);
	long child_size = branch_size >> 1;

	if(child_idx & 0x1) child_loc.x += child_size;
	if(child_idx & 0x2) child_loc.y += child_size;
	if(child_idx & 0x4) child_loc.z += child_size;

	return child_loc;
}

// OCTREE MEMBERS

template<typename T>
octree<T>::octree(int log_2_size, const location& first_loc, const T& default_value)
	: log_2_size_(log_2_size), first_loc_(first_loc), root_node_(default_value)
{ }

template<typename T>
octree<T>::octree(const octree<T>& tree)
	: log_2_size_(tree.log_2_size_), first_loc_(tree.first_loc_), root_node_(tree.root_node_)
{ }

template<typename T>
octree<T>::~octree()
{ }

template<typename T>
const octree<T>& octree<T>::operator = (const octree<T>& tree)
{
	log_2_size_ = tree.log_2_size_;
	first_loc_ = tree.first_loc_;
	root_node_ = tree.root_node_;
	return *this;
}

template<typename T>
const typename octree<T>::branch_or_leaf_node_t* octree<T>::get_leaf(
		const location& loc, location& leaf_loc, long& leaf_size) const
{
	assert(node_contains(loc, first_loc_, size()));

	const branch_or_leaf_node_t* leaf_node(&root_node_);
	leaf_loc = first_loc_;
	leaf_size = size();

	while(const branch_node_t* p_branch_node = boost::get<branch_node_t>(leaf_node))
	{
		size_t child_idx = index_of_child_containing(loc, leaf_loc, leaf_size);
		leaf_node = &(p_branch_node->children[child_idx]);
		leaf_loc = location_of_child(child_idx, leaf_loc, leaf_size);
		leaf_size >>= 1;
	}

	assert(boost::get<T>(leaf_node) != NULL);

	return leaf_node;
}

template<typename T>
typename octree<T>::branch_or_leaf_node_t* octree<T>::get_leaf(const location& loc, location& leaf_loc, long& leaf_size)
{
	// do const-casting tricks to re-use the const implementation above
	return const_cast<branch_or_leaf_node_t*>(const_cast<const octree<T>*>(this)->get_leaf(loc, leaf_loc, leaf_size));
}

template<typename T>
const T& octree<T>::get(const location& loc) const
{
	location leaf_loc(first_loc_);
	long leaf_size(size());
	const branch_or_leaf_node_t* leaf_node(get_leaf(loc, leaf_loc, leaf_size));
	assert(leaf_node != NULL);
	return *boost::get<T>(leaf_node);
}

template<typename T>
void octree<T>::set(const location& loc, const T& val)
{
	// find the leaf node for this location
	location leaf_loc(first_loc_);
	long leaf_size(size());
	branch_or_leaf_node_t* leaf_node(get_leaf(loc, leaf_loc, leaf_size));

	assert(leaf_node != NULL);

	// if already the correct value, skip
	if(boost::get<T>(*leaf_node) == val)
		return;

	// if already the lowest leaf size, set directly
	if(leaf_size == 1)
	{
		*leaf_node = val;
		return;
	}

	// otherwise we need to turn this leaf into a branch node and recurse down to the leaf
	T previous_value = boost::get<T>(*leaf_node);
	while(leaf_size > 1)
	{
		*leaf_node = branch_node_t(previous_value);

		// move one level down the tree
		branch_node_t* p_branch = boost::get<branch_node_t>(leaf_node);
		assert(p_branch != NULL);

		size_t child_idx = index_of_child_containing(loc, leaf_loc, leaf_size);
		leaf_node = &(p_branch->children[child_idx]);
		leaf_loc = location_of_child(child_idx, leaf_loc, leaf_size);
		leaf_size >>= 1;
	}

	assert(leaf_size == 1);
	assert(boost::get<T>(leaf_node) != NULL);
	*boost::get<T>(leaf_node) = val;
}

template<typename T>
long octree<T>::node_count() const
{
	std::deque<const branch_or_leaf_node_t*> nodes_to_examine;
	nodes_to_examine.push_back(&root_node_);

	// visit nodes in a depth-first way
	long count = 0;
	while(nodes_to_examine.size() > 0)
	{
		++count;

		const branch_or_leaf_node_t* next_node = nodes_to_examine.back();
		nodes_to_examine.pop_back();

		if(const branch_node_t* branch_p = boost::get<branch_node_t>(next_node))
		{
			for(int i=0; i<8; ++i)
			{
				nodes_to_examine.push_back(&branch_p->children[i]);
			}
		}
	}

	return count;
}

template<typename T>
void octree<T>::compact_node(branch_or_leaf_node_t& n)
{
	if(NULL != boost::get<T>(&n))
		return;

	branch_node_t* branch_p = boost::get<branch_node_t>(&n);
	assert(branch_p != NULL);

	for(int i=0; i<8; ++i)
	{
		compact_node(branch_p->children[i]);
	}

	bool all_leaf_children = true;
	for(int i=0; (i<8) && all_leaf_children; ++i)
	{
		if(NULL != boost::get<branch_node_t>(&branch_p->children[i]))
			all_leaf_children = false;
	}

	// carry on if this has non-leaf children
	if(!all_leaf_children)
		return;

	// see if all the leaves have the same value
	const T& match_val = boost::get<T>(branch_p->children[0]);
	bool all_children_identical = true;
	for(int i=1; (i<8) && all_children_identical; ++i)
	{
		if(match_val != boost::get<T>(branch_p->children[i]))
			all_children_identical = false;
	}

	// carry on if children are non-identical
	if(!all_children_identical)
		return;

	// we can merge this into a single leaf
	n = match_val;
}

template<typename T>
void octree<T>::compact()
{
	compact_node(root_node_);
}

struct intersection_result
{
	float    distance;
	location loc;
	long     size;

	intersection_result()
		: distance(0.f), loc(0,0,0), size(0)
	{ }

	intersection_result(float distance, const location& loc, long size)
		: distance(distance), loc(loc), size(size)
	{ }

	bool operator < (const intersection_result& rhs) const { return distance < rhs.distance; }
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

template<typename T>
struct child_intersection
{
	aabox    box;
	float    distance;
	location loc;
	long     size;
	const typename octree<T>::branch_or_leaf_node_t* child_p;

	child_intersection()
		: loc(0,0,0)
	{ }

	bool operator < (const child_intersection& rhs) const { return distance < rhs.distance; }
};

template<typename T>
bool octree<T>::ray_intersect(const ray& r, float& out_distance, location& out_loc, long& out_size) const
{
	typedef std::pair<extent, const branch_or_leaf_node_t*> node_record;

	// convert the ray to one in this tree's frame
	ray transformed_ray;
	make_ray(r.x-first_loc_.x, r.y-first_loc_.y, r.z-first_loc_.z, r.i, r.j, r.k, &transformed_ray);

	intersection_result result;

	std::deque<node_record> nodes;
	nodes.push_back(node_record(extent(first_loc_, size()), &root_node_));

	while(nodes.size() > 0)
	{
		const node_record& record = nodes.back();
		nodes.pop_back();

		const extent& node_ext = record.first;
		const branch_or_leaf_node_t* node_p = record.second;

		aabox node_box;
		::make_aabox(node_ext.loc.x, node_ext.loc.y, node_ext.loc.z,
				node_ext.loc.x+node_ext.size, node_ext.loc.y+node_ext.size, node_ext.loc.z+node_ext.size,
				&node_box);

		if(const leaf_node_t* leaf = boost::get<leaf_node_t>(node_p))
		{
			// do we intersect this node?
			if(!slopeint_mul(&transformed_ray, &node_box, &out_distance))
				continue;

			// is it transparent?
			if(*leaf == mc::Air)
				continue;

			out_loc = node_ext.loc;
			out_size = node_ext.size;
			return true;
		}
		else if(const branch_node_t* branch = boost::get<branch_node_t>(node_p))
		{
			// do we intersect this node?
			if(!slope(&transformed_ray, &node_box))
				continue;

			// we do intersect somewhere, create a deque of children whose bounding boxes we intersect
			std::deque< child_intersection<T> > intersections;

			for(size_t i=0; i<8; ++i)
			{
				child_intersection<T> c_int;

				c_int.child_p = &branch->children[i];
				c_int.loc = location_of_child(i, node_ext.loc, node_ext.size);
				c_int.size = node_ext.size >> 1;
				::make_aabox(c_int.loc.x, c_int.loc.y, c_int.loc.z,
						c_int.loc.x+c_int.size, c_int.loc.y+c_int.size, c_int.loc.z+c_int.size,
						&c_int.box);

				if(slopeint_mul(&transformed_ray, &c_int.box, &c_int.distance))
					intersections.push_back(c_int);
			}

			// sort children by intersection distance
			std::sort(intersections.begin(), intersections.end());

			// try to intersect recursively, push nearest on the back
			BOOST_FOREACH(const child_intersection<T>& c_int,
					std::make_pair(intersections.rbegin(), intersections.rend()))
			{
				nodes.push_back(node_record(extent(c_int.loc, c_int.size), c_int.child_p));
			}
		}
		else
		{
			assert(false);
		}
	}

	return false;
}

}
