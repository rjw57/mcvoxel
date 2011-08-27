#ifndef INSIDE_OCTREE_HPP
#  error "This file must only be included by octree.hpp"
#endif

#include <algorithm>
#include <cassert>
#include <iostream>
#include <stack>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include <rayslope/aabox.h>
#include <rayslope/ray.h>
#include <rayslope/slope.h>
#include <rayslope/slopeint_mul.h>

#include "io.hpp"

namespace octree
{

// EXTENT METHODS

inline ::aabox extent::make_aabox() const
{
	::aabox box;
	::make_aabox(loc.x, loc.y, loc.z, loc.x+size, loc.y+size, loc.z+size, &box);
	return box;
}

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

inline bool node_contains(const location& loc, const extent& branch_ext)
{
	if((loc.x < branch_ext.loc.x) || (loc.x >= branch_ext.loc.x + branch_ext.size)) return false;
	if((loc.y < branch_ext.loc.y) || (loc.y >= branch_ext.loc.y + branch_ext.size)) return false;
	if((loc.z < branch_ext.loc.z) || (loc.z >= branch_ext.loc.z + branch_ext.size)) return false;
	return true;
}

inline size_t index_of_child_containing(const location& loc, const extent& branch_ext)
{
	assert(branch_ext.size > 1);
	assert(node_contains(loc, branch_ext));

	long child_size = branch_ext.size >> 1;
	int dx = (loc.x >= branch_ext.loc.x + child_size) ? 1 : 0;
	int dy = (loc.y >= branch_ext.loc.y + child_size) ? 1 : 0;
	int dz = (loc.z >= branch_ext.loc.z + child_size) ? 1 : 0;

	return dx + (dy<<1) + (dz<<2);
}

inline location location_of_child(size_t child_idx, const extent& branch_ext)
{
	assert(branch_ext.size > 1);

	location child_loc(branch_ext.loc);
	long child_size = branch_ext.size >> 1;

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
octree<T>::octree(std::istream& is)
{
	is >> log_2_size_ >> first_loc_.x >> first_loc_.y >> first_loc_.z >> root_node_;
}

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
	assert(node_contains(loc, extent(first_loc_, size())));

	const branch_or_leaf_node_t* leaf_node(&root_node_);
	leaf_loc = first_loc_;
	leaf_size = size();

	while(const branch_node_t* p_branch_node = boost::get<branch_node_t>(leaf_node))
	{
		size_t child_idx = index_of_child_containing(loc, extent(leaf_loc, leaf_size));
		leaf_node = &(p_branch_node->children[child_idx]);
		leaf_loc = location_of_child(child_idx, extent(leaf_loc, leaf_size));
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

		size_t child_idx = index_of_child_containing(loc, extent(leaf_loc, leaf_size));
		leaf_node = &(p_branch->children[child_idx]);
		leaf_loc = location_of_child(child_idx, extent(leaf_loc, leaf_size));
		leaf_size >>= 1;
	}

	assert(leaf_size == 1);
	assert(boost::get<T>(leaf_node) != NULL);
	*boost::get<T>(leaf_node) = val;
}

template<typename T>
long octree<T>::node_count() const
{
	std::stack<const branch_or_leaf_node_t*> nodes_to_examine;
	nodes_to_examine.push(&root_node_);

	// visit nodes in a depth-first way
	long count = 0;
	while(nodes_to_examine.size() > 0)
	{
		++count;

		const branch_or_leaf_node_t* next_node = nodes_to_examine.top();
		nodes_to_examine.pop();

		if(const branch_node_t* branch_p = boost::get<branch_node_t>(next_node))
		{
			for(int i=0; i<8; ++i)
			{
				nodes_to_examine.push(&branch_p->children[i]);
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

template<typename T>
std::ostream& octree<T>::serialise(std::ostream& os) const
{
	os << log_2_size_ << '\n';
	os << first_loc_.x << ' ' << first_loc_.y << ' ' << first_loc_.z << '\n';

	std::stack<const branch_or_leaf_node_t*> node_stack;
	node_stack.push(&root_node_);

	while(!node_stack.empty())
	{
		const branch_or_leaf_node_t* node_p = node_stack.top();
		node_stack.pop();

		if(const leaf_node_t* leaf_p = boost::get<leaf_node_t>(node_p))
		{
			os << 'L';
			os.write(reinterpret_cast<char*>(const_cast<leaf_node_t*>(leaf_p)), sizeof(leaf_node_t));
		}
		else if(const branch_node_t* branch_p = boost::get<branch_node_t>(node_p))
		{
			os << 'B';
			for(int i=0; i<8; ++i)
			{
				node_stack.push(&branch_p->children[7-i]);
			}
		}
		else
		{
			assert(false);
		}
	}

	os << std::endl;

	return os;
}

template<typename T>
std::istream& deserialise_node(std::istream& is, typename octree<T>::branch_or_leaf_node_t& node)
{
	char type;

	is >> type;

	switch(type)
	{
		case 'L':
			{
				typename octree<T>::leaf_node_t leaf;
				is.read(reinterpret_cast<char*>(&leaf), sizeof(leaf));
				node = leaf;
			}
			break;
		case 'B':
			{
				typename octree<T>::branch_node_t branch;
				for(int i=0; i<8; ++i)
				{
					deserialise_node<T>(is, branch.children[i]);
				}
				node = branch;
			}
			break;
		default:
			throw std::runtime_error(std::string("Unknown node type: ") + type);
	}

	return is;
}

template<typename T>
std::istream& octree<T>::deserialise(std::istream& is)
{
	is >> log_2_size_;
	is >> first_loc_.x >> first_loc_.y >> first_loc_.z;
	deserialise_node<T>(is, root_node_);

	return is;
}

// CRYSTALISED OCTREE

template<typename T>
crystalised_octree::crystalised_octree(const octree<T>& octree)
	: log_2_size_(0), min_loc_(0,0,0)
{
	assign_from(octree);
}

template<typename T>
void crystalised_octree::assign_from(const octree<T>& ot)
{
	typedef typename octree<T>::leaf_node_t leaf_node_t;
	typedef typename octree<T>::branch_node_t branch_node_t;
	typedef typename octree<T>::branch_or_leaf_node_t branch_or_leaf_node_t;

	log_2_size_ = ot.log_2_size();
	min_loc_ = ot.first_loc();

	// we will build up the representation as a deque of uint32_t's.
	std::deque<uint32_t> repr;

	// create a stack of nodes to process
	std::stack<const branch_or_leaf_node_t*> node_stack;
	node_stack.push(&ot.root());

	while(!node_stack.empty())
	{
		const branch_or_leaf_node_t* node_p = node_stack.top();
		node_stack.pop();

		if(const leaf_node_t* leaf_p = boost::get<leaf_node_t>(node_p))
		{
			// convert the leaf to a int32_t and check that the high-bit is unset
			int32_t data = static_cast<int32_t>(*leaf_p);
			assert(data >= 0);

			// promote the int32 to a uint32 and check that the integer representation on this
			// platform is sane!
			uint32_t u32_data = data;
			assert((u32_data & 0x80000000u) == 0);

			repr.push_back(u32_data);
		}
		else if(const branch_node_t* branch_p = boost::get<branch_node_t>(node_p))
		{
			// how many descendant nodes in total does this node have?
			int32_t n_descendants = 0;
			std::stack<const branch_or_leaf_node_t*> nodes_to_examine;
			nodes_to_examine.push(node_p);

			// visit nodes in a depth-first way
			while(nodes_to_examine.size() > 0)
			{
				++n_descendants;

				const branch_or_leaf_node_t* next_node = nodes_to_examine.top();
				nodes_to_examine.pop();

				if(const branch_node_t* branch_p = boost::get<branch_node_t>(next_node))
				{
					for(int i=0; i<8; ++i)
					{
						nodes_to_examine.push(&branch_p->children[i]);
					}
				}
			}

			// don't count this node as one of its descendants
			assert(n_descendants > 0);
			--n_descendants;

			// a branch is encoded by specifying the number of descendant nodes and setting the
			// high bit. Therefore one can skip this node by skipping the n_descendants following
			// nodes.
			uint32_t data = 0x80000000u | static_cast<uint32_t>(n_descendants);

			// signal this is a branch by pushing a uint32 with all bits set
			repr.push_back(data);

			// push children onto the stack, last child first
			for(int i=0; i<8; ++i)
			{
				node_stack.push(&branch_p->children[7-i]);
			}
		}
		else
		{
			assert(false);
		}
	}

	// copy the representation into a new vector
	data_ = boost::shared_ptr<std::vector<uint32_t> >(new std::vector<uint32_t>(repr.begin(), repr.end()));

	std::cout << "crystalised octree with " << ot.node_count() << " nodes into "
		<< data_->size() * sizeof(uint32_t) << " bytes." << std::endl;
}

// this is used by crystalised_octree::ray_intersect. Looking forward to C++11's lambdas yet?
// we use >= because we want to sort in the order we push onto the stack (i.e. farthest first)
inline bool record_cmp_(const boost::tuple<extent, size_t, float>& a, const boost::tuple<extent, size_t, float>& b)
{
	return boost::get<2>(a) >= boost::get<2>(b);
}

template<typename T>
bool crystalised_octree::ray_intersect(const ray& eye_ray_, sub_location& out_sub_loc) const
{
	// have to make a copy because of const correctness and I don't like const_casting...
	ray eye_ray(eye_ray_);

	// a node record is the extent of the node, it's index and the distance it's intersection is at
	typedef boost::tuple<extent, size_t, float> node_record;

	// cache the ray's origin
	location eye_ray_origin(eye_ray.x, eye_ray.y, eye_ray.z);

	// do we intersect this tree at all?
	aabox root_box(extent_().make_aabox());

	float distance;
	if(!slopeint_mul(&eye_ray, &root_box, &distance))
		return false;

	intersection_result result;

	// optimisation: only nodes which definitely intersect the ray are
	// pushed on this stack. The worst case usage for this stack is 8 nodes
	// right from root to leaf: log_2_size * 8.
	node_record stack[log_2_size_ << 3]; int stack_top = 0; //< where to insert the next node
	stack[stack_top++] = node_record(extent_(), 0, distance);

	while(stack_top > 0)
	{
		// pop the top-most record from the stack
		const node_record& record = stack[--stack_top];
		const extent node_ext = boost::get<0>(record);
		size_t node_idx = boost::get<1>(record);
		float distance = boost::get<2>(record);

		// NB: record gets clobbered below.

		if(is_branch(node_idx))
		{
			// we'll store child intersections on the stack
			node_record* intersections = stack + stack_top;
			size_t n_intersections = 0;

			// iterate over all children
			for(size_t i=0, child_idx = node_idx + 1; i<8; ++i)
			{
				// create a speculative record
				extent child_ext(location_of_child(i, node_ext), node_ext.size >> 1);
				size_t saved_child_idx(child_idx);

				// advance to next child
				if(is_branch(child_idx))
				{
					child_idx += 1 + (data_->at(child_idx) & 0x7fffffffu);
				}
				else
				{
					child_idx += 1;
				}

				// NB: DO NOT USE child_idx BELOW HERE, USE saved_child_idx

				// optimisation: if this child is a leaf node and is transparent, skip it
				if(!is_branch(saved_child_idx))
				{
					//if(child_ext.contains(eye_ray_origin))
					//	continue;

					T child_data(static_cast<int32_t>(data_->at(saved_child_idx)));

					// skip transparent leaf nodes
					if(child_data.is_transparent())
						continue;
				}

				// create a bounding box for the child
				aabox child_box(child_ext.make_aabox());

				// attempt to intersect ray with child
				float distance;
				if(slopeint_mul(&eye_ray, &child_box, &distance))
				{
					intersections[n_intersections] = node_record(child_ext, saved_child_idx, distance);
					++n_intersections;
				}
			}

			// if no intersections, bail
			if(n_intersections == 0)
				continue;

			// sort children by intersection distance
			std::sort(intersections, intersections + n_intersections, record_cmp_);

			// update stack pointer
			stack_top += n_intersections;
			assert(stack_top <= (log_2_size_ << 3));
		}
		else
		{
			if(distance <= 0.f)
				continue;

			// extract the leaf data
			T node_data(static_cast<int32_t>(data_->at(node_idx)));

			// is it transparent?
			if(node_data.is_transparent())
				continue;

			out_sub_loc.coords[0] = eye_ray.x + eye_ray.i * distance;
			out_sub_loc.coords[1] = eye_ray.y + eye_ray.j * distance;
			out_sub_loc.coords[2] = eye_ray.z + eye_ray.k * distance;
			out_sub_loc.node_extent = node_ext;

			return true;
		}
	}

	return false;
}

}

// IO OPERATORS

template<typename T>
std::istream& operator >> (std::istream& is, typename octree::octree<T>& tree)
{
	return tree.deserialise(is);
}

template<typename T>
std::ostream& operator << (std::ostream& os, const typename octree::octree<T>& tree)
{
	return tree.serialise(os);
}
