#ifndef _OCTREE_HPP
#define _OCTREE_HPP

#include <iostream>
#include <stdint.h>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/variant.hpp>

#include <rayslope/ray.h>
#include <rayslope/aabox.h>

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
	branch_or_leaf_node_t dummy; // required with GCC 4.6.0 FSR

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

	bool contains(const location& tloc) const
	{
		if((tloc.x < loc.x) || (tloc.x >= loc.x + size)) return false;
		if((tloc.y < loc.y) || (tloc.y >= loc.y + size)) return false;
		if((tloc.z < loc.z) || (tloc.z >= loc.z + size)) return false;
		return true;
	}

	::aabox make_aabox() const;
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

// an octree than can *only* store +ve int32 values (i.e. high bit must be 0)
// since crystalised_octrees are immutable, their data is implicitly shared.
class crystalised_octree
{
	public:
		crystalised_octree(int32_t log_2_size, const location& min_loc = location(0,0,0));
		~crystalised_octree();

		// copy and assignment constructors
		crystalised_octree(const crystalised_octree& octree);
		const crystalised_octree& operator = (const crystalised_octree& octree);

		template<typename T>
		explicit crystalised_octree(const octree<T>& octree);

		// simple inline functions
		int32_t size() const { return static_cast<int32_t>(1) << log_2_size_; }

		// get the extent of this octree
		extent extent_() const { return extent(min_loc_, size()); }

		// retrieve data
		const int32_t get(const location& loc) const;
		const int32_t get(long x, long y, long z) const { return get(location(x,y,z)); }

		// intersect ray (treating lead data as T). This requires T::is_transparent().
		template<typename T>
		bool ray_intersect(const ray& r, sub_location& out_sub_loc) const;

		// serialisation
		std::ostream& serialise(std::ostream& os) const;
		std::istream& deserialise(std::istream& is);

	protected:
		bool is_branch(size_t idx) const { return (data_->at(idx) & 0x80000000u); }

		size_t child_containing(const location& loc, size_t node_idx, const extent& node_ext,
				extent& r_child_ext) const;

		template<typename T>
		void assign_from(const octree<T>& octree);

	private:
		int32_t                                   log_2_size_;
		location                                  min_loc_;
		boost::shared_ptr<std::vector<uint32_t> > data_;
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
