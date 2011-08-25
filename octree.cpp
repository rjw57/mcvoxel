#include "io.hpp"
#include "octree.hpp"

namespace octree
{

location::location(long x, long y, long z)
	: x(x), y(y), z(z)
{ }

location::location(const location& loc)
	: x(loc.x), y(loc.y), z(loc.z)
{ }

const location& location::operator = (const location& loc)
{
	x = loc.x; y = loc.y; z = loc.z;
	return *this;
}

location location::operator + (const location& rhs) const
{
	return location(x+rhs.x, y+rhs.y, z+rhs.z);
}

location location::operator - (const location& rhs) const
{
	return location(x-rhs.x, y-rhs.y, z-rhs.z);
}

// CRYSTALISED OCTREE

crystalised_octree::crystalised_octree(int32_t log_2_size_, const location& min_loc)
	: log_2_size_(log_2_size_), min_loc_(min_loc_), data_(new std::vector<uint32_t>())
{ }

crystalised_octree::~crystalised_octree()
{ }

crystalised_octree::crystalised_octree(const crystalised_octree& ot)
	: log_2_size_(ot.log_2_size_), min_loc_(ot.min_loc_), data_(ot.data_)
{ }

const crystalised_octree& crystalised_octree::operator = (const crystalised_octree& ot)
{
	log_2_size_ = ot.log_2_size_;
	min_loc_ = ot.min_loc_;
	data_ = ot.data_;
	return *this;
}

size_t crystalised_octree::child_containing(const location& loc, size_t node_idx, const extent& node_ext,
		extent& r_child_ext) const
{
	assert(is_branch(node_idx));

	// sanity check size
	assert(node_ext.size > 1);

	// sanity check location
	assert(node_ext.contains(loc));

	// which child _should_ contain the point?
	size_t which_child = index_of_child_containing(loc, node_ext);

	// move to that child
	size_t child_idx = node_idx + 1;
	for(size_t child=0; child != which_child; ++child)
	{
		assert(child_idx < data_->size());

		// advance to next child
		if(is_branch(child_idx))
		{
			child_idx += 1 + (data_->at(child_idx) & 0x7fffffffu);
		}
		else
		{
			child_idx += 1;
		}
	}

	// set the child extent and return the child index
	r_child_ext = extent(location_of_child(which_child, node_ext), node_ext.size >> 1);
	return child_idx;
}

const int32_t crystalised_octree::get(const location& loc) const
{
	assert(extent_().contains(loc));
	assert(data_->size() > 0);

	// start at the root
	size_t node_idx = 0;
	extent node_ext(extent_());

	// recurse to the leaf
	while(is_branch(node_idx))
	{
		node_idx = child_containing(loc, node_idx, node_ext, node_ext);
	}

	assert(!is_branch(node_idx));
	return static_cast<int32_t>(data_->at(node_idx));
}

std::ostream& crystalised_octree::serialise(std::ostream& out) const
{
	// check - write 0xbeefface
	io::nbo::write(out, uint32_t(0xbeeffaceu));

	// write the log-2 size
	io::nbo::write(out, static_cast<uint32_t>(log_2_size_));

	// write the location
	io::nbo::write(out, static_cast<int32_t>(min_loc_.x));
	io::nbo::write(out, static_cast<int32_t>(min_loc_.y));
	io::nbo::write(out, static_cast<int32_t>(min_loc_.z));

	// write the data
	uint32_t n_data = data_->size();
	io::nbo::write(out, n_data);
	for(uint32_t idx=0; idx<n_data; ++idx)
	{
		io::nbo::write(out, data_->at(idx));
	}

	// check - write 0xdeadbeef
	io::nbo::write(out, uint32_t(0xdeadbeefu));

	return out;
}

std::istream& crystalised_octree::deserialise(std::istream& in)
{
	uint32_t check, l2s;
	int32_t mlx, mly, mlz;

	io::nbo::read(in, check);
	assert(check == 0xbeeffaceu);

	// read the log-2 size
	io::nbo::read(in, l2s);
	log_2_size_ = l2s;

	// read the location
	io::nbo::read(in, mlx);
	io::nbo::read(in, mly);
	io::nbo::read(in, mlz);
	min_loc_ = location(mlx, mly, mlz);

	// read the data
	uint32_t n_data = 0;
	io::nbo::read(in, n_data);

	data_ = boost::shared_ptr<std::vector<uint32_t> >(new std::vector<uint32_t>());
	data_->reserve(n_data);

	for(uint32_t idx=0; idx<n_data; ++idx)
	{
		uint32_t v;
		io::nbo::read(in, v);
		data_->push_back(v);
	}

	io::nbo::read(in, check);
	assert(check == 0xdeadbeefu);

	return in;
}

}
