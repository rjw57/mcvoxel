#ifndef INSIDE_IO_HPP
#  error "This file must only be included by io.hpp"
#endif

#include <boost/foreach.hpp>
#include "octree.hpp"

#if 0
template< typename T >
std::ostream& operator << (std::ostream& os, const octree::tree<T>& tree)
{
	const octree::node<T>& root(tree.root());
	os.write(reinterpret_cast<const char*>(&root.min_loc), sizeof(root.min_loc));
	os.write(reinterpret_cast<const char*>(&root.log_2_size), sizeof(root.log_2_size));

	std::deque<const octree::node<T>*> serialise_list;
	serialise_list.push_back(&root);

	while(!serialise_list.empty())
	{
		const octree::node<T>* next_node = serialise_list.back();
		serialise_list.pop_back();
		char type;

		if(next_node == NULL)
		{
			type = 'N';
			os.write(&type, 1);
		}
		else if(next_node->is_leaf())
		{
			type = 'L';
			os.write(&type, 1);
			os.write(reinterpret_cast<const char*>(&next_node->value), sizeof(next_node->value));
		}
		else
		{
			type = 'B';
			os.write(&type, 1);
			assert(next_node->children.size() == 8);
			BOOST_FOREACH(const octree::node<T>* child_node, next_node->children)
			{
				serialise_list.push_back(child_node);
			}
		}
	}

	return os;
}
#endif
