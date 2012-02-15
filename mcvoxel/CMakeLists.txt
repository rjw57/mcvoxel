include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_library(mcvoxel
	camera.cpp	camera.hpp
			colour.hpp
			datamodel.hpp
	io.cpp		io.hpp
	octree.cpp	octree.hpp
			ray.hpp
	sampling.cpp	sampling.hpp
	sky.cpp		sky.hpp
	world.cpp	world.hpp
)

target_link_libraries(mcvoxel
	c10t-mc c10t-image rayslope hdrloader util
	${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_IOSTREAMS_LIBRARY}
	${ZLIB_LIBRARY} ${GLIB_LIBRARIES} ${GTHREAD_LIBRARIES} ${PNG_LIBRARIES}
)

