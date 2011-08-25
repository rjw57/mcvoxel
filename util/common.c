/* common.c */

#include <stdlib.h>
#include <stdio.h>

#include "common.h"

void*
pt_aligned_malloc(size_t size)
{
#if defined(__APPLE__) 
	/* OSX's malloc is already aligned to 16 bytes. */
	return malloc(size);
#else
	void* ptr = NULL;
	if(0 != posix_memalign(&ptr, 16, size)) {
		fprintf(stderr, "Warning: Memory allocation failed.\n");
	}
	return ptr;
#endif
}

void
pt_aligned_free(size_t PT_UNUSED(size), void* ptr)
{
	free(ptr);
}

static void
free_local_rand(gpointer rand)
{
	g_rand_free((GRand*)rand);
}

float
pt_random_float()
{
	/* we use per-thread storage here to use a different GRand per thread. */
	static GStaticPrivate this_thread_rand = G_STATIC_PRIVATE_INIT;
	GRand* rand = g_static_private_get(&this_thread_rand);
	if(NULL == rand) {
		rand = g_rand_new();
		g_static_private_set(&this_thread_rand, rand, free_local_rand);
	}
	return (float)(g_rand_double(rand));
}

/* vim:sw=4:ts=4:cindent
 */
