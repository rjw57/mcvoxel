/* sampling.h */

#ifndef _PT_SAMPLING_H
#define _PT_SAMPLING_H

#include <glib.h>

#include "vector.h"
#include "common.h"

G_BEGIN_DECLS

/* Draw a sample uniformly over the sphere */
Vector 
pt_sampling_sphere();

/* Draw a sample uniformly from the +ve hemisphere centred
 * on normal. */
Vector 
pt_sampling_hemisphere(Vector normal);

/* Draw a sample from a cosine weighted distribution over
 * the hemisphere centred on normal. */
Vector 
pt_sampling_cosine(Vector normal);

/* Draw a sample from a exponentiated cosine weighted 
 * distribution over the hemisphere centred on normal. */
Vector 
pt_sampling_exponentiated_cosine(Vector normal, guint power);

/* Return the likelihood of having drawn dir from the
 * exponentiated cosine weighted distribution over the 
 * hemisphere centred on normal. */
static __inline__ float __attribute__((__always_inline__))
pt_sampling_exponentiated_cosine_likelihood(Vector normal, guint power, Vector dir)
{
	float cos = MAX(0.f, pt_vector_get_w(pt_vector_dot3(dir, normal)));
	return pt_power_fi(cos, power);
}

/* Return the likelihood of having drawn dir from the
 * cosine weighted distribution over the hemisphere centred 
 * on normal. */
static __inline__ float __attribute__((__always_inline__))
pt_sampling_cosine_likelihood(Vector normal, Vector dir)
{
	return MAX(0.f, pt_vector_get_w(pt_vector_dot3(dir, normal)));
}

G_END_DECLS

#endif /* _PT_SAMPLING_H */

/* vim:sw=4:ts=4:cindent
 */
