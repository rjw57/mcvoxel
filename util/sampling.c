/* sampling.c */

#include <math.h>

#include "sampling.h"
#include "common.h"

/*========================================================================*/
/* Draw a sample uniformly over the sphere */
Vector 
pt_sampling_sphere()
{
	float z = 2.f * (pt_random_float() - 0.5f);
	float t = 2.0 * G_PI * pt_random_float();
	float s = sqrtf(1.f - (z*z));
	float x = s * cosf(t);
	float y = s * sinf(t);

	return pt_vector_make(x,y,z,0.f);
}

/*========================================================================*/
/* Draw a sample uniformly from the +ve hemisphere centred
 * on normal. */
Vector 
pt_sampling_hemisphere(Vector normal)
{
	/* find direction uniform over +ve z hemisphere */
	float z = pt_random_float();
	float t = 2.0 * G_PI * pt_random_float();
	float s = sqrtf(1.f - (z*z));
	float x = s * cosf(t);
	float y = s * sinf(t);

	Vector e1 = pt_vector_perpendicular(normal);
	Vector e2 = pt_vector_cross(e1, normal);

	return pt_vector_scale(x,e1) + pt_vector_scale(y,e2) + pt_vector_scale(z,normal);
}

/*========================================================================*/
/* Draw a sample from a cosine weighted distribution over
 * the hemisphere centred on normal. */
Vector 
pt_sampling_cosine(Vector normal)
{
	float two_pi_r1 = pt_random_float() * 2.f * G_PI;
	float sr2 = sqrtf(pt_random_float());

	float x = cosf(two_pi_r1) * sr2;
	float y = sinf(two_pi_r1) * sr2;
	float z = sqrtf(1.f - (sr2*sr2));

	Vector e1 = pt_vector_perpendicular(normal);
	Vector e2 = pt_vector_cross(e1, normal);

	return pt_vector_scale(x,e1) + pt_vector_scale(y,e2) + pt_vector_scale(z,normal);
}

/*========================================================================*/
/* Draw a sample from a exponentiated cosine weighted 
 * distribution over the hemisphere centred on normal. */
Vector 
pt_sampling_exponentiated_cosine(Vector normal, guint power)
{
	/* see: http://http.developer.nvidia.com/GPUGems3/gpugems3_ch20.html */

	float z1 = pt_random_float();
	float z2 = pt_random_float();

	float z = powf(z1, 1.f / (power+1.f));
	float theta = acosf(z);
	float phi = 2.f * 3.14159f * z2;

	float s = sinf(theta);
	float x = cosf(phi) * s;
	float y = sinf(phi) * s;
	z = fabs(z);

	Vector e1 = pt_vector_perpendicular(normal);
	Vector e2 = pt_vector_cross(e1, normal);

	return pt_vector_scale(x,e1) + pt_vector_scale(y,e2) + pt_vector_scale(z,normal);
}

/* vim:sw=4:ts=4:cindent
 */
