/* vector.c */

#include <math.h>

#include "vector.h"

/*========================================================================*/
Triangle 
pt_triangle_make(Vector a, Vector b, Vector c)
{
	Triangle t;

	t.origin = pt_vector_set_w(a, 1.f);
	t.delta[0] = pt_vector_set_w(b-a, 0.f);
	t.delta[1] = pt_vector_set_w(c-a, 0.f);

	t.normal = pt_vector_normalise3(pt_vector_cross(t.delta[0], t.delta[1]));

	t.origin = pt_vector_set_w(t.origin,
			pt_vector_get_x(pt_vector_dot3(t.delta[0], t.delta[1])));

	t.delta[0] = pt_vector_set_w(t.delta[0],
			pt_vector_get_x(pt_vector_dot3(t.delta[0], t.delta[0])));

	t.delta[1] = pt_vector_set_w(t.delta[1],
			pt_vector_get_x(pt_vector_dot3(t.delta[1], t.delta[1])));

	t.normal = pt_vector_set_w(t.normal,
			1.f / pt_vector_get_w(
				pt_vector_w_sub(
					pt_vector_w_mul(t.origin, t.origin),
					pt_vector_w_mul(t.delta[0], t.delta[1])
					)));
	return t;
}

static __inline__ float __attribute__((__always_inline__))
dot(Vector a, Vector b)
{
	return pt_vector_get_w(pt_vector_dot4(a,pt_vector_zero_w(b)));
}

/*========================================================================*/
float 
pt_ray_intersect_triangle_max_dist(const Ray ray, const Triangle triangle,
		float max_dist, float* p_bary_coords)
{
    Vector    w;          /* ray vectors */
    float     r, a, b;        /* params to calc ray-plane intersect */
	Vector    I;

    a = -dot(triangle.normal,ray.origin - triangle.origin);
    b = dot(triangle.normal,ray.direction);
    if (b == 0.f) {
		/* ray is parallel to triangle plane */
		return -1.f;
    }

    /* get intersect point of ray with triangle plane */
    r = a / b;
    if ((r < 0.f) || (r > max_dist))    /* ray goes away from triangle */
        return -1.f;                    /* => no intersect */

	I = pt_ray_point_along(ray, r);

    /* is I inside T? */
    w = pt_vector_zero_w(I - triangle.origin);
	
	Vector wdotu = pt_vector_dot4(w, triangle.delta[0]);
	Vector wdotv = pt_vector_dot4(w, triangle.delta[1]);

    /* get and test parametric coords */
    float s, t;

	s = pt_vector_get_w(pt_vector_w_mul(triangle.normal,
				pt_vector_w_sub(
					pt_vector_w_mul(triangle.origin, wdotv),
					pt_vector_w_mul(triangle.delta[1], wdotu)
					)
			));

	t = pt_vector_get_w(pt_vector_w_mul(triangle.normal,
				pt_vector_w_sub(
					pt_vector_w_mul(triangle.origin, wdotu),
					pt_vector_w_mul(triangle.delta[0], wdotv)
					)
			));

    if (s < 0.f || t < 0.f || (s + t) > 1.f)  /* I is outside T */
        return -1.f;

	if(p_bary_coords != NULL) {
		p_bary_coords[0] = s;
		p_bary_coords[1] = t;
	}

	return r;
}

/*========================================================================*/
float 
pt_ray_intersect_sphere(const Ray ray, const Sphere sphere)
{
	Vector dst, B, C, D, lambda;
	gboolean is_valid;

	/* Algorithm from http://devmaster.net/wiki/Ray-sphere_intersection */

	dst = ray.origin - sphere; /* w component is now garbage */

	/* The w=0.0 of direction will nobble the w component of dst */
	B = pt_vector_dot4(dst, ray.direction); 

	/* We do most of the math in the w-component (0th) now */
	C = pt_vector_w_sub(pt_vector_dot3(dst,dst), pt_vector_w_mul(sphere, sphere));
	D = pt_vector_w_sub(pt_vector_w_mul(B,B), C);

	is_valid = _mm_comigt_ss(D, pt_vector_make_zero());
	if(!is_valid) {
		return -1;
	}

	lambda = pt_vector_make_zero();
	lambda = pt_vector_w_sub(lambda, B);
	lambda = pt_vector_w_sub(lambda, _mm_sqrt_ss(D));

	return pt_vector_get_w(lambda);
}


/* vim:sw=4:ts=4:cindent
 */
