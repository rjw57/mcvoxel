/* vector.h */

#ifndef _PT_VECTOR_H
#define _PT_VECTOR_H

#include <glib.h>
#include <xmmintrin.h>

G_BEGIN_DECLS

/* This is all highly GCC/ICC specific. */

/*========================================================================*/
/* Note that we store vectors in memory as {w, x, y, z} but the API exposes
 * a notional x,y,z,w order. */
typedef float Vector __attribute__ ((vector_size(16)));

#define pt_vector_make(x, y, z, w) (_mm_set_ps((z), (y), (x), (w)))

#define pt_vector_make_zero() (_mm_setzero_ps())

/* Create a vector with all values equal to f. */
#define pt_vector_make_broadcast(f) (_mm_set_ps1(f))

/* W-component only operations */
#define pt_vector_w_add(a,b) (_mm_add_ss((a),(b)))
#define pt_vector_w_sub(a,b) (_mm_sub_ss((a),(b)))
#define pt_vector_w_mul(a,b) (_mm_mul_ss((a),(b)))
#define pt_vector_w_div(a,b) (_mm_div_ss((a),(b)))
#define pt_vector_w_sqrt(vec) (_mm_sqrt_ss(vec))
#define pt_vector_w_rsqrt(vec) (_mm_rsqrt_ss(vec))
#define pt_vector_w_reciprocal(v) (_mm_rcp_ss(v))

#define pt_vector_w_max(a,b) (_mm_max_ss((a),(b)))
#define pt_vector_w_min(a,b) (_mm_min_ss((a),(b)))

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_w_abs(Vector vec)
{
	return pt_vector_w_max(
			vec, 
			pt_vector_w_sub(pt_vector_make_zero(), vec));
}

/* Element-wise 4-way operations */
#define pt_vector_add(a,b) (_mm_add_ps((a),(b)))
#define pt_vector_sub(a,b) (_mm_sub_ps((a),(b)))
#define pt_vector_mul(a,b) (_mm_mul_ps((a),(b)))
#define pt_vector_div(a,b) (_mm_div_ps((a),(b)))
#define pt_vector_sqrt(vec) (_mm_sqrt_ps(vec))
#define pt_vector_rsqrt(vec) (_mm_rsqrt_ps(vec))
#define pt_vector_reciprocal(v) (_mm_rcp_ps(v))

#define pt_vector_max(a,b) (_mm_max_ps((a),(b)))
#define pt_vector_min(a,b) (_mm_min_ps((a),(b)))

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_abs(Vector vec)
{
	return pt_vector_max(
			vec, 
			pt_vector_sub(pt_vector_make_zero(), vec));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_neg(Vector vec)
{
	return pt_vector_sub(pt_vector_make_zero(), vec);
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_scale(float scale, Vector vec)
{
	return pt_vector_mul(pt_vector_make_broadcast(scale), vec);
}

#define pt_vector_get_x(vec) (__builtin_ia32_vec_ext_v4sf((vec),1))
#define pt_vector_get_y(vec) (__builtin_ia32_vec_ext_v4sf((vec),2))
#define pt_vector_get_z(vec) (__builtin_ia32_vec_ext_v4sf((vec),3))
#define pt_vector_get_w(vec) (__builtin_ia32_vec_ext_v4sf((vec),0))

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_set_x(Vector vec, float val)
{
	return pt_vector_make(
			val,					pt_vector_get_y(vec),
			pt_vector_get_z(vec),	pt_vector_get_w(vec) );
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_set_y(Vector vec, float val)
{
	return pt_vector_make(
			pt_vector_get_x(vec),	val,
			pt_vector_get_z(vec),	pt_vector_get_w(vec) );
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_set_z(Vector vec, float val)
{
	return pt_vector_make(
			pt_vector_get_x(vec),	pt_vector_get_y(vec),
			val,					pt_vector_get_w(vec) );
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_set_w(Vector vec, float val)
{
	return pt_vector_make(
			pt_vector_get_x(vec),	pt_vector_get_y(vec),
			pt_vector_get_z(vec),	val);
}

#if 0
/* Return a vector where the w component is the logical
 * OR of all components. */
static __inline__ Vector __attribute__((__always_inline__))
pt_vector_any(Vector vec)
{
	Vector s1 = _mm_or_ps( vec,
			_mm_shuffle_ps(vec,vec,_MM_SHUFFLE(0,1,2,3)) );
	Vector s2 = _mm_shuffle_ps(s1,s1,_MM_SHUFFLE(0,0,0,2));
	return _mm_or_ps_ss(s1,s2);
}
#endif

#define pt_vector_zero_x(vec) (pt_vector_mul(vec, pt_vector_make(0,1,1,1)))
#define pt_vector_zero_y(vec) (pt_vector_mul(vec, pt_vector_make(1,0,1,1)))
#define pt_vector_zero_z(vec) (pt_vector_mul(vec, pt_vector_make(1,1,0,1)))
#define pt_vector_zero_w(vec) (pt_vector_mul(vec, pt_vector_make(1,1,1,0)))

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_broadcast_x(Vector vec)
{
	return _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(1,1,1,1));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_broadcast_y(Vector vec)
{
	return _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(2,2,2,2));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_broadcast_z(Vector vec)
{
	return _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(3,3,3,3));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_broadcast_w(Vector vec)
{
	return _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(0,0,0,0));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_horiz_sum(Vector vec)
{
	Vector abs1, abs2, abs3;
	abs1 = vec + _mm_shuffle_ps(vec, vec, _MM_SHUFFLE(0,1,2,3));
	abs2 = _mm_shuffle_ps(abs1, abs1, _MM_SHUFFLE(0,0,0,2));
	abs3 = _mm_add_ss(abs2, abs1);
	return _mm_shuffle_ps(abs3, abs3, _MM_SHUFFLE(0,0,0,0));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_dot2(Vector a, Vector b)
{
	Vector ab = a * b * pt_vector_make(1,1,0,0);
	return pt_vector_horiz_sum(ab);
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_dot3(Vector a, Vector b)
{
	Vector ab = a * b * pt_vector_make(1,1,1,0);
	return pt_vector_horiz_sum(ab);
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_dot4(Vector a, Vector b)
{
	Vector ab = a * b;
	return pt_vector_horiz_sum(ab);
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_normalise2(Vector vec)
{
	return vec * pt_vector_broadcast_w(pt_vector_w_rsqrt(pt_vector_dot2(vec,vec)));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_normalise3(Vector vec)
{
	return vec * pt_vector_broadcast_w(pt_vector_w_rsqrt(pt_vector_dot3(vec,vec)));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_normalise4(Vector vec)
{
	return vec * pt_vector_broadcast_w(pt_vector_w_rsqrt(pt_vector_dot4(vec,vec)));
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_cross(Vector a, Vector b)
{
	/* R.x = A.y * B.z - A.z * B.y */
	/* R.y = A.z * B.x - A.x * B.z */
	/* R.z = A.x * B.y - A.y * B.x */

	Vector a1, a2, b1, b2;

	a1 = _mm_shuffle_ps(a, a, _MM_SHUFFLE(1,3,2,0));
	a2 = _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,1,3,0));
	b1 = _mm_shuffle_ps(b, b, _MM_SHUFFLE(2,1,3,0));
	b2 = _mm_shuffle_ps(b, b, _MM_SHUFFLE(1,3,2,0));

	/* Zero the w-component */
	return pt_vector_zero_w(a1*b1 - a2*b2);
}

/* Return some vector perpendicular to vec. Vec should be normalised. */
static __inline__ Vector __attribute__((__always_inline__))
pt_vector_perpendicular(Vector vec)
{
	Vector abs_vec = pt_vector_abs(vec);
	Vector min_vec = pt_vector_make(0,0,1,0);
	if(pt_vector_get_x(abs_vec) < pt_vector_get_y(abs_vec)) {
		if(pt_vector_get_z(abs_vec) > pt_vector_get_x(abs_vec)) {
			min_vec = pt_vector_make(1,0,0,0);
		}
	} else {
		if(pt_vector_get_z(abs_vec) > pt_vector_get_y(abs_vec)) {
			min_vec = pt_vector_make(0,1,0,0);
		}
	}
	return pt_vector_cross(min_vec, vec);
}

static __inline__ Vector __attribute__((__always_inline__))
pt_vector_clamp(Vector vec, Vector min, Vector max)
{
	return pt_vector_max(min, pt_vector_min(max, vec));
}

/*========================================================================*/
/* A sphere stores its origin in x,y,z and it's radius in w. */
typedef Vector Sphere;

/* A sphere has a center and a radius */
#define pt_sphere_make(x,y,z,r) (pt_vector_make(x,y,z,r))

static __inline__ Vector __attribute__((__always_inline__))
pt_sphere_normal_at_point(Vector sphere, Vector point)
{
	return pt_vector_normalise3(point - sphere);
}

/*========================================================================*/
typedef struct {
	Vector	origin;		/* w = delta[0].delta[1] */
	Vector  normal;		/* w = 1 / (origin.w*origin.w - delta[1].w*delta[0].w) */
	Vector	delta[2];	/* w = delta.delta */
} Triangle;

Triangle 
pt_triangle_make(Vector a, Vector b, Vector c);

/*========================================================================*/
static __inline__ Vector __attribute__((__always_inline__))
pt_triangle_centroid(Triangle t)
{
	return t.origin + 
		pt_vector_make_broadcast(1.f/3.f) * (t.delta[0]+t.delta[1]);
}

/*========================================================================*/
typedef struct {
	Vector	origin;		/* w = 1 */
	Vector	direction;	/* w = 0 */
	Vector	one_over_direction;	/* w = 0 */
} Ray;

/*========================================================================*/
static __inline__ Ray __attribute__((__always_inline__))
pt_ray_make(float ox, float oy, float oz, float dx, float dy, float dz)
{
	Vector dir = pt_vector_normalise4(pt_vector_make(dx,dy,dz,0.f));
	Ray output;
	output.origin = pt_vector_make(ox,oy,oz,1.f);
	output.direction = dir;
	output.one_over_direction = pt_vector_reciprocal(dir);
	return output;
}

/*========================================================================*/
static __inline__ Ray __attribute__((__always_inline__))
pt_ray_make_vector(Vector origin, Vector direction)
{
	Vector dir = pt_vector_normalise3(direction);
	Ray output;
	output.origin = pt_vector_mul(origin, pt_vector_make(1,1,1,0));
	output.direction = dir;
	output.one_over_direction = pt_vector_reciprocal(dir);
	return output;
}

/*========================================================================*/
static __inline__ Vector __attribute__((__always_inline__))
pt_ray_point_along(const Ray ray, float lambda)
{
	return ray.origin + ray.direction * pt_vector_make_broadcast(lambda);
}

/*========================================================================*/
/* lambda must have all elements equal to the distance */
static __inline__ Vector __attribute__((__always_inline__))
pt_ray_point_along_vector(const Ray ray, Vector lambda)
{
	return ray.origin + ray.direction * lambda;
}

/* The intersect functions return a -ve distance if there is no
 * intersection. */

/*========================================================================*/
float 
pt_ray_intersect_sphere(const Ray ray, const Sphere sphere);

/*========================================================================*/
float 
pt_ray_intersect_triangle_max_dist(const Ray ray, const Triangle triangle,
		float max_dist, float* p_bary_coords);

/*========================================================================*/
static __inline__ float __attribute__((__always_inline__))
pt_ray_intersect_triangle(const Ray ray, const Triangle triangle, 
		float* p_bary_coords)
{
	return pt_ray_intersect_triangle_max_dist(ray, triangle, G_MAXFLOAT, p_bary_coords);
}

G_END_DECLS

#endif /* _PT_VECTOR_H */

/* vim:sw=4:ts=4:cindent
 */
