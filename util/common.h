/* common.h */

#ifndef _PT_COMMON_H
#define _PT_COMMON_H

#include <glib.h>

G_BEGIN_DECLS

/**
 * PT_UNUSED:
 * 
 * Mark a function argument as unused to quieten a compiler warning.
 */
#ifdef PT_UNUSED
#elif defined(__GNUC__)
# define PT_UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define PT_UNUSED(x) /*@unused@*/ x
#else
# define PT_UNUSED(x) x
#endif

void*
pt_aligned_malloc(size_t size);

void
pt_aligned_free(size_t size, void* ptr);

float
pt_random_float();

/**
 * pt_power_fi:
 *
 * Raise a float to an integer power.
 *
 * See: NVIDIA copyright below.
 */
static __inline__ float __attribute__((__always_inline__))
pt_power_fi(float a, int b)
{
	guint e = MAX(b,-b);
	float r = 1.f;

  	while (1) {
	  	if ((e & 1) != 0) {
			r = r * a;
		}
		e = e >> 1;
		if (e == 0) {
			return b < 0 ? 1.0f/r : r;
		}
		a = a * a;
	}
}

/**
 * Some portions of this file (identified as suce) are under the following
 * license:
 *
 * Copyright 1993-2008 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code are
 * hereby granted a nonexclusive, royalty-free license to use this code in
 * individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE CODE FOR
 * ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
 * ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOURCE CODE,
 * INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND
 * FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT SHALL NVIDIA BE LIABLE FOR
 * ANY SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,  WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR
 * IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.   This source code is a "commercial item" as that
 * term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of "commercial
 * computer  software"  and "commercial computer software documentation" as
 * such terms are  used in 48 C.F.R. 12.212 (SEPT 1995) and is provided to the
 * U.S. Government only as a commercial end item.  Consistent with 48
 * C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 227.7202-4 (JUNE 1995), all
 * U.S. Government End Users acquire the source code with only those rights set
 * forth herein. 
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code, the
 * above Disclaimer and U.S. Government End Users Notice.
 */

G_END_DECLS

#endif /* _PT_COMMON_H */

/* vim:sw=4:ts=4:cindent
 */
