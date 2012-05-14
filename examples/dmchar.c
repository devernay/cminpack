#include <math.h>
#include "cminpack.h"
#define real __cminpack_real__

extern void dmchar_(int *ibeta, int *it, int *irnd, 
                   int *ngrd, int *machep, int *negep, int *iexp, 
                   int *minexp, int *maxexp, real *eps, real *epsneg,
                   real *xmin, real *xmax);

void dmchar_(int *ibeta, int *it, int *irnd, 
             int *ngrd, int *machep, int *negep, int *iexp, 
             int *minexp, int *maxexp, real *eps, real *epsneg,
             real *xmin, real *xmax)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static real a, b;
    static int i__, j, k;
    static real y, z__;
    static int iz, mx;
    static real one, beta, zero, betam1, betain;



/*     this subroutine is intended to determine the characteristics */
/*     of the floating-point arithmetic system that are specified */
/*     below.  the first three are determined according to an */
/*     algorithm due to m. malcolm, cacm 15 (1972), pp. 949-951, */
/*     incorporating some, but not all, of the improvements */
/*     suggested by m. gentleman and s. marovich, cacm 17 (1974), */
/*     pp. 276-277. */


/*       ibeta   - the radix of the floating-point representation */
/*       it      - the number of base ibeta digits in the floating-point */
/*                 significand */
/*       irnd    - 0 if floating-point addition chops, */
/*                 1 if floating-point addition rounds */
/*       ngrd    - the number of guard digits for multiplication.  it is */
/*                 0 if  irnd=1, or if  irnd=0  and only  it  base  ibeta */
/*                   digits participate in the post normalization shift */
/*                   of the floating-point significand in multiplication */
/*                 1 if  irnd=0  and more than  it  base  ibeta  digits */
/*                   participate in the post normalization shift of the */
/*                   floating-point significand in multiplication */
/*       machep  - the largest negative integer such that */
/*                 1.0+float(ibeta)**machep .ne. 1.0, except that */
/*                 machep is bounded below by  -(it+3) */
/*       negeps  - the largest negative integer such that */
/*                 1.0-float(ibeta)**negeps .ne. 1.0, except that */
/*                 negeps is bounded below by  -(it+3) */
/*       iexp    - the number of bits (decimal places if ibeta = 10) */
/*                 reserved for the representation of the exponent */
/*                 (including the bias or sign) of a floating-point */
/*                 number */
/*       minexp  - the largest in magnitude negative integer such that */
/*                 float(ibeta)**minexp is a positive floating-point */
/*                 number */
/*       maxexp  - the largest positive integer exponent for a finite */
/*                 floating-point number */
/*       eps     - the smallest positive floating-point number such */
/*                 that  1.0+eps .ne. 1.0. in particular, if either */
/*                 ibeta = 2  or  irnd = 0, eps = float(ibeta)**machep. */
/*                 otherwise,  eps = (float(ibeta)**machep)/2 */
/*       epsneg  - a small positive floating-point number such that */
/*                 1.0-epsneg .ne. 1.0. in particular, if ibeta = 2 */
/*                 or  irnd = 0, epsneg = float(ibeta)**negeps. */
/*                 otherwise,  epsneg = (ibeta**negeps)/2.  because */
/*                 negeps is bounded below by -(it+3), epsneg may not */
/*                 be the smallest number which can alter 1.0 by */
/*                 subtraction. */
/*       xmin    - the smallest non-vanishing floating-point power of the */
/*                 radix.  in particular,  xmin = float(ibeta)**minexp */
/*       xmax    - the largest finite floating-point number.  in */
/*                 particular   xmax = (1.0-epsneg)*float(ibeta)**maxexp */
/*                 note - on some machines  xmax  will be only the */
/*                 second, or perhaps third, largest number, being */
/*                 too small by 1 or 2 units in the last digit of */
/*                 the significand. */

/*     latest revision - october 22, 1979 */

/*     author - w. j. cody */
/*              argonne national laboratory */

/* ----------------------------------------------------------------- */
    one = 1.;
    zero = 0.;
/* ----------------------------------------------------------------- */
/*     determine ibeta,beta ala malcolm */
/* ----------------------------------------------------------------- */
    a = one;
L10:
    a += a;
    if (a + one - a - one == zero) {
	goto L10;
    }
    b = one;
L20:
    b += b;
    if (a + b - a == zero) {
	goto L20;
    }
    *ibeta = (int) ((real) (a + b - a));
    beta = (real) ((real) (*ibeta));
/* ----------------------------------------------------------------- */
/*     determine it, irnd */
/* ----------------------------------------------------------------- */
    *it = 0;
    b = one;
L100:
    ++(*it);
    b *= beta;
    if (b + one - b - one == zero) {
	goto L100;
    }
    *irnd = 0;
    betam1 = beta - one;
    if (a + betam1 - a != zero) {
	*irnd = 1;
    }
/* ----------------------------------------------------------------- */
/*     determine negep, epsneg */
/* ----------------------------------------------------------------- */
    *negep = *it + 3;
    betain = one / beta;
    a = one;

    i__1 = *negep;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a *= betain;
/* L200: */
    }

    b = a;
L210:
    if (one - a - one != zero) {
	goto L220;
    }
    a *= beta;
    --(*negep);
    goto L210;
L220:
    *negep = -(*negep);
    *epsneg = a;
    if (*ibeta == 2 || *irnd == 0) {
	goto L300;
    }
    a = a * (one + a) / (one + one);
    if (one - a - one != zero) {
	*epsneg = a;
    }
/* ----------------------------------------------------------------- */
/*     determine machep, eps */
/* ----------------------------------------------------------------- */
L300:
    *machep = -(*it) - 3;
    a = b;
L310:
    if (one + a - one != zero) {
	goto L320;
    }
    a *= beta;
    ++(*machep);
    goto L310;
L320:
    *eps = a;
    if (*ibeta == 2 || *irnd == 0) {
	goto L350;
    }
    a = a * (one + a) / (one + one);
    if (one + a - one != zero) {
	*eps = a;
    }
/* ----------------------------------------------------------------- */
/*     determine ngrd */
/* ----------------------------------------------------------------- */
L350:
    *ngrd = 0;
    if (*irnd == 0 && (one + *eps) * one - one != zero) {
	*ngrd = 1;
    }
/* ----------------------------------------------------------------- */
/*     determine iexp, minexp, xmin */

/*     loop to determine largest i and k = 2**i such that */
/*         (1/beta) ** (2**(i)) */
/*     does not underflow */
/*     exit from loop is signaled by an underflow. */
/* ----------------------------------------------------------------- */
    i__ = 0;
    k = 1;
    z__ = betain;
L400:
    y = z__;
    z__ = y * y;
/* ----------------------------------------------------------------- */
/*        check for underflow here */
/* ----------------------------------------------------------------- */
    a = z__ * one;
    if (a + a == zero || fabs(z__) >= y) {
	goto L410;
    }
    ++i__;
    k += k;
    goto L400;
L410:
    if (*ibeta == 10) {
	goto L420;
    }
    *iexp = i__ + 1;
    mx = k + k;
    goto L450;
/* ----------------------------------------------------------------- */
/*     for decimal machines only */
/* ----------------------------------------------------------------- */
L420:
    *iexp = 2;
    iz = *ibeta;
L430:
    if (k < iz) {
	goto L440;
    }
    iz *= *ibeta;
    ++(*iexp);
    goto L430;
L440:
    mx = iz + iz - 1;
/* ----------------------------------------------------------------- */
/*     loop to determine minexp, xmin */
/*     exit from loop is signaled by an underflow. */
/* ----------------------------------------------------------------- */
L450:
    *xmin = y;
    y *= betain;
/* ----------------------------------------------------------------- */
/*        check for underflow here */
/* ----------------------------------------------------------------- */
    a = y * one;
    if (a + a == zero || fabs(y) >= *xmin) {
	goto L460;
    }
    ++k;
    goto L450;
L460:
    *minexp = -k;
/* ----------------------------------------------------------------- */
/*     determine maxexp, xmax */
/* ----------------------------------------------------------------- */
    if (mx > k + k - 3 || *ibeta == 10) {
	goto L500;
    }
    mx += mx;
    ++(*iexp);
L500:
    *maxexp = mx + *minexp;
/* ----------------------------------------------------------------- */
/*     adjust for machines with implicit leading */
/*     bit in binary significand and machines with */
/*     radix point at extreme right of significand */
/* ----------------------------------------------------------------- */
    i__ = *maxexp + *minexp;
    if (*ibeta == 2 && i__ == 0) {
	--(*maxexp);
    }
    if (i__ > 20) {
	--(*maxexp);
    }
    if (a != y) {
	*maxexp += -2;
    }
    *xmax = one - *epsneg;
    if (*xmax * one != *xmax) {
	*xmax = one - beta * *epsneg;
    }
    *xmax /= beta * beta * beta * *xmin;
    i__ = *maxexp + *minexp + 3;
    if (i__ <= 0) {
	goto L520;
    }

    i__1 = i__;
    for (j = 1; j <= i__1; ++j) {
	if (*ibeta == 2) {
	    *xmax += *xmax;
	}
	if (*ibeta != 2) {
	    *xmax *= beta;
	}
/* L510: */
    }

L520:
    return;
/*     ---------- last card of dmchar ---------- */
} /* dmchar_ */

