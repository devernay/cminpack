#include <stdio.h>
#include <math.h>
#include "cminpack.h"
#define real __cminpack_real__

extern void dmchar_(int *ibeta, int *it, int *irnd, 
                    int *ngrd, int *machep, int *negep, int *iexp, 
                    int *minexp, int *maxexp, real *eps, real *epsneg,
                    real *xmin, real *xmax);

/*     ********** */

/*     this program checks the constants of machine precision and */
/*     smallest and largest machine representable numbers specified in */
/*     function dpmpar, against the corresponding hardware-determined */
/*     machine constants obtained by dmchar, a subroutine due to */
/*     w. j. cody. */

/*     data statements in dpmpar corresponding to the machine used must */
/*     be activated by removing c in column 1. */

/*     the printed output consists of the machine constants obtained by */
/*     dmchar and comparisons of the dpmpar constants with their */
/*     dmchar counterparts. descriptions of the machine constants are */
/*     given in the prologue comments of dmchar. */

/*     subprograms called */

/*       minpack-supplied ... dmchar,dpmpar */

/*     argonne national laboratory. minpack project. march 1980. */
/*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

/*     ********** */
/* Main program */
int main(int argc, char **argv)
{
    /* Local variables */
    int it;
    real eps;
    int ngrd, irnd, iexp;
    real rerr[3], xmin, xmax;
    int ibeta, negep;
    real giant, dwarf;
    int machep;
    real epsmch, epsneg;
    int minexp, maxexp;

/*     determine the machine constants dynamically from dmchar. */

    dmchar_(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &
	    maxexp, &eps, &epsneg, &xmin, &xmax);

/*     compare the dpmpar constants with their dmchar counterparts and */
/*     store the relative differences in rerr. */

    epsmch = dpmpar(1);
    dwarf = dpmpar(2);
    giant = dpmpar(3);
    rerr[0] = (epsmch - eps) / epsmch;
    rerr[1] = (dwarf - xmin) / dwarf;
    rerr[2] = (xmax - giant) / giant;

/*     write the dmchar constants. */

    printf("\f dmchar constants\n\n\n");
    printf(" ibeta =%6i\n\n", ibeta);
    printf(" it    =%6i\n\n", it);
    printf(" irnd  =%6i\n\n", irnd);
    printf(" ngrd  =%6i\n\n", ngrd);
    printf(" machep =%6i\n\n", machep);
    printf(" negep =%6i\n\n", negep);
    printf(" iexp =%6i\n\n", iexp);
    printf(" minexp =%6i\n\n", minexp);
    printf(" maxexp =%6i\n\n", maxexp);
    printf(" eps =%15.7e\n\n", (double)eps);
    printf(" epsneg =%15.7e\n\n", (double)epsneg);
    printf(" xmin =%15.7e\n\n", (double)xmin);
    printf(" xmax =%15.7e\n\n", (double)xmax);


/*     write the dpmpar constants and the relative differences. */

    printf("\n\n dpmpar constants and relative differences\n\n\n");
    printf(" epsmch =%15.7e\n", (double)epsmch);
    printf(" rerr(1) =%15.7e\n\n", (double)rerr[0]);
    printf(" dwarf =%15.7e\n", (double)dwarf);
    printf(" rerr(2) =%15.7e\n\n", (double)rerr[1]);
    printf(" giant =%15.7e\n", (double)giant);
    printf(" rerr(3) =%15.7e\n", (double)rerr[2]);

    return 0;
} /* MAIN__ */

