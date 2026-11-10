/*      driver for lmder with a parallel finite-difference Jacobian. */

/* lmdif approximates the Jacobian by forward differences, calling fcn n times
   one after the other. When fcn is expensive and thread-safe, these n calls
   can run in parallel: fcn_fdjac below turns an lmdif function fcn into an
   lmder function that computes the same forward differences, one column per
   thread, and lmder is called instead of lmdif. This example does it with
   OpenMP (compile with -fopenmp or equivalent; without OpenMP the loop simply
   runs serially). The application decides at runtime whether to use threads,
   and how many (e.g. with OMP_NUM_THREADS).

   The problem is the one of tlmdifc.c. The program checks that lmder finds
   the same solution as lmdif. */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <cminpack.h>
#define real __cminpack_real__

/* the following struct defines the data points */
typedef struct  {
    int m;
    real *y;
} fcndata_t;

int fcn(void *p, int m, int n, const real *x, real *fvec, int iflag);
int fcn_fdjac(void *p, int m, int n, const real *x, real *fvec, real *fjac,
              int ldfjac, int iflag);

int main()
{
  int j, ldfjac, maxfev, mode, nprint, info, nfev, njev, info_lmdif, nfev_lmdif;
  int ipvt[3];
  real ftol, xtol, gtol, factor, fnorm;
  real x[3], x_lmdif[3], fvec[15], fjac[15*3], diag[3], qtf[3],
    wa1[3], wa2[3], wa3[3], wa4[15];
  int same;
  const int m = 15;
  const int n = 3;
  /* auxiliary data (e.g. measurements) */
  real y[15] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
                  3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};
  fcndata_t data;
  data.m = m;
  data.y = y;

/*      the following starting values provide a rough fit. */

  x[0] = 1.;
  x[1] = 1.;
  x[2] = 1.;
  memcpy(x_lmdif, x, sizeof(x));

  ldfjac = 15;

  /*      set ftol and xtol to the square root of the machine */
  /*      and gtol to zero. unless high solutions are */
  /*      required, these are the recommended settings. */

  ftol = sqrt(__cminpack_func__(dpmpar)(1));
  xtol = sqrt(__cminpack_func__(dpmpar)(1));
  gtol = 0.;

  maxfev = 400;
  mode = 1;
  factor = 1.e2;
  nprint = 0;

  info = __cminpack_func__(lmder)(fcn_fdjac, &data, m, n, x, fvec, fjac, ldfjac,
        ftol, xtol, gtol, maxfev, diag, mode, factor, nprint, &nfev, &njev,
        ipvt, qtf, wa1, wa2, wa3, wa4);
  fnorm = __cminpack_func__(enorm)(m, fvec);
  printf("      final l2 norm of the residuals%15.7g\n\n", (double)fnorm);
  printf("      number of function evaluations%10i\n\n", nfev);
  printf("      number of Jacobian evaluations%10i\n\n", njev);
  printf("      exit parameter                %10i\n\n", info);
  printf("      final approximate solution\n");
  for (j=0; j<n; ++j) {
    printf("%s%15.7g", j%3==0?"\n     ":"", (double)x[j]);
  }
  printf("\n\n");

  /*      solve the same problem with lmdif, which computes the same */
  /*      forward differences serially, and compare the solutions. */
  /*      lmdif also counts the n evaluations of each Jacobian in */
  /*      maxfev, so it is given more of them. Both stop when the */
  /*      relative step is below xtol, so their solutions agree to */
  /*      a few times xtol (exactly, if they round the same way). */

  info_lmdif = __cminpack_func__(lmdif)(fcn, &data, m, n, x_lmdif, fvec,
        ftol, xtol, gtol, 2*maxfev, 0., diag, mode, factor, nprint, &nfev_lmdif,
        fjac, ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4);
  same = (info_lmdif == info);
  for (j=0; j<n; ++j) {
    same = same &&
      fabs((double)(x_lmdif[j] - x[j])) <= 10 * xtol * fabs((double)x[j]);
  }
  printf("      lmdif: exit parameter %d, %s\n", info_lmdif,
         same ? "same solution" : "DIFFERENT SOLUTION");
  return same ? 0 : 1;
}

int fcn(void *p, int m, int n, const real *x, real *fvec, int iflag)
{

/*      subroutine fcn for lmdif example. it only reads the shared */
/*      data, so it can be called by several threads at the same time. */

  int i;
  real tmp1, tmp2, tmp3;
  const real *y = ((fcndata_t*)p)->y;
  assert(m == 15 && n == 3);

  if (iflag == 0) {
    /*      insert print statements here when nprint is positive. */
    return 0;
  }

  /* compute residuals */
  for (i = 0; i < 15; ++i) {
    tmp1 = i + 1;
    tmp2 = 15 - i;
    tmp3 = (i > 7) ? tmp2 : tmp1;
    fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
  }
  return 0;
}

int fcn_fdjac(void *p, int m, int n, const real *x, real *fvec, real *fjac,
              int ldfjac, int iflag)
{

/*      lmder function made from the lmdif function fcn: with iflag = 2, */
/*      compute the Jacobian by forward differences, as lmdif does, */
/*      one column per thread. any other iflag is passed on to fcn. */

  /* the step of lmdif, whose epsfcn argument is 0 here */
  const real eps = sqrt(__cminpack_func__(dpmpar)(1));
  int j, ret = 0;

  if (iflag != 2) {
    return fcn(p, m, n, x, fvec, iflag);
  }
  assert(m <= 15 && n <= 3);
  /* fvec holds the residuals at x */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (j = 0; j < n; ++j) {
    /* each thread perturbs its own copy of x */
    real xh[3], fvech[15], h;
    int i, r;

    h = eps * fabs((double)x[j]);
    if (h == 0.) {
      h = eps;
    }
    memcpy(xh, x, n * sizeof(real));
    xh[j] = x[j] + h;
    /* lmdif also calls fcn with iflag = 2 for the Jacobian */
    r = fcn(p, m, n, xh, fvech, 2);
    if (r < 0) {
      /* fcn asks to stop: lmder will return this value */
#ifdef _OPENMP
#pragma omp critical
#endif
      ret = r;
    }
    for (i = 0; i < m; ++i) {
      fjac[i + j * ldfjac] = (fvech[i] - fvec[i]) / h;
    }
  }
  return ret;
}
