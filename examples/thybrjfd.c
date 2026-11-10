/*      driver for hybrj with a parallel finite-difference Jacobian. */

/* hybrd approximates the Jacobian by forward differences, calling fcn n times
   one after the other. When fcn is expensive and thread-safe, these n calls
   can run in parallel: fcn_fdjac below turns a hybrd function fcn into a
   hybrj function that computes the same forward differences, one column per
   thread, and hybrj is called instead of hybrd. This example does it with
   OpenMP (compile with -fopenmp or equivalent; without OpenMP the loop simply
   runs serially). The application decides at runtime whether to use threads,
   and how many (e.g. with OMP_NUM_THREADS). See tlmderfd.c for the same with
   lmdif and lmder.

   The problem is the one of thybrdc.c. Its Jacobian is tridiagonal, and with
   ml = mu = 1 hybrd itself would need only 3 evaluations of fcn per Jacobian
   instead of n (see thybrdc.c): for a banded Jacobian, use that first. Here
   the Jacobian is treated as dense (ml = mu = n-1), as for most problems.
   The program checks that hybrj finds the same solution as hybrd. */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <cminpack.h>
#define real __cminpack_real__

int fcn(void *p, int n, const real *x, real *fvec, int iflag);
int fcn_fdjac(void *p, int n, const real *x, real *fvec, real *fjac,
              int ldfjac, int iflag);

int main()
{
  int j, n, ldfjac, maxfev, mode, nprint, info, nfev, njev, lr;
  int info_hybrd, nfev_hybrd, same;
  real xtol, factor, fnorm;
  real x[9], x_hybrd[9], fvec[9], fjac[9*9], diag[9], r[45], qtf[9],
    wa1[9], wa2[9], wa3[9], wa4[9];

  n = 9;

/*      the following starting values provide a rough solution. */

  for (j=1; j<=9; j++) {
    x[j-1] = -1.;
  }
  memcpy(x_hybrd, x, sizeof(x));

  ldfjac = 9;
  lr = 45;

/*      set xtol to the square root of the machine precision. */
/*      unless high solutions are required, */
/*      this is the recommended setting. */

  xtol = sqrt(__cminpack_func__(dpmpar)(1));

  maxfev = 1000;
  mode = 2;
  for (j=1; j<=9; j++) {
    diag[j-1] = 1.;
  }
  factor = 1.e2;
  nprint = 0;

  info = __cminpack_func__(hybrj)(fcn_fdjac, NULL, n, x, fvec, fjac, ldfjac,
        xtol, maxfev, diag, mode, factor, nprint, &nfev, &njev, r, lr, qtf,
        wa1, wa2, wa3, wa4);
  fnorm = __cminpack_func__(enorm)(n, fvec);
  printf("     final l2 norm of the residuals%15.7g\n\n", (double)fnorm);
  printf("     number of function evaluations%10i\n\n", nfev);
  printf("     number of jacobian evaluations%10i\n\n", njev);
  printf("     exit parameter                %10i\n\n", info);
  printf("     final approximate solution\n\n");
  for (j=1; j<=n; j++) {
    printf("%s%15.7g", j%3==1?"\n     ":"", (double)x[j-1]);
  }
  printf("\n\n");

  /*      solve the same problem with hybrd, which computes the same */
  /*      forward differences serially (dense Jacobian: ml = mu = n-1), */
  /*      and compare the solutions. hybrd also counts the n */
  /*      evaluations of each Jacobian in maxfev, so it is given more */
  /*      of them. Both stop when the relative step is below xtol, so */
  /*      their solutions agree to a few times xtol (exactly, if they */
  /*      round the same way). */

  for (j=1; j<=9; j++) {
    diag[j-1] = 1.;
  }
  info_hybrd = __cminpack_func__(hybrd)(fcn, NULL, n, x_hybrd, fvec, xtol,
        2*maxfev, n-1, n-1, 0., diag, mode, factor, nprint, &nfev_hybrd,
        fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4);
  same = (info_hybrd == info);
  for (j=0; j<n; ++j) {
    same = same &&
      fabs((double)(x_hybrd[j] - x[j])) <= 10 * xtol * fabs((double)x[j]);
  }
  printf("     hybrd: exit parameter %d, %s\n", info_hybrd,
         same ? "same solution" : "DIFFERENT SOLUTION");
  return same ? 0 : 1;
}

int fcn(void *p, int n, const real *x, real *fvec, int iflag)
{

/*      subroutine fcn for hybrd example. it does not modify any shared */
/*      data, so it can be called by several threads at the same time. */

  int k;
  real temp, temp1, temp2;
  (void)p;
  assert(n == 9);

  if (iflag == 0) {
    /*      insert print statements here when nprint is positive. */
    return 0;
  }

  /* compute residuals */
  for (k=1; k<=n; k++) {
    temp = (3 - 2*x[k-1])*x[k-1];
    temp1 = 0;
    if (k != 1) {
      temp1 = x[k-1-1];
    }
    temp2 = 0;
    if (k != n) {
      temp2 = x[k+1-1];
    }
    fvec[k-1] = temp - temp1 - 2*temp2 + 1;
  }
  return 0;
}

int fcn_fdjac(void *p, int n, const real *x, real *fvec, real *fjac,
              int ldfjac, int iflag)
{

/*      hybrj function made from the hybrd function fcn: with iflag = 2, */
/*      compute the Jacobian by forward differences, as hybrd does, */
/*      one column per thread. any other iflag is passed on to fcn. */

  /* the step of hybrd, whose epsfcn argument is 0 here */
  const real eps = sqrt(__cminpack_func__(dpmpar)(1));
  int j, ret = 0;

  if (iflag != 2) {
    return fcn(p, n, x, fvec, iflag);
  }
  assert(n <= 9);
  /* fvec holds the residuals at x */
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (j = 0; j < n; ++j) {
    /* each thread perturbs its own copy of x */
    real xh[9], fvech[9], h;
    int i, r;

    h = eps * fabs((double)x[j]);
    if (h == 0.) {
      h = eps;
    }
    memcpy(xh, x, n * sizeof(real));
    xh[j] = x[j] + h;
    /* hybrd also calls fcn with iflag = 2 for the Jacobian */
    r = fcn(p, n, xh, fvech, 2);
    if (r < 0) {
      /* fcn asks to stop: hybrj will return this value */
#ifdef _OPENMP
#pragma omp critical
#endif
      ret = r;
    }
    for (i = 0; i < n; ++i) {
      fjac[i + j * ldfjac] = (fvech[i] - fvec[i]) / h;
    }
  }
  return ret;
}
