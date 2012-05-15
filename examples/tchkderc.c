/*      driver for chkder example. */

#include <stdio.h>
#include <math.h>
#include <cminpack.h>
#define real __cminpack_real__

/* the following struct defines the data points */
typedef struct  {
    int m;
    real *y;
} fcndata_t;

int fcn(void *p, int m, int n, const real *x, real *fvec,
	 real *fjac, int ldfjac, int iflag);

int main()
{
  int i, ldfjac;
  real x[3], fvec[15], fjac[15*3], xp[3], fvecp[15], 
    err[15];
  const int m = 15;
  const int n = 3;
  /* auxiliary data (e.g. measurements) */
  real y[15] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
                  3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};
  fcndata_t data;
  data.m = m;
  data.y = y;

  /*      the following values should be suitable for */
  /*      checking the jacobian matrix. */

  x[0] = 9.2e-1;
  x[1] = 1.3e-1;
  x[2] = 5.4e-1;

  ldfjac = 15;

  /* compute xp from x */
  __cminpack_func__(chkder)(m, n, x, NULL, NULL, ldfjac, xp, NULL, 1, NULL);
  /* compute fvec at x (all components of fvec should be != 0).*/
  fcn(&data, m, n, x, fvec, NULL, ldfjac, 1);
  /* compute fjac at x */
  fcn(&data, m, n, x, NULL, fjac, ldfjac, 2);
  /* compute fvecp at xp (all components of fvecp should be != 0)*/
  fcn(&data, m, n, xp, fvecp, NULL, ldfjac, 1);
  /* check Jacobian, put the result in err */
  __cminpack_func__(chkder)(m, n, x, fvec, fjac, ldfjac, NULL, fvecp, 2, err);
  /* Output values:
     err[i] = 1.: i-th gradient is correct
     err[i] = 0.: i-th gradient is incorrect
     err[I] > 0.5: i-th gradient is probably correct
  */

  for (i=0; i<m; ++i)
    {
      fvecp[i] = fvecp[i] - fvec[i];
    }
  printf("\n      fvec\n");  
  for (i=0; i<m; ++i) printf("%s%15.7g",i%3==0?"\n     ":"", (double)fvec[i]);
  printf("\n      fvecp - fvec\n");  
  for (i=0; i<m; ++i) printf("%s%15.7g",i%3==0?"\n     ":"", (double)fvecp[i]);
  printf("\n      err\n");  
  for (i=0; i<m; ++i) printf("%s%15.7g",i%3==0?"\n     ":"", (double)err[i]);
  printf("\n");
  return 0;
}

int fcn(void *p, int m, int n, const real *x, real *fvec,
	 real *fjac, int ldfjac, int iflag)
{
  /*      subroutine fcn for chkder example. */

  int i;
  real tmp1, tmp2, tmp3, tmp4;
  const real *y = ((fcndata_t*)p)->y;

  if (iflag == 0) 
    {
      /*      insert print statements here when nprint is positive. */
      return 0;
    }

  if (iflag != 2) 
    {
      for (i=0; i < 15; ++i)
	{
	  tmp1 = i + 1;
	  tmp2 = 15 - i;
	  tmp3 = (i > 7) ? tmp2 : tmp1;
	  fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
	}
    }
  else
    {
      for (i=0; i<15; ++i)
	{
	  tmp1 = i + 1;
	  tmp2 = 15 - i;
	  /* error introduced into next statement for illustration. */
	  /* corrected statement should read    tmp3 = (i > 7) ? tmp2 : tmp1 . */
	  tmp3 = (i > 7) ? tmp2 : tmp2;
	  tmp4 = (x[1]*tmp2 + x[2]*tmp3); tmp4 = tmp4*tmp4;
	  fjac[i + ldfjac*0] = -1.;
	  fjac[i + ldfjac*1] = tmp1*tmp2/tmp4;
	  fjac[i + ldfjac*2] = tmp1*tmp3/tmp4;
	};
    }
  return 0;
}
