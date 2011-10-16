/*      driver for chkder example. */

#include <stdio.h>
#include <math.h>
#include <cminpack.h>

int fcn(void *p, int m, int n, const double *x, double *fvec,
	 double *fjac, int ldfjac, int iflag);

int main()
{
  int i, m, n, ldfjac;
  double x[3], fvec[15], fjac[15*3], xp[3], fvecp[15], 
    err[15];

  m = 15;
  n = 3;

  /*      the following values should be suitable for */
  /*      checking the jacobian matrix. */

  x[1-1] = 9.2e-1;
  x[2-1] = 1.3e-1;
  x[3-1] = 5.4e-1;

  ldfjac = 15;

  /* compute xp from x */
  chkder(m, n, x, NULL, NULL, ldfjac, xp, NULL, 1, NULL);
  /* compute fvec at x (all components of fvec should be != 0).*/
  fcn(0, m, n, x, fvec, NULL, ldfjac, 1);
  /* compute fjac at x */
  fcn(0, m, n, x, NULL, fjac, ldfjac, 2);
  /* compute fvecp at xp (all components of fvecp should be != 0)*/
  fcn(0, m, n, xp, fvecp, NULL, ldfjac, 1);
  /* check Jacobian, put the result in err */
  chkder(m, n, x, fvec, fjac, ldfjac, NULL, fvecp, 2, err);
  /* Output values:
     err[i] = 1.: i-th gradient is correct
     err[i] = 0.: i-th gradient is incorrect
     err[I] > 0.5: i-th gradient is probably correct
  */

  for (i=1; i<=m; i++)
    {
      fvecp[i-1] = fvecp[i-1] - fvec[i-1];
    }
  printf("\n      fvec\n");  
  for (i=1; i<=m; i++) printf("%s%15.7g",i%3==1?"\n     ":"", fvec[i-1]);
  printf("\n      fvecp - fvec\n");  
  for (i=1; i<=m; i++) printf("%s%15.7g",i%3==1?"\n     ":"", fvecp[i-1]);
  printf("\n      err\n");  
  for (i=1; i<=m; i++) printf("%s%15.7g",i%3==1?"\n     ":"", err[i-1]);
  printf("\n");
  return 0;
}

int fcn(void *p, int m, int n, const double *x, double *fvec,
	 double *fjac, int ldfjac, int iflag)
{
  /*      subroutine fcn for chkder example. */

  int i;
  double tmp1, tmp2, tmp3, tmp4;
  double y[15]={1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
		3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};

  
  if (iflag == 0) 
    {
      /*      insert print statements here when nprint is positive. */
      return 0;
    }

  if (iflag != 2) 

    for (i=1; i<=15; i++)
      {
	tmp1 = i;
	tmp2 = 16 - i;
	tmp3 = tmp1;
	if (i > 8) tmp3 = tmp2;
	fvec[i-1] = y[i-1] - (x[1-1] + tmp1/(x[2-1]*tmp2 + x[3-1]*tmp3));
      }
  else
    {
      for (i = 1; i <= 15; i++)
	{
	  tmp1 = i;
	  tmp2 = 16 - i;
	  
	  /* error introduced into next statement for illustration. */
	  /* corrected statement should read    tmp3 = tmp1 . */
	  
	  tmp3 = tmp2;
	  if (i > 8) tmp3 = tmp2;
	  tmp4 = (x[2-1]*tmp2 + x[3-1]*tmp3); tmp4=tmp4*tmp4;
	  fjac[i-1+ ldfjac*(1-1)] = -1.;
	  fjac[i-1+ ldfjac*(2-1)] = tmp1*tmp2/tmp4;
	  fjac[i-1+ ldfjac*(3-1)] = tmp1*tmp3/tmp4;
	}
    }
  return 0;
}
