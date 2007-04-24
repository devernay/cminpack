/*      driver for lmder1 example. */


#include <stdio.h>
#include <math.h>
#include <cminpack.h>

int fcn(void *p, int m, int n, const double *x, double *fvec, double *fjac, 
	 int ldfjac, int iflag);

int main()
{
  int j, m, n, ldfjac, info, lwa;
  int ipvt[3];
  double tol, fnorm;
  double x[3], fvec[15], fjac[15*3], wa[30];

  m = 15;
  n = 3;

/*      the following starting values provide a rough fit. */

  x[1-1] = 1.;
  x[2-1] = 1.;
  x[3-1] = 1.;

  ldfjac = 15;
  lwa = 30;

/*      set tol to the square root of the machine precision. */
/*      unless high solutions are required, */
/*      this is the recommended setting. */

  tol = sqrt(dpmpar(1));

  info = lmder1(fcn, 0, m, n, x, fvec, fjac, ldfjac, tol, 
	  ipvt, wa, lwa);
  fnorm = enorm(m, fvec);
  printf("      final l2 norm of the residuals%15.7g\n\n", fnorm);
  printf("      exit parameter                %10i\n\n", info);
  printf("      final approximate solution\n");
  for (j=1; j<=n; j++) printf("%s%15.7g", j%3==1?"\n     ":"", x[j-1]);
  printf("\n");

  return 0;
}

int fcn(void *p, int m, int n, const double *x, double *fvec, double *fjac, 
	 int ldfjac, int iflag)
{

/*      subroutine fcn for lmder1 example. */

  int i;
  double tmp1, tmp2, tmp3, tmp4;
  double y[15] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
		  3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};

  if (iflag != 2)
    {
      for (i = 1; i <= 15; i++)
	{
	  tmp1 = i;
	  tmp2 = 16 - i;
	  tmp3 = tmp1;
	  if (i > 8) tmp3 = tmp2;
	  fvec[i-1] = y[i-1] - (x[1-1] + tmp1/(x[2-1]*tmp2 + x[3-1]*tmp3));
	}
    }
  else
    {
      for ( i = 1; i <= 15; i++)
	{
	  tmp1 = i;
	  tmp2 = 16 - i;
	  tmp3 = tmp1;
	  if (i > 8) tmp3 = tmp2;
	  tmp4 = (x[2-1]*tmp2 + x[3-1]*tmp3); tmp4 = tmp4*tmp4;
	  fjac[i-1 + ldfjac*(1-1)] = -1.;
	  fjac[i-1 + ldfjac*(2-1)] = tmp1*tmp2/tmp4;
	  fjac[i-1 + ldfjac*(3-1)] = tmp1*tmp3/tmp4;
	}
    }
  return 0;
}
