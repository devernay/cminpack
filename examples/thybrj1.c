/*      driver for hybrj1 example. */

#include <stdio.h>
#include <math.h>
#include <cminpack.h>

int fcn(void *p, int n, const double *x, double *fvec, double *fjac, int ldfjac, 
	 int iflag);

int main()
{
  int j, n, ldfjac, info, lwa;
  double tol, fnorm;
  double x[9], fvec[9], fjac[9*9], wa[99];

  n = 9;

/*      the following starting values provide a rough solution. */

  for (j=1; j<=9; j++)
    {
      x[j-1] = -1.;
    }

  ldfjac = 9;
  lwa = 99;

/*      set tol to the square root of the machine precision. */
/*      unless high solutions are required, */
/*      this is the recommended setting. */

  tol = sqrt(dpmpar(1));

  info = hybrj1(fcn, 0, n, x, fvec, fjac, ldfjac, tol, wa, lwa);

  fnorm = enorm(n, fvec);

  printf("      final l2 norm of the residuals%15.7g\n\n", fnorm);
  printf("      exit parameter                %10i\n\n", info);
  printf("      final approximate solution\n");
  for (j=1; j<=n; j++) printf("%s%15.7g", j%3==1?"\n     ":"", x[j-1]);
  printf("\n");

  return 0;
}

int fcn(void *p, int n, const double *x, double *fvec, double *fjac, int ldfjac, 
	 int iflag)
{
  /*      subroutine fcn for hybrj1 example. */

  int j, k;
  double one=1, temp, temp1, temp2, three=3, two=2, zero=0, four=4;

  if (iflag != 2)
    {
      for (k = 1; k <= n; k++)
	{
	  temp = (three - two*x[k-1])*x[k-1];
	  temp1 = zero;
	  if (k != 1) temp1 = x[k-1-1];
	  temp2 = zero;
	  if (k != n) temp2 = x[k+1-1];
	  fvec[k-1] = temp - temp1 - two*temp2 + one;
	}
    }
  else
    {
     for (k = 1; k <= n; k++)
       {
	 for (j = 1; j <= n; j++)
	   {
	     fjac[k-1 + ldfjac*(j-1)] = zero;
	   }
         fjac[k-1 + ldfjac*(k-1)] = three - four*x[k-1];
         if (k != 1) fjac[k-1 + ldfjac*(k-1-1)] = -one;
         if (k != n) fjac[k-1 + ldfjac*(k+1-1)] = -two;
       }
    }
  return 0;
}
