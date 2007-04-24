/*      driver for hybrd example. */

#include <stdio.h>
#include <math.h>
#include <minpack.h>

void fcn(const int *n, const double *x, double *fvec, int *iflag);

int main()
{
  int j, n, maxfev, ml, mu, mode, nprint, info, nfev, ldfjac, lr;
  double xtol, epsfcn, factor, fnorm;
  double x[9], fvec[9], diag[9], fjac[9*9], r[45], qtf[9],
    wa1[9], wa2[9], wa3[9], wa4[9];
  int one=1;

  n = 9;

/*      the following starting values provide a rough solution. */

  for (j=1; j<=9; j++)
    {
      x[j-1] = -1.;
    }

  ldfjac = 9;
  lr = 45;

/*      set xtol to the square root of the machine precision. */
/*      unless high solutions are required, */
/*      this is the recommended setting. */

  xtol = sqrt(dpmpar_(&one));

  maxfev = 2000;
  ml = 1;
  mu = 1;
  epsfcn = 0.;
  mode = 2;
  for (j=1; j<=9; j++)
    {
      diag[j-1] = 1.;
    }

  factor = 1.e2;
  nprint = 0;

  hybrd_(&fcn, &n, x, fvec, &xtol, &maxfev, &ml, &mu, &epsfcn,
	 diag, &mode, &factor, &nprint, &info, &nfev,
	 fjac, &ldfjac, r, &lr, qtf, wa1, wa2, wa3, wa4);
  fnorm = enorm_(&n, fvec);
  printf("     final l2 norm of the residuals %15.7g\n\n", fnorm);
  printf("     number of function evaluations  %10i\n\n", nfev);
  printf("     exit parameter                  %10i\n\n", info);
  printf("     final approximate solution\n");
  for (j=1; j<=n; j++) printf("%s%15.7g", j%3==1?"\n     ":"", x[j-1]);
  printf("\n");
  return 0;
}


void fcn(const int *n, const double *x, double *fvec, int *iflag)
{
  /*      subroutine fcn for hybrd example. */

  int k;
  double one=1, temp, temp1, temp2, three=3, two=2, zero=0;

  if (iflag == 0)
    {
      /*      insert print statements here when nprint is positive. */
      return;
    }
  for (k=1; k<=*n; k++)
    {
      
      temp = (three - two*x[k-1])*x[k-1];
      temp1 = zero;
      if (k != 1) temp1 = x[k-1-1];
      temp2 = zero;
      if (k != *n) temp2 = x[k+1-1];
      fvec[k-1] = temp - temp1 - two*temp2 + one;
    }
  return;
}

