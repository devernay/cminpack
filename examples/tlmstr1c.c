/*    driver for lmstr1 example. */

#include <stdio.h>
#include <math.h>
#include <cminpack.h>

/* the following struct defines the data points */
typedef struct  {
    int m;
    double *y;
} fcndata_t;

int  fcn(void *p, int m, int n, const double *x, double *fvec, double *fjrow, int iflag);

int main()
{
  int ldfjac, info, lwa, ipvt[3];
  double tol, fnorm;
  double x[3], fvec[15], fjac[9], wa[30];
  const int m = 15;
  const int n = 3;
  /* auxiliary data (e.g. measurements) */
  double y[15] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
                  3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};
  fcndata_t data;
  data.m = m;
  data.y = y;

  /*     the following starting values provide a rough fit. */

  x[0] = 1.;
  x[1] = 1.;
  x[2] = 1.;

  ldfjac = 3;
  lwa = 30;

  /*     set tol to the square root of the machine precision.
     unless high precision solutions are required,
     this is the recommended setting. */

  tol = sqrt(dpmpar(1));

  info = lmstr1(fcn, &data, m, n, 
	  x, fvec, fjac, ldfjac, 
	  tol, ipvt, wa, lwa);

  fnorm = enorm(m, fvec);

  printf("      FINAL L2 NORM OF THE RESIDUALS%15.7g\n\n", fnorm);
  printf("      EXIT PARAMETER                %10i\n\n", info);
  printf("      FINAL APPROXIMATE SOLUTION\n\n%15.7g%15.7g%15.7g\n",
	 x[0], x[1], x[2]);

  return 0;
}

int  fcn(void *p, int m, int n, const double *x, double *fvec, double *fjrow, int iflag)
{
  /*  subroutine fcn for lmstr1 example. */
  int i;
  double tmp1, tmp2, tmp3, tmp4;
  const double *y = ((fcndata_t*)p)->y;

  if (iflag < 2)
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
      i = iflag - 2;
      tmp1 = i + 1;
      tmp2 = 15 - i;
      tmp3 = (i > 7) ? tmp2 : tmp1;
      tmp4 = (x[1]*tmp2 + x[2]*tmp3); tmp4 = tmp4*tmp4;
      fjrow[0] = -1.;
      fjrow[1] = tmp1*tmp2/tmp4;
      fjrow[2] = tmp1*tmp3/tmp4;
    }
  return 0;
}
