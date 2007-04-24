

/*     driver for lmdif1 example. */

#include <stdio.h>
#include <math.h>
#include <minpack.h>

void fcn(const int *m, const int *n, const double *x, double *fvec, int *iflag);

int main()
{
  int m, n, info, lwa, iwa[3], one=1;
  double tol, fnorm, x[3], fvec[15], wa[75];

  m = 15;
  n = 3;

  /* the following starting values provide a rough fit. */

  x[0] = 1.e0;
  x[1] = 1.e0;
  x[2] = 1.e0;

  lwa = 75;

  /* set tol to the square root of the machine precision.  unless high
     precision solutions are required, this is the recommended
     setting. */

  tol = sqrt(dpmpar_(&one));

  lmdif1_(&fcn, &m, &n, x, fvec, &tol, &info, iwa, wa, &lwa);

  fnorm = enorm_(&m, fvec);

  printf("      FINAL L2 NORM OF THE RESIDUALS%15.7f\n\n",fnorm);
  printf("      EXIT PARAMETER                %10i\n\n", info);
  printf("      FINAL APPROXIMATE SOLUTION\n\n %15.7f%15.7f%15.7f\n",
	 x[0], x[1], x[2]);
  return 0;
}

void fcn(const int *m, const int *n, const double *x, double *fvec, int *iflag)
{
  /* function fcn for lmdif1 example */

  int i;
  double tmp1,tmp2,tmp3;
  double y[15]={1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
		3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0};

  for (i=0; i<15; i++)
    {
      tmp1 = i+1;
      tmp2 = 15 - i;
      tmp3 = tmp1;
      
      if (i >= 8) tmp3 = tmp2;
      fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
    }
}
