/*      driver for fdjac2 example. */
/*      The test works by running chkder both on the Jacobian computed
        by forward-differences and on the real Jacobian */

#include <stdio.h>
#include <math.h>
#include <cminpack.h>

int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag);
void fcnjac(int m, int n, const double *x, double *fjac, int ldfjac);

int main()
{
  int i, m, n, ldfjac;
  double epsfcn;
  double x[3], fvec[15], fjac[15*3], fdjac[15*3], xp[3], fvecp[15], 
      err[15], errd[15], wa[15];

  m = 15;
  n = 3;

  /*      the following values should be suitable for */
  /*      checking the jacobian matrix. */

  x[1-1] = 9.2e-1;
  x[2-1] = 1.3e-1;
  x[3-1] = 5.4e-1;

  ldfjac = 15;

  epsfcn = 0.;

  chkder(m, n, x, fvec, fjac, ldfjac, xp, fvecp, 1, err);
  fcn(0, m, n, x, fvec, 1);
  fdjac2(fcn, 0, m, n, x, fvec, fdjac, ldfjac, epsfcn, wa);
  fcnjac(m, n, x, fjac, ldfjac);
  fcn(0, m, n, xp, fvecp, 1);
  chkder(m, n, x, fvec, fdjac, ldfjac, xp, fvecp, 2, errd);
  chkder(m, n, x, fvec, fjac, ldfjac, xp, fvecp, 2, err);

  for (i=1; i<=m; i++)
    {
      fvecp[i-1] = fvecp[i-1] - fvec[i-1];
    }
  printf("\n      fvec\n");  
  for (i=1; i<=m; i++) printf("%s%15.7g",i%3==1?"\n     ":"", fvec[i-1]);
  printf("\n      fvecp - fvec\n");  
  for (i=1; i<=m; i++) printf("%s%15.7g",i%3==1?"\n     ":"", fvecp[i-1]);
  printf("\n      errd\n");  
  for (i=1; i<=m; i++) printf("%s%15.7g",i%3==1?"\n     ":"", errd[i-1]);
  printf("\n      err\n");  
  for (i=1; i<=m; i++) printf("%s%15.7g",i%3==1?"\n     ":"", err[i-1]);
  printf("\n");
  return 0;
}

int fcn(void *p, int m, int n, const double *x, double *fvec, int iflag)
{

/*      subroutine fcn for fdjac2 example. */

  int i;
  double tmp1, tmp2, tmp3;
  double y[15]={1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
		3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};

  if (iflag == 0)
    {
      /*      insert print statements here when nprint is positive. */
      return 0;
    }
  for (i = 1; i <= 15; i++)
    {
      tmp1 = i;
      tmp2 = 16 - i;
      tmp3 = tmp1;
      if (i > 8) tmp3 = tmp2;
      fvec[i-1] = y[i-1] - (x[1-1] + tmp1/(x[2-1]*tmp2 + x[3-1]*tmp3));
    }
  return 0;
}

void fcnjac(int m, int n, const double *x,
            double *fjac, int ldfjac)
{
  /*      Jacobian of fcn (corrected version from tchkder). */

  int i;
  double tmp1, tmp2, tmp3, tmp4;

  for (i = 1; i <= 15; i++)
    {
      tmp1 = i;
      tmp2 = 16 - i;
	  
      tmp3 = tmp1;
      if (i > 8) tmp3 = tmp2;
      tmp4 = (x[2-1]*tmp2 + x[3-1]*tmp3); tmp4=tmp4*tmp4;
      fjac[i-1+ ldfjac*(1-1)] = -1.;
      fjac[i-1+ ldfjac*(2-1)] = tmp1*tmp2/tmp4;
      fjac[i-1+ ldfjac*(3-1)] = tmp1*tmp3/tmp4;
    }
}
