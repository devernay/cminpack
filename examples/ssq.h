#ifndef __CMINPACK_SSQ_H__
#define __CMINPACK_SSQ_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void initpt(int n, double *x, int nprob, double factor);

void ssqfcn(int m, int n, const double *x, double *fvec, int nprob);

void ssqjac(int m, int n, const double *x, double *fjac, int ldfjac, int nprob);

#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif /* __CMINPACK_SSQ_H__ */
