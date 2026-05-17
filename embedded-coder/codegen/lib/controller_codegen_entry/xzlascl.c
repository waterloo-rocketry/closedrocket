/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzlascl.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

/* Include Files */
#include "xzlascl.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : double cfrom
 *                double cto
 *                double A[11]
 * Return Type  : void
 */
void b_xzlascl(double cfrom, double cto, double A[11])
{
  double cfromc;
  double ctoc;
  int i;
  boolean_T notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }
    for (i = 0; i <= 8; i += 2) {
      __m128d b_r;
      b_r = _mm_loadu_pd(&A[i]);
      _mm_storeu_pd(&A[i], _mm_mul_pd(b_r, _mm_set1_pd(mul)));
    }
    A[10] *= mul;
  }
}

/*
 * Arguments    : double cfrom
 *                double cto
 *                double A[121]
 * Return Type  : void
 */
void xzlascl(double cfrom, double cto, double A[121])
{
  double cfromc;
  double ctoc;
  int i;
  boolean_T notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }
    for (i = 0; i <= 118; i += 2) {
      __m128d b_r;
      b_r = _mm_loadu_pd(&A[i]);
      b_r = _mm_mul_pd(b_r, _mm_set1_pd(mul));
      _mm_storeu_pd(&A[i], b_r);
    }
    A[120] *= mul;
  }
}

/*
 * File trailer for xzlascl.c
 *
 * [EOF]
 */
