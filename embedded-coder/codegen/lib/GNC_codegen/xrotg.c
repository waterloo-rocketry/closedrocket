/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xrotg.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

/* Include Files */
#include "xrotg.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : double *a
 *                double *b
 *                double *s
 * Return Type  : double
 */
double xrotg(double *a, double *b, double *s)
{
  double absa;
  double absb;
  double b_c;
  double roe;
  double scale;
  roe = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *a;
  }
  scale = absa + absb;
  if (scale == 0.0) {
    scale = 0.0;
    b_c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    double bds;
    b_c = absa / scale;
    bds = absb / scale;
    bds = scale * sqrt(b_c * b_c + bds * bds);
    if (roe < 0.0) {
      bds = -bds;
    }
    b_c = *a / bds;
    scale = *b / bds;
    if (absa > absb) {
      *b = scale;
    } else if (b_c != 0.0) {
      *b = 1.0 / b_c;
    } else {
      *b = 1.0;
    }
    *a = bds;
  }
  *s = scale;
  return b_c;
}

/*
 * File trailer for xrotg.c
 *
 * [EOF]
 */
