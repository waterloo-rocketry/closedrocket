/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xnrm2.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

/* Include Files */
#include "xnrm2.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : int n
 *                const double x[11]
 *                int ix0
 * Return Type  : double
 */
double b_xnrm2(int n, const double x[11], int ix0)
{
  double scale;
  double y;
  int k;
  int kend;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = ix0 + n;
  for (k = ix0; k < kend; k++) {
    double absxk;
    absxk = fabs(x[k - 1]);
    if (absxk > scale) {
      double b_t;
      b_t = scale / absxk;
      y = y * b_t * b_t + 1.0;
      scale = absxk;
    } else {
      double b_t;
      b_t = absxk / scale;
      y += b_t * b_t;
    }
  }
  return scale * sqrt(y);
}

/*
 * Arguments    : int n
 *                const double x[121]
 *                int ix0
 * Return Type  : double
 */
double xnrm2(int n, const double x[121], int ix0)
{
  double scale;
  double y;
  int k;
  int kend;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = ix0 + n;
  for (k = ix0; k < kend; k++) {
    double absxk;
    absxk = fabs(x[k - 1]);
    if (absxk > scale) {
      double b_t;
      b_t = scale / absxk;
      y = y * b_t * b_t + 1.0;
      scale = absxk;
    } else {
      double b_t;
      b_t = absxk / scale;
      y += b_t * b_t;
    }
  }
  return scale * sqrt(y);
}

/*
 * File trailer for xnrm2.c
 *
 * [EOF]
 */
