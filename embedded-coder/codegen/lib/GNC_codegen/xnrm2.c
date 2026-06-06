/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xnrm2.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 05-Jun-2026 20:31:45
 */

/* Include Files */
#include "xnrm2.h"
#include "rt_nonfinite.h"
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
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        double absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
      if (rtIsNaN(y)) {
        int b_k;
        b_k = ix0;
        int exitg1;
        do {
          exitg1 = 0;
          if (b_k <= kend) {
            if (rtIsNaN(x[b_k - 1])) {
              exitg1 = 1;
            } else {
              b_k++;
            }
          } else {
            y = rtInf;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }
    }
  }
  return y;
}

/*
 * Arguments    : int n
 *                const double x[121]
 *                int ix0
 * Return Type  : double
 */
double xnrm2(int n, const double x[121], int ix0)
{
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        double absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
      if (rtIsNaN(y)) {
        int b_k;
        b_k = ix0;
        int exitg1;
        do {
          exitg1 = 0;
          if (b_k <= kend) {
            if (rtIsNaN(x[b_k - 1])) {
              exitg1 = 1;
            } else {
              b_k++;
            }
          } else {
            y = rtInf;
            exitg1 = 1;
          }
        } while (exitg1 == 0);
      }
    }
  }
  return y;
}

/*
 * File trailer for xnrm2.c
 *
 * [EOF]
 */
