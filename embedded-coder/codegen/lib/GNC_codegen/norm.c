/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: norm.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 05-Jun-2026 20:31:45
 */

/* Include Files */
#include "norm.h"
#include "rt_nonfinite.h"
#include "svd.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : const double x[3]
 * Return Type  : double
 */
double b_norm(const double x[3])
{
  double absxk;
  double scale;
  double t;
  double y;
  bool b;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  y = scale * sqrt(y);
  b = rtIsNaN(y);
  if (b) {
    int k;
    k = 0;
    int exitg1;
    do {
      exitg1 = 0;
      if (k < 3) {
        if (rtIsNaN(x[k])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        y = rtInf;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return y;
}

/*
 * Arguments    : const double x[4]
 * Return Type  : double
 */
double c_norm(const double x[4])
{
  double absxk;
  double scale;
  double t;
  double y;
  bool b;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[3]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  y = scale * sqrt(y);
  b = rtIsNaN(y);
  if (b) {
    int k;
    k = 0;
    int exitg1;
    do {
      exitg1 = 0;
      if (k < 4) {
        if (rtIsNaN(x[k])) {
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        y = rtInf;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return y;
}

/*
 * Arguments    : const double x[121]
 * Return Type  : double
 */
double d_norm(const double x[121])
{
  double y;
  int i;
  y = 0.0;
  for (i = 0; i < 121; i++) {
    double absx;
    absx = fabs(x[i]);
    if (rtIsNaN(absx) || (absx > y)) {
      y = absx;
    }
  }
  if ((!rtIsInf(y)) && (!rtIsNaN(y))) {
    double dv[11];
    svd(x, dv);
    y = dv[0];
  }
  return y;
}

/*
 * File trailer for norm.c
 *
 * [EOF]
 */
