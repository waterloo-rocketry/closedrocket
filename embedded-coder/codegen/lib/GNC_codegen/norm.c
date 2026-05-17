/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: norm.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
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
  double b_t;
  double scale;
  double y;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    b_t = absxk / 3.3121686421112381E-170;
    y = b_t * b_t;
  }
  absxk = fabs(x[1]);
  if (absxk > scale) {
    b_t = scale / absxk;
    y = y * b_t * b_t + 1.0;
    scale = absxk;
  } else {
    b_t = absxk / scale;
    y += b_t * b_t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    b_t = scale / absxk;
    y = y * b_t * b_t + 1.0;
    scale = absxk;
  } else {
    b_t = absxk / scale;
    y += b_t * b_t;
  }
  return scale * sqrt(y);
}

/*
 * Arguments    : const double x[4]
 * Return Type  : double
 */
double c_norm(const double x[4])
{
  double absxk;
  double b_t;
  double scale;
  double y;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    b_t = absxk / 3.3121686421112381E-170;
    y = b_t * b_t;
  }
  absxk = fabs(x[1]);
  if (absxk > scale) {
    b_t = scale / absxk;
    y = y * b_t * b_t + 1.0;
    scale = absxk;
  } else {
    b_t = absxk / scale;
    y += b_t * b_t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    b_t = scale / absxk;
    y = y * b_t * b_t + 1.0;
    scale = absxk;
  } else {
    b_t = absxk / scale;
    y += b_t * b_t;
  }
  absxk = fabs(x[3]);
  if (absxk > scale) {
    b_t = scale / absxk;
    y = y * b_t * b_t + 1.0;
    scale = absxk;
  } else {
    b_t = absxk / scale;
    y += b_t * b_t;
  }
  return scale * sqrt(y);
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
