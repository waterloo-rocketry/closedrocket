/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mpower.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 15:50:44
 */

/* Include Files */
#include "mpower.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : double a
 *                double b
 * Return Type  : double
 */
double mpower(double a, double b)
{
  double c;
  if (rtIsNaN(a) || rtIsNaN(b)) {
    c = rtNaN;
  } else {
    double d;
    c = fabs(a);
    d = fabs(b);
    if (rtIsInf(b)) {
      if (c == 1.0) {
        c = 1.0;
      } else if (c > 1.0) {
        if (b > 0.0) {
          c = rtInf;
        } else {
          c = 0.0;
        }
      } else if (b > 0.0) {
        c = 0.0;
      } else {
        c = rtInf;
      }
    } else if (d == 0.0) {
      c = 1.0;
    } else if (d == 1.0) {
      if (b > 0.0) {
        c = a;
      } else {
        c = 1.0 / a;
      }
    } else if (b == 2.0) {
      c = a * a;
    } else if ((b == 0.5) && (a >= 0.0)) {
      c = sqrt(a);
    } else if ((a < 0.0) && (b > floor(b))) {
      c = rtNaN;
    } else {
      c = pow(a, b);
    }
  }
  return c;
}

/*
 * File trailer for mpower.c
 *
 * [EOF]
 */
