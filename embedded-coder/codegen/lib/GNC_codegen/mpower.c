/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mpower.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
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
  double b_c;
  if (rtIsNaN(a) || rtIsNaN(b)) {
    b_c = rtNaN;
  } else {
    double d;
    b_c = fabs(a);
    d = fabs(b);
    if (rtIsInf(b)) {
      if (b_c == 1.0) {
        b_c = 1.0;
      } else if (b_c > 1.0) {
        if (b > 0.0) {
          b_c = rtInf;
        } else {
          b_c = 0.0;
        }
      } else if (b > 0.0) {
        b_c = 0.0;
      } else {
        b_c = rtInf;
      }
    } else if (d == 0.0) {
      b_c = 1.0;
    } else if (d == 1.0) {
      if (b > 0.0) {
        b_c = a;
      } else {
        b_c = 1.0 / a;
      }
    } else if (b == 2.0) {
      b_c = a * a;
    } else if ((b == 0.5) && (a >= 0.0)) {
      b_c = sqrt(a);
    } else if ((a < 0.0) && (b > floor(b))) {
      b_c = rtNaN;
    } else {
      b_c = pow(a, b);
    }
  }
  return b_c;
}

/*
 * File trailer for mpower.c
 *
 * [EOF]
 */
