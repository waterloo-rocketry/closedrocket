/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: xzlangeM.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 05-Jun-2026 20:31:45
 */

/* Include Files */
#include "xzlangeM.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * Arguments    : const double x[121]
 * Return Type  : double
 */
double xzlangeM(const double x[121])
{
  double y;
  int k;
  bool exitg1;
  y = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 121)) {
    double absxk;
    absxk = fabs(x[k]);
    if (rtIsNaN(absxk)) {
      y = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > y) {
        y = absxk;
      }
      k++;
    }
  }
  return y;
}

/*
 * File trailer for xzlangeM.c
 *
 * [EOF]
 */
