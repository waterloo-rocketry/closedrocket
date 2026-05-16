/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzlangeM.c
 *
 * Code generation for function 'xzlangeM'
 *
 */

/* Include files */
#include "xzlangeM.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
double xzlangeM(const double x[121])
{
  double y;
  int k;
  boolean_T exitg1;
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

/* End of code generation (xzlangeM.c) */
