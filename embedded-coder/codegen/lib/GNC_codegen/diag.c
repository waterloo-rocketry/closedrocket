/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: diag.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 15:50:44
 */

/* Include Files */
#include "diag.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : double d[4]
 * Return Type  : void
 */
void diag(double d[4])
{
  d[1] = 0.0;
  d[2] = 0.0;
  d[0] = 1.0E-5;
  d[3] = 1.0E-9;
}

/*
 * File trailer for diag.c
 *
 * [EOF]
 */
