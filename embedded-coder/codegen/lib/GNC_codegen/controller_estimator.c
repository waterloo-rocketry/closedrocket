/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_estimator.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

/* Include Files */
#include "controller_estimator.h"
#include "GNC_codegen_data.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void controller_estimator_init(void)
{
  w_old_not_empty = false;
  P_minus_not_empty = false;
  t = -0.01;
  /*  for /(time - t) */
  c[0] = 2.0;
  c[1] = 0.0;
  /*  initial coefficient guess */
  d_old = 0.0;
  w_dot_old = 0.0;
}

/*
 * File trailer for controller_estimator.c
 *
 * [EOF]
 */
