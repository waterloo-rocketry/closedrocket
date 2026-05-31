/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: GNC_codegen_initialize.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 15:50:44
 */

/* Include Files */
#include "GNC_codegen_initialize.h"
#include "GNC_codegen_data.h"
#include "dynamics.h"
#include "dynamics_jacobian.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void GNC_codegen_initialize(void)
{
  dynamics_init();
  dynamics_jacobian_init();
  isInitialized_GNC_codegen = true;
}

/*
 * File trailer for GNC_codegen_initialize.c
 *
 * [EOF]
 */
