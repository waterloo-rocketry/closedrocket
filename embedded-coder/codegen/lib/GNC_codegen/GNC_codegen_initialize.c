/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: GNC_codegen_initialize.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

/* Include Files */
#include "GNC_codegen_initialize.h"
#include "GNC_codegen_data.h"
#include "controller_estimator.h"
#include "dynamics.h"
#include "dynamics_jacobian.h"
#include "navigation_codegen_entry.h"
#include "pad_filter.h"
#include "rt_nonfinite.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void GNC_codegen_initialize(void)
{
  controller_estimator_init();
  navigation_codegen_entry_init();
  pad_filter_init();
  dynamics_init();
  dynamics_jacobian_init();
  isInitialized_GNC_codegen = true;
}

/*
 * File trailer for GNC_codegen_initialize.c
 *
 * [EOF]
 */
