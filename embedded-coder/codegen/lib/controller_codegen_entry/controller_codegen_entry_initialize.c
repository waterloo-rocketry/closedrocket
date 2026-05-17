/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry_initialize.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

/* Include Files */
#include "controller_codegen_entry_initialize.h"
#include "controller_codegen_entry_data.h"
#include "controller_estimator.h"
#include "dynamics.h"
#include "dynamics_jacobian.h"
#include "navigation_codegen_entry.h"
#include "pad_filter.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void controller_codegen_entry_initialize(void)
{
  omp_init_nest_lock(&controller_codegen_entry_nestLockGlobal);
  controller_estimator_init();
  navigation_codegen_entry_init();
  pad_filter_init();
  dynamics_init();
  dynamics_jacobian_init();
  isInitialized_controller_codegen_entry = true;
}

/*
 * File trailer for controller_codegen_entry_initialize.c
 *
 * [EOF]
 */
