/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry_terminate.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

/* Include Files */
#include "controller_codegen_entry_terminate.h"
#include "controller_codegen_entry_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
/*
 * Arguments    : void
 * Return Type  : void
 */
void controller_codegen_entry_terminate(void)
{
  omp_destroy_nest_lock(&controller_codegen_entry_nestLockGlobal);
  isInitialized_controller_codegen_entry = false;
}

/*
 * File trailer for controller_codegen_entry_terminate.c
 *
 * [EOF]
 */
