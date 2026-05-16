/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * navigation_codegen_entry_initialize.c
 *
 * Code generation for function 'navigation_codegen_entry_initialize'
 *
 */

/* Include files */
#include "navigation_codegen_entry_initialize.h"
#include "dynamics.h"
#include "dynamics_jacobian.h"
#include "navigation_codegen_entry.h"
#include "navigation_codegen_entry_data.h"
#include "pad_filter.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
void navigation_codegen_entry_initialize(void)
{
  omp_init_nest_lock(&navigation_codegen_entry_nestLockGlobal);
  navigation_codegen_entry_init();
  pad_filter_init();
  dynamics_init();
  dynamics_jacobian_init();
  isInitialized_navigation_codegen_entry = true;
}

/* End of code generation (navigation_codegen_entry_initialize.c) */
