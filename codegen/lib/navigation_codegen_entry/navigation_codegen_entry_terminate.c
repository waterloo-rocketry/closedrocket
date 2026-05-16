/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * navigation_codegen_entry_terminate.c
 *
 * Code generation for function 'navigation_codegen_entry_terminate'
 *
 */

/* Include files */
#include "navigation_codegen_entry_terminate.h"
#include "navigation_codegen_entry_data.h"
#include "rt_nonfinite.h"
#include "omp.h"

/* Function Definitions */
void navigation_codegen_entry_terminate(void)
{
  omp_destroy_nest_lock(&navigation_codegen_entry_nestLockGlobal);
  isInitialized_navigation_codegen_entry = false;
}

/* End of code generation (navigation_codegen_entry_terminate.c) */
