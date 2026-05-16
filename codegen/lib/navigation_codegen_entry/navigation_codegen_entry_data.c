/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * navigation_codegen_entry_data.c
 *
 * Code generation for function 'navigation_codegen_entry_data'
 *
 */

/* Include files */
#include "navigation_codegen_entry_data.h"
#include "navigation_codegen_entry_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
boolean_T board_accel_f_not_empty;

b_struct_T param;

b_struct_T b_param;

omp_nest_lock_t navigation_codegen_entry_nestLockGlobal;

const b_struct_T r = {
    {0.46, 0.0, 0.0, 0.0, 49.5, 0.0, 0.0, 0.0, 49.5}, /* J */
    {2.1739130434782608, 0.0, 0.0, 0.0, 0.020202020202020204, 0.0, 0.0, 0.0,
     0.020202020202020204}, /* Jinv */
    {-9.81, 0.0, 0.0}       /* g */
};

boolean_T isInitialized_navigation_codegen_entry = false;

/* End of code generation (navigation_codegen_entry_data.c) */
