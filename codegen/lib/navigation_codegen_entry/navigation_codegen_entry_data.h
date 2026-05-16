/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * navigation_codegen_entry_data.h
 *
 * Code generation for function 'navigation_codegen_entry_data'
 *
 */

#ifndef NAVIGATION_CODEGEN_ENTRY_DATA_H
#define NAVIGATION_CODEGEN_ENTRY_DATA_H

/* Include files */
#include "navigation_codegen_entry_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

/* Variable Declarations */
extern boolean_T board_accel_f_not_empty;
extern b_struct_T param;
extern b_struct_T b_param;
extern omp_nest_lock_t navigation_codegen_entry_nestLockGlobal;
extern const b_struct_T r;
extern boolean_T isInitialized_navigation_codegen_entry;

#endif
/* End of code generation (navigation_codegen_entry_data.h) */
