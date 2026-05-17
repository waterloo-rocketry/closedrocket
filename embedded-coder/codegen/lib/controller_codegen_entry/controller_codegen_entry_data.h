/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry_data.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

#ifndef CONTROLLER_CODEGEN_ENTRY_DATA_H
#define CONTROLLER_CODEGEN_ENTRY_DATA_H

/* Include Files */
#include "controller_codegen_entry_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

/* Variable Declarations */
extern double t;
extern double c[2];
extern boolean_T P_minus_not_empty;
extern boolean_T w_old_not_empty;
extern double d_old;
extern double w_dot_old;
extern boolean_T board_accel_f_not_empty;
extern b_struct_T param;
extern b_struct_T b_param;
extern omp_nest_lock_t controller_codegen_entry_nestLockGlobal;
extern const b_struct_T r;
extern boolean_T isInitialized_controller_codegen_entry;

#endif
/*
 * File trailer for controller_codegen_entry_data.h
 *
 * [EOF]
 */
