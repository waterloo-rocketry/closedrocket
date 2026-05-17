/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry_data.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

/* Include Files */
#include "controller_codegen_entry_data.h"
#include "controller_codegen_entry_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
double t;

double c[2];

boolean_T P_minus_not_empty;

boolean_T w_old_not_empty;

double d_old;

double w_dot_old;

boolean_T board_accel_f_not_empty;

b_struct_T param;

b_struct_T b_param;

omp_nest_lock_t controller_codegen_entry_nestLockGlobal;

const b_struct_T r = {
    {0.46, 0.0, 0.0, 0.0, 49.5, 0.0, 0.0, 0.0, 49.5}, /* J */
    {2.1739130434782608, 0.0, 0.0, 0.0, 0.020202020202020204, 0.0, 0.0, 0.0,
     0.020202020202020204}, /* Jinv */
    {-9.81, 0.0, 0.0}       /* g */
};

boolean_T isInitialized_controller_codegen_entry = false;

/*
 * File trailer for controller_codegen_entry_data.c
 *
 * [EOF]
 */
