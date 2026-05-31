/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 15:50:44
 */

#ifndef CONTROLLER_CODEGEN_ENTRY_H
#define CONTROLLER_CODEGEN_ENTRY_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void controller_codegen_entry(
    double b_time, double dt_ctrl, const double xR[2], double pdyn,
    double delta, double w_old, const double coeffs[2], const double P_minus[4],
    double d_old, double w_dot_old, double *u, double *b_r,
    double coeffs_ret[2], double *w_old_ret, double P_minus_ret[4],
    double *d_old_ret, double *w_dot_old_ret);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for controller_codegen_entry.h
 *
 * [EOF]
 */
