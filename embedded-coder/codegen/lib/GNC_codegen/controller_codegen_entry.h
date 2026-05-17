/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
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
extern void controller_codegen_entry(double b_time, const double xR[2],
                                     double pdyn, double delta, double *u,
                                     double *b_r, double *C_l_delta);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for controller_codegen_entry.h
 *
 * [EOF]
 */
