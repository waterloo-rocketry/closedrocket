/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 05-Jun-2026 20:31:45
 */

#ifndef CONTROLLER_CODEGEN_ENTRY_H
#define CONTROLLER_CODEGEN_ENTRY_H

/* Include Files */
#include "GNC_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void controller_codegen_entry(double b_time, double dt_ctrl,
                                     const double xR[2], double pdyn,
                                     double delta, const struct0_T *ctrl_mem_in,
                                     double *u, double *b_r,
                                     struct0_T *ctrl_mem_out);

void controller_codegen_entry_init(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for controller_codegen_entry.h
 *
 * [EOF]
 */
