/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: navigation_codegen_entry.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 23:24:33
 */

#ifndef NAVIGATION_CODEGEN_ENTRY_H
#define NAVIGATION_CODEGEN_ENTRY_H

/* Include Files */
#include "GNC_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void navigation_codegen_entry(
    double dt, bool flight_phase, const double x[11], const double P[121],
    struct1_T *bias, const struct2_T *sens_filt, const struct3_T *sens_input,
    double x_ret[11], double P_ret[121], struct1_T *bias_ret,
    struct2_T *sens_filt_ret);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for navigation_codegen_entry.h
 *
 * [EOF]
 */
