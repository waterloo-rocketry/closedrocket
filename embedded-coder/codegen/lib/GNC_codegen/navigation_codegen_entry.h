/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: navigation_codegen_entry.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 01-Jun-2026 00:25:13
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
    const struct0_T *b, const struct1_T *sf, const struct2_T *board_accel,
    const struct2_T *board_gyro, const struct2_T *mti_accel,
    const struct2_T *mti_gyro, const struct2_T *ad_accel,
    const struct2_T *ad_gyro, const struct3_T *board_baro,
    const struct2_T *board_mag, const struct3_T *mti_baro,
    const struct2_T *mti_mag, double x_ret[11], double P_ret[121],
    struct0_T *b_ret, struct1_T *sf_ret);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for navigation_codegen_entry.h
 *
 * [EOF]
 */
