/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: navigation_codegen_entry.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
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
    double dt, bool flight_phase, const struct0_T *board_accel,
    const struct0_T *board_gyro, const struct0_T *mti_accel,
    const struct0_T *mti_gyro, const struct0_T *ad_accel,
    const struct0_T *ad_gyro, const struct1_T *board_baro,
    const struct0_T *board_mag, const struct1_T *mti_baro,
    const struct0_T *mti_mag, struct2_T *state, double *cov_norm,
    struct3_T *airdata, double roll_state[2]);

void navigation_codegen_entry_init(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for navigation_codegen_entry.h
 *
 * [EOF]
 */
