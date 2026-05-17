/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_GNC_codegen_api.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

#ifndef _CODER_GNC_CODEGEN_API_H
#define _CODER_GNC_CODEGEN_API_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  real_T meas[3];
  boolean_T status;
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  real_T meas;
  boolean_T status;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
typedef struct {
  real_T q[4];
  real_T w[3];
  real_T v[3];
  real_T alt;
  real_T x[11];
} struct2_T;
#endif /* typedef_struct2_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
typedef struct {
  real_T pressure;
  real_T temperature;
  real_T density;
  real_T sonic_speed;
  real_T mach;
  real_T dynamic_pressure;
} struct3_T;
#endif /* typedef_struct3_T */

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void GNC_codegen_atexit(void);

void GNC_codegen_initialize(void);

void GNC_codegen_terminate(void);

void GNC_codegen_xil_shutdown(void);

void GNC_codegen_xil_terminate(void);

void controller_codegen_entry(real_T b_time, real_T xR[2], real_T pdyn,
                              real_T delta, real_T *u, real_T *r,
                              real_T *C_l_delta);

void controller_codegen_entry_api(const mxArray *const prhs[4], int32_T nlhs,
                                  const mxArray *plhs[3]);

void navigation_codegen_entry(real_T dt, boolean_T flight_phase,
                              struct0_T *board_accel, struct0_T *board_gyro,
                              struct0_T *mti_accel, struct0_T *mti_gyro,
                              struct0_T *ad_accel, struct0_T *ad_gyro,
                              struct1_T *board_baro, struct0_T *board_mag,
                              struct1_T *mti_baro, struct0_T *mti_mag,
                              struct2_T *state, real_T *cov_norm,
                              struct3_T *airdata, real_T roll_state[2]);

void navigation_codegen_entry_api(const mxArray *const prhs[12], int32_T nlhs,
                                  const mxArray *plhs[4]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_GNC_codegen_api.h
 *
 * [EOF]
 */
