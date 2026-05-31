/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_GNC_codegen_api.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 14:39:19
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
  real_T board_gyro[3];
  real_T mti_gyro[3];
  real_T ad_gyro[3];
  real_T board_mag_earth[3];
  real_T mti_mag_earth[3];
  real_T board_baro;
  real_T mti_baro;
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  real_T board_accel_f[3];
  real_T board_gyro_f[3];
  real_T mti_accel_f[3];
  real_T mti_gyro_f[3];
  real_T ad_accel_f[3];
  real_T ad_gyro_f[3];
  real_T board_baro_f;
  real_T board_mag_f[3];
  real_T mti_baro_f;
  real_T mti_mag_f[3];
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
typedef struct {
  real_T meas[3];
  boolean_T status;
} struct2_T;
#endif /* typedef_struct2_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
typedef struct {
  real_T meas[3];
  real_T status;
} struct3_T;
#endif /* typedef_struct3_T */

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  real_T meas;
  boolean_T status;
} struct4_T;
#endif /* typedef_struct4_T */

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

void controller_codegen_entry(real_T b_time, real_T dt_ctrl, real_T xR[2],
                              real_T pdyn, real_T delta, real_T w_old,
                              real_T coeffs[2], real_T P_minus[4], real_T d_old,
                              real_T w_dot_old, real_T *u, real_T *r,
                              real_T coeffs_ret[2], real_T *w_old_ret,
                              real_T P_minus_ret[4], real_T *d_old_ret,
                              real_T *w_dot_old_ret);

void controller_codegen_entry_api(const mxArray *const prhs[10], int32_T nlhs,
                                  const mxArray *plhs[7]);

void navigation_codegen_entry(real_T dt, boolean_T flight_phase, real_T x[11],
                              real_T P[121], struct0_T *b, struct1_T *sf,
                              struct2_T *board_accel, struct2_T *board_gyro,
                              struct2_T *mti_accel, struct2_T *mti_gyro,
                              struct2_T *ad_accel, struct3_T *ad_gyro,
                              struct4_T *board_baro, struct2_T *board_mag,
                              struct4_T *mti_baro, struct2_T *mti_mag,
                              real_T x_ret[11], real_T P_ret[121],
                              struct0_T *b_ret, struct1_T *sf_ret);

void navigation_codegen_entry_api(const mxArray *const prhs[16], int32_T nlhs,
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
