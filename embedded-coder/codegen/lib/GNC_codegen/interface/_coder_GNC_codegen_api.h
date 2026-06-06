/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_GNC_codegen_api.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 05-Jun-2026 20:31:45
 */

#ifndef _CODER_GNC_CODEGEN_API_H
#define _CODER_GNC_CODEGEN_API_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  real_T board_gyro[3];
  real_T mti_gyro[3];
  real_T ad_gyro[3];
  real_T board_mag_earth[3];
  real_T mti_mag_earth[3];
  real_T board_baro;
  real_T mti_baro;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
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
} struct2_T;
#endif /* typedef_struct2_T */

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  real_T meas[3];
  boolean_T status;
} struct4_T;
#endif /* typedef_struct4_T */

#ifndef typedef_struct5_T
#define typedef_struct5_T
typedef struct {
  real_T meas;
  boolean_T status;
} struct5_T;
#endif /* typedef_struct5_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
typedef struct {
  struct4_T board_accel;
  struct4_T board_gyro;
  struct4_T mti_accel;
  struct4_T mti_gyro;
  struct4_T ad_accel;
  struct4_T ad_gyro;
  struct5_T board_baro;
  struct4_T board_mag;
  struct5_T mti_baro;
  struct4_T mti_mag;
} struct3_T;
#endif /* typedef_struct3_T */

#ifndef typedef_struct6_T
#define typedef_struct6_T
typedef struct {
  real_T pressure;
  real_T temperature;
  real_T density;
  real_T sonic_speed;
  real_T mach;
  real_T dynamic_pressure;
} struct6_T;
#endif /* typedef_struct6_T */

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  real_T coeffs[2];
  real_T w_old;
  real_T P_minus[4];
  real_T d_old;
  real_T w_dot_old;
} struct0_T;
#endif /* typedef_struct0_T */

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
                              real_T pdyn, real_T delta, struct0_T *ctrl_mem_in,
                              real_T *u, real_T *r, struct0_T *ctrl_mem_out);

void controller_codegen_entry_api(const mxArray *const prhs[6], int32_T nlhs,
                                  const mxArray *plhs[3]);

void navigation_codegen_entry(real_T dt, boolean_T flight_phase, real_T x[11],
                              real_T P[121], struct1_T *bias,
                              struct2_T *sens_filt, struct3_T *sens_input,
                              real_T x_ret[11], real_T P_ret[121],
                              struct1_T *bias_ret, struct2_T *sens_filt_ret,
                              real_T *cov_norm, struct6_T *airdata,
                              real_T roll_state[2]);

void navigation_codegen_entry_api(const mxArray *const prhs[7], int32_T nlhs,
                                  const mxArray *plhs[7]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_GNC_codegen_api.h
 *
 * [EOF]
 */
