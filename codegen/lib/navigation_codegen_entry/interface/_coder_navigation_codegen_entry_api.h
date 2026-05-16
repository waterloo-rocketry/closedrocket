/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_navigation_codegen_entry_api.h
 *
 * Code generation for function 'navigation_codegen_entry'
 *
 */

#ifndef _CODER_NAVIGATION_CODEGEN_ENTRY_API_H
#define _CODER_NAVIGATION_CODEGEN_ENTRY_API_H

/* Include files */
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

void navigation_codegen_entry_atexit(void);

void navigation_codegen_entry_initialize(void);

void navigation_codegen_entry_terminate(void);

void navigation_codegen_entry_xil_shutdown(void);

void navigation_codegen_entry_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_navigation_codegen_entry_api.h) */
