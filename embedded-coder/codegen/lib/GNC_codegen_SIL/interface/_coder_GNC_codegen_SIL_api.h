#ifndef _CODER_GNC_CODEGEN_SIL_API_H
#define _CODER_GNC_CODEGEN_SIL_API_H

#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

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
#endif

#ifndef typedef_struct2_T
#define typedef_struct2_T
typedef struct {
  real_T board_accel[3];
  real_T board_gyro[3];
  real_T mti_accel[3];
  real_T mti_gyro[3];
  real_T ad_accel[3];
  real_T ad_gyro[3];
  real_T board_baro;
  real_T board_mag[3];
  real_T mti_baro;
  real_T mti_mag[3];
} struct2_T;
#endif

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  real_T meas[3];
  boolean_T status;
} struct4_T;
#endif

#ifndef typedef_struct5_T
#define typedef_struct5_T
typedef struct {
  real_T meas;
  boolean_T status;
} struct5_T;
#endif

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
#endif

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  real_T coeffs[2];
  real_T w;
  real_T P[4];
  real_T delta_lp;
  real_T w_dot_lp;
} struct0_T;
#endif

extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

void GNC_codegen_SIL_atexit(void);

void GNC_codegen_SIL_initialize(void);

void GNC_codegen_SIL_terminate(void);

void GNC_codegen_SIL_xil_shutdown(void);

void GNC_codegen_SIL_xil_terminate(void);

void controller_codegen_entryXilWrapper(real_T b_time, real_T dt_ctrl,
                                        const real_T where_it_is[2],
                                        real_T pdyn, real_T delta_encoder,
                                        struct0_T *ctrl_mem, real_T *u_motor,
                                        real_T where_it_isnt[2],
                                        boolean_T *w_status_ctrl);

void controller_codegen_entry_api(const mxArray *const prhs[6], int32_T nlhs,
                                  const mxArray *plhs[4]);

void navigation_codegen_entryXilWrapper(real_T dt, boolean_T flight_phase,
                                        real_T x[11], real_T P[121],
                                        struct1_T *bias, struct2_T *sens_filt,
                                        const struct3_T *sens_in,
                                        real_T *cov_norm, real_T roll_state[2],
                                        real_T *pdyn, boolean_T *w_status_nav);

void navigation_codegen_entry_api(const mxArray *const prhs[7], int32_T nlhs,
                                  const mxArray *plhs[8]);

void xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
