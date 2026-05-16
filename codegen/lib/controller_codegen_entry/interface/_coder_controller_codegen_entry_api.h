/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_controller_codegen_entry_api.h
 *
 * Code generation for function 'controller_codegen_entry'
 *
 */

#ifndef _CODER_CONTROLLER_CODEGEN_ENTRY_API_H
#define _CODER_CONTROLLER_CODEGEN_ENTRY_API_H

/* Include files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <string.h>

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void controller_codegen_entry(real_T b_time, real_T xR[2], real_T pdyn,
                              real_T delta, real_T *u, real_T *r,
                              real_T *C_l_delta);

void controller_codegen_entry_api(const mxArray *const prhs[4], int32_T nlhs,
                                  const mxArray *plhs[3]);

void controller_codegen_entry_atexit(void);

void controller_codegen_entry_initialize(void);

void controller_codegen_entry_terminate(void);

void controller_codegen_entry_xil_shutdown(void);

void controller_codegen_entry_xil_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (_coder_controller_codegen_entry_api.h) */
