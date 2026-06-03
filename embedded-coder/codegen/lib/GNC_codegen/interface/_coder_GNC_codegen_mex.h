/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_GNC_codegen_mex.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 22:53:09
 */

#ifndef _CODER_GNC_CODEGEN_MEX_H
#define _CODER_GNC_CODEGEN_MEX_H

/* Include Files */
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS(void);

void unsafe_controller_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[6]);

void unsafe_navigation_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[4],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[16]);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for _coder_GNC_codegen_mex.h
 *
 * [EOF]
 */
