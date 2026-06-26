#ifndef _CODER_GNC_CODEGEN_SIL_MEX_H
#define _CODER_GNC_CODEGEN_SIL_MEX_H

#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"

#ifdef __cplusplus
extern "C" {
#endif

MEXFUNCTION_LINKAGE void mexFunction(int32_T nlhs, mxArray *plhs[],
                                     int32_T nrhs, const mxArray *prhs[]);

emlrtCTX mexFunctionCreateRootTLS(void);

void unsafe_controller_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[4],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[6]);

void unsafe_navigation_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[8],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[7]);

#ifdef __cplusplus
}
#endif

#endif
