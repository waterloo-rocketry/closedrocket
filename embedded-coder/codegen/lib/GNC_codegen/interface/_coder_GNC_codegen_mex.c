/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_GNC_codegen_mex.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 22:53:09
 */

/* Include Files */
#include "_coder_GNC_codegen_mex.h"
#include "_coder_GNC_codegen_api.h"

/* Function Definitions */
/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[]
 *                int32_T nrhs
 *                const mxArray *prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *c_prhs[16];
  const mxArray *b_prhs[6];
  int32_T i;
  int32_T i1;
  const char_T *entryPointTemplateNames[2] = {"controller_codegen_entry",
                                              "navigation_codegen_entry"};
  mexAtExit(&GNC_codegen_atexit);
  GNC_codegen_initialize();
  st.tls = emlrtRootTLSGlobal;
  switch (emlrtGetEntryPointIndexR2016a(
      &st, nrhs, &prhs[0], (const char_T **)&entryPointTemplateNames[0], 2)) {
  case 0:
    for (i = 0; i < 6; i++) {
      b_prhs[i] = prhs[i + 1];
    }
    unsafe_controller_codegen_entry_mexFunction(nlhs, plhs, nrhs - 1, b_prhs);
    break;
  case 1:
    for (i1 = 0; i1 < 16; i1++) {
      c_prhs[i1] = prhs[i1 + 1];
    }
    unsafe_navigation_codegen_entry_mexFunction(nlhs, plhs, nrhs - 1, c_prhs);
    break;
  }
  GNC_codegen_terminate();
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[3]
 *                int32_T nrhs
 *                const mxArray *prhs[6]
 * Return Type  : void
 */
void unsafe_controller_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[6])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *b_prhs[6];
  const mxArray *outputs[3];
  int32_T i;
  int32_T i1;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 6, 4,
                        24, "controller_codegen_entry");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 24,
                        "controller_codegen_entry");
  }
  /* Call the function. */
  for (i = 0; i < 6; i++) {
    b_prhs[i] = prhs[i];
  }
  controller_codegen_entry_api(b_prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i1 = 1;
  } else {
    i1 = nlhs;
  }
  emlrtReturnArrays(i1, &plhs[0], &outputs[0]);
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[4]
 *                int32_T nrhs
 *                const mxArray *prhs[16]
 * Return Type  : void
 */
void unsafe_navigation_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[4],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[16])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *b_prhs[16];
  const mxArray *outputs[4];
  int32_T i;
  int32_T i1;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 16) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 16, 4,
                        24, "navigation_codegen_entry");
  }
  if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 24,
                        "navigation_codegen_entry");
  }
  /* Call the function. */
  for (i = 0; i < 16; i++) {
    b_prhs[i] = prhs[i];
  }
  navigation_codegen_entry_api(b_prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i1 = 1;
  } else {
    i1 = nlhs;
  }
  emlrtReturnArrays(i1, &plhs[0], &outputs[0]);
}

/*
 * File trailer for _coder_GNC_codegen_mex.c
 *
 * [EOF]
 */
