#include "_coder_GNC_codegen_mex.h"
#include "_coder_GNC_codegen_api.h"

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[]) {
  emlrtStack st = {NULL, NULL, NULL};
  const mxArray *c_prhs[7];
  const mxArray *b_prhs[6];
  int32_T i;
  int32_T i1;
  const char_T *entryPointTemplateNames[2] = {"controller_codegen_entry",
                                              "navigation_codegen_entry"};
  mexAtExit(&GNC_codegen_atexit);
  GNC_codegen_initialize();
  st.tls = emlrtRootTLSGlobal;
  switch (emlrtGetEntryPointIndexR2016a(
      &st, nrhs, &prhs[0], (const char_T **)(&entryPointTemplateNames[0]), 2)) {
  case 0:
    for (i = 0; i < 6; i++) {
      b_prhs[i] = prhs[i + 1];
    }
    unsafe_controller_codegen_entry_mexFunction(nlhs, plhs, nrhs - 1, b_prhs);
    break;
  case 1:
    for (i1 = 0; i1 < 7; i1++) {
      c_prhs[i1] = prhs[i1 + 1];
    }
    unsafe_navigation_codegen_entry_mexFunction(nlhs, plhs, nrhs - 1, c_prhs);
    break;
  default:
    /* no actions */
    break;
  }
  GNC_codegen_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void) {
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "windows-1252", true);
  return emlrtRootTLSGlobal;
}

void unsafe_controller_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[3],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[6]) {
  emlrtStack st = {NULL, NULL, NULL};
  const mxArray *b_prhs[6];
  const mxArray *outputs[3];
  int32_T i;
  int32_T i1;
  st.tls = emlrtRootTLSGlobal;

  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 6, 4,
                        24, "controller_codegen_entry");
  }
  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 24,
                        "controller_codegen_entry");
  }

  for (i = 0; i < 6; i++) {
    b_prhs[i] = prhs[i];
  }
  controller_codegen_entry_api(b_prhs, nlhs, outputs);

  if (nlhs < 1) {
    i1 = 1;
  } else {
    i1 = nlhs;
  }
  emlrtReturnArrays(i1, &plhs[0], &outputs[0]);
}

void unsafe_navigation_codegen_entry_mexFunction(int32_T nlhs, mxArray *plhs[6],
                                                 int32_T nrhs,
                                                 const mxArray *prhs[7]) {
  emlrtStack st = {NULL, NULL, NULL};
  const mxArray *b_prhs[7];
  const mxArray *outputs[6];
  int32_T i;
  int32_T i1;
  st.tls = emlrtRootTLSGlobal;

  if (nrhs != 7) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 7, 4,
                        24, "navigation_codegen_entry");
  }
  if (nlhs > 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 24,
                        "navigation_codegen_entry");
  }

  for (i = 0; i < 7; i++) {
    b_prhs[i] = prhs[i];
  }
  navigation_codegen_entry_api(b_prhs, nlhs, outputs);

  if (nlhs < 1) {
    i1 = 1;
  } else {
    i1 = nlhs;
  }
  emlrtReturnArrays(i1, &plhs[0], &outputs[0]);
}
