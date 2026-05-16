/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_navigation_codegen_entry_api.c
 *
 * Code generation for function 'navigation_codegen_entry'
 *
 */

/* Include files */
#include "_coder_navigation_codegen_entry_api.h"
#include "_coder_navigation_codegen_entry_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131675U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "navigation_codegen_entry",                           /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

static const char_T *sv[2] = {"meas", "status"};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static const mxArray *b_emlrt_marshallOut(const real_T u);

static boolean_T c_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static const mxArray *c_emlrt_marshallOut(const struct3_T u);

static boolean_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static const mxArray *d_emlrt_marshallOut(real_T u[2]);

static struct0_T e_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static void emlrtExitTimeCleanupDtorFcn(const void *r);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier);

static const mxArray *emlrt_marshallOut(const struct2_T *u);

static struct0_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3]);

static struct1_T h_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static struct1_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static boolean_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId);

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[3]);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static boolean_T c_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  boolean_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *c_emlrt_marshallOut(const struct3_T u)
{
  static const char_T *b_sv[6] = {"pressure", "temperature",
                                  "density",  "sonic_speed",
                                  "mach",     "dynamic_pressure"};
  const mxArray *y;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 6, (const char_T **)&b_sv[0]));
  emlrtSetFieldR2017b(y, 0, "pressure", b_emlrt_marshallOut(u.pressure), 0);
  emlrtSetFieldR2017b(y, 0, "temperature", b_emlrt_marshallOut(u.temperature),
                      1);
  emlrtSetFieldR2017b(y, 0, "density", b_emlrt_marshallOut(u.density), 2);
  emlrtSetFieldR2017b(y, 0, "sonic_speed", b_emlrt_marshallOut(u.sonic_speed),
                      3);
  emlrtSetFieldR2017b(y, 0, "mach", b_emlrt_marshallOut(u.mach), 4);
  emlrtSetFieldR2017b(y, 0, "dynamic_pressure",
                      b_emlrt_marshallOut(u.dynamic_pressure), 5);
  return y;
}

static boolean_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId)
{
  boolean_T y;
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *d_emlrt_marshallOut(real_T u[2])
{
  static const int32_T i = 0;
  static const int32_T i1 = 2;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

static struct0_T e_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  struct0_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static void emlrtExitTimeCleanupDtorFcn(const void *r)
{
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *emlrt_marshallOut(const struct2_T *u)
{
  static const int32_T i = 4;
  static const int32_T i1 = 3;
  static const int32_T i2 = 3;
  static const int32_T i3 = 11;
  static const char_T *b_sv[5] = {"q", "w", "v", "alt", "x"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *m;
  const mxArray *y;
  real_T *pData;
  int32_T b_i;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 5, (const char_T **)&b_sv[0]));
  b_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  pData[0] = u->q[0];
  pData[1] = u->q[1];
  pData[2] = u->q[2];
  pData[3] = u->q[3];
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, "q", b_y, 0);
  c_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i1, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  pData[0] = u->w[0];
  pData[1] = u->w[1];
  pData[2] = u->w[2];
  emlrtAssign(&c_y, m);
  emlrtSetFieldR2017b(y, 0, "w", c_y, 1);
  d_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i2, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  pData[0] = u->v[0];
  pData[1] = u->v[1];
  pData[2] = u->v[2];
  emlrtAssign(&d_y, m);
  emlrtSetFieldR2017b(y, 0, "v", d_y, 2);
  emlrtSetFieldR2017b(y, 0, "alt", b_emlrt_marshallOut(u->alt), 3);
  e_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i3, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  for (b_i = 0; b_i < 11; b_i++) {
    pData[b_i] = u->x[b_i];
  }
  emlrtAssign(&e_y, m);
  emlrtSetFieldR2017b(y, 0, "x", e_y, 4);
  return y;
}

static struct0_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId)
{
  static const int32_T dims = 0;
  emlrtMsgIdentifier thisId;
  struct0_T y;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 2,
                         (const char_T **)&sv[0], 0U, (const void *)&dims);
  thisId.fIdentifier = "meas";
  g_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "meas")),
      &thisId, y.meas);
  thisId.fIdentifier = "status";
  y.status = d_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "status")),
      &thisId);
  emlrtDestroyArray(&u);
  return y;
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static struct1_T h_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  struct1_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = i_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static struct1_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId)
{
  static const int32_T dims = 0;
  emlrtMsgIdentifier thisId;
  struct1_T y;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 2,
                         (const char_T **)&sv[0], 0U, (const void *)&dims);
  thisId.fIdentifier = "meas";
  y.meas = b_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "meas")),
      &thisId);
  thisId.fIdentifier = "status";
  y.status = d_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "status")),
      &thisId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static boolean_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  boolean_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "logical", false, 0U,
                          (const void *)&dims);
  ret = *emlrtMxGetLogicals(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[3])
{
  static const int32_T dims = 3;
  real_T(*r)[3];
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                          (const void *)&dims);
  r = (real_T(*)[3])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  ret[2] = (*r)[2];
  emlrtDestroyArray(&src);
}

void navigation_codegen_entry_api(const mxArray *const prhs[12], int32_T nlhs,
                                  const mxArray *plhs[4])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  struct0_T ad_accel;
  struct0_T ad_gyro;
  struct0_T board_accel;
  struct0_T board_gyro;
  struct0_T board_mag;
  struct0_T mti_accel;
  struct0_T mti_gyro;
  struct0_T mti_mag;
  struct1_T board_baro;
  struct1_T mti_baro;
  struct2_T state;
  struct3_T airdata;
  real_T(*roll_state)[2];
  real_T cov_norm;
  real_T dt;
  boolean_T flight_phase;
  st.tls = emlrtRootTLSGlobal;
  roll_state = (real_T(*)[2])mxMalloc(sizeof(real_T[2]));
  /* Marshall function inputs */
  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "dt");
  flight_phase = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "flight_phase");
  board_accel = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "board_accel");
  board_gyro = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "board_gyro");
  mti_accel = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "mti_accel");
  mti_gyro = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "mti_gyro");
  ad_accel = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "ad_accel");
  ad_gyro = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "ad_gyro");
  board_baro = h_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "board_baro");
  board_mag = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "board_mag");
  mti_baro = h_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "mti_baro");
  mti_mag = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "mti_mag");
  /* Invoke the target function */
  navigation_codegen_entry(dt, flight_phase, &board_accel, &board_gyro,
                           &mti_accel, &mti_gyro, &ad_accel, &ad_gyro,
                           &board_baro, &board_mag, &mti_baro, &mti_mag, &state,
                           &cov_norm, &airdata, *roll_state);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(&state);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(cov_norm);
  }
  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(airdata);
  }
  if (nlhs > 3) {
    plhs[3] = d_emlrt_marshallOut(*roll_state);
  }
}

void navigation_codegen_entry_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtPushHeapReferenceStackR2021a(
      &st, false, NULL, (void *)&emlrtExitTimeCleanupDtorFcn, NULL, NULL, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  navigation_codegen_entry_xil_terminate();
  navigation_codegen_entry_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void navigation_codegen_entry_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

void navigation_codegen_entry_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (_coder_navigation_codegen_entry_api.c) */
