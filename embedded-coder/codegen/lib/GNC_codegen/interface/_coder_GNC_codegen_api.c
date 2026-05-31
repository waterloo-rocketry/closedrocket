/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_GNC_codegen_api.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 10:28:26
 */

/* Include Files */
#include "_coder_GNC_codegen_api.h"
#include "_coder_GNC_codegen_mex.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;

emlrtContext emlrtContextGlobal = {
    true,                                                 /* bFirstTime */
    false,                                                /* bInitialized */
    131675U,                                              /* fVersionInfo */
    NULL,                                                 /* fErrorFunction */
    "GNC_codegen",                                        /* fFunctionName */
    NULL,                                                 /* fRTCallStack */
    false,                                                /* bDebugMode */
    {2045744189U, 2170104910U, 2743257031U, 4284093946U}, /* fSigWrd */
    NULL                                                  /* fSigMem */
};

static const char_T *sv[2] = {"meas", "status"};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static const mxArray *b_emlrt_marshallOut(real_T u[2]);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[2];

static const mxArray *c_emlrt_marshallOut(real_T u[4]);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2];

static const mxArray *d_emlrt_marshallOut(const struct2_T *u);

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[4];

static const mxArray *e_emlrt_marshallOut(const struct3_T u);

static void emlrtExitTimeCleanupDtorFcn(const void *r);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier);

static const mxArray *emlrt_marshallOut(const real_T u);

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[4];

static boolean_T g_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static boolean_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static struct0_T i_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static struct0_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3]);

static struct1_T l_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static struct1_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2];

static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[4];

static boolean_T q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId);

static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId, real_T ret[3]);

/* Function Definitions */
/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T
 */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : real_T u[2]
 * Return Type  : const mxArray *
 */
static const mxArray *b_emlrt_marshallOut(real_T u[2])
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

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real_T (*)[2]
 */
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[2]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[2];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : real_T u[4]
 * Return Type  : const mxArray *
 */
static const mxArray *c_emlrt_marshallOut(real_T u[4])
{
  static const int32_T iv[2] = {0, 0};
  static const int32_T iv1[2] = {2, 2};
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &iv1[0], 2);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[2]
 */
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2]
{
  real_T(*y)[2];
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const struct2_T *u
 * Return Type  : const mxArray *
 */
static const mxArray *d_emlrt_marshallOut(const struct2_T *u)
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
  emlrtSetFieldR2017b(y, 0, "alt", emlrt_marshallOut(u->alt), 3);
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

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real_T (*)[4]
 */
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[4]
{
  emlrtMsgIdentifier thisId;
  real_T(*y)[4];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const struct3_T u
 * Return Type  : const mxArray *
 */
static const mxArray *e_emlrt_marshallOut(const struct3_T u)
{
  static const char_T *b_sv[6] = {"pressure", "temperature",
                                  "density",  "sonic_speed",
                                  "mach",     "dynamic_pressure"};
  const mxArray *y;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 6, (const char_T **)&b_sv[0]));
  emlrtSetFieldR2017b(y, 0, "pressure", emlrt_marshallOut(u.pressure), 0);
  emlrtSetFieldR2017b(y, 0, "temperature", emlrt_marshallOut(u.temperature), 1);
  emlrtSetFieldR2017b(y, 0, "density", emlrt_marshallOut(u.density), 2);
  emlrtSetFieldR2017b(y, 0, "sonic_speed", emlrt_marshallOut(u.sonic_speed), 3);
  emlrtSetFieldR2017b(y, 0, "mach", emlrt_marshallOut(u.mach), 4);
  emlrtSetFieldR2017b(y, 0, "dynamic_pressure",
                      emlrt_marshallOut(u.dynamic_pressure), 5);
  return y;
}

/*
 * Arguments    : const void *r
 * Return Type  : void
 */
static void emlrtExitTimeCleanupDtorFcn(const void *r)
{
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : real_T
 */
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

/*
 * Arguments    : const real_T u
 * Return Type  : const mxArray *
 */
static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : real_T (*)[4]
 */
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[4]
{
  real_T(*y)[4];
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : boolean_T
 */
static boolean_T g_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  boolean_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = h_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : boolean_T
 */
static boolean_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId)
{
  boolean_T y;
  y = q_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : struct0_T
 */
static struct0_T i_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  struct0_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : struct0_T
 */
static struct0_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
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
  k_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "meas")),
      &thisId, y.meas);
  thisId.fIdentifier = "status";
  y.status = h_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "status")),
      &thisId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 *                real_T y[3]
 * Return Type  : void
 */
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3])
{
  r_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *nullptr
 *                const char_T *identifier
 * Return Type  : struct1_T
 */
static struct1_T l_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  struct1_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = m_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *u
 *                const emlrtMsgIdentifier *parentId
 * Return Type  : struct1_T
 */
static struct1_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
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
  y.status = h_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "status")),
      &thisId);
  emlrtDestroyArray(&u);
  return y;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T
 */
static real_T n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[2]
 */
static real_T (*o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2]
{
  static const int32_T dims = 2;
  real_T(*ret)[2];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[2])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : real_T (*)[4]
 */
static real_T (*p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[4]
{
  static const int32_T dims[2] = {2, 2};
  real_T(*ret)[4];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[4])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 * Return Type  : boolean_T
 */
static boolean_T q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

/*
 * Arguments    : const emlrtStack *sp
 *                const mxArray *src
 *                const emlrtMsgIdentifier *msgId
 *                real_T ret[3]
 * Return Type  : void
 */
static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

/*
 * Arguments    : void
 * Return Type  : void
 */
void GNC_codegen_atexit(void)
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
  GNC_codegen_xil_terminate();
  GNC_codegen_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void GNC_codegen_initialize(void)
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

/*
 * Arguments    : void
 * Return Type  : void
 */
void GNC_codegen_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/*
 * Arguments    : const mxArray * const prhs[10]
 *                int32_T nlhs
 *                const mxArray *plhs[7]
 * Return Type  : void
 */
void controller_codegen_entry_api(const mxArray *const prhs[10], int32_T nlhs,
                                  const mxArray *plhs[7])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  real_T(*P_minus)[4];
  real_T(*P_minus_ret)[4];
  real_T(*coeffs)[2];
  real_T(*coeffs_ret)[2];
  real_T(*xR)[2];
  real_T b_time;
  real_T d_old;
  real_T d_old_ret;
  real_T delta;
  real_T dt_ctrl;
  real_T pdyn;
  real_T r;
  real_T u;
  real_T w_dot_old;
  real_T w_dot_old_ret;
  real_T w_old;
  real_T w_old_ret;
  st.tls = emlrtRootTLSGlobal;
  coeffs_ret = (real_T(*)[2])mxMalloc(sizeof(real_T[2]));
  P_minus_ret = (real_T(*)[4])mxMalloc(sizeof(real_T[4]));
  /* Marshall function inputs */
  b_time = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "time");
  dt_ctrl = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "dt_ctrl");
  xR = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "xR");
  pdyn = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "pdyn");
  delta = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "delta");
  w_old = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "w_old");
  coeffs = c_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "coeffs");
  P_minus = e_emlrt_marshallIn(&st, emlrtAlias(prhs[7]), "P_minus");
  d_old = emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "d_old");
  w_dot_old = emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "w_dot_old");
  /* Invoke the target function */
  controller_codegen_entry(b_time, dt_ctrl, *xR, pdyn, delta, w_old, *coeffs,
                           *P_minus, d_old, w_dot_old, &u, &r, *coeffs_ret,
                           &w_old_ret, *P_minus_ret, &d_old_ret,
                           &w_dot_old_ret);
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(u);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(r);
  }
  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(*coeffs_ret);
  }
  if (nlhs > 3) {
    plhs[3] = emlrt_marshallOut(w_old_ret);
  }
  if (nlhs > 4) {
    plhs[4] = c_emlrt_marshallOut(*P_minus_ret);
  }
  if (nlhs > 5) {
    plhs[5] = emlrt_marshallOut(d_old_ret);
  }
  if (nlhs > 6) {
    plhs[6] = emlrt_marshallOut(w_dot_old_ret);
  }
}

/*
 * Arguments    : const mxArray * const prhs[12]
 *                int32_T nlhs
 *                const mxArray *plhs[4]
 * Return Type  : void
 */
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
  flight_phase = g_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "flight_phase");
  board_accel = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "board_accel");
  board_gyro = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "board_gyro");
  mti_accel = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "mti_accel");
  mti_gyro = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "mti_gyro");
  ad_accel = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "ad_accel");
  ad_gyro = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "ad_gyro");
  board_baro = l_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "board_baro");
  board_mag = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "board_mag");
  mti_baro = l_emlrt_marshallIn(&st, emlrtAliasP(prhs[10]), "mti_baro");
  mti_mag = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[11]), "mti_mag");
  /* Invoke the target function */
  navigation_codegen_entry(dt, flight_phase, &board_accel, &board_gyro,
                           &mti_accel, &mti_gyro, &ad_accel, &ad_gyro,
                           &board_baro, &board_mag, &mti_baro, &mti_mag, &state,
                           &cov_norm, &airdata, *roll_state);
  /* Marshall function outputs */
  plhs[0] = d_emlrt_marshallOut(&state);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(cov_norm);
  }
  if (nlhs > 2) {
    plhs[2] = e_emlrt_marshallOut(airdata);
  }
  if (nlhs > 3) {
    plhs[3] = b_emlrt_marshallOut(*roll_state);
  }
}

/*
 * File trailer for _coder_GNC_codegen_api.c
 *
 * [EOF]
 */
