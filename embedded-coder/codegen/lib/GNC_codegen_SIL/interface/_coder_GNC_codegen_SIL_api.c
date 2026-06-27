#include "_coder_GNC_codegen_SIL_api.h"
#include "_coder_GNC_codegen_SIL_mex.h"
#include "xil_host_interface.h"

emlrtCTX emlrtRootTLSGlobal = NULL;

static boolean_T xilAlreadyInited;

emlrtContext emlrtContextGlobal = {
    true,
    false,
    131675U,
    NULL,
    "GNC_codegen_SIL",
    NULL,
    false,
    {2045744189U, 2170104910U, 2743257031U, 4284093946U},
    NULL};

static const char_T *sv[10] = {
    "board_accel", "board_gyro", "mti_accel", "mti_gyro", "ad_accel",
    "ad_gyro",     "board_baro", "board_mag", "mti_baro", "mti_mag"};

static const char_T *sv1[2] = {"meas", "status"};

static void GNC_codegen_SIL_once(void);

static void ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[2]);

static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static const mxArray *b_emlrt_marshallOut(real_T u[2]);

static void b_xilHostDeserializer(real_T r[2]);

static void b_xilHostSerializer(const real_T r[2]);

static void bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[4]);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[2];

static const mxArray *c_emlrt_marshallOut(const struct0_T *u);

static void c_xilHostDeserializer(struct0_T *r);

static void c_xilHostSerializer(const struct0_T *r);

static boolean_T cb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2];

static const mxArray *d_emlrt_marshallOut(const boolean_T u);

static void d_xilHostDeserializer(real_T r[4]);

static void d_xilHostSerializer(const real_T r[4]);

static real_T (*db_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[11];

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct0_T *y);

static void e_emlrt_marshallOut(real_T u[11], const mxArray *y);

static void e_xilHostDeserializer(boolean_T *r);

static void e_xilHostSerializer(const boolean_T *r);

static real_T (*eb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[121];

static void emlrtExitTimeCleanupDtorFcn(const void *r);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier);

static const mxArray *emlrt_marshallOut(const real_T u);

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct0_T *y);

static void f_emlrt_marshallOut(real_T u[121], const mxArray *y);

static void f_xilHostDeserializer(real_T r[11]);

static void f_xilHostSerializer(const real_T r[11]);

static void fb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[3]);

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[2]);

static const mxArray *g_emlrt_marshallOut(const struct1_T *u);

static void g_xilHostDeserializer(real_T r[121]);

static void g_xilHostSerializer(const real_T r[121]);

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[4]);

static const mxArray *h_emlrt_marshallOut(const struct2_T *u);

static void h_xilHostDeserializer(struct1_T *r);

static void h_xilHostSerializer(const struct1_T *r);

static boolean_T i_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static void i_xilHostDeserializer(real_T r[3]);

static void i_xilHostSerializer(const real_T r[3]);

static boolean_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static void j_xilHostDeserializer(struct2_T *r);

static void j_xilHostSerializer(const struct2_T *r);

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[11];

static void k_xilHostSerializer(const struct3_T *r);

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[11];

static void l_xilHostSerializer(const struct4_T *r);

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[121];

static void m_xilHostSerializer(const struct5_T *r);

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[121];

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct1_T *y);

static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct1_T *y);

static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[3]);

static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct2_T *y);

static void s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct2_T *y);

static void t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct3_T *y);

static void u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct3_T *y);

static struct4_T v_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static struct5_T w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static real_T x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static void xilHostDeserializer(real_T *r);

static void xilHostSerializer(const real_T *r);

static real_T (*y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2];

static void GNC_codegen_SIL_once(void) { xilAlreadyInited = false; }

static void ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId,
                                real_T ret[2]) {
  static const int32_T dims = 2;
  real_T(*r)[2];
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                          (const void *)&dims);
  r = (real_T(*)[2])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  emlrtDestroyArray(&src);
}

static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId) {
  real_T y;
  y = x_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *b_emlrt_marshallOut(real_T u[2]) {
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

static void b_xilHostDeserializer(real_T r[2]) {
  xilReadData((uint8_T *)&r[0], (size_t)2, MEM_UNIT_DOUBLE_TYPE);
}

static void b_xilHostSerializer(const real_T r[2]) {
  xilWriteData((uint8_T *)&r[0], (size_t)2, MEM_UNIT_DOUBLE_TYPE);
}

static void bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId,
                                real_T ret[4]) {
  static const int32_T dims[2] = {2, 2};
  real_T(*r)[4];
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                          (const void *)&dims[0]);
  r = (real_T(*)[4])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  ret[2] = (*r)[2];
  ret[3] = (*r)[3];
  emlrtDestroyArray(&src);
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[2] {
  emlrtMsgIdentifier thisId;
  real_T(*y)[2];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *c_emlrt_marshallOut(const struct0_T *u) {
  static const int32_T iv[2] = {2, 2};
  static const int32_T i = 2;
  static const char_T *b_sv[5] = {"coeffs", "w", "P", "delta_lp", "w_dot_lp"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  const mxArray *m1;
  const mxArray *y;
  real_T *b_pData;
  real_T *pData;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 5, (const char_T **)&b_sv[0]));
  b_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  pData[0] = u->coeffs[0];
  pData[1] = u->coeffs[1];
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, "coeffs", b_y, 0);
  emlrtSetFieldR2017b(y, 0, "w", emlrt_marshallOut(u->w), 1);
  c_y = NULL;
  m1 = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  b_pData = emlrtMxGetPr(m1);
  b_pData[0] = u->P[0];
  b_pData[1] = u->P[1];
  b_pData[2] = u->P[2];
  b_pData[3] = u->P[3];
  emlrtAssign(&c_y, m1);
  emlrtSetFieldR2017b(y, 0, "P", c_y, 2);
  emlrtSetFieldR2017b(y, 0, "delta_lp", emlrt_marshallOut(u->delta_lp), 3);
  emlrtSetFieldR2017b(y, 0, "w_dot_lp", emlrt_marshallOut(u->w_dot_lp), 4);
  return y;
}

static void c_xilHostDeserializer(struct0_T *r) {
  b_xilHostDeserializer(r->coeffs);
  xilHostDeserializer(&r->w);
  d_xilHostDeserializer(r->P);
  xilHostDeserializer(&r->delta_lp);
  xilHostDeserializer(&r->w_dot_lp);
}

static void c_xilHostSerializer(const struct0_T *r) {
  b_xilHostSerializer(r->coeffs);
  xilHostSerializer(&r->w);
  d_xilHostSerializer(r->P);
  xilHostSerializer(&r->delta_lp);
  xilHostSerializer(&r->w_dot_lp);
}

static boolean_T cb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId) {
  static const int32_T dims = 0;
  boolean_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "logical", false, 0U,
                          (const void *)&dims);
  ret = *emlrtMxGetLogicals(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2] {
  real_T(*y)[2];
  y = y_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static const mxArray *d_emlrt_marshallOut(const boolean_T u) {
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateLogicalScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static void d_xilHostDeserializer(real_T r[4]) {
  xilReadData((uint8_T *)&r[0], (size_t)4, MEM_UNIT_DOUBLE_TYPE);
}

static void d_xilHostSerializer(const real_T r[4]) {
  xilWriteData((uint8_T *)&r[0], (size_t)4, MEM_UNIT_DOUBLE_TYPE);
}

static real_T (*db_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[11] {
  static const int32_T dims = 11;
  real_T(*ret)[11];
  int32_T i;
  boolean_T b = false;
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 1U,
                            (const void *)&dims, &b, &i);
  ret = (real_T(*)[11])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct0_T *y) {
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static void e_emlrt_marshallOut(real_T u[11], const mxArray *y) {
  static const int32_T i = 11;
  emlrtMxSetData((mxArray *)y, &u[0]);
  emlrtSetDimensions((mxArray *)y, &i, 1);
}

static void e_xilHostDeserializer(boolean_T *r) {
  xilReadData((uint8_T *)r, (size_t)1, MEM_UNIT_BOOLEAN_TYPE);
}

static void e_xilHostSerializer(const boolean_T *r) {
  xilWriteData((uint8_T *)r, (size_t)1, MEM_UNIT_BOOLEAN_TYPE);
}

static real_T (*eb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[121] {
  static const int32_T dims[2] = {11, 11};
  real_T(*ret)[121];
  int32_T iv[2];
  boolean_T bv[2] = {false, false};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret = (real_T(*)[121])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void emlrtExitTimeCleanupDtorFcn(const void *r) {
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier) {
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static const mxArray *emlrt_marshallOut(const real_T u) {
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct0_T *y) {
  static const int32_T dims = 0;
  static const char_T *fieldNames[5] = {"coeffs", "w", "P", "delta_lp",
                                        "w_dot_lp"};
  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 5,
                         (const char_T **)&fieldNames[0], 0U,
                         (const void *)&dims);
  thisId.fIdentifier = "coeffs";
  g_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "coeffs")),
      &thisId, y->coeffs);
  thisId.fIdentifier = "w";
  y->w = b_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "w")),
      &thisId);
  thisId.fIdentifier = "P";
  h_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 2, "P")),
      &thisId, y->P);
  thisId.fIdentifier = "delta_lp";
  y->delta_lp = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 3, "delta_lp")),
      &thisId);
  thisId.fIdentifier = "w_dot_lp";
  y->w_dot_lp = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 4, "w_dot_lp")),
      &thisId);
  emlrtDestroyArray(&u);
}

static void f_emlrt_marshallOut(real_T u[121], const mxArray *y) {
  static const int32_T iv[2] = {11, 11};
  emlrtMxSetData((mxArray *)y, &u[0]);
  emlrtSetDimensions((mxArray *)y, &iv[0], 2);
}

static void f_xilHostDeserializer(real_T r[11]) {
  xilReadData((uint8_T *)&r[0], (size_t)11, MEM_UNIT_DOUBLE_TYPE);
}

static void f_xilHostSerializer(const real_T r[11]) {
  xilWriteData((uint8_T *)&r[0], (size_t)11, MEM_UNIT_DOUBLE_TYPE);
}

static void fb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId,
                                real_T ret[3]) {
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

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[2]) {
  ab_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *g_emlrt_marshallOut(const struct1_T *u) {
  static const int32_T i = 3;
  static const int32_T i1 = 3;
  static const int32_T i2 = 3;
  static const int32_T i3 = 3;
  static const int32_T i4 = 3;
  static const char_T *b_sv[7] = {
      "board_gyro",    "mti_gyro",   "ad_gyro", "board_mag_earth",
      "mti_mag_earth", "board_baro", "mti_baro"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *m;
  const mxArray *m1;
  const mxArray *m2;
  const mxArray *m3;
  const mxArray *m4;
  const mxArray *y;
  real_T *b_pData;
  real_T *c_pData;
  real_T *d_pData;
  real_T *e_pData;
  real_T *pData;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 7, (const char_T **)&b_sv[0]));
  b_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  pData[0] = u->board_gyro[0];
  pData[1] = u->board_gyro[1];
  pData[2] = u->board_gyro[2];
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, "board_gyro", b_y, 0);
  c_y = NULL;
  m1 = emlrtCreateNumericArray(1, (const void *)&i1, mxDOUBLE_CLASS, mxREAL);
  b_pData = emlrtMxGetPr(m1);
  b_pData[0] = u->mti_gyro[0];
  b_pData[1] = u->mti_gyro[1];
  b_pData[2] = u->mti_gyro[2];
  emlrtAssign(&c_y, m1);
  emlrtSetFieldR2017b(y, 0, "mti_gyro", c_y, 1);
  d_y = NULL;
  m2 = emlrtCreateNumericArray(1, (const void *)&i2, mxDOUBLE_CLASS, mxREAL);
  c_pData = emlrtMxGetPr(m2);
  c_pData[0] = u->ad_gyro[0];
  c_pData[1] = u->ad_gyro[1];
  c_pData[2] = u->ad_gyro[2];
  emlrtAssign(&d_y, m2);
  emlrtSetFieldR2017b(y, 0, "ad_gyro", d_y, 2);
  e_y = NULL;
  m3 = emlrtCreateNumericArray(1, (const void *)&i3, mxDOUBLE_CLASS, mxREAL);
  d_pData = emlrtMxGetPr(m3);
  d_pData[0] = u->board_mag_earth[0];
  d_pData[1] = u->board_mag_earth[1];
  d_pData[2] = u->board_mag_earth[2];
  emlrtAssign(&e_y, m3);
  emlrtSetFieldR2017b(y, 0, "board_mag_earth", e_y, 3);
  f_y = NULL;
  m4 = emlrtCreateNumericArray(1, (const void *)&i4, mxDOUBLE_CLASS, mxREAL);
  e_pData = emlrtMxGetPr(m4);
  e_pData[0] = u->mti_mag_earth[0];
  e_pData[1] = u->mti_mag_earth[1];
  e_pData[2] = u->mti_mag_earth[2];
  emlrtAssign(&f_y, m4);
  emlrtSetFieldR2017b(y, 0, "mti_mag_earth", f_y, 4);
  emlrtSetFieldR2017b(y, 0, "board_baro", emlrt_marshallOut(u->board_baro), 5);
  emlrtSetFieldR2017b(y, 0, "mti_baro", emlrt_marshallOut(u->mti_baro), 6);
  return y;
}

static void g_xilHostDeserializer(real_T r[121]) {
  xilReadData((uint8_T *)&r[0], (size_t)121, MEM_UNIT_DOUBLE_TYPE);
}

static void g_xilHostSerializer(const real_T r[121]) {
  xilWriteData((uint8_T *)&r[0], (size_t)121, MEM_UNIT_DOUBLE_TYPE);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[4]) {
  bb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *h_emlrt_marshallOut(const struct2_T *u) {
  static const int32_T i = 3;
  static const int32_T i1 = 3;
  static const int32_T i2 = 3;
  static const int32_T i3 = 3;
  static const int32_T i4 = 3;
  static const int32_T i5 = 3;
  static const int32_T i6 = 3;
  static const int32_T i7 = 3;
  static const char_T *b_sv[10] = {
      "board_accel", "board_gyro", "mti_accel", "mti_gyro", "ad_accel",
      "ad_gyro",     "board_baro", "board_mag", "mti_baro", "mti_mag"};
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  const mxArray *m;
  const mxArray *m1;
  const mxArray *m2;
  const mxArray *m3;
  const mxArray *m4;
  const mxArray *m5;
  const mxArray *m6;
  const mxArray *m7;
  const mxArray *y;
  real_T *b_pData;
  real_T *c_pData;
  real_T *d_pData;
  real_T *e_pData;
  real_T *f_pData;
  real_T *g_pData;
  real_T *h_pData;
  real_T *pData;
  y = NULL;
  emlrtAssign(&y, emlrtCreateStructMatrix(1, 1, 10, (const char_T **)&b_sv[0]));
  b_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  pData = emlrtMxGetPr(m);
  pData[0] = u->board_accel[0];
  pData[1] = u->board_accel[1];
  pData[2] = u->board_accel[2];
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, "board_accel", b_y, 0);
  c_y = NULL;
  m1 = emlrtCreateNumericArray(1, (const void *)&i1, mxDOUBLE_CLASS, mxREAL);
  b_pData = emlrtMxGetPr(m1);
  b_pData[0] = u->board_gyro[0];
  b_pData[1] = u->board_gyro[1];
  b_pData[2] = u->board_gyro[2];
  emlrtAssign(&c_y, m1);
  emlrtSetFieldR2017b(y, 0, "board_gyro", c_y, 1);
  d_y = NULL;
  m2 = emlrtCreateNumericArray(1, (const void *)&i2, mxDOUBLE_CLASS, mxREAL);
  c_pData = emlrtMxGetPr(m2);
  c_pData[0] = u->mti_accel[0];
  c_pData[1] = u->mti_accel[1];
  c_pData[2] = u->mti_accel[2];
  emlrtAssign(&d_y, m2);
  emlrtSetFieldR2017b(y, 0, "mti_accel", d_y, 2);
  e_y = NULL;
  m3 = emlrtCreateNumericArray(1, (const void *)&i3, mxDOUBLE_CLASS, mxREAL);
  d_pData = emlrtMxGetPr(m3);
  d_pData[0] = u->mti_gyro[0];
  d_pData[1] = u->mti_gyro[1];
  d_pData[2] = u->mti_gyro[2];
  emlrtAssign(&e_y, m3);
  emlrtSetFieldR2017b(y, 0, "mti_gyro", e_y, 3);
  f_y = NULL;
  m4 = emlrtCreateNumericArray(1, (const void *)&i4, mxDOUBLE_CLASS, mxREAL);
  e_pData = emlrtMxGetPr(m4);
  e_pData[0] = u->ad_accel[0];
  e_pData[1] = u->ad_accel[1];
  e_pData[2] = u->ad_accel[2];
  emlrtAssign(&f_y, m4);
  emlrtSetFieldR2017b(y, 0, "ad_accel", f_y, 4);
  g_y = NULL;
  m5 = emlrtCreateNumericArray(1, (const void *)&i5, mxDOUBLE_CLASS, mxREAL);
  f_pData = emlrtMxGetPr(m5);
  f_pData[0] = u->ad_gyro[0];
  f_pData[1] = u->ad_gyro[1];
  f_pData[2] = u->ad_gyro[2];
  emlrtAssign(&g_y, m5);
  emlrtSetFieldR2017b(y, 0, "ad_gyro", g_y, 5);
  emlrtSetFieldR2017b(y, 0, "board_baro", emlrt_marshallOut(u->board_baro), 6);
  h_y = NULL;
  m6 = emlrtCreateNumericArray(1, (const void *)&i6, mxDOUBLE_CLASS, mxREAL);
  g_pData = emlrtMxGetPr(m6);
  g_pData[0] = u->board_mag[0];
  g_pData[1] = u->board_mag[1];
  g_pData[2] = u->board_mag[2];
  emlrtAssign(&h_y, m6);
  emlrtSetFieldR2017b(y, 0, "board_mag", h_y, 7);
  emlrtSetFieldR2017b(y, 0, "mti_baro", emlrt_marshallOut(u->mti_baro), 8);
  i_y = NULL;
  m7 = emlrtCreateNumericArray(1, (const void *)&i7, mxDOUBLE_CLASS, mxREAL);
  h_pData = emlrtMxGetPr(m7);
  h_pData[0] = u->mti_mag[0];
  h_pData[1] = u->mti_mag[1];
  h_pData[2] = u->mti_mag[2];
  emlrtAssign(&i_y, m7);
  emlrtSetFieldR2017b(y, 0, "mti_mag", i_y, 9);
  return y;
}

static void h_xilHostDeserializer(struct1_T *r) {
  i_xilHostDeserializer(r->board_gyro);
  i_xilHostDeserializer(r->mti_gyro);
  i_xilHostDeserializer(r->ad_gyro);
  i_xilHostDeserializer(r->board_mag_earth);
  i_xilHostDeserializer(r->mti_mag_earth);
  xilHostDeserializer(&r->board_baro);
  xilHostDeserializer(&r->mti_baro);
}

static void h_xilHostSerializer(const struct1_T *r) {
  i_xilHostSerializer(r->board_gyro);
  i_xilHostSerializer(r->mti_gyro);
  i_xilHostSerializer(r->ad_gyro);
  i_xilHostSerializer(r->board_mag_earth);
  i_xilHostSerializer(r->mti_mag_earth);
  xilHostSerializer(&r->board_baro);
  xilHostSerializer(&r->mti_baro);
}

static boolean_T i_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier) {
  emlrtMsgIdentifier thisId;
  boolean_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = j_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static void i_xilHostDeserializer(real_T r[3]) {
  xilReadData((uint8_T *)&r[0], (size_t)3, MEM_UNIT_DOUBLE_TYPE);
}

static void i_xilHostSerializer(const real_T r[3]) {
  xilWriteData((uint8_T *)&r[0], (size_t)3, MEM_UNIT_DOUBLE_TYPE);
}

static boolean_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId) {
  boolean_T y;
  y = cb_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void j_xilHostDeserializer(struct2_T *r) {
  i_xilHostDeserializer(r->board_accel);
  i_xilHostDeserializer(r->board_gyro);
  i_xilHostDeserializer(r->mti_accel);
  i_xilHostDeserializer(r->mti_gyro);
  i_xilHostDeserializer(r->ad_accel);
  i_xilHostDeserializer(r->ad_gyro);
  xilHostDeserializer(&r->board_baro);
  i_xilHostDeserializer(r->board_mag);
  xilHostDeserializer(&r->mti_baro);
  i_xilHostDeserializer(r->mti_mag);
}

static void j_xilHostSerializer(const struct2_T *r) {
  i_xilHostSerializer(r->board_accel);
  i_xilHostSerializer(r->board_gyro);
  i_xilHostSerializer(r->mti_accel);
  i_xilHostSerializer(r->mti_gyro);
  i_xilHostSerializer(r->ad_accel);
  i_xilHostSerializer(r->ad_gyro);
  xilHostSerializer(&r->board_baro);
  i_xilHostSerializer(r->board_mag);
  xilHostSerializer(&r->mti_baro);
  i_xilHostSerializer(r->mti_mag);
}

static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[11] {
  emlrtMsgIdentifier thisId;
  real_T(*y)[11];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = l_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static void k_xilHostSerializer(const struct3_T *r) {
  l_xilHostSerializer(&r->board_accel);
  l_xilHostSerializer(&r->board_gyro);
  l_xilHostSerializer(&r->mti_accel);
  l_xilHostSerializer(&r->mti_gyro);
  l_xilHostSerializer(&r->ad_accel);
  l_xilHostSerializer(&r->ad_gyro);
  m_xilHostSerializer(&r->board_baro);
  l_xilHostSerializer(&r->board_mag);
  m_xilHostSerializer(&r->mti_baro);
  l_xilHostSerializer(&r->mti_mag);
}

static real_T (*l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[11] {
  real_T(*y)[11];
  y = db_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void l_xilHostSerializer(const struct4_T *r) {
  i_xilHostSerializer(r->meas);
  e_xilHostSerializer(&r->status);
}

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[121] {
  emlrtMsgIdentifier thisId;
  real_T(*y)[121];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = n_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static void m_xilHostSerializer(const struct5_T *r) {
  xilHostSerializer(&r->meas);
  e_xilHostSerializer(&r->status);
}

static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[121] {
  real_T(*y)[121];
  y = eb_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct1_T *y) {
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  p_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct1_T *y) {
  static const int32_T dims = 0;
  static const char_T *fieldNames[7] = {
      "board_gyro",    "mti_gyro",   "ad_gyro", "board_mag_earth",
      "mti_mag_earth", "board_baro", "mti_baro"};
  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 7,
                         (const char_T **)&fieldNames[0], 0U,
                         (const void *)&dims);
  thisId.fIdentifier = "board_gyro";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "board_gyro")),
      &thisId, y->board_gyro);
  thisId.fIdentifier = "mti_gyro";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "mti_gyro")),
      &thisId, y->mti_gyro);
  thisId.fIdentifier = "ad_gyro";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 2, "ad_gyro")),
      &thisId, y->ad_gyro);
  thisId.fIdentifier = "board_mag_earth";
  q_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 3,
                                                    "board_mag_earth")),
                     &thisId, y->board_mag_earth);
  thisId.fIdentifier = "mti_mag_earth";
  q_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 4,
                                                    "mti_mag_earth")),
                     &thisId, y->mti_mag_earth);
  thisId.fIdentifier = "board_baro";
  y->board_baro = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 5, "board_baro")),
      &thisId);
  thisId.fIdentifier = "mti_baro";
  y->mti_baro = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 6, "mti_baro")),
      &thisId);
  emlrtDestroyArray(&u);
}

static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               real_T y[3]) {
  fb_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct2_T *y) {
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  s_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static void s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct2_T *y) {
  static const int32_T dims = 0;
  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 10,
                         (const char_T **)&sv[0], 0U, (const void *)&dims);
  thisId.fIdentifier = "board_accel";
  q_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0,
                                                    "board_accel")),
                     &thisId, y->board_accel);
  thisId.fIdentifier = "board_gyro";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "board_gyro")),
      &thisId, y->board_gyro);
  thisId.fIdentifier = "mti_accel";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 2, "mti_accel")),
      &thisId, y->mti_accel);
  thisId.fIdentifier = "mti_gyro";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 3, "mti_gyro")),
      &thisId, y->mti_gyro);
  thisId.fIdentifier = "ad_accel";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 4, "ad_accel")),
      &thisId, y->ad_accel);
  thisId.fIdentifier = "ad_gyro";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 5, "ad_gyro")),
      &thisId, y->ad_gyro);
  thisId.fIdentifier = "board_baro";
  y->board_baro = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 6, "board_baro")),
      &thisId);
  thisId.fIdentifier = "board_mag";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 7, "board_mag")),
      &thisId, y->board_mag);
  thisId.fIdentifier = "mti_baro";
  y->mti_baro = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 8, "mti_baro")),
      &thisId);
  thisId.fIdentifier = "mti_mag";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 9, "mti_mag")),
      &thisId, y->mti_mag);
  emlrtDestroyArray(&u);
}

static void t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct3_T *y) {
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  u_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static void u_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               struct3_T *y) {
  static const int32_T dims = 0;
  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 10,
                         (const char_T **)&sv[0], 0U, (const void *)&dims);
  thisId.fIdentifier = "board_accel";
  y->board_accel =
      v_emlrt_marshallIn(sp,
                         emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0,
                                                        0, "board_accel")),
                         &thisId);
  thisId.fIdentifier = "board_gyro";
  y->board_gyro = v_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "board_gyro")),
      &thisId);
  thisId.fIdentifier = "mti_accel";
  y->mti_accel = v_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 2, "mti_accel")),
      &thisId);
  thisId.fIdentifier = "mti_gyro";
  y->mti_gyro = v_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 3, "mti_gyro")),
      &thisId);
  thisId.fIdentifier = "ad_accel";
  y->ad_accel = v_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 4, "ad_accel")),
      &thisId);
  thisId.fIdentifier = "ad_gyro";
  y->ad_gyro = v_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 5, "ad_gyro")),
      &thisId);
  thisId.fIdentifier = "board_baro";
  y->board_baro = w_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 6, "board_baro")),
      &thisId);
  thisId.fIdentifier = "board_mag";
  y->board_mag = v_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 7, "board_mag")),
      &thisId);
  thisId.fIdentifier = "mti_baro";
  y->mti_baro = w_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 8, "mti_baro")),
      &thisId);
  thisId.fIdentifier = "mti_mag";
  y->mti_mag = v_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 9, "mti_mag")),
      &thisId);
  emlrtDestroyArray(&u);
}

static struct4_T v_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId) {
  static const int32_T dims = 0;
  emlrtMsgIdentifier thisId;
  struct4_T y;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 2,
                         (const char_T **)&sv1[0], 0U, (const void *)&dims);
  thisId.fIdentifier = "meas";
  q_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "meas")),
      &thisId, y.meas);
  thisId.fIdentifier = "status";
  y.status = j_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "status")),
      &thisId);
  emlrtDestroyArray(&u);
  return y;
}

static struct5_T w_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId) {
  static const int32_T dims = 0;
  emlrtMsgIdentifier thisId;
  struct5_T y;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 2,
                         (const char_T **)&sv1[0], 0U, (const void *)&dims);
  thisId.fIdentifier = "meas";
  y.meas = b_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0, "meas")),
      &thisId);
  thisId.fIdentifier = "status";
  y.status = j_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "status")),
      &thisId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T x_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId) {
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void xilHostDeserializer(real_T *r) {
  xilReadData((uint8_T *)r, (size_t)1, MEM_UNIT_DOUBLE_TYPE);
}

static void xilHostSerializer(const real_T *r) {
  xilWriteData((uint8_T *)r, (size_t)1, MEM_UNIT_DOUBLE_TYPE);
}

static real_T (*y_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                   const emlrtMsgIdentifier *msgId))[2] {
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

void GNC_codegen_SIL_atexit(void) {
  emlrtStack st = {NULL, NULL, NULL};
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtPushHeapReferenceStackR2021a(
      &st, false, NULL, (void *)&emlrtExitTimeCleanupDtorFcn, NULL, NULL, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  GNC_codegen_SIL_xil_terminate();
  GNC_codegen_SIL_xil_shutdown();
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void GNC_codegen_SIL_initialize(void) {
  emlrtStack st = {NULL, NULL, NULL};
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    GNC_codegen_SIL_once();
  }
  if (!xilAlreadyInited) {
    xilInitializeHost(&xil_terminate);
    xilAlreadyInited = true;
  }
}

void GNC_codegen_SIL_terminate(void) {
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void controller_codegen_entryXilWrapper(real_T b_time, real_T dt_ctrl,
                                        const real_T xR[2], real_T pdyn,
                                        real_T delta_encoder,
                                        struct0_T *ctrl_mem, real_T *u_motor,
                                        real_T r[2], boolean_T *w_status_ctrl) {

  xilHostSerializer(&b_time);

  xilHostSerializer(&dt_ctrl);

  b_xilHostSerializer(xR);

  xilHostSerializer(&pdyn);

  xilHostSerializer(&delta_encoder);

  c_xilHostSerializer(ctrl_mem);

  xilEntryPointHost(1U);

  xilHostDeserializer(u_motor);

  b_xilHostDeserializer(r);

  c_xilHostDeserializer(ctrl_mem);

  e_xilHostDeserializer(w_status_ctrl);
}

void controller_codegen_entry_api(const mxArray *const prhs[6], int32_T nlhs,
                                  const mxArray *plhs[4]) {
  emlrtStack st = {NULL, NULL, NULL};
  struct0_T ctrl_mem;
  real_T(*r)[2];
  real_T(*xR)[2];
  real_T b_time;
  real_T delta_encoder;
  real_T dt_ctrl;
  real_T pdyn;
  real_T u_motor;
  boolean_T w_status_ctrl;
  st.tls = emlrtRootTLSGlobal;
  r = (real_T(*)[2])mxMalloc(sizeof(real_T[2]));

  b_time = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "time");
  dt_ctrl = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "dt_ctrl");
  xR = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "xR");
  pdyn = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "pdyn");
  delta_encoder = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "delta_encoder");
  e_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "ctrl_mem", &ctrl_mem);

  xilPreEntryPointHost(1U);

  controller_codegen_entryXilWrapper(b_time, dt_ctrl, *xR, pdyn, delta_encoder,
                                     &ctrl_mem, &u_motor, *r, &w_status_ctrl);

  xilPostEntryPointHost(1U);

  plhs[0] = emlrt_marshallOut(u_motor);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut(*r);
  }
  if (nlhs > 2) {
    plhs[2] = c_emlrt_marshallOut(&ctrl_mem);
  }
  if (nlhs > 3) {
    plhs[3] = d_emlrt_marshallOut(w_status_ctrl);
  }
}

void navigation_codegen_entryXilWrapper(real_T dt, boolean_T flight_phase,
                                        real_T x[11], real_T P[121],
                                        struct1_T *bias, struct2_T *sens_filt,
                                        const struct3_T *sens_in,
                                        real_T *cov_norm, real_T roll_state[2],
                                        real_T *pdyn, boolean_T *w_status_nav) {

  xilHostSerializer(&dt);

  e_xilHostSerializer(&flight_phase);

  f_xilHostSerializer(x);

  g_xilHostSerializer(P);

  h_xilHostSerializer(bias);

  j_xilHostSerializer(sens_filt);

  k_xilHostSerializer(sens_in);

  xilEntryPointHost(2U);

  f_xilHostDeserializer(x);

  g_xilHostDeserializer(P);

  h_xilHostDeserializer(bias);

  j_xilHostDeserializer(sens_filt);

  xilHostDeserializer(cov_norm);

  b_xilHostDeserializer(roll_state);

  xilHostDeserializer(pdyn);

  e_xilHostDeserializer(w_status_nav);
}

void navigation_codegen_entry_api(const mxArray *const prhs[7], int32_T nlhs,
                                  const mxArray *plhs[8]) {
  emlrtStack st = {NULL, NULL, NULL};
  const mxArray *prhs_copy_idx_2;
  const mxArray *prhs_copy_idx_3;
  struct1_T bias;
  struct2_T sens_filt;
  struct3_T sens_in;
  real_T(*P)[121];
  real_T(*x)[11];
  real_T(*roll_state)[2];
  real_T cov_norm;
  real_T dt;
  real_T pdyn;
  boolean_T flight_phase;
  boolean_T w_status_nav;
  st.tls = emlrtRootTLSGlobal;
  roll_state = (real_T(*)[2])mxMalloc(sizeof(real_T[2]));
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, true, -1);
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, true, -1);

  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "dt");
  flight_phase = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "flight_phase");
  x = k_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "x");
  P = m_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_3), "P");
  o_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "bias", &bias);
  r_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "sens_filt", &sens_filt);
  t_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "sens_in", &sens_in);

  xilPreEntryPointHost(2U);

  navigation_codegen_entryXilWrapper(dt, flight_phase, *x, *P, &bias,
                                     &sens_filt, &sens_in, &cov_norm,
                                     *roll_state, &pdyn, &w_status_nav);

  xilPostEntryPointHost(2U);

  e_emlrt_marshallOut(*x, prhs_copy_idx_2);
  plhs[0] = prhs_copy_idx_2;
  if (nlhs > 1) {
    f_emlrt_marshallOut(*P, prhs_copy_idx_3);
    plhs[1] = prhs_copy_idx_3;
  }
  if (nlhs > 2) {
    plhs[2] = g_emlrt_marshallOut(&bias);
  }
  if (nlhs > 3) {
    plhs[3] = h_emlrt_marshallOut(&sens_filt);
  }
  if (nlhs > 4) {
    plhs[4] = emlrt_marshallOut(cov_norm);
  }
  if (nlhs > 5) {
    plhs[5] = b_emlrt_marshallOut(*roll_state);
  }
  if (nlhs > 6) {
    plhs[6] = emlrt_marshallOut(pdyn);
  }
  if (nlhs > 7) {
    plhs[7] = d_emlrt_marshallOut(w_status_nav);
  }
}

void xil_terminate(void) {
  xilAlreadyInited = false;
  GNC_codegen_SIL_terminate();
}
