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

static const char_T *sv[2] = {"meas", "status"};

static void GNC_codegen_SIL_once(void);

static void ab_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[2]);

static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static const mxArray *b_emlrt_marshallOut(const struct0_T *u);

static void b_xilHostDeserializer(struct0_T *r);

static void b_xilHostSerializer(const real_T r[2]);

static void bb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[4]);

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                   const char_T *identifier))[2];

static const mxArray *c_emlrt_marshallOut(real_T u[11]);

static void c_xilHostDeserializer(real_T r[2]);

static void c_xilHostSerializer(const struct0_T *r);

static boolean_T cb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                     const emlrtMsgIdentifier *msgId);

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                   const emlrtMsgIdentifier *parentId))[2];

static const mxArray *d_emlrt_marshallOut(real_T u[121]);

static void d_xilHostDeserializer(real_T r[4]);

static void d_xilHostSerializer(const real_T r[4]);

static real_T (*db_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                    const emlrtMsgIdentifier *msgId))[11];

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier, struct0_T *y);

static const mxArray *e_emlrt_marshallOut(const struct1_T *u);

static void e_xilHostDeserializer(real_T r[11]);

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

static const mxArray *f_emlrt_marshallOut(const struct2_T *u);

static void f_xilHostDeserializer(real_T r[121]);

static void f_xilHostSerializer(const real_T r[11]);

static void fb_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                const emlrtMsgIdentifier *msgId, real_T ret[3]);

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[2]);

static const mxArray *g_emlrt_marshallOut(const struct6_T u);

static void g_xilHostDeserializer(struct1_T *r);

static void g_xilHostSerializer(const real_T r[121]);

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId, real_T y[4]);

static const mxArray *h_emlrt_marshallOut(real_T u[2]);

static void h_xilHostDeserializer(real_T r[3]);

static void h_xilHostSerializer(const struct1_T *r);

static boolean_T i_emlrt_marshallIn(const emlrtStack *sp,
                                    const mxArray *nullptr,
                                    const char_T *identifier);

static void i_xilHostDeserializer(struct2_T *r);

static void i_xilHostSerializer(const real_T r[3]);

static boolean_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                    const emlrtMsgIdentifier *parentId);

static void j_xilHostDeserializer(struct6_T *r);

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

static const mxArray *b_emlrt_marshallOut(const struct0_T *u) {
  static const int32_T iv[2] = {2, 2};
  static const int32_T i = 2;
  static const char_T *b_sv[5] = {"coeffs", "w_old", "P_minus", "d_old",
                                  "w_dot_old"};
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
  emlrtSetFieldR2017b(y, 0, "w_old", emlrt_marshallOut(u->w_old), 1);
  c_y = NULL;
  m1 = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  b_pData = emlrtMxGetPr(m1);
  b_pData[0] = u->P_minus[0];
  b_pData[1] = u->P_minus[1];
  b_pData[2] = u->P_minus[2];
  b_pData[3] = u->P_minus[3];
  emlrtAssign(&c_y, m1);
  emlrtSetFieldR2017b(y, 0, "P_minus", c_y, 2);
  emlrtSetFieldR2017b(y, 0, "d_old", emlrt_marshallOut(u->d_old), 3);
  emlrtSetFieldR2017b(y, 0, "w_dot_old", emlrt_marshallOut(u->w_dot_old), 4);
  return y;
}

static void b_xilHostDeserializer(struct0_T *r) {
  c_xilHostDeserializer(r->coeffs);
  xilHostDeserializer(&r->w_old);
  d_xilHostDeserializer(r->P_minus);
  xilHostDeserializer(&r->d_old);
  xilHostDeserializer(&r->w_dot_old);
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

static const mxArray *c_emlrt_marshallOut(real_T u[11]) {
  static const int32_T i = 0;
  static const int32_T i1 = 11;
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&i, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &i1, 1);
  emlrtAssign(&y, m);
  return y;
}

static void c_xilHostDeserializer(real_T r[2]) {
  xilReadData((uint8_T *)&r[0], (size_t)2, MEM_UNIT_DOUBLE_TYPE);
}

static void c_xilHostSerializer(const struct0_T *r) {
  b_xilHostSerializer(r->coeffs);
  xilHostSerializer(&r->w_old);
  d_xilHostSerializer(r->P_minus);
  xilHostSerializer(&r->d_old);
  xilHostSerializer(&r->w_dot_old);
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

static const mxArray *d_emlrt_marshallOut(real_T u[121]) {
  static const int32_T iv[2] = {0, 0};
  static const int32_T iv1[2] = {11, 11};
  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, &u[0]);
  emlrtSetDimensions((mxArray *)m, &iv1[0], 2);
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

static const mxArray *e_emlrt_marshallOut(const struct1_T *u) {
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

static void e_xilHostDeserializer(real_T r[11]) {
  xilReadData((uint8_T *)&r[0], (size_t)11, MEM_UNIT_DOUBLE_TYPE);
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
  static const char_T *fieldNames[5] = {"coeffs", "w_old", "P_minus", "d_old",
                                        "w_dot_old"};
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
  thisId.fIdentifier = "w_old";
  y->w_old = b_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1, "w_old")),
      &thisId);
  thisId.fIdentifier = "P_minus";
  h_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 2, "P_minus")),
      &thisId, y->P_minus);
  thisId.fIdentifier = "d_old";
  y->d_old = b_emlrt_marshallIn(
      sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 3, "d_old")),
      &thisId);
  thisId.fIdentifier = "w_dot_old";
  y->w_dot_old = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 4, "w_dot_old")),
      &thisId);
  emlrtDestroyArray(&u);
}

static const mxArray *f_emlrt_marshallOut(const struct2_T *u) {
  static const int32_T i = 3;
  static const int32_T i1 = 3;
  static const int32_T i2 = 3;
  static const int32_T i3 = 3;
  static const int32_T i4 = 3;
  static const int32_T i5 = 3;
  static const int32_T i6 = 3;
  static const int32_T i7 = 3;
  static const char_T *b_sv[10] = {
      "board_accel_f", "board_gyro_f", "mti_accel_f",  "mti_gyro_f",
      "ad_accel_f",    "ad_gyro_f",    "board_baro_f", "board_mag_f",
      "mti_baro_f",    "mti_mag_f"};
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
  pData[0] = u->board_accel_f[0];
  pData[1] = u->board_accel_f[1];
  pData[2] = u->board_accel_f[2];
  emlrtAssign(&b_y, m);
  emlrtSetFieldR2017b(y, 0, "board_accel_f", b_y, 0);
  c_y = NULL;
  m1 = emlrtCreateNumericArray(1, (const void *)&i1, mxDOUBLE_CLASS, mxREAL);
  b_pData = emlrtMxGetPr(m1);
  b_pData[0] = u->board_gyro_f[0];
  b_pData[1] = u->board_gyro_f[1];
  b_pData[2] = u->board_gyro_f[2];
  emlrtAssign(&c_y, m1);
  emlrtSetFieldR2017b(y, 0, "board_gyro_f", c_y, 1);
  d_y = NULL;
  m2 = emlrtCreateNumericArray(1, (const void *)&i2, mxDOUBLE_CLASS, mxREAL);
  c_pData = emlrtMxGetPr(m2);
  c_pData[0] = u->mti_accel_f[0];
  c_pData[1] = u->mti_accel_f[1];
  c_pData[2] = u->mti_accel_f[2];
  emlrtAssign(&d_y, m2);
  emlrtSetFieldR2017b(y, 0, "mti_accel_f", d_y, 2);
  e_y = NULL;
  m3 = emlrtCreateNumericArray(1, (const void *)&i3, mxDOUBLE_CLASS, mxREAL);
  d_pData = emlrtMxGetPr(m3);
  d_pData[0] = u->mti_gyro_f[0];
  d_pData[1] = u->mti_gyro_f[1];
  d_pData[2] = u->mti_gyro_f[2];
  emlrtAssign(&e_y, m3);
  emlrtSetFieldR2017b(y, 0, "mti_gyro_f", e_y, 3);
  f_y = NULL;
  m4 = emlrtCreateNumericArray(1, (const void *)&i4, mxDOUBLE_CLASS, mxREAL);
  e_pData = emlrtMxGetPr(m4);
  e_pData[0] = u->ad_accel_f[0];
  e_pData[1] = u->ad_accel_f[1];
  e_pData[2] = u->ad_accel_f[2];
  emlrtAssign(&f_y, m4);
  emlrtSetFieldR2017b(y, 0, "ad_accel_f", f_y, 4);
  g_y = NULL;
  m5 = emlrtCreateNumericArray(1, (const void *)&i5, mxDOUBLE_CLASS, mxREAL);
  f_pData = emlrtMxGetPr(m5);
  f_pData[0] = u->ad_gyro_f[0];
  f_pData[1] = u->ad_gyro_f[1];
  f_pData[2] = u->ad_gyro_f[2];
  emlrtAssign(&g_y, m5);
  emlrtSetFieldR2017b(y, 0, "ad_gyro_f", g_y, 5);
  emlrtSetFieldR2017b(y, 0, "board_baro_f", emlrt_marshallOut(u->board_baro_f),
                      6);
  h_y = NULL;
  m6 = emlrtCreateNumericArray(1, (const void *)&i6, mxDOUBLE_CLASS, mxREAL);
  g_pData = emlrtMxGetPr(m6);
  g_pData[0] = u->board_mag_f[0];
  g_pData[1] = u->board_mag_f[1];
  g_pData[2] = u->board_mag_f[2];
  emlrtAssign(&h_y, m6);
  emlrtSetFieldR2017b(y, 0, "board_mag_f", h_y, 7);
  emlrtSetFieldR2017b(y, 0, "mti_baro_f", emlrt_marshallOut(u->mti_baro_f), 8);
  i_y = NULL;
  m7 = emlrtCreateNumericArray(1, (const void *)&i7, mxDOUBLE_CLASS, mxREAL);
  h_pData = emlrtMxGetPr(m7);
  h_pData[0] = u->mti_mag_f[0];
  h_pData[1] = u->mti_mag_f[1];
  h_pData[2] = u->mti_mag_f[2];
  emlrtAssign(&i_y, m7);
  emlrtSetFieldR2017b(y, 0, "mti_mag_f", i_y, 9);
  return y;
}

static void f_xilHostDeserializer(real_T r[121]) {
  xilReadData((uint8_T *)&r[0], (size_t)121, MEM_UNIT_DOUBLE_TYPE);
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

static const mxArray *g_emlrt_marshallOut(const struct6_T u) {
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

static void g_xilHostDeserializer(struct1_T *r) {
  h_xilHostDeserializer(r->board_gyro);
  h_xilHostDeserializer(r->mti_gyro);
  h_xilHostDeserializer(r->ad_gyro);
  h_xilHostDeserializer(r->board_mag_earth);
  h_xilHostDeserializer(r->mti_mag_earth);
  xilHostDeserializer(&r->board_baro);
  xilHostDeserializer(&r->mti_baro);
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

static const mxArray *h_emlrt_marshallOut(real_T u[2]) {
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

static void h_xilHostDeserializer(real_T r[3]) {
  xilReadData((uint8_T *)&r[0], (size_t)3, MEM_UNIT_DOUBLE_TYPE);
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

static void i_xilHostDeserializer(struct2_T *r) {
  h_xilHostDeserializer(r->board_accel_f);
  h_xilHostDeserializer(r->board_gyro_f);
  h_xilHostDeserializer(r->mti_accel_f);
  h_xilHostDeserializer(r->mti_gyro_f);
  h_xilHostDeserializer(r->ad_accel_f);
  h_xilHostDeserializer(r->ad_gyro_f);
  xilHostDeserializer(&r->board_baro_f);
  h_xilHostDeserializer(r->board_mag_f);
  xilHostDeserializer(&r->mti_baro_f);
  h_xilHostDeserializer(r->mti_mag_f);
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

static void j_xilHostDeserializer(struct6_T *r) {
  xilHostDeserializer(&r->pressure);
  xilHostDeserializer(&r->temperature);
  xilHostDeserializer(&r->density);
  xilHostDeserializer(&r->sonic_speed);
  xilHostDeserializer(&r->mach);
  xilHostDeserializer(&r->dynamic_pressure);
}

static void j_xilHostSerializer(const struct2_T *r) {
  i_xilHostSerializer(r->board_accel_f);
  i_xilHostSerializer(r->board_gyro_f);
  i_xilHostSerializer(r->mti_accel_f);
  i_xilHostSerializer(r->mti_gyro_f);
  i_xilHostSerializer(r->ad_accel_f);
  i_xilHostSerializer(r->ad_gyro_f);
  xilHostSerializer(&r->board_baro_f);
  i_xilHostSerializer(r->board_mag_f);
  xilHostSerializer(&r->mti_baro_f);
  i_xilHostSerializer(r->mti_mag_f);
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
  static const char_T *fieldNames[10] = {
      "board_accel_f", "board_gyro_f", "mti_accel_f",  "mti_gyro_f",
      "ad_accel_f",    "ad_gyro_f",    "board_baro_f", "board_mag_f",
      "mti_baro_f",    "mti_mag_f"};
  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 10,
                         (const char_T **)&fieldNames[0], 0U,
                         (const void *)&dims);
  thisId.fIdentifier = "board_accel_f";
  q_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 0,
                                                    "board_accel_f")),
                     &thisId, y->board_accel_f);
  thisId.fIdentifier = "board_gyro_f";
  q_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 1,
                                                    "board_gyro_f")),
                     &thisId, y->board_gyro_f);
  thisId.fIdentifier = "mti_accel_f";
  q_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 2,
                                                    "mti_accel_f")),
                     &thisId, y->mti_accel_f);
  thisId.fIdentifier = "mti_gyro_f";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 3, "mti_gyro_f")),
      &thisId, y->mti_gyro_f);
  thisId.fIdentifier = "ad_accel_f";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 4, "ad_accel_f")),
      &thisId, y->ad_accel_f);
  thisId.fIdentifier = "ad_gyro_f";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 5, "ad_gyro_f")),
      &thisId, y->ad_gyro_f);
  thisId.fIdentifier = "board_baro_f";
  y->board_baro_f =
      b_emlrt_marshallIn(sp,
                         emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0,
                                                        6, "board_baro_f")),
                         &thisId);
  thisId.fIdentifier = "board_mag_f";
  q_emlrt_marshallIn(sp,
                     emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 7,
                                                    "board_mag_f")),
                     &thisId, y->board_mag_f);
  thisId.fIdentifier = "mti_baro_f";
  y->mti_baro_f = b_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 8, "mti_baro_f")),
      &thisId);
  thisId.fIdentifier = "mti_mag_f";
  q_emlrt_marshallIn(
      sp,
      emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 9, "mti_mag_f")),
      &thisId, y->mti_mag_f);
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
  static const char_T *fieldNames[10] = {
      "board_accel", "board_gyro", "mti_accel", "mti_gyro", "ad_accel",
      "ad_gyro",     "board_baro", "board_mag", "mti_baro", "mti_mag"};
  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 10,
                         (const char_T **)&fieldNames[0], 0U,
                         (const void *)&dims);
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
                         (const char_T **)&sv[0], 0U, (const void *)&dims);
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
                         (const char_T **)&sv[0], 0U, (const void *)&dims);
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
                                        real_T delta,
                                        const struct0_T *ctrl_mem_in, real_T *u,
                                        real_T *r, struct0_T *ctrl_mem_out) {

  xilHostSerializer(&b_time);

  xilHostSerializer(&dt_ctrl);

  b_xilHostSerializer(xR);

  xilHostSerializer(&pdyn);

  xilHostSerializer(&delta);

  c_xilHostSerializer(ctrl_mem_in);

  xilEntryPointHost(1U);

  xilHostDeserializer(u);

  xilHostDeserializer(r);

  b_xilHostDeserializer(ctrl_mem_out);
}

void controller_codegen_entry_api(const mxArray *const prhs[6], int32_T nlhs,
                                  const mxArray *plhs[3]) {
  emlrtStack st = {NULL, NULL, NULL};
  struct0_T ctrl_mem_in;
  struct0_T ctrl_mem_out;
  real_T(*xR)[2];
  real_T b_time;
  real_T delta;
  real_T dt_ctrl;
  real_T pdyn;
  real_T r;
  real_T u;
  st.tls = emlrtRootTLSGlobal;

  b_time = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "time");
  dt_ctrl = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "dt_ctrl");
  xR = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "xR");
  pdyn = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "pdyn");
  delta = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "delta");
  e_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "ctrl_mem_in", &ctrl_mem_in);

  xilPreEntryPointHost(1U);

  controller_codegen_entryXilWrapper(b_time, dt_ctrl, *xR, pdyn, delta,
                                     &ctrl_mem_in, &u, &r, &ctrl_mem_out);

  xilPostEntryPointHost(1U);

  plhs[0] = emlrt_marshallOut(u);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(r);
  }
  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(&ctrl_mem_out);
  }
}

void navigation_codegen_entryXilWrapper(
    real_T dt, boolean_T flight_phase, const real_T x[11], const real_T P[121],
    const struct1_T *bias, const struct2_T *sens_filt,
    const struct3_T *sens_input, real_T x_ret[11], real_T P_ret[121],
    struct1_T *bias_ret, struct2_T *sens_filt_ret, struct6_T *airdata,
    real_T roll_state[2]) {

  xilHostSerializer(&dt);

  e_xilHostSerializer(&flight_phase);

  f_xilHostSerializer(x);

  g_xilHostSerializer(P);

  h_xilHostSerializer(bias);

  j_xilHostSerializer(sens_filt);

  k_xilHostSerializer(sens_input);

  xilEntryPointHost(2U);

  e_xilHostDeserializer(x_ret);

  f_xilHostDeserializer(P_ret);

  g_xilHostDeserializer(bias_ret);

  i_xilHostDeserializer(sens_filt_ret);

  j_xilHostDeserializer(airdata);

  c_xilHostDeserializer(roll_state);
}

void navigation_codegen_entry_api(const mxArray *const prhs[7], int32_T nlhs,
                                  const mxArray *plhs[6]) {
  emlrtStack st = {NULL, NULL, NULL};
  struct1_T bias;
  struct1_T bias_ret;
  struct2_T sens_filt;
  struct2_T sens_filt_ret;
  struct3_T sens_input;
  struct6_T airdata;
  real_T(*P)[121];
  real_T(*P_ret)[121];
  real_T(*x)[11];
  real_T(*x_ret)[11];
  real_T(*roll_state)[2];
  real_T dt;
  boolean_T flight_phase;
  st.tls = emlrtRootTLSGlobal;
  x_ret = (real_T(*)[11])mxMalloc(sizeof(real_T[11]));
  P_ret = (real_T(*)[121])mxMalloc(sizeof(real_T[121]));
  roll_state = (real_T(*)[2])mxMalloc(sizeof(real_T[2]));

  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "dt");
  flight_phase = i_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "flight_phase");
  x = k_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "x");
  P = m_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "P");
  o_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "bias", &bias);
  r_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "sens_filt", &sens_filt);
  t_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "sens_input", &sens_input);

  xilPreEntryPointHost(2U);

  navigation_codegen_entryXilWrapper(
      dt, flight_phase, *x, *P, &bias, &sens_filt, &sens_input, *x_ret, *P_ret,
      &bias_ret, &sens_filt_ret, &airdata, *roll_state);

  xilPostEntryPointHost(2U);

  plhs[0] = c_emlrt_marshallOut(*x_ret);
  if (nlhs > 1) {
    plhs[1] = d_emlrt_marshallOut(*P_ret);
  }
  if (nlhs > 2) {
    plhs[2] = e_emlrt_marshallOut(&bias_ret);
  }
  if (nlhs > 3) {
    plhs[3] = f_emlrt_marshallOut(&sens_filt_ret);
  }
  if (nlhs > 4) {
    plhs[4] = g_emlrt_marshallOut(airdata);
  }
  if (nlhs > 5) {
    plhs[5] = h_emlrt_marshallOut(*roll_state);
  }
}

void xil_terminate(void) {
  xilAlreadyInited = false;
  GNC_codegen_SIL_terminate();
}
