#include "_coder_GNC_codegen_SIL_target.h"
#include "GNC_codegen_SIL.h"
#include "GNC_codegen_SIL_types.h"

static GNC_codegen_SILStackData *SD;

static void b_xilTargetDeserializer(double r[2]);

static void b_xilTargetSerializer(const struct0_T *r);

static void c_xilTargetDeserializer(struct0_T *r);

static void c_xilTargetSerializer(const double r[2]);

static void d_xilTargetDeserializer(double r[4]);

static void d_xilTargetSerializer(const double r[4]);

static void e_xilTargetDeserializer(boolean_T *r);

static void e_xilTargetSerializer(const boolean_T *r);

static void f_xilTargetDeserializer(double r[11]);

static void f_xilTargetSerializer(const double r[11]);

static void g_xilTargetDeserializer(double r[121]);

static void g_xilTargetSerializer(const double r[121]);

static GNC_codegen_SILStackData *getGNC_codegen_SILStackData(void);

static void h_xilTargetDeserializer(struct1_T *r);

static void h_xilTargetSerializer(const struct1_T *r);

static void i_xilTargetDeserializer(double r[3]);

static void i_xilTargetSerializer(const double r[3]);

static void j_xilTargetDeserializer(struct2_T *r);

static void j_xilTargetSerializer(const struct2_T *r);

static void k_xilTargetDeserializer(struct3_T *r);

static void l_xilTargetDeserializer(struct4_T *r);

static void m_xilTargetDeserializer(struct5_T *r);

static void xilTargetDeserializer(double *r);

static void xilTargetSerializer(const double *r);

static void b_xilTargetDeserializer(double r[2]) {
  xilReadData((MemUnit_T *)&r[0], sizeof(double) << 1);
}

static void b_xilTargetSerializer(const struct0_T *r) {
  c_xilTargetSerializer(r->coeffs);
  xilTargetSerializer(&r->w_old);
  d_xilTargetSerializer(r->P_minus);
  xilTargetSerializer(&r->d_old);
  xilTargetSerializer(&r->w_dot_old);
  xilTargetSerializer(&r->delta_lp);
  xilTargetSerializer(&r->w_dot_lp);
  xilTargetSerializer(&r->w);
  d_xilTargetSerializer(r->P);
}

static void c_xilTargetDeserializer(struct0_T *r) {
  b_xilTargetDeserializer(r->coeffs);
  xilTargetDeserializer(&r->w_old);
  d_xilTargetDeserializer(r->P_minus);
  xilTargetDeserializer(&r->d_old);
  xilTargetDeserializer(&r->w_dot_old);
  xilTargetDeserializer(&r->delta_lp);
  xilTargetDeserializer(&r->w_dot_lp);
  xilTargetDeserializer(&r->w);
  d_xilTargetDeserializer(r->P);
}

static void c_xilTargetSerializer(const double r[2]) {
  xilWriteData((const MemUnit_T *)&r[0], sizeof(double) << 1);
}

static void d_xilTargetDeserializer(double r[4]) {
  xilReadData((MemUnit_T *)&r[0], sizeof(double) << 2);
}

static void d_xilTargetSerializer(const double r[4]) {
  xilWriteData((const MemUnit_T *)&r[0], sizeof(double) << 2);
}

static void e_xilTargetDeserializer(boolean_T *r) {
  xilReadData((MemUnit_T *)r, sizeof(boolean_T));
}

static void e_xilTargetSerializer(const boolean_T *r) {
  xilWriteData((const MemUnit_T *)r, sizeof(boolean_T));
}

static void f_xilTargetDeserializer(double r[11]) {
  xilReadData((MemUnit_T *)&r[0], (unsigned int)(11 * (int)sizeof(double)));
}

static void f_xilTargetSerializer(const double r[11]) {
  xilWriteData((const MemUnit_T *)&r[0],
               (unsigned int)(11 * (int)sizeof(double)));
}

static void g_xilTargetDeserializer(double r[121]) {
  xilReadData((MemUnit_T *)&r[0], (unsigned int)(121 * (int)sizeof(double)));
}

static void g_xilTargetSerializer(const double r[121]) {
  xilWriteData((const MemUnit_T *)&r[0],
               (unsigned int)(121 * (int)sizeof(double)));
}

static GNC_codegen_SILStackData *getGNC_codegen_SILStackData(void) {
  static GNC_codegen_SILPersistentData PersistentData;
  static GNC_codegen_SILStackData GlobalStackData;
  GlobalStackData.pd = &PersistentData;
  return &GlobalStackData;
}

static void h_xilTargetDeserializer(struct1_T *r) {
  i_xilTargetDeserializer(r->board_gyro);
  i_xilTargetDeserializer(r->mti_gyro);
  i_xilTargetDeserializer(r->ad_gyro);
  i_xilTargetDeserializer(r->board_mag_earth);
  i_xilTargetDeserializer(r->mti_mag_earth);
  xilTargetDeserializer(&r->board_baro);
  xilTargetDeserializer(&r->mti_baro);
}

static void h_xilTargetSerializer(const struct1_T *r) {
  i_xilTargetSerializer(r->board_gyro);
  i_xilTargetSerializer(r->mti_gyro);
  i_xilTargetSerializer(r->ad_gyro);
  i_xilTargetSerializer(r->board_mag_earth);
  i_xilTargetSerializer(r->mti_mag_earth);
  xilTargetSerializer(&r->board_baro);
  xilTargetSerializer(&r->mti_baro);
}

static void i_xilTargetDeserializer(double r[3]) {
  xilReadData((MemUnit_T *)&r[0], (unsigned int)(3 * (int)sizeof(double)));
}

static void i_xilTargetSerializer(const double r[3]) {
  xilWriteData((const MemUnit_T *)&r[0],
               (unsigned int)(3 * (int)sizeof(double)));
}

static void j_xilTargetDeserializer(struct2_T *r) {
  i_xilTargetDeserializer(r->board_accel_f);
  i_xilTargetDeserializer(r->board_gyro_f);
  i_xilTargetDeserializer(r->mti_accel_f);
  i_xilTargetDeserializer(r->mti_gyro_f);
  i_xilTargetDeserializer(r->ad_accel_f);
  i_xilTargetDeserializer(r->ad_gyro_f);
  xilTargetDeserializer(&r->board_baro_f);
  i_xilTargetDeserializer(r->board_mag_f);
  xilTargetDeserializer(&r->mti_baro_f);
  i_xilTargetDeserializer(r->mti_mag_f);
  i_xilTargetDeserializer(r->board_accel);
  i_xilTargetDeserializer(r->board_gyro);
  i_xilTargetDeserializer(r->mti_accel);
  i_xilTargetDeserializer(r->mti_gyro);
  i_xilTargetDeserializer(r->ad_accel);
  i_xilTargetDeserializer(r->ad_gyro);
  xilTargetDeserializer(&r->board_baro);
  i_xilTargetDeserializer(r->board_mag);
  xilTargetDeserializer(&r->mti_baro);
  i_xilTargetDeserializer(r->mti_mag);
}

static void j_xilTargetSerializer(const struct2_T *r) {
  i_xilTargetSerializer(r->board_accel_f);
  i_xilTargetSerializer(r->board_gyro_f);
  i_xilTargetSerializer(r->mti_accel_f);
  i_xilTargetSerializer(r->mti_gyro_f);
  i_xilTargetSerializer(r->ad_accel_f);
  i_xilTargetSerializer(r->ad_gyro_f);
  xilTargetSerializer(&r->board_baro_f);
  i_xilTargetSerializer(r->board_mag_f);
  xilTargetSerializer(&r->mti_baro_f);
  i_xilTargetSerializer(r->mti_mag_f);
  i_xilTargetSerializer(r->board_accel);
  i_xilTargetSerializer(r->board_gyro);
  i_xilTargetSerializer(r->mti_accel);
  i_xilTargetSerializer(r->mti_gyro);
  i_xilTargetSerializer(r->ad_accel);
  i_xilTargetSerializer(r->ad_gyro);
  xilTargetSerializer(&r->board_baro);
  i_xilTargetSerializer(r->board_mag);
  xilTargetSerializer(&r->mti_baro);
  i_xilTargetSerializer(r->mti_mag);
}

static void k_xilTargetDeserializer(struct3_T *r) {
  l_xilTargetDeserializer(&r->board_accel);
  l_xilTargetDeserializer(&r->board_gyro);
  l_xilTargetDeserializer(&r->mti_accel);
  l_xilTargetDeserializer(&r->mti_gyro);
  l_xilTargetDeserializer(&r->ad_accel);
  l_xilTargetDeserializer(&r->ad_gyro);
  m_xilTargetDeserializer(&r->board_baro);
  l_xilTargetDeserializer(&r->board_mag);
  m_xilTargetDeserializer(&r->mti_baro);
  l_xilTargetDeserializer(&r->mti_mag);
}

static void l_xilTargetDeserializer(struct4_T *r) {
  i_xilTargetDeserializer(r->meas);
  e_xilTargetDeserializer(&r->status);
}

static void m_xilTargetDeserializer(struct5_T *r) {
  xilTargetDeserializer(&r->meas);
  e_xilTargetDeserializer(&r->status);
}

static void xilTargetDeserializer(double *r) {
  xilReadData((MemUnit_T *)r, sizeof(double));
}

static void xilTargetSerializer(const double *r) {
  xilWriteData((const MemUnit_T *)r, sizeof(double));
}

void XILTarget_initialize(unsigned int fcnId) {
  SD = getGNC_codegen_SILStackData();

  xilPreEntryPoint(fcnId);

  GNC_codegen_SIL_initialize(SD);

  xilPostEntryPoint(fcnId);
}

void XILTarget_terminate(unsigned int fcnId) {

  xilPreEntryPoint(fcnId);

  GNC_codegen_SIL_terminate();

  xilPostTerminatePoint(fcnId);
}

XIL_PROCESSDATA_ERROR_CODE
xilTarget_controller_codegen_entry(unsigned int fcnId) {
  struct0_T ctrl_mem;
  double xR[2];
  double b_time;
  double delta_encoder;
  double dt_ctrl;
  double pdyn;
  double r;
  double u_motor;
  boolean_T w_status_ctrl;
  SD = getGNC_codegen_SILStackData();

  xilTargetDeserializer(&b_time);

  xilTargetDeserializer(&dt_ctrl);

  b_xilTargetDeserializer(xR);

  xilTargetDeserializer(&pdyn);

  xilTargetDeserializer(&delta_encoder);

  c_xilTargetDeserializer(&ctrl_mem);

  xilPreEntryPoint(fcnId);

  controller_codegen_entry(SD, b_time, dt_ctrl, xR, pdyn, delta_encoder,
                           &ctrl_mem, &u_motor, &r, &w_status_ctrl);

  xilPostEntryPoint(fcnId);

  xilTargetSerializer(&u_motor);

  xilTargetSerializer(&r);

  b_xilTargetSerializer(&ctrl_mem);

  e_xilTargetSerializer(&w_status_ctrl);
  return XIL_PROCESSDATA_SUCCESS;
}

XIL_PROCESSDATA_ERROR_CODE
xilTarget_navigation_codegen_entry(unsigned int fcnId) {
  struct1_T bias;
  struct2_T sens_filt;
  struct3_T sens_in;
  double P[121];
  double x[11];
  double roll_state[2];
  double dt;
  double pdyn;
  boolean_T flight_phase;
  boolean_T w_status_nav;
  SD = getGNC_codegen_SILStackData();

  xilTargetDeserializer(&dt);

  e_xilTargetDeserializer(&flight_phase);

  f_xilTargetDeserializer(x);

  g_xilTargetDeserializer(P);

  h_xilTargetDeserializer(&bias);

  j_xilTargetDeserializer(&sens_filt);

  k_xilTargetDeserializer(&sens_in);

  xilPreEntryPoint(fcnId);

  navigation_codegen_entry(SD, dt, flight_phase, x, P, &bias, &sens_filt,
                           &sens_in, roll_state, &pdyn, &w_status_nav);

  xilPostEntryPoint(fcnId);

  f_xilTargetSerializer(x);

  g_xilTargetSerializer(P);

  h_xilTargetSerializer(&bias);

  j_xilTargetSerializer(&sens_filt);

  c_xilTargetSerializer(roll_state);

  xilTargetSerializer(&pdyn);

  e_xilTargetSerializer(&w_status_nav);
  return XIL_PROCESSDATA_SUCCESS;
}
