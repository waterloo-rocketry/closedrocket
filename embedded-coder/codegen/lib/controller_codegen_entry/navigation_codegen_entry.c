/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: navigation_codegen_entry.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

/* Include Files */
#include "navigation_codegen_entry.h"
#include "airdata_atmos.h"
#include "atan2.h"
#include "controller_codegen_entry_data.h"
#include "controller_codegen_entry_initialize.h"
#include "controller_codegen_entry_types.h"
#include "ekf_correct.h"
#include "eye.h"
#include "inv.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include "omp.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  double board_gyro[3];
  double mti_gyro[3];
  double ad_gyro[3];
  double board_mag_earth[3];
  double mti_mag_earth[3];
  double board_baro;
  double mti_baro;
} struct_T;
#endif /* typedef_struct_T */

/* Variable Definitions */
static boolean_T x_not_empty;

static double P[121];

static boolean_T b_not_empty;

/* Function Definitions */
/*
 * Calls the pad and flight filters.
 *
 * Arguments    : double dt
 *                boolean_T flight_phase
 *                const struct0_T *board_accel
 *                const struct0_T *board_gyro
 *                const struct0_T *mti_accel
 *                const struct0_T *mti_gyro
 *                const struct0_T *ad_accel
 *                const struct0_T *ad_gyro
 *                const struct1_T *board_baro
 *                const struct0_T *board_mag
 *                const struct1_T *mti_baro
 *                const struct0_T *mti_mag
 *                struct2_T *state
 *                double *cov_norm
 *                struct3_T *airdata
 *                double roll_state[2]
 * Return Type  : void
 */
void navigation_codegen_entry(
    double dt, boolean_T flight_phase, const struct0_T *board_accel,
    const struct0_T *board_gyro, const struct0_T *mti_accel,
    const struct0_T *mti_gyro, const struct0_T *ad_accel,
    const struct0_T *ad_gyro, const struct1_T *board_baro,
    const struct0_T *board_mag, const struct1_T *mti_baro,
    const struct0_T *mti_mag, struct2_T *state, double *cov_norm,
    struct3_T *airdata, double roll_state[2])
{
  static struct_T c_b;
  static const double Q[121] = {
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.01,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.01,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.01,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.001};
  static const double R[9] = {1.0E-9, 0.0, 0.0, 0.0,   1.0E-9,
                              0.0,    0.0, 0.0, 1.0E-9};
  static const double b_b[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  static double b_x[11];
  static double ad_accel_f[3];
  static double ad_gyro_f[3];
  static double board_accel_f[3];
  static double board_gyro_f[3];
  static double board_mag_f[3];
  static double mti_accel_f[3];
  static double mti_gyro_f[3];
  static double mti_mag_f[3];
  static double board_baro_f;
  static double mti_baro_f;
  __m128d b_r;
  __m128d r1;
  double F[121];
  double P_pred[121];
  double b_F[121];
  double K[33];
  double b_K[33];
  double W_dt[16];
  double aBuffer[16];
  double b_W_dt[16];
  double b_q[16];
  double c_W_dt[16];
  double dv2[16];
  double c_dt[12];
  double xhat[11];
  double ST[9];
  double b_a[9];
  double b_w_exp_tilde[9];
  double w_exp_tilde[9];
  double w_exp_tilde_tmp[9];
  double c_q[4];
  double q[4];
  double A[3];
  double a[3];
  double b_dt[3];
  double dn[3];
  double dv1[3];
  double n[3];
  double b_expl_temp;
  double d6;
  double expl_temp;
  double qw;
  double qy;
  double qz;
  double t1_density;
  int P_pred_tmp;
  int i;
  int i1;
  int i3;
  int i4;
  int i5;
  int i6;
  if (!isInitialized_controller_codegen_entry) {
    controller_codegen_entry_initialize();
  }
  /*     %% initialize at beginning */
  memset(&xhat[0], 0, 11U * sizeof(double));
  xhat[0] = 1.0;
  if (!x_not_empty) {
    memcpy(&b_x[0], &xhat[0], 11U * sizeof(double));
    x_not_empty = true;
  }
  /*     %% Pad filter iteration */
  if ((!flight_phase) || (!b_not_empty)) {
    /*  only before ignition (or if not run before) */
    /*  Computes on pad: inital state for flight filter, and bias values for the
     * sensors */
    /*  Outputs: initial state x, sensor biases bias */
    /*     %% Settings */
    /*  [s], low pass time constant */
    /*     %% parameters */
    /*     %% Initialization */
    /* %% remember filtered values from last iteration */
    if (!board_accel_f_not_empty) {
      board_accel_f_not_empty = true;
      board_baro_f = board_baro->meas;
      mti_baro_f = mti_baro->meas;
      board_accel_f[0] = board_accel->meas[0];
      board_gyro_f[0] = board_gyro->meas[0];
      mti_accel_f[0] = mti_accel->meas[0];
      mti_gyro_f[0] = mti_gyro->meas[0];
      ad_accel_f[0] = ad_accel->meas[0];
      ad_gyro_f[0] = ad_gyro->meas[0];
      board_mag_f[0] = board_mag->meas[0];
      mti_mag_f[0] = mti_mag->meas[0];
      board_accel_f[1] = board_accel->meas[1];
      board_gyro_f[1] = board_gyro->meas[1];
      mti_accel_f[1] = mti_accel->meas[1];
      mti_gyro_f[1] = mti_gyro->meas[1];
      ad_accel_f[1] = ad_accel->meas[1];
      ad_gyro_f[1] = ad_gyro->meas[1];
      board_mag_f[1] = board_mag->meas[1];
      mti_mag_f[1] = mti_mag->meas[1];
      board_accel_f[2] = board_accel->meas[2];
      board_gyro_f[2] = board_gyro->meas[2];
      mti_accel_f[2] = mti_accel->meas[2];
      mti_gyro_f[2] = mti_gyro->meas[2];
      ad_accel_f[2] = ad_accel->meas[2];
      ad_gyro_f[2] = ad_gyro->meas[2];
      board_mag_f[2] = board_mag->meas[2];
      mti_mag_f[2] = mti_mag->meas[2];
    }
    /*     %% lowpass filter */
    /*  filtered = filtered + alpha*(measured-filtered); */
    /* %% lowpass filter function used in pad filter */
    if (board_accel->status) {
      b_r = _mm_loadu_pd(&board_accel_f[0]);
      _mm_storeu_pd(&board_accel_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&board_accel->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      board_accel_f[2] =
          0.0005 * board_accel->meas[2] + 0.9995 * board_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (board_gyro->status) {
      b_r = _mm_loadu_pd(&board_gyro_f[0]);
      _mm_storeu_pd(&board_gyro_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&board_gyro->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      board_gyro_f[2] = 0.0005 * board_gyro->meas[2] + 0.9995 * board_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_accel->status) {
      b_r = _mm_loadu_pd(&mti_accel_f[0]);
      _mm_storeu_pd(&mti_accel_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&mti_accel->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      mti_accel_f[2] = 0.0005 * mti_accel->meas[2] + 0.9995 * mti_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_gyro->status) {
      b_r = _mm_loadu_pd(&mti_gyro_f[0]);
      _mm_storeu_pd(&mti_gyro_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&mti_gyro->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      mti_gyro_f[2] = 0.0005 * mti_gyro->meas[2] + 0.9995 * mti_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (ad_accel->status) {
      b_r = _mm_loadu_pd(&ad_accel_f[0]);
      _mm_storeu_pd(&ad_accel_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&ad_accel->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      ad_accel_f[2] = 0.0005 * ad_accel->meas[2] + 0.9995 * ad_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (ad_gyro->status) {
      b_r = _mm_loadu_pd(&ad_gyro_f[0]);
      _mm_storeu_pd(&ad_gyro_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&ad_gyro->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      ad_gyro_f[2] = 0.0005 * ad_gyro->meas[2] + 0.9995 * ad_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (board_baro->status) {
      board_baro_f = 0.0005 * board_baro->meas + 0.9995 * board_baro_f;
    }
    /* %% lowpass filter function used in pad filter */
    if (board_mag->status) {
      b_r = _mm_loadu_pd(&board_mag_f[0]);
      _mm_storeu_pd(&board_mag_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&board_mag->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      board_mag_f[2] = 0.0005 * board_mag->meas[2] + 0.9995 * board_mag_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_baro->status) {
      mti_baro_f = 0.0005 * mti_baro->meas + 0.9995 * mti_baro_f;
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_mag->status) {
      b_r = _mm_loadu_pd(&mti_mag_f[0]);
      _mm_storeu_pd(&mti_mag_f[0],
                    _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.0005),
                                          _mm_loadu_pd(&mti_mag->meas[0])),
                               _mm_mul_pd(_mm_set1_pd(0.9995), b_r)));
      mti_mag_f[2] = 0.0005 * mti_mag->meas[2] + 0.9995 * mti_mag_f[2];
    }
    /*     %% Initial state determination     */
    /* %% Orientation */
    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;
    /*  acceleration  */
    if (board_accel->status) {
      /*  only add alive IMUs to average */
      a[0] = board_accel_f[0];
      a[1] = board_accel_f[1];
      a[2] = board_accel_f[2];
    }
    if (mti_accel->status) {
      b_r = _mm_loadu_pd(&a[0]);
      r1 = _mm_loadu_pd(&mti_accel_f[0]);
      _mm_storeu_pd(&a[0], _mm_add_pd(b_r, r1));
      a[2] += mti_accel_f[2];
    }
    if (ad_accel->status) {
      b_r = _mm_loadu_pd(&a[0]);
      r1 = _mm_loadu_pd(&ad_accel_f[0]);
      _mm_storeu_pd(&a[0], _mm_add_pd(b_r, r1));
      a[2] += ad_accel_f[2];
    }
    /* %% computes initial orientation of stationary body from gravity
     * acceleration  */
    /* %% Input: 3D acceleration vector */
    /* %% Output: Orientation quaternion */
    /* %% normed gravity vector in body-fixed frame */
    qy = b_norm(a) + 1.0E-6;
    b_r = _mm_loadu_pd(&a[0]);
    _mm_storeu_pd(&A[0], _mm_div_pd(b_r, _mm_set1_pd(qy)));
    /*  unit vector of gravity direction */
    /* %% determine initial orientation quaternion */
    qw = sqrt(0.5 * A[0] + 0.5);
    if (qw == 0.0) {
      /*  exact upside down case */
      qy = 1.0;
      /*  either qy = 1 or qz = 1, this is arbitrary  */
      qz = 0.0;
    } else {
      qy = 0.5 * (a[2] / qy) / qw;
      qz = -0.5 * A[1] / qw;
    }
    q[0] = qw;
    q[1] = 0.0;
    q[2] = qy;
    q[3] = qz;
    /*  a gets normed inside function */
    /* %% launch altitude */
    /* %% set constant initials */
    /*  stationary on rail */
    /* %% conconct state vector */
    b_r = _mm_loadu_pd(&q[0]);
    r1 = _mm_set1_pd(c_norm(q));
    b_r = _mm_div_pd(b_r, r1);
    _mm_storeu_pd(&q[0], b_r);
    _mm_storeu_pd(&b_x[0], b_r);
    b_r = _mm_loadu_pd(&q[2]);
    b_r = _mm_div_pd(b_r, r1);
    _mm_storeu_pd(&q[2], b_r);
    _mm_storeu_pd(&b_x[2], b_r);
    b_x[10] = 420.0;
    /*     %% Bias determination */
    /* %% gyroscope */
    /* %% earth magnetic field */
    /*  computes rotation matrix from quaternion */
    /* %% norm quaternions */
    /*  q = quaternion_norm(q); */
    /* %% quaternion definition */
    /* %% skew symetric quaternion matrix */
    /* %% rotation matrix */
    b_x[4] = 0.0;
    b_x[7] = 0.0;
    c_b.board_gyro[0] = board_gyro_f[0];
    c_b.mti_gyro[0] = mti_gyro_f[0];
    c_b.ad_gyro[0] = ad_gyro_f[0];
    b_x[5] = 0.0;
    b_x[8] = 0.0;
    c_b.board_gyro[1] = board_gyro_f[1];
    c_b.mti_gyro[1] = mti_gyro_f[1];
    c_b.ad_gyro[1] = ad_gyro_f[1];
    b_x[6] = 0.0;
    b_x[9] = 0.0;
    c_b.board_gyro[2] = board_gyro_f[2];
    c_b.mti_gyro[2] = mti_gyro_f[2];
    c_b.ad_gyro[2] = ad_gyro_f[2];
    qy = q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + q[3] * q[3]);
    qw = 2.0 * q[0];
    /*  Stevens */
    /* %% for hardcoding:  */
    /*  [qw^2 + qx^2 - qy^2 - qz^2,         2*qw*qz + 2*qx*qy,         2*qx*qz -
     * 2*qw*qy] */
    /*  [        2*qx*qy - 2*qw*qz, qw^2 - qx^2 + qy^2 - qz^2,         2*qw*qx +
     * 2*qy*qz] */
    /*  [        2*qw*qy + 2*qx*qz,         2*qy*qz - 2*qw*qx, qw^2 - qx^2 -
     * qy^2 + qz^2] */
    for (i = 0; i < 3; i++) {
      qz = 2.0 * q[i + 1];
      ST[3 * i] = qy * b_b[i] + qz * q[1];
      ST[3 * i + 1] = qy * b_b[i + 3] + qz * q[2];
      ST[3 * i + 2] = qy * b_b[i + 6] + qz * q[3];
    }
    qy = qw * 0.0;
    b_a[0] = qy;
    b_a[1] = qw * -q[3];
    b_a[2] = qw * q[2];
    b_a[3] = qw * q[3];
    b_a[4] = qy;
    b_a[5] = qw * -q[1];
    b_a[6] = qw * -q[2];
    b_a[7] = qw * q[1];
    b_a[8] = qy;
    b_r = _mm_loadu_pd(&ST[0]);
    r1 = _mm_loadu_pd(&b_a[0]);
    _mm_storeu_pd(&ST[0], _mm_sub_pd(b_r, r1));
    b_r = _mm_loadu_pd(&ST[2]);
    r1 = _mm_loadu_pd(&b_a[2]);
    _mm_storeu_pd(&ST[2], _mm_sub_pd(b_r, r1));
    b_r = _mm_loadu_pd(&ST[4]);
    r1 = _mm_loadu_pd(&b_a[4]);
    _mm_storeu_pd(&ST[4], _mm_sub_pd(b_r, r1));
    b_r = _mm_loadu_pd(&ST[6]);
    r1 = _mm_loadu_pd(&b_a[6]);
    _mm_storeu_pd(&ST[6], _mm_sub_pd(b_r, r1));
    ST[8] -= qy;
    /*  launch attitude */
    c_b.board_mag_earth[0] = 0.0;
    c_b.board_mag_earth[1] = 0.0;
    c_b.board_mag_earth[2] = 0.0;
    b_r = _mm_loadu_pd(&ST[0]);
    r1 = _mm_loadu_pd(&c_b.board_mag_earth[0]);
    _mm_storeu_pd(&c_b.board_mag_earth[0],
                  _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(board_mag_f[0]))));
    c_b.board_mag_earth[2] += board_mag_f[0] * ST[2];
    c_b.mti_mag_earth[0] = 0.0;
    b_r = _mm_loadu_pd(&ST[3]);
    r1 = _mm_loadu_pd(&c_b.board_mag_earth[0]);
    _mm_storeu_pd(&c_b.board_mag_earth[0],
                  _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(board_mag_f[1]))));
    c_b.board_mag_earth[2] += board_mag_f[1] * ST[5];
    c_b.mti_mag_earth[1] = 0.0;
    b_r = _mm_loadu_pd(&ST[6]);
    r1 = _mm_loadu_pd(&c_b.board_mag_earth[0]);
    _mm_storeu_pd(&c_b.board_mag_earth[0],
                  _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(board_mag_f[2]))));
    c_b.board_mag_earth[2] += board_mag_f[2] * ST[8];
    c_b.mti_mag_earth[2] = 0.0;
    b_r = _mm_loadu_pd(&ST[0]);
    r1 = _mm_loadu_pd(&c_b.mti_mag_earth[0]);
    _mm_storeu_pd(&c_b.mti_mag_earth[0],
                  _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(mti_mag_f[0]))));
    c_b.mti_mag_earth[2] += mti_mag_f[0] * ST[2];
    b_r = _mm_loadu_pd(&ST[3]);
    r1 = _mm_loadu_pd(&c_b.mti_mag_earth[0]);
    _mm_storeu_pd(&c_b.mti_mag_earth[0],
                  _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(mti_mag_f[1]))));
    c_b.mti_mag_earth[2] += mti_mag_f[1] * ST[5];
    b_r = _mm_loadu_pd(&ST[6]);
    r1 = _mm_loadu_pd(&c_b.mti_mag_earth[0]);
    _mm_storeu_pd(&c_b.mti_mag_earth[0],
                  _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(mti_mag_f[2]))));
    c_b.mti_mag_earth[2] += mti_mag_f[2] * ST[8];
    /* %% barometer */
    qy = airdata_atmos(420.0, &qy, &t1_density, &qz, &qw, &expl_temp);
    /*  what the pressure should be at launch elevation */
    c_b.board_baro = board_baro_f - qy;
    c_b.mti_baro = mti_baro_f - qy;
    b_not_empty = true;
  }
  /*     %% Flight filter iteration */
  if (flight_phase) {
    __m128d r2;
    __m128d r3;
    __m128d r4;
    __m128d r5;
    __m128d r6;
    __m128d r7;
    __m128d r8;
    __m128d r9;
    double dv[4];
    double A_tmp_tmp;
    double C_total_w_idx_0;
    double C_total_w_idx_1;
    double b;
    double b_A_tmp_tmp;
    double c_A_tmp_tmp;
    double c_w_exp_tilde_tmp;
    double d;
    double d1;
    double d2;
    double d3;
    double d4;
    double d5;
    double dphi;
    double dphi_tmp;
    double x;
    int ST_tmp;
    int b_w_exp_tilde_tmp;
    int d_w_exp_tilde_tmp;
    /*  in flight */
    /*  Computes state in flight */
    /*  Input variables: time step dt, old state x, old covariance P,  */
    /*  Input parameters: sensor biases bias; */
    /*  Input measurements: sensorgroup_sensortype */
    /*  Outputs: new state x, new covariance P */
    /*     %% IMU Prediction + Correction steps */
    /* %% x = [   q(4),           w(3),         v(3),    alt(1)] */
    /* %% Q is a square 11 matrix, tuning for expected dynamics noise magnitude
     * E(noise) */
    /* %% R is a square 3 matrix, tuning for expected measurement noise
     * magnitude E(noise) of the gyroscope */
    /* %% computes average acceleration and angular rates from multiple IMUs. */
    /* %% includes correction of gyroscope bias and centrifugal acceleration. */
    /*     %% confidences */
    /* %% base confidences (tune per sensor) */
    /*  use accelerometer bias standard deviation */
    /*  use gyroscope noise standard deviation */
    /*     %% confidence calculations */
    /*  sensor status */
    /*  normalize (Hadamard division) */
    qy = 9.9999999999999981E+9 * (double)ad_gyro->status;
    A_tmp_tmp = 1.0000000000000002E+14 * (double)board_accel->status;
    b_A_tmp_tmp = 1.0000000000000002E+14 * (double)mti_accel->status;
    c_A_tmp_tmp = 1.0000000000000002E+14 * (double)ad_accel->status;
    qz = (A_tmp_tmp + b_A_tmp_tmp) + c_A_tmp_tmp;
    A[0] = qz;
    d = 9.9999999999999981E+9 * (double)board_gyro->status;
    d1 = 9.9999999999999981E+9 * (double)mti_gyro->status;
    C_total_w_idx_1 = d + d1;
    C_total_w_idx_0 = C_total_w_idx_1 + qy;
    qy /= C_total_w_idx_0;
    a[0] = qy;
    A[1] = qz;
    qy = 0.0 / C_total_w_idx_1;
    a[1] = qy;
    A[2] = qz;
    a[2] = qy;
    /*     %% parameters */
    /*     %% averaging */
    /*  gyro bias correction */
    /*  weighted angular rates */
    /*  [rad/s] */
    /*  centrifugal acceleration correction */
    /*  w_tilde = math_tilde(w); */
    /*  w_tilde_sq = w_tilde * w_tilde; */
    /*  - w_tilde_sq * param.d_board; */
    /*  - w_tilde_sq * param.d_mti; */
    /*  - w_tilde_sq * param.d_ad; */
    /*  weighted acceleration */
    /*  [m/s^2] */
    /*  Computes EKF prediction+correction step for IMU data. */
    /*  Input variables: time step dt, old state x, old covariance P;
     * acceleration a, angular rate w;  */
    /*  Input parameters: dynamics weighting Q, gyroscope weighting R;   */
    /*  Outputs: new state x, new covariance P */
    /*     %% Prediction */
    /*  computes a-priori state and covariance estimates */
    /*  Uses discrete-time dynamics and analytical Jacobian */
    /* %% discrete dynamics update */
    /*  Computes state update with dynamics model and time integration */
    /*     %% decomp */
    /*  decompose state vector: [q(4); w(3); v(3); alt] */
    /*     %% load parameters */
    /*     %% time updates */
    /*  quaternion update */
    /*  computes new quaternion from old quaternion and body rates */
    /* %% norm quaternions */
    /*  norms quaternion */
    /* %% inverse quaternion  */
    b_r = _mm_loadu_pd(&b_x[0]);
    r1 = _mm_set1_pd(c_norm(&b_x[0]));
    _mm_storeu_pd(&q[0], _mm_div_pd(b_r, r1));
    b_r = _mm_loadu_pd(&b_x[2]);
    _mm_storeu_pd(&q[2], _mm_div_pd(b_r, r1));
    /* %% incremental quaternion */
    dphi_tmp = b_norm(&b_x[4]);
    qz = dphi_tmp * dt;
    dphi = qz / 2.0;
    if (dphi_tmp == 0.0) {
      dn[0] = 0.0;
      dn[1] = 0.0;
      dn[2] = 0.0;
      n[0] = 0.0;
      n[1] = 0.0;
      n[2] = 0.0;
    } else {
      b_r = _mm_loadu_pd(&b_x[4]);
      r1 = _mm_set1_pd(dphi_tmp);
      _mm_storeu_pd(&dn[0], _mm_div_pd(b_r, r1));
      qy = b_x[6] / dphi_tmp;
      dn[2] = qy;
      b_r = _mm_loadu_pd(&b_x[4]);
      _mm_storeu_pd(&n[0], _mm_div_pd(b_r, r1));
      n[2] = qy;
    }
    b = sin(dphi);
    /* %% quaternion update */
    /*  quaternion multiplication */
    /*  Quaternion product matrix */
    /* %% product */
    /*  rate update */
    /*  computes matrix exponential of rotation */
    /* %% incremental angle */
    /* %% normed skew-symmetric matrix */
    /*  skew symmetric matrix / cross-product jacobian */
    ST[0] = 0.0;
    ST[3] = -n[2];
    ST[6] = n[1];
    ST[1] = n[2];
    ST[4] = 0.0;
    ST[7] = -n[0];
    ST[2] = -n[1];
    ST[5] = n[0];
    ST[8] = 0.0;
    /* %% Rodrigues rotation formula */
    qw = sin(qz);
    x = cos(qz);
    b_eye(w_exp_tilde_tmp);
    memset(&w_exp_tilde[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qy = ST[3 * i];
      b_r = _mm_loadu_pd(&ST[0]);
      r1 = _mm_loadu_pd(&w_exp_tilde[3 * i]);
      _mm_storeu_pd(&w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_w_exp_tilde_tmp = 3 * i + 2;
      w_exp_tilde[b_w_exp_tilde_tmp] += ST[2] * qy;
      qy = ST[3 * i + 1];
      b_r = _mm_loadu_pd(&ST[3]);
      r1 = _mm_loadu_pd(&w_exp_tilde[3 * i]);
      _mm_storeu_pd(&w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      w_exp_tilde[b_w_exp_tilde_tmp] += ST[5] * qy;
      qy = ST[b_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&ST[6]);
      r1 = _mm_loadu_pd(&w_exp_tilde[3 * i]);
      _mm_storeu_pd(&w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      w_exp_tilde[b_w_exp_tilde_tmp] += 0.0 * qy;
    }
    b_r = _mm_loadu_pd(&ST[0]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[0]);
    r2 = _mm_loadu_pd(&w_exp_tilde[0]);
    r3 = _mm_set1_pd(qw);
    r4 = _mm_set1_pd(1.0 - x);
    _mm_storeu_pd(
        &w_exp_tilde[0],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    b_r = _mm_loadu_pd(&ST[2]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[2]);
    r2 = _mm_loadu_pd(&w_exp_tilde[2]);
    _mm_storeu_pd(
        &w_exp_tilde[2],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    b_r = _mm_loadu_pd(&ST[4]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[4]);
    r2 = _mm_loadu_pd(&w_exp_tilde[4]);
    _mm_storeu_pd(
        &w_exp_tilde[4],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    b_r = _mm_loadu_pd(&ST[6]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[6]);
    r2 = _mm_loadu_pd(&w_exp_tilde[6]);
    _mm_storeu_pd(
        &w_exp_tilde[6],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    c_w_exp_tilde_tmp = w_exp_tilde_tmp[8] - qw * 0.0;
    w_exp_tilde[8] = c_w_exp_tilde_tmp + (1.0 - x) * w_exp_tilde[8];
    /*  aerodynamics model */
    /* %% air data   */
    /*  appends airadata struct with dynamic air data from velocity vector */
    /*  dynamic air data: dynamic pressure, mach number */
    qy = b_norm(&b_x[7]);
    /*  return values */
    /*  [-], Mach number */
    /*  [Pa] */
    airdata_atmos(b_x[10], &qz, &t1_density, &qw, &expl_temp, &b_expl_temp);
    qz = 0.5 * t1_density * (qy * qy);
    /* %% angle of attack/sideslip */
    /* %% torques */
    /* torque_canards = Cl *  delta * param.c_canard * p_dyn *[1;0;0]; */
    /* + param.Cn_omega*[0; w(2); w(3)] ) * param.c_aero; % commented */
    /*  out because timeline */
    /* torque = torque_canards + torque_aero; */
    /*  w_new = w + dt * param.Jinv * (torque - cross(w, param.J*w)); % old
     * version */
    /* %% possibly more accurate: for Jx < Jy = Jz : u = (Jy-Jx)/Jy * wx, and */
    /* %% wx_new = wx, [wy, wz]_new = Sx(u*dt)*[wy, wz] with Sx = [c(u), s(u);
     * -s(u), c(u)] */
    /*  velocity update  */
    /*  computes rotation matrix from quaternion */
    /* %% norm quaternions */
    /*  q = quaternion_norm(q); */
    /* %% quaternion definition */
    /* %% skew symetric quaternion matrix */
    /* %% rotation matrix */
    qw = b_x[0] * b_x[0] -
         ((b_x[1] * b_x[1] + b_x[2] * b_x[2]) + b_x[3] * b_x[3]);
    qy = 2.0 * b_x[0];
    expl_temp = qy * 0.0;
    b_a[0] = expl_temp;
    b_a[3] = qy * -b_x[3];
    b_a[6] = qy * b_x[2];
    b_a[1] = qy * b_x[3];
    b_a[4] = expl_temp;
    b_a[7] = qy * -b_x[1];
    b_a[2] = qy * -b_x[2];
    b_a[5] = qy * b_x[1];
    b_a[8] = expl_temp;
    r2 = _mm_loadu_pd(&b_a[0]);
    /*  Stevens */
    /* %% for hardcoding:  */
    /*  [qw^2 + qx^2 - qy^2 - qz^2,         2*qw*qz + 2*qx*qy,         2*qx*qz -
     * 2*qw*qy] */
    /*  [        2*qx*qy - 2*qw*qz, qw^2 - qx^2 + qy^2 - qz^2,         2*qw*qx +
     * 2*qy*qz] */
    /*  [        2*qw*qy + 2*qx*qz,         2*qy*qz - 2*qw*qx, qw^2 - qx^2 -
     * qy^2 + qz^2] */
    /*  inertial to body frame */
    /*  v_new = v + dt * (a - cross(w,v) + S*param.g); % old version */
    /*  altitude update */
    /*     %% concoct state update vector */
    b_q[0] = q[0];
    b_q[4] = -q[1];
    b_q[8] = -q[2];
    b_q[12] = -q[3];
    b_q[1] = q[1];
    b_q[5] = q[0];
    b_q[9] = -q[3];
    b_q[13] = q[2];
    b_q[2] = q[2];
    b_q[6] = q[3];
    b_q[10] = q[0];
    b_q[14] = -q[1];
    b_q[3] = q[3];
    b_q[7] = -q[2];
    b_q[11] = q[1];
    b_q[15] = q[0];
    dv[0] = cos(dphi);
    memset(&b_w_exp_tilde[0], 0, 9U * sizeof(double));
    memset(&n[0], 0, 3U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qy = b_x[i + 1];
      ST[3 * i] = qw * b_b[3 * i] + 2.0 * b_x[1] * qy;
      b_w_exp_tilde_tmp = 3 * i + 1;
      ST[b_w_exp_tilde_tmp] = qw * b_b[b_w_exp_tilde_tmp] + 2.0 * b_x[2] * qy;
      ST_tmp = 3 * i + 2;
      ST[ST_tmp] = qw * b_b[ST_tmp] + 2.0 * b_x[3] * qy;
      dv[i + 1] = dn[i] * b;
      qy = param.J[3 * i];
      b_r = _mm_loadu_pd(&w_exp_tilde[0]);
      r1 = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      _mm_storeu_pd(&b_w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_w_exp_tilde[ST_tmp] += w_exp_tilde[2] * qy;
      qy = param.J[b_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&w_exp_tilde[3]);
      r1 = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      _mm_storeu_pd(&b_w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_w_exp_tilde[ST_tmp] += w_exp_tilde[5] * qy;
      qy = param.J[ST_tmp];
      b_r = _mm_loadu_pd(&w_exp_tilde[6]);
      r1 = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      _mm_storeu_pd(&b_w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_w_exp_tilde[ST_tmp] += w_exp_tilde[8] * qy;
      b_r = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      r1 = _mm_loadu_pd(&n[0]);
      qy = b_x[i + 4];
      _mm_storeu_pd(&n[0], _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      n[2] += b_w_exp_tilde[ST_tmp] * qy;
    }
    b_r = _mm_loadu_pd(&ST[0]);
    _mm_storeu_pd(&ST[0], _mm_sub_pd(b_r, r2));
    b_r = _mm_loadu_pd(&ST[2]);
    r2 = _mm_loadu_pd(&b_a[2]);
    _mm_storeu_pd(&ST[2], _mm_sub_pd(b_r, r2));
    b_r = _mm_loadu_pd(&ST[4]);
    r2 = _mm_loadu_pd(&b_a[4]);
    _mm_storeu_pd(&ST[4], _mm_sub_pd(b_r, r2));
    b_r = _mm_loadu_pd(&ST[6]);
    r2 = _mm_loadu_pd(&b_a[6]);
    _mm_storeu_pd(&ST[6], _mm_sub_pd(b_r, r2));
    ST[8] -= expl_temp;
    dn[0] = qz * -0.0;
    dn[1] = qz * (-0.16182736457722724 * sin(b_atan2(b_x[9], b_x[7])));
    dn[2] = qz * (-0.16182736457722724 * -sin(b_atan2(b_x[8], b_x[7])));
    memset(&dv1[0], 0, 3U * sizeof(double));
    memset(&b_dt[0], 0, 3U * sizeof(double));
    expl_temp = dv1[0];
    dphi = dv1[1];
    d2 = dv1[2];
    b = b_dt[0];
    b_expl_temp = b_dt[1];
    d3 = b_dt[2];
    d4 = b_x[7];
    t1_density = b_x[8];
    d5 = b_x[9];
    for (i = 0; i < 3; i++) {
      qy = param.Jinv[3 * i];
      qz = n[i];
      expl_temp += qy * qz;
      qw = dn[i];
      b += dt * qy * qw;
      qy = param.Jinv[3 * i + 1];
      dphi += qy * qz;
      b_expl_temp += dt * qy * qw;
      qy = param.Jinv[3 * i + 2];
      d2 += qy * qz;
      d3 += dt * qy * qw;
      qy = A[i];
      n[i] = ((w_exp_tilde[i] * d4 + w_exp_tilde[i + 3] * t1_density) +
              w_exp_tilde[i + 6] * d5) +
             dt * ((A_tmp_tmp / qy * board_accel->meas[i] +
                    b_A_tmp_tmp / qy * mti_accel->meas[i]) +
                   c_A_tmp_tmp / qy * ad_accel->meas[i]);
    }
    b_dt[2] = d3;
    b_dt[1] = b_expl_temp;
    b_dt[0] = b;
    dv1[2] = d2;
    dv1[1] = dphi;
    dv1[0] = expl_temp;
    memset(&dn[0], 0, 3U * sizeof(double));
    qy = dn[0];
    qz = dn[1];
    qw = dn[2];
    expl_temp = b_x[7];
    dphi = b_x[8];
    b = b_x[9];
    for (i = 0; i < 3; i++) {
      b_expl_temp = ST[3 * i];
      t1_density = param.g[i];
      qy += dt * b_expl_temp * t1_density;
      d4 = b_expl_temp * expl_temp;
      b_expl_temp = ST[3 * i + 1];
      qz += dt * b_expl_temp * t1_density;
      d4 += b_expl_temp * dphi;
      b_expl_temp = ST[3 * i + 2];
      qw += dt * b_expl_temp * t1_density;
      d4 += b_expl_temp * b;
      A[i] = d4;
    }
    dn[2] = qw;
    dn[1] = qz;
    dn[0] = qy;
    memset(&c_q[0], 0, sizeof(double) << 2);
    for (i = 0; i < 4; i++) {
      b_w_exp_tilde_tmp = i << 2;
      b_r = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp]);
      r1 = _mm_loadu_pd(&c_q[0]);
      r2 = _mm_set1_pd(dv[i]);
      _mm_storeu_pd(&c_q[0], _mm_add_pd(r1, _mm_mul_pd(b_r, r2)));
      b_r = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp + 2]);
      r1 = _mm_loadu_pd(&c_q[2]);
      _mm_storeu_pd(&c_q[2], _mm_add_pd(r1, _mm_mul_pd(b_r, r2)));
    }
    xhat[0] = c_q[0];
    xhat[1] = c_q[1];
    xhat[2] = c_q[2];
    xhat[3] = c_q[3];
    b_r = _mm_loadu_pd(&dv1[0]);
    r1 = _mm_loadu_pd(&b_dt[0]);
    _mm_storeu_pd(&xhat[4], _mm_add_pd(b_r, r1));
    b_r = _mm_loadu_pd(&n[0]);
    r1 = _mm_loadu_pd(&dn[0]);
    _mm_storeu_pd(&xhat[7], _mm_add_pd(b_r, r1));
    xhat[6] = d2 + d3;
    xhat[9] = n[2] + qw;
    xhat[10] = b_x[10] + dt * A[0];
    /* %% discrete Jacobian: F = df/dx */
    /*  Computes state derivative with predictive model. Use ODE solver to
     * compute next state. */
    /*     %% decomp */
    /*  decompose state vector: [q(4); w(3); v(3)] */
    /*     %% load parameters */
    /*     %% create empty Jacobian  */
    memset(&F[0], 0, 121U * sizeof(double));
    /*  could also initialize as identity eye(length(x)), as right now all */
    /*  sub-Jacobians have identity on the main diagonal */
    /*     %% quaternion attitude rows (q, 1:4) */
    /*  computes quaternion derivative from quaternion and body rates */
    /*  norm quaternions */
    /*  norms quaternion */
    /* %% inverse quaternion  */
    /*  Quaternion product matrix */
    /*  Big Omega matrix */
    /*  quaternion derivative Jacobians */
    t1_density = 0.5 * dt;
    qz = t1_density * 0.0;
    W_dt[0] = qz;
    qw = t1_density * -b_x[4];
    W_dt[4] = qw;
    qy = t1_density * -b_x[5];
    W_dt[8] = qy;
    expl_temp = t1_density * -b_x[6];
    W_dt[12] = expl_temp;
    b_expl_temp = t1_density * b_x[4];
    W_dt[1] = b_expl_temp;
    W_dt[5] = qz;
    dphi = t1_density * b_x[6];
    W_dt[9] = dphi;
    W_dt[13] = qy;
    qy = t1_density * b_x[5];
    W_dt[2] = qy;
    W_dt[6] = expl_temp;
    W_dt[10] = qz;
    W_dt[14] = b_expl_temp;
    W_dt[3] = dphi;
    W_dt[7] = qy;
    W_dt[11] = qw;
    W_dt[15] = qz;
    memset(&aBuffer[0], 0, sizeof(double) << 4);
    c_eye(dv2);
    memset(&b_W_dt[0], 0, sizeof(double) << 4);
    memset(&b_q[0], 0, sizeof(double) << 4);
    memset(&c_W_dt[0], 0, sizeof(double) << 4);
    for (i1 = 0; i1 < 4; i1++) {
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        ST_tmp = i << 2;
        b_r = _mm_loadu_pd(&W_dt[ST_tmp]);
        r2 = _mm_loadu_pd(&aBuffer[b_w_exp_tilde_tmp]);
        r8 = _mm_set1_pd(W_dt[i + b_w_exp_tilde_tmp]);
        r1 = _mm_mul_pd(b_r, r8);
        _mm_storeu_pd(&aBuffer[b_w_exp_tilde_tmp], _mm_add_pd(r2, r1));
        b_r = _mm_loadu_pd(&b_W_dt[b_w_exp_tilde_tmp]);
        _mm_storeu_pd(&b_W_dt[b_w_exp_tilde_tmp], _mm_add_pd(b_r, r1));
        b_r = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp]);
        _mm_storeu_pd(&b_q[b_w_exp_tilde_tmp], _mm_add_pd(b_r, r1));
        b_r = _mm_loadu_pd(&W_dt[ST_tmp + 2]);
        r2 = _mm_loadu_pd(&aBuffer[b_w_exp_tilde_tmp + 2]);
        r1 = _mm_mul_pd(b_r, r8);
        _mm_storeu_pd(&aBuffer[b_w_exp_tilde_tmp + 2], _mm_add_pd(r2, r1));
        b_r = _mm_loadu_pd(&b_W_dt[b_w_exp_tilde_tmp + 2]);
        _mm_storeu_pd(&b_W_dt[b_w_exp_tilde_tmp + 2], _mm_add_pd(b_r, r1));
        b_r = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp + 2]);
        _mm_storeu_pd(&b_q[b_w_exp_tilde_tmp + 2], _mm_add_pd(b_r, r1));
      }
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        ST_tmp = i << 2;
        b_r = _mm_loadu_pd(&W_dt[ST_tmp]);
        r1 = _mm_loadu_pd(&c_W_dt[b_w_exp_tilde_tmp]);
        r2 = _mm_set1_pd(b_q[i + b_w_exp_tilde_tmp]);
        _mm_storeu_pd(&c_W_dt[b_w_exp_tilde_tmp],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, r2)));
        b_r = _mm_loadu_pd(&W_dt[ST_tmp + 2]);
        r1 = _mm_loadu_pd(&c_W_dt[b_w_exp_tilde_tmp + 2]);
        _mm_storeu_pd(&c_W_dt[b_w_exp_tilde_tmp + 2],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, r2)));
      }
    }
    memset(&b_q[0], 0, sizeof(double) << 4);
    r5 = _mm_set1_pd(0.5);
    r6 = _mm_set1_pd(0.16666666666666666);
    r7 = _mm_set1_pd(0.041666666666666664);
    for (i1 = 0; i1 < 4; i1++) {
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        ST_tmp = i << 2;
        b_r = _mm_loadu_pd(&aBuffer[ST_tmp]);
        r1 = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp]);
        r2 = _mm_set1_pd(aBuffer[i + b_w_exp_tilde_tmp]);
        _mm_storeu_pd(&b_q[b_w_exp_tilde_tmp],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, r2)));
        b_r = _mm_loadu_pd(&aBuffer[ST_tmp + 2]);
        r1 = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp + 2]);
        _mm_storeu_pd(&b_q[b_w_exp_tilde_tmp + 2],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, r2)));
      }
      b_w_exp_tilde_tmp = i1 << 2;
      b_r = _mm_loadu_pd(&dv2[b_w_exp_tilde_tmp]);
      r1 = _mm_loadu_pd(&W_dt[b_w_exp_tilde_tmp]);
      r2 = _mm_loadu_pd(&b_W_dt[b_w_exp_tilde_tmp]);
      r8 = _mm_loadu_pd(&c_W_dt[b_w_exp_tilde_tmp]);
      r9 = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp]);
      _mm_storeu_pd(&F[11 * i1],
                    _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(b_r, r1),
                                                     _mm_mul_pd(r5, r2)),
                                          _mm_mul_pd(r6, r8)),
                               _mm_mul_pd(r7, r9)));
      b_r = _mm_loadu_pd(&dv2[b_w_exp_tilde_tmp + 2]);
      r1 = _mm_loadu_pd(&W_dt[b_w_exp_tilde_tmp + 2]);
      r2 = _mm_loadu_pd(&b_W_dt[b_w_exp_tilde_tmp + 2]);
      r8 = _mm_loadu_pd(&c_W_dt[b_w_exp_tilde_tmp + 2]);
      r9 = _mm_loadu_pd(&b_q[b_w_exp_tilde_tmp + 2]);
      _mm_storeu_pd(&F[11 * i1 + 2],
                    _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(b_r, r1),
                                                     _mm_mul_pd(r5, r2)),
                                          _mm_mul_pd(r6, r8)),
                               _mm_mul_pd(r7, r9)));
    }
    /*  + 1/2 * Q_dt^2 + 1/6 * Q_dt^3 + 1/24 * Q_dt^4; */
    qz = t1_density * q[0];
    b_q[0] = qz;
    qw = t1_density * -q[1];
    b_q[4] = qw;
    expl_temp = t1_density * -q[2];
    b_q[8] = expl_temp;
    qy = t1_density * -q[3];
    b_q[12] = qy;
    b = t1_density * q[1];
    b_q[1] = b;
    b_q[5] = qz;
    b_q[9] = qy;
    qy = t1_density * q[2];
    b_q[13] = qy;
    b_q[2] = qy;
    qy = t1_density * q[3];
    b_q[6] = qy;
    b_q[10] = qz;
    b_q[14] = qw;
    b_q[3] = qy;
    b_q[7] = expl_temp;
    b_q[11] = b;
    b_q[15] = qz;
    for (i = 0; i < 3; i++) {
      b_w_exp_tilde_tmp = (i + 1) << 2;
      ST_tmp = 11 * (i + 4);
      F[ST_tmp] = b_q[b_w_exp_tilde_tmp];
      F[ST_tmp + 1] = b_q[b_w_exp_tilde_tmp + 1];
      F[ST_tmp + 2] = b_q[b_w_exp_tilde_tmp + 2];
      F[ST_tmp + 3] = b_q[b_w_exp_tilde_tmp + 3];
    }
    /*  column q (attitude) */
    /*  column w (rates) */
    /*     %% angular rate rows (w, 5:7) */
    /*  aerodynamics partial derivatives */
    /* %% air data  */
    /* torque_vx = Cl * delta * param.c_canard * [v(1), v(2), v(3); 0, 0, 0; 0,
     * 0, 0]; */
    /* torque_v =  airdata.density * (torque_vx + torque_vyz); */
    airdata_atmos(b_x[10], &qy, &t1_density, &qz, &qw, &expl_temp);
    /*  computes matrix exponential of rotation */
    /* %% incremental angle */
    /* %% normed skew-symmetric matrix */
    if (dphi_tmp == 0.0) {
      n[0] = 0.0;
      n[1] = 0.0;
      n[2] = 0.0;
    } else {
      b_r = _mm_loadu_pd(&b_x[4]);
      _mm_storeu_pd(&n[0], _mm_div_pd(b_r, _mm_set1_pd(dphi_tmp)));
      n[2] = b_x[6] / dphi_tmp;
    }
    /*  skew symmetric matrix / cross-product jacobian */
    ST[0] = 0.0;
    ST[3] = -n[2];
    ST[6] = n[1];
    ST[1] = n[2];
    ST[4] = 0.0;
    ST[7] = -n[0];
    ST[2] = -n[1];
    ST[5] = n[0];
    ST[8] = 0.0;
    /* %% Rodrigues rotation formula */
    memset(&w_exp_tilde[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qy = ST[3 * i];
      b_r = _mm_loadu_pd(&ST[0]);
      r1 = _mm_loadu_pd(&w_exp_tilde[3 * i]);
      _mm_storeu_pd(&w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_w_exp_tilde_tmp = 3 * i + 2;
      w_exp_tilde[b_w_exp_tilde_tmp] += ST[2] * qy;
      qy = ST[3 * i + 1];
      b_r = _mm_loadu_pd(&ST[3]);
      r1 = _mm_loadu_pd(&w_exp_tilde[3 * i]);
      _mm_storeu_pd(&w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      w_exp_tilde[b_w_exp_tilde_tmp] += ST[5] * qy;
      qy = ST[b_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&ST[6]);
      r1 = _mm_loadu_pd(&w_exp_tilde[3 * i]);
      _mm_storeu_pd(&w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      w_exp_tilde[b_w_exp_tilde_tmp] += 0.0 * qy;
    }
    b_r = _mm_loadu_pd(&ST[0]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[0]);
    r2 = _mm_loadu_pd(&w_exp_tilde[0]);
    _mm_storeu_pd(
        &w_exp_tilde[0],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    b_r = _mm_loadu_pd(&ST[2]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[2]);
    r2 = _mm_loadu_pd(&w_exp_tilde[2]);
    _mm_storeu_pd(
        &w_exp_tilde[2],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    b_r = _mm_loadu_pd(&ST[4]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[4]);
    r2 = _mm_loadu_pd(&w_exp_tilde[4]);
    _mm_storeu_pd(
        &w_exp_tilde[4],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    b_r = _mm_loadu_pd(&ST[6]);
    r1 = _mm_loadu_pd(&w_exp_tilde_tmp[6]);
    r2 = _mm_loadu_pd(&w_exp_tilde[6]);
    _mm_storeu_pd(
        &w_exp_tilde[6],
        _mm_add_pd(_mm_sub_pd(r1, _mm_mul_pd(r3, b_r)), _mm_mul_pd(r4, r2)));
    w_exp_tilde[8] = c_w_exp_tilde_tmp + (1.0 - x) * w_exp_tilde[8];
    memset(&b_a[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qy = w_exp_tilde[3 * i];
      b_r = _mm_loadu_pd(&b_param.Jinv[0]);
      r1 = _mm_loadu_pd(&b_a[3 * i]);
      _mm_storeu_pd(&b_a[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_w_exp_tilde_tmp = 3 * i + 2;
      b_a[b_w_exp_tilde_tmp] += b_param.Jinv[2] * qy;
      ST_tmp = 11 * (i + 4);
      F[ST_tmp + 4] = 0.0;
      qy = w_exp_tilde[3 * i + 1];
      b_r = _mm_loadu_pd(&b_param.Jinv[3]);
      r1 = _mm_loadu_pd(&b_a[3 * i]);
      _mm_storeu_pd(&b_a[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_a[b_w_exp_tilde_tmp] += b_param.Jinv[5] * qy;
      F[ST_tmp + 5] = 0.0;
      qy = w_exp_tilde[b_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&b_param.Jinv[6]);
      r1 = _mm_loadu_pd(&b_a[3 * i]);
      _mm_storeu_pd(&b_a[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_a[b_w_exp_tilde_tmp] += b_param.Jinv[8] * qy;
      F[ST_tmp + 6] = 0.0;
    }
    b_w_exp_tilde[1] = t1_density * (-0.08091368228861362 * b_x[9]);
    qz = t1_density * -0.0;
    b_w_exp_tilde[4] = qz;
    b_w_exp_tilde[7] = t1_density * (-0.08091368228861362 * b_x[7]);
    b_w_exp_tilde[2] = t1_density * (-0.08091368228861362 * -b_x[8]);
    b_w_exp_tilde[5] = t1_density * (-0.08091368228861362 * -b_x[7]);
    b_w_exp_tilde[8] = qz;
    /*  column w */
    /*  column v */
    /*     %% velocity rows (v, 8:10) */
    /*  computes rotation matrix from quaternion */
    /* %% norm quaternions */
    /*  q = quaternion_norm(q); */
    /* %% quaternion definition */
    /* %% rotation matrix */
    /*  S = (qw^2-qv'*qv)*eye(3) + 2*qv*qv' - 2*qw*q_tilde; % Stevens */
    /* %% Jacobian of rotation wrt quaternion */
    qw = 0.0;
    r2 = _mm_set1_pd(dt);
    for (i = 0; i < 3; i++) {
      int i2;
      b_w_exp_tilde[3 * i] = qz;
      qy = b_param.J[3 * i];
      b_r = _mm_loadu_pd(&b_a[0]);
      b_w_exp_tilde_tmp = 11 * (i + 4);
      r1 = _mm_loadu_pd(&F[b_w_exp_tilde_tmp + 4]);
      _mm_storeu_pd(&F[b_w_exp_tilde_tmp + 4],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      F[b_w_exp_tilde_tmp + 6] += b_a[2] * qy;
      ST_tmp = 11 * (i + 7);
      F[ST_tmp + 4] = 0.0;
      d_w_exp_tilde_tmp = 3 * i + 1;
      qy = b_param.J[d_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&b_a[3]);
      r1 = _mm_loadu_pd(&F[b_w_exp_tilde_tmp + 4]);
      _mm_storeu_pd(&F[b_w_exp_tilde_tmp + 4],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      F[b_w_exp_tilde_tmp + 6] += b_a[5] * qy;
      F[ST_tmp + 5] = 0.0;
      i2 = 3 * i + 2;
      qy = b_param.J[i2];
      b_r = _mm_loadu_pd(&b_a[6]);
      r1 = _mm_loadu_pd(&F[b_w_exp_tilde_tmp + 4]);
      _mm_storeu_pd(&F[b_w_exp_tilde_tmp + 4],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      F[b_w_exp_tilde_tmp + 6] += b_a[8] * qy;
      F[ST_tmp + 6] = 0.0;
      b_r = _mm_loadu_pd(&b_param.Jinv[0]);
      r1 = _mm_loadu_pd(&F[ST_tmp + 4]);
      _mm_storeu_pd(
          &F[ST_tmp + 4],
          _mm_add_pd(r1, _mm_mul_pd(_mm_mul_pd(r2, b_r), _mm_set1_pd(qz))));
      F[ST_tmp + 6] += dt * b_param.Jinv[2] * qz;
      qy = b_w_exp_tilde[d_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&b_param.Jinv[3]);
      r1 = _mm_loadu_pd(&F[ST_tmp + 4]);
      _mm_storeu_pd(
          &F[ST_tmp + 4],
          _mm_add_pd(r1, _mm_mul_pd(_mm_mul_pd(r2, b_r), _mm_set1_pd(qy))));
      F[ST_tmp + 6] += dt * b_param.Jinv[5] * qy;
      qy = b_w_exp_tilde[i2];
      b_r = _mm_loadu_pd(&b_param.Jinv[6]);
      r1 = _mm_loadu_pd(&F[ST_tmp + 4]);
      _mm_storeu_pd(
          &F[ST_tmp + 4],
          _mm_add_pd(r1, _mm_mul_pd(_mm_mul_pd(r2, b_r), _mm_set1_pd(qy))));
      F[ST_tmp + 6] += dt * b_param.Jinv[8] * qy;
      qw += b_x[i + 1] * b_param.g[i];
    }
    /*  Sola */
    /* %% for hardcoding: */
    /*  [2*qw*vx - 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz, 2*qx*vy -
     * 2*qw*vz - 2*qy*vx, 2*qw*vy + 2*qx*vz - 2*qz*vx] */
    /*  [2*qw*vy + 2*qx*vz - 2*qz*vx, 2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qx*vx +
     * 2*qy*vy + 2*qz*vz, 2*qy*vz - 2*qw*vx - 2*qz*vy] */
    /*  [2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qz*vx - 2*qx*vz - 2*qw*vy, 2*qw*vx -
     * 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz] */
    dv1[0] = b_x[2] * b_param.g[2] - b_param.g[1] * b_x[3];
    dv1[1] = b_param.g[0] * b_x[3] - b_x[1] * b_param.g[2];
    dv1[2] = b_x[1] * b_param.g[1] - b_param.g[0] * b_x[2];
    qy = b_x[0] * 0.0;
    b_w_exp_tilde[0] = qy;
    b_w_exp_tilde[3] = b_x[0] * -b_param.g[2];
    b_w_exp_tilde[6] = b_x[0] * b_param.g[1];
    b_w_exp_tilde[1] = b_x[0] * b_param.g[2];
    b_w_exp_tilde[4] = qy;
    b_w_exp_tilde[7] = b_x[0] * -b_param.g[0];
    b_w_exp_tilde[2] = b_x[0] * -b_param.g[1];
    b_w_exp_tilde[5] = b_x[0] * b_param.g[0];
    b_w_exp_tilde[8] = qy;
    for (i = 0; i < 3; i++) {
      F[i + 7] = dt * (2.0 * (b_x[0] * b_param.g[i] - dv1[i]));
      qy = b_x[i + 1];
      ST_tmp = 11 * (i + 1);
      F[ST_tmp + 7] = dt * (2.0 * (((qw * b_b[3 * i] + b_x[1] * b_param.g[i]) -
                                    b_param.g[0] * qy) +
                                   b_w_exp_tilde[3 * i]));
      b_w_exp_tilde_tmp = 3 * i + 1;
      F[ST_tmp + 8] =
          dt * (2.0 * (((qw * b_b[b_w_exp_tilde_tmp] + b_x[2] * b_param.g[i]) -
                        b_param.g[1] * qy) +
                       b_w_exp_tilde[b_w_exp_tilde_tmp]));
      b_w_exp_tilde_tmp = 3 * i + 2;
      F[ST_tmp + 9] =
          dt * (2.0 * (((qw * b_b[b_w_exp_tilde_tmp] + b_x[3] * b_param.g[i]) -
                        b_param.g[2] * qy) +
                       b_w_exp_tilde[b_w_exp_tilde_tmp]));
    }
    /*  jacobian of {exp(-tilde(w)*dt)*v} wrt w    */
    /*  skew symmetric matrix / cross-product jacobian */
    ST[0] = 0.0;
    ST[3] = -b_x[9];
    ST[6] = b_x[8];
    ST[1] = b_x[9];
    ST[4] = 0.0;
    ST[7] = -b_x[7];
    ST[2] = -b_x[8];
    ST[5] = b_x[7];
    ST[8] = 0.0;
    /*  skew symmetric matrix / cross-product jacobian */
    w_exp_tilde_tmp[0] = 0.0;
    w_exp_tilde_tmp[3] = -b_x[6];
    w_exp_tilde_tmp[6] = b_x[5];
    w_exp_tilde_tmp[1] = b_x[6];
    w_exp_tilde_tmp[4] = 0.0;
    w_exp_tilde_tmp[7] = -b_x[4];
    w_exp_tilde_tmp[2] = -b_x[5];
    w_exp_tilde_tmp[5] = b_x[4];
    w_exp_tilde_tmp[8] = 0.0;
    qw = 0.5 * (dt * dt);
    memset(&b_w_exp_tilde[0], 0, 9U * sizeof(double));
    memset(&b_a[0], 0, 9U * sizeof(double));
    /*  column q */
    /*  column w */
    /*  column v */
    /*     %% altitude row (alt, 11) */
    /*  inverse quaternion */
    /* %% quaternion definition */
    /* %% inverse quaternion  */
    q[0] = b_x[0];
    /*  computes rotation matrix from quaternion */
    /* %% norm quaternions */
    /*  q = quaternion_norm(q); */
    /* %% quaternion definition */
    /* %% rotation matrix */
    /*  S = (qw^2-qv'*qv)*eye(3) + 2*qv*qv' - 2*qw*q_tilde; % Stevens */
    /* %% Jacobian of rotation wrt quaternion */
    b = 0.0;
    r8 = _mm_set1_pd(2.0);
    for (i = 0; i < 3; i++) {
      qy = w_exp_tilde_tmp[3 * i];
      qz = ST[3 * i];
      b_r = _mm_loadu_pd(&ST[0]);
      r1 = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      _mm_storeu_pd(&b_w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_r = _mm_loadu_pd(&w_exp_tilde_tmp[0]);
      r1 = _mm_loadu_pd(&b_a[3 * i]);
      _mm_storeu_pd(&b_a[3 * i], _mm_add_pd(r1, _mm_mul_pd(_mm_mul_pd(r8, b_r),
                                                           _mm_set1_pd(qz))));
      d_w_exp_tilde_tmp = 3 * i + 2;
      b_w_exp_tilde[d_w_exp_tilde_tmp] += ST[2] * qy;
      b_a[d_w_exp_tilde_tmp] += 2.0 * w_exp_tilde_tmp[2] * qz;
      b_w_exp_tilde_tmp = 3 * i + 1;
      qy = w_exp_tilde_tmp[b_w_exp_tilde_tmp];
      qz = ST[b_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&ST[3]);
      r1 = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      _mm_storeu_pd(&b_w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_r = _mm_loadu_pd(&w_exp_tilde_tmp[3]);
      r1 = _mm_loadu_pd(&b_a[3 * i]);
      _mm_storeu_pd(&b_a[3 * i], _mm_add_pd(r1, _mm_mul_pd(_mm_mul_pd(r8, b_r),
                                                           _mm_set1_pd(qz))));
      b_w_exp_tilde[d_w_exp_tilde_tmp] += ST[5] * qy;
      b_a[d_w_exp_tilde_tmp] += 2.0 * w_exp_tilde_tmp[5] * qz;
      qy = w_exp_tilde_tmp[d_w_exp_tilde_tmp];
      qz = ST[d_w_exp_tilde_tmp];
      b_r = _mm_loadu_pd(&ST[6]);
      r1 = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      _mm_storeu_pd(&b_w_exp_tilde[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
      b_r = _mm_loadu_pd(&w_exp_tilde_tmp[6]);
      r1 = _mm_loadu_pd(&b_a[3 * i]);
      _mm_storeu_pd(&b_a[3 * i], _mm_add_pd(r1, _mm_mul_pd(_mm_mul_pd(r8, b_r),
                                                           _mm_set1_pd(qz))));
      b_w_exp_tilde[d_w_exp_tilde_tmp] += 0.0 * qy;
      b_a[d_w_exp_tilde_tmp] += 0.0 * qz;
      b_r = _mm_loadu_pd(&ST[3 * i]);
      r1 = _mm_loadu_pd(&b_w_exp_tilde[3 * i]);
      r2 = _mm_loadu_pd(&b_a[3 * i]);
      b_w_exp_tilde_tmp = 11 * (i + 4);
      _mm_storeu_pd(
          &F[b_w_exp_tilde_tmp + 7],
          _mm_add_pd(_mm_mul_pd(_mm_set1_pd(dt), b_r),
                     _mm_mul_pd(_mm_set1_pd(qw), _mm_sub_pd(r1, r2))));
      b_r = _mm_loadu_pd(&w_exp_tilde[3 * i]);
      ST_tmp = 11 * (i + 7);
      _mm_storeu_pd(&F[ST_tmp + 7], b_r);
      F[b_w_exp_tilde_tmp + 9] =
          dt * qz +
          qw * (b_w_exp_tilde[d_w_exp_tilde_tmp] - b_a[d_w_exp_tilde_tmp]);
      F[ST_tmp + 9] = w_exp_tilde[d_w_exp_tilde_tmp];
      qy = -b_x[i + 1];
      q[i + 1] = qy;
      b += qy * b_x[i + 7];
    }
    __m128d r10;
    /*  Sola */
    /* %% for hardcoding: */
    /*  [2*qw*vx - 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz, 2*qx*vy -
     * 2*qw*vz - 2*qy*vx, 2*qw*vy + 2*qx*vz - 2*qz*vx] */
    /*  [2*qw*vy + 2*qx*vz - 2*qz*vx, 2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qx*vx +
     * 2*qy*vy + 2*qz*vz, 2*qy*vz - 2*qw*vx - 2*qz*vy] */
    /*  [2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qz*vx - 2*qx*vz - 2*qw*vy, 2*qw*vx -
     * 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz] */
    qz = q[3] * b_x[8];
    qw = q[2] * b_x[9];
    expl_temp = q[3] * b_x[7];
    b_expl_temp = q[1] * b_x[9];
    dphi = q[0] * b_x[7];
    c_dt[0] = dt * (2.0 * (dphi - (qw - qz)));
    b_r = _mm_loadu_pd(&q[1]);
    r1 = _mm_loadu_pd(&b_x[7]);
    r2 = _mm_loadu_pd(&ST[0]);
    r8 = _mm_set1_pd(dt);
    r6 = _mm_set1_pd(2.0);
    r9 = _mm_set1_pd(b);
    r7 = _mm_set1_pd(0.0);
    r5 = _mm_set1_pd(q[0]);
    r3 = _mm_set1_pd(q[1]);
    _mm_storeu_pd(
        &c_dt[3],
        _mm_mul_pd(
            r8,
            _mm_mul_pd(
                r6, _mm_add_pd(
                        _mm_sub_pd(
                            _mm_add_pd(
                                _mm_mul_pd(r9, _mm_loadu_pd(&b_b[0])),
                                _mm_add_pd(
                                    r7, _mm_mul_pd(b_r, _mm_set1_pd(b_x[7])))),
                            _mm_mul_pd(r1, r3)),
                        _mm_mul_pd(r5, r2)))));
    qy = b * 0.0;
    c_dt[5] = dt * (2.0 * (((qy + expl_temp) - b_expl_temp) + q[0] * -b_x[8]));
    c_dt[1] = dt * (2.0 * (q[0] * b_x[8] - (expl_temp - b_expl_temp)));
    b_r = _mm_loadu_pd(&q[1]);
    r1 = _mm_loadu_pd(&b_x[7]);
    r2 = _mm_loadu_pd(&ST[3]);
    r4 = _mm_set1_pd(q[2]);
    _mm_storeu_pd(
        &c_dt[6],
        _mm_mul_pd(
            r8,
            _mm_mul_pd(
                r6, _mm_add_pd(
                        _mm_sub_pd(
                            _mm_add_pd(
                                _mm_mul_pd(r9, _mm_loadu_pd(&b_b[3])),
                                _mm_add_pd(
                                    r7, _mm_mul_pd(b_r, _mm_set1_pd(b_x[8])))),
                            _mm_mul_pd(r1, r4)),
                        _mm_mul_pd(r5, r2)))));
    c_dt[8] = dt * (2.0 * (((qy + qz) - qw) + dphi));
    c_dt[2] = dt * (2.0 * (q[0] * b_x[9] - (q[1] * b_x[8] - q[2] * b_x[7])));
    b_r = _mm_loadu_pd(&q[1]);
    r1 = _mm_loadu_pd(&b_x[7]);
    r2 = _mm_loadu_pd(&ST[6]);
    r10 = _mm_set1_pd(q[3]);
    _mm_storeu_pd(
        &c_dt[9],
        _mm_mul_pd(
            r8,
            _mm_mul_pd(
                r6, _mm_add_pd(
                        _mm_sub_pd(
                            _mm_add_pd(
                                _mm_mul_pd(r9, _mm_loadu_pd(&b_b[6])),
                                _mm_add_pd(
                                    r7, _mm_mul_pd(b_r, _mm_set1_pd(b_x[9])))),
                            _mm_mul_pd(r1, r10)),
                        _mm_mul_pd(r5, r2)))));
    F[10] = c_dt[0];
    F[21] = c_dt[3];
    F[32] = c_dt[6];
    F[43] = c_dt[9];
    /*  only use altitude from position vector */
    /*  computes rotation matrix from quaternion */
    /* %% norm quaternions */
    /*  q = quaternion_norm(q); */
    /* %% quaternion definition */
    /* %% skew symetric quaternion matrix */
    /* %% rotation matrix */
    qy = q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + q[3] * q[3]);
    qw = 2.0 * q[0];
    /*  Stevens */
    /* %% for hardcoding:  */
    /*  [qw^2 + qx^2 - qy^2 - qz^2,         2*qw*qz + 2*qx*qy,         2*qx*qz -
     * 2*qw*qy] */
    /*  [        2*qx*qy - 2*qw*qz, qw^2 - qx^2 + qy^2 - qz^2,         2*qw*qx +
     * 2*qy*qz] */
    /*  [        2*qw*qy + 2*qx*qz,         2*qy*qz - 2*qw*qx, qw^2 - qx^2 -
     * qy^2 + qz^2] */
    b_r = _mm_loadu_pd(&q[1]);
    r1 = _mm_set1_pd(qy);
    _mm_storeu_pd(
        &b_a[0],
        _mm_add_pd(_mm_mul_pd(r1, _mm_loadu_pd(&b_b[0])),
                   _mm_add_pd(r7, _mm_mul_pd(_mm_mul_pd(r6, b_r), r3))));
    qy *= 0.0;
    qz = 2.0 * q[3];
    b_a[2] = qy + qz * q[1];
    F[87] = dt * (b_a[0] - qw * 0.0);
    b_r = _mm_loadu_pd(&q[1]);
    _mm_storeu_pd(
        &b_a[3],
        _mm_add_pd(_mm_mul_pd(r1, _mm_loadu_pd(&b_b[3])),
                   _mm_add_pd(r7, _mm_mul_pd(_mm_mul_pd(r6, b_r), r4))));
    b_a[5] = qy + qz * q[2];
    F[98] = dt * (b_a[3] - qw * -q[3]);
    b_r = _mm_loadu_pd(&q[1]);
    _mm_storeu_pd(
        &b_a[6],
        _mm_add_pd(_mm_mul_pd(r1, _mm_loadu_pd(&b_b[6])),
                   _mm_add_pd(r7, _mm_mul_pd(_mm_mul_pd(r6, b_r), r10))));
    F[109] = dt * (b_a[6] - qw * q[2]);
    /*  r_r = eye(3); */
    /*  column q */
    /*  column v */
    F[120] = 1.0;
    /*  column alt */
    /* %% discrete covariance */
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      b_w_exp_tilde_tmp = 11 * i + 10;
      for (i1 = 0; i1 < 11; i1++) {
        qy = P[i1 + 11 * i];
        for (i4 = 0; i4 <= 8; i4 += 2) {
          b_r = _mm_loadu_pd(&F[i4 + 11 * i1]);
          ST_tmp = i4 + 11 * i;
          r1 = _mm_loadu_pd(&b_F[ST_tmp]);
          _mm_storeu_pd(&b_F[ST_tmp],
                        _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
        }
        b_F[b_w_exp_tilde_tmp] += F[11 * i1 + 10] * qy;
      }
    }
#pragma omp parallel for num_threads(omp_get_max_threads()) private(           \
        d6, i5, i6, P_pred_tmp)

    for (i3 = 0; i3 < 11; i3++) {
      for (i5 = 0; i5 < 11; i5++) {
        d6 = 0.0;
        for (i6 = 0; i6 < 11; i6++) {
          d6 += b_F[i3 + 11 * i6] * F[i5 + 11 * i6];
        }
        P_pred_tmp = i3 + 11 * i5;
        P_pred[P_pred_tmp] = d6 + Q[P_pred_tmp];
      }
    }
    /*     %% Correction */
    /*  computes a-posteriori state and covariance estimates */
    /*  Uses discrete-time measurement model and analytical Jacobian */
    /* %% compute expected measurement and difference to measured values */
    /*  hardcoded measurement model, state vector: [q(4); w(3); v(3); alt] */
    /* %% compute Jacobian: H = dh/dx */
    /*  H = zeros(3, 11); % is hardcoded */
    /*  H(:, 5:7) = eye(3); % hardcoded measurement jacobian */
    /* %% compute Kalman gain (and helper matrices) */
    /*  hardcoded H*P*H'  */
    b_r = _mm_loadu_pd(&P_pred[48]);
    _mm_storeu_pd(&b_w_exp_tilde[0], _mm_add_pd(b_r, _mm_loadu_pd(&R[0])));
    b_w_exp_tilde[2] = P_pred[50];
    b_r = _mm_loadu_pd(&P_pred[59]);
    _mm_storeu_pd(&b_w_exp_tilde[3], _mm_add_pd(b_r, _mm_loadu_pd(&R[3])));
    b_w_exp_tilde[5] = P_pred[61];
    b_r = _mm_loadu_pd(&P_pred[70]);
    _mm_storeu_pd(&b_w_exp_tilde[6], _mm_add_pd(b_r, _mm_loadu_pd(&R[6])));
    b_w_exp_tilde[8] = P_pred[72] + 1.0E-9;
    inv(b_w_exp_tilde, ST);
    memset(&K[0], 0, 33U * sizeof(double));
    for (i = 0; i < 3; i++) {
      b_w_exp_tilde_tmp = 11 * i + 10;
      for (i1 = 0; i1 < 3; i1++) {
        qy = ST[i1 + 3 * i];
        for (i4 = 0; i4 <= 8; i4 += 2) {
          b_r = _mm_loadu_pd(&P_pred[i4 + 11 * (i1 + 4)]);
          ST_tmp = i4 + 11 * i;
          r1 = _mm_loadu_pd(&K[ST_tmp]);
          _mm_storeu_pd(&K[ST_tmp],
                        _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
        }
        K[b_w_exp_tilde_tmp] += P_pred[11 * (i1 + 4) + 10] * qy;
      }
    }
    /*  hardcoded P*H */
    d_eye(F);
    memcpy(&P[0], &F[0], 44U * sizeof(double));
    for (i1 = 0; i1 < 3; i1++) {
      for (i = 0; i <= 8; i += 2) {
        b_w_exp_tilde_tmp = i + 11 * (i1 + 4);
        b_r = _mm_loadu_pd(&F[b_w_exp_tilde_tmp]);
        r1 = _mm_loadu_pd(&K[i + 11 * i1]);
        _mm_storeu_pd(&P[b_w_exp_tilde_tmp], _mm_sub_pd(b_r, r1));
      }
      b_w_exp_tilde_tmp = 11 * (i1 + 4) + 10;
      P[b_w_exp_tilde_tmp] = F[b_w_exp_tilde_tmp] - K[11 * i1 + 10];
    }
    memcpy(&P[77], &F[77], 44U * sizeof(double));
    /*  hardcoded K*H */
    /* %% correct covariance estimate */
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      b_w_exp_tilde_tmp = 11 * i + 10;
      for (i1 = 0; i1 < 11; i1++) {
        qy = P_pred[i1 + 11 * i];
        for (i4 = 0; i4 <= 8; i4 += 2) {
          b_r = _mm_loadu_pd(&P[i4 + 11 * i1]);
          ST_tmp = i4 + 11 * i;
          r1 = _mm_loadu_pd(&b_F[ST_tmp]);
          _mm_storeu_pd(&b_F[ST_tmp],
                        _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
        }
        b_F[b_w_exp_tilde_tmp] += P[11 * i1 + 10] * qy;
      }
    }
    memset(&b_K[0], 0, 33U * sizeof(double));
    for (i = 0; i < 3; i++) {
      b_w_exp_tilde_tmp = 11 * i + 10;
      for (i1 = 0; i1 < 3; i1++) {
        qy = R[i1 + 3 * i];
        for (i4 = 0; i4 <= 8; i4 += 2) {
          b_r = _mm_loadu_pd(&K[i4 + 11 * i1]);
          ST_tmp = i4 + 11 * i;
          r1 = _mm_loadu_pd(&b_K[ST_tmp]);
          _mm_storeu_pd(&b_K[ST_tmp],
                        _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
        }
        b_K[b_w_exp_tilde_tmp] += K[11 * i1 + 10] * qy;
      }
    }
    memset(&F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      b_w_exp_tilde_tmp = 11 * i + 10;
      for (i1 = 0; i1 < 11; i1++) {
        qy = P[i + 11 * i1];
        for (i4 = 0; i4 <= 8; i4 += 2) {
          b_r = _mm_loadu_pd(&b_F[i4 + 11 * i1]);
          ST_tmp = i4 + 11 * i;
          r1 = _mm_loadu_pd(&F[ST_tmp]);
          _mm_storeu_pd(&F[ST_tmp],
                        _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
        }
        F[b_w_exp_tilde_tmp] += b_F[11 * i1 + 10] * qy;
      }
    }
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      b_w_exp_tilde_tmp = 11 * i + 10;
      for (i1 = 0; i1 < 3; i1++) {
        qy = K[i + 11 * i1];
        for (i4 = 0; i4 <= 8; i4 += 2) {
          b_r = _mm_loadu_pd(&b_K[i4 + 11 * i1]);
          ST_tmp = i4 + 11 * i;
          r1 = _mm_loadu_pd(&b_F[ST_tmp]);
          _mm_storeu_pd(&b_F[ST_tmp],
                        _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(qy))));
        }
        b_F[b_w_exp_tilde_tmp] += b_K[11 * i1 + 10] * qy;
      }
    }
    for (i = 0; i <= 118; i += 2) {
      b_r = _mm_loadu_pd(&F[i]);
      r1 = _mm_loadu_pd(&b_F[i]);
      _mm_storeu_pd(&P[i], _mm_add_pd(b_r, r1));
    }
    P[120] = F[120] + b_F[120];
    /*  joseph stabilized */
    /* %% correct state estimate */
    C_total_w_idx_0 =
        ((d / C_total_w_idx_0 * (board_gyro->meas[0] - c_b.board_gyro[0]) +
          d1 / C_total_w_idx_0 * (mti_gyro->meas[0] - c_b.mti_gyro[0])) +
         a[0] * (ad_gyro->meas[0] - c_b.ad_gyro[0])) -
        xhat[4];
    qz = d / C_total_w_idx_1;
    qy = d1 / C_total_w_idx_1;
    C_total_w_idx_1 = ((qz * (board_gyro->meas[1] - c_b.board_gyro[1]) +
                        qy * (mti_gyro->meas[1] - c_b.mti_gyro[1])) +
                       a[1] * (ad_gyro->meas[1] - c_b.ad_gyro[1])) -
                      xhat[5];
    qy = ((qz * (board_gyro->meas[2] - c_b.board_gyro[2]) +
           qy * (mti_gyro->meas[2] - c_b.mti_gyro[2])) +
          a[2] * (ad_gyro->meas[2] - c_b.ad_gyro[2])) -
         xhat[6];
    for (i = 0; i <= 8; i += 2) {
      b_r = _mm_loadu_pd(&K[i]);
      r1 = _mm_mul_pd(b_r, _mm_set1_pd(C_total_w_idx_0));
      b_r = _mm_loadu_pd(&K[i + 11]);
      b_r = _mm_mul_pd(b_r, _mm_set1_pd(C_total_w_idx_1));
      r1 = _mm_add_pd(r1, b_r);
      b_r = _mm_loadu_pd(&K[i + 22]);
      b_r = _mm_mul_pd(b_r, _mm_set1_pd(qy));
      b_r = _mm_add_pd(r1, b_r);
      r1 = _mm_loadu_pd(&xhat[i]);
      b_r = _mm_add_pd(r1, b_r);
      _mm_storeu_pd(&b_x[i], b_r);
    }
    b_x[10] = xhat[10] + ((K[10] * C_total_w_idx_0 + K[21] * C_total_w_idx_1) +
                          K[32] * qy);
    /*  norms quaternion */
    /* %% inverse quaternion  */
    b_r = _mm_loadu_pd(&b_x[0]);
    r1 = _mm_set1_pd(c_norm(&b_x[0]));
    _mm_storeu_pd(&b_x[0], _mm_div_pd(b_r, r1));
    b_r = _mm_loadu_pd(&b_x[2]);
    _mm_storeu_pd(&b_x[2], _mm_div_pd(b_r, r1));
    /*  norm quaternions */
    /*     %% Correction steps, sequential for each additional sensor */
    /* %% Correction with barometer, magnetometer */
    /* %% R is a square matrix (size length of sensor vector), tuning for
     * expected measurement noise magnitude E(noise) */
    /* %% Barometer */
    if (board_baro->status) {
      /*  only correct with alive IMUs */
      /* %% y = [ P(1) ] */
      memcpy(&xhat[0], &b_x[0], 11U * sizeof(double));
      memcpy(&P_pred[0], &P[0], 121U * sizeof(double));
      b_ekf_correct(xhat, P_pred, board_baro->meas, c_b.board_baro, b_x, P);
    }
    if (mti_baro->status) {
      /*  only correct with alive IMUs */
      /* %% y = [ P(1) ] */
      memcpy(&xhat[0], &b_x[0], 11U * sizeof(double));
      memcpy(&P_pred[0], &P[0], 121U * sizeof(double));
      b_ekf_correct(xhat, P_pred, mti_baro->meas, c_b.mti_baro, b_x, P);
    }
    /* %% Magnetometer */
    if (board_mag->status) {
      /* %% y = [  Mag(3) ] */
      memcpy(&xhat[0], &b_x[0], 11U * sizeof(double));
      memcpy(&P_pred[0], &P[0], 121U * sizeof(double));
      ekf_correct(xhat, P_pred, board_mag->meas, c_b.board_mag_earth, b_b, b_x,
                  P);
    }
    if (mti_mag->status) {
      /* %% y = [  Mag(3) ] */
      memcpy(&xhat[0], &b_x[0], 11U * sizeof(double));
      memcpy(&P_pred[0], &P[0], 121U * sizeof(double));
      ekf_correct(xhat, P_pred, mti_mag->meas, c_b.mti_mag_earth, b_b, b_x, P);
    }
  }
  /*     %% Pack state as struct */
  /* %% use union in C or smth */
  state->q[0] = b_x[0];
  state->q[1] = b_x[1];
  state->q[2] = b_x[2];
  state->q[3] = b_x[3];
  /*  [1], attitude quaternion */
  /*  [rad/s], angular rate  */
  state->w[0] = b_x[4];
  state->v[0] = b_x[7];
  state->w[1] = b_x[5];
  state->v[1] = b_x[8];
  state->w[2] = b_x[6];
  state->v[2] = b_x[9];
  /*  [m/s], velocity */
  state->alt = b_x[10];
  /*  [m], altitude */
  memcpy(&state->x[0], &b_x[0], 11U * sizeof(double));
  /*  also full state as vector is needed in simulink */
  /*     %% Compute variance norm  */
  /* %% for evaluating EKF quality */
  *cov_norm = d_norm(P);
  /*  scalar, 2-norm of the covariance matrix */
  /*     %% Compute air data */
  airdata->pressure = airdata_atmos(b_x[10], &airdata->temperature,
                                    &airdata->density, &airdata->sonic_speed,
                                    &airdata->mach, &airdata->dynamic_pressure);
  /*  appends airadata struct with dynamic air data from velocity vector */
  /*  dynamic air data: dynamic pressure, mach number */
  qy = b_norm(&b_x[7]);
  /*  return values */
  airdata->mach = qy / airdata->sonic_speed;
  /*  [-], Mach number */
  airdata->dynamic_pressure = 0.5 * airdata->density * (qy * qy);
  /*  [Pa] */
  /*     %% controller input vector */
  /*  computes euler angles (yaw-pitch-roll sequence) from quaternion */
  /* %% norm quaternions */
  /*  q = quaternion_norm(q); */
  /* %% quaternion definition */
  /* %% Euler angles, after Zipfel */
  /*  yaw = atan2( 2*(qx*qy + qw*qz), (qw^2+qx^2-qy^2-qz^2) ); */
  /*  pitch = asin( - 2*(qx*qz - qw*qy) ); */
  roll_state[0] =
      b_atan2(2.0 * (b_x[2] * b_x[3] + b_x[0] * b_x[1]),
              ((b_x[0] * b_x[0] - b_x[1] * b_x[1]) - b_x[2] * b_x[2]) +
                  b_x[3] * b_x[3]);
  roll_state[1] = b_x[4];
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void navigation_codegen_entry_init(void)
{
  b_not_empty = false;
  x_not_empty = false;
  memset(&P[0], 0, 121U * sizeof(double));
}

/*
 * File trailer for navigation_codegen_entry.c
 *
 * [EOF]
 */
