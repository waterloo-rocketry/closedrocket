/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: navigation_codegen_entry.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

/* Include Files */
#include "navigation_codegen_entry.h"
#include "GNC_codegen_data.h"
#include "GNC_codegen_initialize.h"
#include "GNC_codegen_types.h"
#include "airdata_atmos.h"
#include "atan2.h"
#include "ekf_correct.h"
#include "eye.h"
#include "inv.h"
#include "norm.h"
#include "rt_nonfinite.h"
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
static bool x_not_empty;

static double P[121];

static bool b_not_empty;

/* Function Definitions */
/*
 * Calls the pad and flight filters.
 *
 * Arguments    : double dt
 *                bool flight_phase
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
    double dt, bool flight_phase, const struct0_T *board_accel,
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
  double F[121];
  double P_pred[121];
  double b_F[121];
  double K[33];
  double b_K[33];
  double aBuffer[16];
  double b_W_dt[16];
  double b_q[16];
  double c_W_dt[16];
  double xhat[11];
  double ST[9];
  double c_a[9];
  double dv3[9];
  double w_exp_tilde[9];
  double c_q[4];
  double q[4];
  double a[3];
  double b_dt[3];
  double b_w_exp_tilde[3];
  double dn[3];
  double dv1[3];
  double b_expl_temp;
  double expl_temp;
  double n_idx_2;
  double qw;
  double qy;
  double qz;
  double t1_density;
  int i;
  int i1;
  int i2;
  if (!isInitialized_GNC_codegen) {
    GNC_codegen_initialize();
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
      board_accel_f[0] =
          0.0005 * board_accel->meas[0] + 0.9995 * board_accel_f[0];
      board_accel_f[1] =
          0.0005 * board_accel->meas[1] + 0.9995 * board_accel_f[1];
      board_accel_f[2] =
          0.0005 * board_accel->meas[2] + 0.9995 * board_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (board_gyro->status) {
      board_gyro_f[0] = 0.0005 * board_gyro->meas[0] + 0.9995 * board_gyro_f[0];
      board_gyro_f[1] = 0.0005 * board_gyro->meas[1] + 0.9995 * board_gyro_f[1];
      board_gyro_f[2] = 0.0005 * board_gyro->meas[2] + 0.9995 * board_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_accel->status) {
      mti_accel_f[0] = 0.0005 * mti_accel->meas[0] + 0.9995 * mti_accel_f[0];
      mti_accel_f[1] = 0.0005 * mti_accel->meas[1] + 0.9995 * mti_accel_f[1];
      mti_accel_f[2] = 0.0005 * mti_accel->meas[2] + 0.9995 * mti_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_gyro->status) {
      mti_gyro_f[0] = 0.0005 * mti_gyro->meas[0] + 0.9995 * mti_gyro_f[0];
      mti_gyro_f[1] = 0.0005 * mti_gyro->meas[1] + 0.9995 * mti_gyro_f[1];
      mti_gyro_f[2] = 0.0005 * mti_gyro->meas[2] + 0.9995 * mti_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (ad_accel->status) {
      ad_accel_f[0] = 0.0005 * ad_accel->meas[0] + 0.9995 * ad_accel_f[0];
      ad_accel_f[1] = 0.0005 * ad_accel->meas[1] + 0.9995 * ad_accel_f[1];
      ad_accel_f[2] = 0.0005 * ad_accel->meas[2] + 0.9995 * ad_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (ad_gyro->status) {
      ad_gyro_f[0] = 0.0005 * ad_gyro->meas[0] + 0.9995 * ad_gyro_f[0];
      ad_gyro_f[1] = 0.0005 * ad_gyro->meas[1] + 0.9995 * ad_gyro_f[1];
      ad_gyro_f[2] = 0.0005 * ad_gyro->meas[2] + 0.9995 * ad_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (board_baro->status) {
      board_baro_f = 0.0005 * board_baro->meas + 0.9995 * board_baro_f;
    }
    /* %% lowpass filter function used in pad filter */
    if (board_mag->status) {
      board_mag_f[0] = 0.0005 * board_mag->meas[0] + 0.9995 * board_mag_f[0];
      board_mag_f[1] = 0.0005 * board_mag->meas[1] + 0.9995 * board_mag_f[1];
      board_mag_f[2] = 0.0005 * board_mag->meas[2] + 0.9995 * board_mag_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_baro->status) {
      mti_baro_f = 0.0005 * mti_baro->meas + 0.9995 * mti_baro_f;
    }
    /* %% lowpass filter function used in pad filter */
    if (mti_mag->status) {
      mti_mag_f[0] = 0.0005 * mti_mag->meas[0] + 0.9995 * mti_mag_f[0];
      mti_mag_f[1] = 0.0005 * mti_mag->meas[1] + 0.9995 * mti_mag_f[1];
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
      a[0] += mti_accel_f[0];
      a[1] += mti_accel_f[1];
      a[2] += mti_accel_f[2];
    }
    if (ad_accel->status) {
      a[0] += ad_accel_f[0];
      a[1] += ad_accel_f[1];
      a[2] += ad_accel_f[2];
    }
    /* %% computes initial orientation of stationary body from gravity
     * acceleration  */
    /* %% Input: 3D acceleration vector */
    /* %% Output: Orientation quaternion */
    /* %% normed gravity vector in body-fixed frame */
    qz = b_norm(a) + 1.0E-6;
    /*  unit vector of gravity direction */
    /* %% determine initial orientation quaternion */
    qw = sqrt(0.5 * (a[0] / qz) + 0.5);
    if (qw == 0.0) {
      /*  exact upside down case */
      qy = 1.0;
      /*  either qy = 1 or qz = 1, this is arbitrary  */
      qz = 0.0;
    } else {
      qy = 0.5 * (a[2] / qz) / qw;
      qz = -0.5 * (a[1] / qz) / qw;
    }
    q[0] = qw;
    q[1] = 0.0;
    q[2] = qy;
    q[3] = qz;
    n_idx_2 = c_norm(q);
    /*  a gets normed inside function */
    /* %% launch altitude */
    /* %% set constant initials */
    /*  stationary on rail */
    /* %% conconct state vector */
    expl_temp = qw / n_idx_2;
    q[0] = expl_temp;
    b_x[0] = expl_temp;
    expl_temp = 0.0 / n_idx_2;
    q[1] = expl_temp;
    b_x[1] = expl_temp;
    expl_temp = qy / n_idx_2;
    q[2] = expl_temp;
    b_x[2] = expl_temp;
    expl_temp = qz / n_idx_2;
    q[3] = expl_temp;
    b_x[3] = expl_temp;
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
    qz = q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + expl_temp * expl_temp);
    qy = 2.0 * q[0];
    /*  Stevens */
    /* %% for hardcoding:  */
    /*  [qw^2 + qx^2 - qy^2 - qz^2,         2*qw*qz + 2*qx*qy,         2*qx*qz -
     * 2*qw*qy] */
    /*  [        2*qx*qy - 2*qw*qz, qw^2 - qx^2 + qy^2 - qz^2,         2*qw*qx +
     * 2*qy*qz] */
    /*  [        2*qw*qy + 2*qx*qz,         2*qy*qz - 2*qw*qx, qw^2 - qx^2 -
     * qy^2 + qz^2] */
    for (i = 0; i < 3; i++) {
      qw = 2.0 * q[i + 1];
      ST[3 * i] = qz * b_b[i] + qw * q[1];
      ST[3 * i + 1] = qz * b_b[i + 3] + qw * q[2];
      ST[3 * i + 2] = qz * b_b[i + 6] + qw * expl_temp;
    }
    qz = qy * 0.0;
    c_a[0] = qz;
    c_a[1] = qy * -expl_temp;
    c_a[2] = qy * q[2];
    c_a[3] = qy * expl_temp;
    c_a[4] = qz;
    c_a[5] = qy * -q[1];
    c_a[6] = qy * -q[2];
    c_a[7] = qy * q[1];
    c_a[8] = qz;
    for (i = 0; i < 9; i++) {
      ST[i] -= c_a[i];
    }
    /*  launch attitude */
    c_b.board_mag_earth[0] = 0.0;
    c_b.board_mag_earth[1] = 0.0;
    c_b.board_mag_earth[2] = 0.0;
    for (i = 0; i < 3; i++) {
      qz = board_mag_f[i];
      c_b.board_mag_earth[0] += ST[3 * i] * qz;
      c_b.board_mag_earth[1] += ST[3 * i + 1] * qz;
      c_b.board_mag_earth[2] += ST[3 * i + 2] * qz;
      c_b.mti_mag_earth[i] = 0.0;
    }
    qz = c_b.mti_mag_earth[0];
    qw = c_b.mti_mag_earth[1];
    n_idx_2 = c_b.mti_mag_earth[2];
    for (i = 0; i < 3; i++) {
      qy = mti_mag_f[i];
      qz += ST[3 * i] * qy;
      qw += ST[3 * i + 1] * qy;
      n_idx_2 += ST[3 * i + 2] * qy;
    }
    c_b.mti_mag_earth[2] = n_idx_2;
    c_b.mti_mag_earth[1] = qw;
    c_b.mti_mag_earth[0] = qz;
    /* %% barometer */
    qz = airdata_atmos(420.0, &qw, &t1_density, &qz, &n_idx_2, &qy);
    /*  what the pressure should be at launch elevation */
    c_b.board_baro = board_baro_f - qz;
    c_b.mti_baro = mti_baro_f - qz;
    b_not_empty = true;
  }
  /*     %% Flight filter iteration */
  if (flight_phase) {
    double W_dt[16];
    double dv2[16];
    double c_dt[12];
    double w_exp_tilde_tmp[9];
    double dv[4];
    double A[3];
    double A_tmp_tmp;
    double b;
    double b_A_tmp_tmp;
    double b_a;
    double c_A_tmp_tmp;
    double d;
    double d1;
    double d10;
    double d2;
    double d3;
    double d4;
    double d5;
    double d6;
    double d7;
    double d8;
    double d9;
    double dphi;
    double dphi_tmp;
    double x;
    int b_w_exp_tilde_tmp;
    int c_w_exp_tilde_tmp;
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
    d = 9.9999999999999981E+9 * (double)ad_gyro->status;
    A_tmp_tmp = 1.0000000000000002E+14 * (double)board_accel->status;
    b_A_tmp_tmp = 1.0000000000000002E+14 * (double)mti_accel->status;
    c_A_tmp_tmp = 1.0000000000000002E+14 * (double)ad_accel->status;
    qw = (A_tmp_tmp + b_A_tmp_tmp) + c_A_tmp_tmp;
    A[0] = qw;
    d1 = 9.9999999999999981E+9 * (double)board_gyro->status;
    d2 = 9.9999999999999981E+9 * (double)mti_gyro->status;
    d3 = d1 + d2;
    d4 = d3 + d;
    d /= d4;
    a[0] = d;
    A[1] = qw;
    d = 0.0 / d3;
    A[2] = qw;
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
    qw = c_norm(&b_x[0]);
    q[0] = b_x[0] / qw;
    q[1] = b_x[1] / qw;
    q[2] = b_x[2] / qw;
    q[3] = b_x[3] / qw;
    /* %% incremental quaternion */
    dphi_tmp = b_norm(&b_x[4]);
    qw = dphi_tmp * dt;
    dphi = qw / 2.0;
    if (dphi_tmp == 0.0) {
      dn[0] = 0.0;
      dn[1] = 0.0;
      dn[2] = 0.0;
      qz = 0.0;
      qy = 0.0;
      n_idx_2 = 0.0;
    } else {
      qz = b_x[4] / dphi_tmp;
      dn[0] = qz;
      qy = b_x[5] / dphi_tmp;
      dn[1] = qy;
      n_idx_2 = b_x[6] / dphi_tmp;
      dn[2] = n_idx_2;
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
    ST[3] = -n_idx_2;
    ST[6] = qy;
    ST[1] = n_idx_2;
    ST[4] = 0.0;
    ST[7] = -qz;
    ST[2] = -qy;
    ST[5] = qz;
    ST[8] = 0.0;
    /* %% Rodrigues rotation formula */
    b_a = sin(qw);
    x = cos(qw);
    b_eye(w_exp_tilde_tmp);
    memset(&w_exp_tilde[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qw = w_exp_tilde[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        qz = ST[i1 + 3 * i];
        qw += ST[3 * i1] * qz;
        w_exp_tilde[b_w_exp_tilde_tmp] += ST[3 * i1 + 1] * qz;
        w_exp_tilde[c_w_exp_tilde_tmp] += ST[3 * i1 + 2] * qz;
      }
      w_exp_tilde[3 * i] = qw;
    }
    for (i = 0; i < 9; i++) {
      w_exp_tilde[i] =
          (w_exp_tilde_tmp[i] - b_a * ST[i]) + (1.0 - x) * w_exp_tilde[i];
    }
    /*  aerodynamics model */
    /* %% air data   */
    /*  appends airadata struct with dynamic air data from velocity vector */
    /*  dynamic air data: dynamic pressure, mach number */
    qw = b_norm(&b_x[7]);
    /*  return values */
    /*  [-], Mach number */
    /*  [Pa] */
    airdata_atmos(b_x[10], &qz, &t1_density, &qy, &expl_temp, &b_expl_temp);
    n_idx_2 = 0.5 * t1_density * (qw * qw);
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
    for (i = 0; i < 3; i++) {
      qz = b_x[i + 1];
      ST[3 * i] = qw * b_b[3 * i] + 2.0 * b_x[1] * qz;
      b_w_exp_tilde_tmp = 3 * i + 1;
      ST[b_w_exp_tilde_tmp] = qw * b_b[b_w_exp_tilde_tmp] + 2.0 * b_x[2] * qz;
      b_w_exp_tilde_tmp = 3 * i + 2;
      ST[b_w_exp_tilde_tmp] = qw * b_b[b_w_exp_tilde_tmp] + 2.0 * b_x[3] * qz;
    }
    qw = qy * 0.0;
    c_a[0] = qw;
    c_a[3] = qy * -b_x[3];
    c_a[6] = qy * b_x[2];
    c_a[1] = qy * b_x[3];
    c_a[4] = qw;
    c_a[7] = qy * -b_x[1];
    c_a[2] = qy * -b_x[2];
    c_a[5] = qy * b_x[1];
    c_a[8] = qw;
    for (i = 0; i < 9; i++) {
      ST[i] -= c_a[i];
    }
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
    memset(&c_a[0], 0, 9U * sizeof(double));
    memset(&b_w_exp_tilde[0], 0, 3U * sizeof(double));
    for (i = 0; i < 3; i++) {
      dv[i + 1] = dn[i] * b;
      qz = c_a[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        qw = param.J[i1 + 3 * i];
        qz += w_exp_tilde[3 * i1] * qw;
        c_a[b_w_exp_tilde_tmp] += w_exp_tilde[3 * i1 + 1] * qw;
        c_a[c_w_exp_tilde_tmp] += w_exp_tilde[3 * i1 + 2] * qw;
      }
      c_a[3 * i] = qz;
      qw = b_x[i + 4];
      b_w_exp_tilde[0] += qz * qw;
      b_w_exp_tilde[1] += c_a[3 * i + 1] * qw;
      b_w_exp_tilde[2] += c_a[3 * i + 2] * qw;
    }
    dn[0] = n_idx_2 * -0.0;
    dn[1] = n_idx_2 * (-0.16182736457722724 * sin(b_atan2(b_x[9], b_x[7])));
    dn[2] = n_idx_2 * (-0.16182736457722724 * -sin(b_atan2(b_x[8], b_x[7])));
    memset(&dv1[0], 0, 3U * sizeof(double));
    memset(&b_dt[0], 0, 3U * sizeof(double));
    d5 = dv1[0];
    d6 = dv1[1];
    d7 = dv1[2];
    d8 = b_dt[0];
    d9 = b_dt[1];
    d10 = b_dt[2];
    n_idx_2 = b_x[7];
    expl_temp = b_x[8];
    b_expl_temp = b_x[9];
    for (i = 0; i < 3; i++) {
      qw = param.Jinv[3 * i];
      qz = b_w_exp_tilde[i];
      d5 += qw * qz;
      qy = dn[i];
      d8 += dt * qw * qy;
      qw = param.Jinv[3 * i + 1];
      d6 += qw * qz;
      d9 += dt * qw * qy;
      qw = param.Jinv[3 * i + 2];
      d7 += qw * qz;
      d10 += dt * qw * qy;
      qw = A[i];
      b_w_exp_tilde[i] =
          ((w_exp_tilde[i] * n_idx_2 + w_exp_tilde[i + 3] * expl_temp) +
           w_exp_tilde[i + 6] * b_expl_temp) +
          dt * ((A_tmp_tmp / qw * board_accel->meas[i] +
                 b_A_tmp_tmp / qw * mti_accel->meas[i]) +
                c_A_tmp_tmp / qw * ad_accel->meas[i]);
    }
    memset(&dn[0], 0, 3U * sizeof(double));
    t1_density = dn[0];
    dphi = dn[1];
    b = dn[2];
    qw = b_x[7];
    qz = b_x[8];
    qy = b_x[9];
    for (i = 0; i < 3; i++) {
      n_idx_2 = ST[3 * i];
      expl_temp = param.g[i];
      t1_density += dt * n_idx_2 * expl_temp;
      b_expl_temp = n_idx_2 * qw;
      n_idx_2 = ST[3 * i + 1];
      dphi += dt * n_idx_2 * expl_temp;
      b_expl_temp += n_idx_2 * qz;
      n_idx_2 = ST[3 * i + 2];
      b += dt * n_idx_2 * expl_temp;
      b_expl_temp += n_idx_2 * qy;
      A[i] = b_expl_temp;
    }
    memset(&c_q[0], 0, sizeof(double) << 2);
    qw = c_q[0];
    qz = c_q[1];
    qy = c_q[2];
    n_idx_2 = c_q[3];
    for (i = 0; i < 4; i++) {
      b_w_exp_tilde_tmp = i << 2;
      expl_temp = dv[i];
      qw += b_q[b_w_exp_tilde_tmp] * expl_temp;
      qz += b_q[b_w_exp_tilde_tmp + 1] * expl_temp;
      qy += b_q[b_w_exp_tilde_tmp + 2] * expl_temp;
      n_idx_2 += b_q[b_w_exp_tilde_tmp + 3] * expl_temp;
    }
    xhat[0] = qw;
    xhat[1] = qz;
    xhat[2] = qy;
    xhat[3] = n_idx_2;
    xhat[4] = d5 + d8;
    xhat[7] = b_w_exp_tilde[0] + t1_density;
    xhat[5] = d6 + d9;
    xhat[8] = b_w_exp_tilde[1] + dphi;
    xhat[6] = d7 + d10;
    xhat[9] = b_w_exp_tilde[2] + b;
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
    qy = t1_density * -b_x[4];
    W_dt[4] = qy;
    qw = t1_density * -b_x[5];
    W_dt[8] = qw;
    n_idx_2 = t1_density * -b_x[6];
    W_dt[12] = n_idx_2;
    expl_temp = t1_density * b_x[4];
    W_dt[1] = expl_temp;
    W_dt[5] = qz;
    b_expl_temp = t1_density * b_x[6];
    W_dt[9] = b_expl_temp;
    W_dt[13] = qw;
    qw = t1_density * b_x[5];
    W_dt[2] = qw;
    W_dt[6] = n_idx_2;
    W_dt[10] = qz;
    W_dt[14] = expl_temp;
    W_dt[3] = b_expl_temp;
    W_dt[7] = qw;
    W_dt[11] = qy;
    W_dt[15] = qz;
    memset(&aBuffer[0], 0, sizeof(double) << 4);
    c_eye(dv2);
    memset(&b_W_dt[0], 0, sizeof(double) << 4);
    memset(&b_q[0], 0, sizeof(double) << 4);
    memset(&c_W_dt[0], 0, sizeof(double) << 4);
    for (i1 = 0; i1 < 4; i1++) {
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        qz = W_dt[i + b_w_exp_tilde_tmp];
        c_w_exp_tilde_tmp = i << 2;
        qw = W_dt[c_w_exp_tilde_tmp] * qz;
        aBuffer[b_w_exp_tilde_tmp] += qw;
        b_W_dt[b_w_exp_tilde_tmp] += qw;
        b_q[b_w_exp_tilde_tmp] += qw;
        qw = W_dt[c_w_exp_tilde_tmp + 1] * qz;
        aBuffer[b_w_exp_tilde_tmp + 1] += qw;
        b_W_dt[b_w_exp_tilde_tmp + 1] += qw;
        b_q[b_w_exp_tilde_tmp + 1] += qw;
        qw = W_dt[c_w_exp_tilde_tmp + 2] * qz;
        aBuffer[b_w_exp_tilde_tmp + 2] += qw;
        b_W_dt[b_w_exp_tilde_tmp + 2] += qw;
        b_q[b_w_exp_tilde_tmp + 2] += qw;
        qw = W_dt[c_w_exp_tilde_tmp + 3] * qz;
        aBuffer[b_w_exp_tilde_tmp + 3] += qw;
        b_W_dt[b_w_exp_tilde_tmp + 3] += qw;
        b_q[b_w_exp_tilde_tmp + 3] += qw;
      }
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        qw = b_q[i + b_w_exp_tilde_tmp];
        c_w_exp_tilde_tmp = i << 2;
        c_W_dt[b_w_exp_tilde_tmp] += W_dt[c_w_exp_tilde_tmp] * qw;
        c_W_dt[b_w_exp_tilde_tmp + 1] += W_dt[c_w_exp_tilde_tmp + 1] * qw;
        c_W_dt[b_w_exp_tilde_tmp + 2] += W_dt[c_w_exp_tilde_tmp + 2] * qw;
        c_W_dt[b_w_exp_tilde_tmp + 3] += W_dt[c_w_exp_tilde_tmp + 3] * qw;
      }
    }
    memset(&b_q[0], 0, sizeof(double) << 4);
    for (i1 = 0; i1 < 4; i1++) {
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        qw = aBuffer[i + b_w_exp_tilde_tmp];
        c_w_exp_tilde_tmp = i << 2;
        b_q[b_w_exp_tilde_tmp] += aBuffer[c_w_exp_tilde_tmp] * qw;
        b_q[b_w_exp_tilde_tmp + 1] += aBuffer[c_w_exp_tilde_tmp + 1] * qw;
        b_q[b_w_exp_tilde_tmp + 2] += aBuffer[c_w_exp_tilde_tmp + 2] * qw;
        b_q[b_w_exp_tilde_tmp + 3] += aBuffer[c_w_exp_tilde_tmp + 3] * qw;
      }
      b_w_exp_tilde_tmp = i1 << 2;
      F[11 * i1] = (((dv2[b_w_exp_tilde_tmp] + W_dt[b_w_exp_tilde_tmp]) +
                     0.5 * b_W_dt[b_w_exp_tilde_tmp]) +
                    0.16666666666666666 * c_W_dt[b_w_exp_tilde_tmp]) +
                   0.041666666666666664 * b_q[b_w_exp_tilde_tmp];
      F[11 * i1 + 1] =
          (((dv2[b_w_exp_tilde_tmp + 1] + W_dt[b_w_exp_tilde_tmp + 1]) +
            0.5 * b_W_dt[b_w_exp_tilde_tmp + 1]) +
           0.16666666666666666 * c_W_dt[b_w_exp_tilde_tmp + 1]) +
          0.041666666666666664 * b_q[b_w_exp_tilde_tmp + 1];
      F[11 * i1 + 2] =
          (((dv2[b_w_exp_tilde_tmp + 2] + W_dt[b_w_exp_tilde_tmp + 2]) +
            0.5 * b_W_dt[b_w_exp_tilde_tmp + 2]) +
           0.16666666666666666 * c_W_dt[b_w_exp_tilde_tmp + 2]) +
          0.041666666666666664 * b_q[b_w_exp_tilde_tmp + 2];
      F[11 * i1 + 3] =
          (((dv2[b_w_exp_tilde_tmp + 3] + W_dt[b_w_exp_tilde_tmp + 3]) +
            0.5 * b_W_dt[b_w_exp_tilde_tmp + 3]) +
           0.16666666666666666 * c_W_dt[b_w_exp_tilde_tmp + 3]) +
          0.041666666666666664 * b_q[b_w_exp_tilde_tmp + 3];
    }
    /*  + 1/2 * Q_dt^2 + 1/6 * Q_dt^3 + 1/24 * Q_dt^4; */
    qz = t1_density * q[0];
    b_q[0] = qz;
    qy = t1_density * -q[1];
    b_q[4] = qy;
    n_idx_2 = t1_density * -q[2];
    b_q[8] = n_idx_2;
    qw = t1_density * -q[3];
    b_q[12] = qw;
    expl_temp = t1_density * q[1];
    b_q[1] = expl_temp;
    b_q[5] = qz;
    b_q[9] = qw;
    qw = t1_density * q[2];
    b_q[13] = qw;
    b_q[2] = qw;
    qw = t1_density * q[3];
    b_q[6] = qw;
    b_q[10] = qz;
    b_q[14] = qy;
    b_q[3] = qw;
    b_q[7] = n_idx_2;
    b_q[11] = expl_temp;
    b_q[15] = qz;
    for (i = 0; i < 3; i++) {
      b_w_exp_tilde_tmp = (i + 1) << 2;
      c_w_exp_tilde_tmp = 11 * (i + 4);
      F[c_w_exp_tilde_tmp] = b_q[b_w_exp_tilde_tmp];
      F[c_w_exp_tilde_tmp + 1] = b_q[b_w_exp_tilde_tmp + 1];
      F[c_w_exp_tilde_tmp + 2] = b_q[b_w_exp_tilde_tmp + 2];
      F[c_w_exp_tilde_tmp + 3] = b_q[b_w_exp_tilde_tmp + 3];
    }
    /*  column q (attitude) */
    /*  column w (rates) */
    /*     %% angular rate rows (w, 5:7) */
    /*  aerodynamics partial derivatives */
    /* %% air data  */
    /* torque_vx = Cl * delta * param.c_canard * [v(1), v(2), v(3); 0, 0, 0; 0,
     * 0, 0]; */
    /* torque_v =  airdata.density * (torque_vx + torque_vyz); */
    airdata_atmos(b_x[10], &qw, &t1_density, &qz, &qy, &n_idx_2);
    /*  computes matrix exponential of rotation */
    /* %% incremental angle */
    /* %% normed skew-symmetric matrix */
    if (dphi_tmp == 0.0) {
      qz = 0.0;
      qy = 0.0;
      n_idx_2 = 0.0;
    } else {
      qz = b_x[4] / dphi_tmp;
      qy = b_x[5] / dphi_tmp;
      n_idx_2 = b_x[6] / dphi_tmp;
    }
    /*  skew symmetric matrix / cross-product jacobian */
    ST[0] = 0.0;
    ST[3] = -n_idx_2;
    ST[6] = qy;
    ST[1] = n_idx_2;
    ST[4] = 0.0;
    ST[7] = -qz;
    ST[2] = -qy;
    ST[5] = qz;
    ST[8] = 0.0;
    /* %% Rodrigues rotation formula */
    memset(&w_exp_tilde[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qw = w_exp_tilde[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        qz = ST[i1 + 3 * i];
        qw += ST[3 * i1] * qz;
        w_exp_tilde[b_w_exp_tilde_tmp] += ST[3 * i1 + 1] * qz;
        w_exp_tilde[c_w_exp_tilde_tmp] += ST[3 * i1 + 2] * qz;
      }
      w_exp_tilde[3 * i] = qw;
    }
    for (i = 0; i < 9; i++) {
      w_exp_tilde[i] =
          (w_exp_tilde_tmp[i] - b_a * ST[i]) + (1.0 - x) * w_exp_tilde[i];
    }
    memset(&dv3[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qw = dv3[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        qz = w_exp_tilde[i1 + 3 * i];
        qw += b_param.Jinv[3 * i1] * qz;
        dv3[b_w_exp_tilde_tmp] += b_param.Jinv[3 * i1 + 1] * qz;
        dv3[c_w_exp_tilde_tmp] += b_param.Jinv[3 * i1 + 2] * qz;
        F[(i1 + 11 * (i + 4)) + 4] = 0.0;
      }
      dv3[3 * i] = qw;
    }
    c_a[1] = t1_density * (-0.08091368228861362 * b_x[9]);
    qz = t1_density * -0.0;
    c_a[4] = qz;
    c_a[7] = t1_density * (-0.08091368228861362 * b_x[7]);
    c_a[2] = t1_density * (-0.08091368228861362 * -b_x[8]);
    c_a[5] = t1_density * (-0.08091368228861362 * -b_x[7]);
    c_a[8] = qz;
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
    qy = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      c_a[3 * i1] = qz;
      b_w_exp_tilde_tmp = 11 * (i1 + 4);
      for (i = 0; i < 3; i++) {
        qw = b_param.J[i + 3 * i1];
        F[b_w_exp_tilde_tmp + 4] += dv3[3 * i] * qw;
        F[b_w_exp_tilde_tmp + 5] += dv3[3 * i + 1] * qw;
        F[b_w_exp_tilde_tmp + 6] += dv3[3 * i + 2] * qw;
        F[(i + 11 * (i1 + 7)) + 4] = 0.0;
      }
      b_w_exp_tilde_tmp = 11 * (i1 + 7);
      for (i = 0; i < 3; i++) {
        qw = c_a[i + 3 * i1];
        F[b_w_exp_tilde_tmp + 4] += dt * b_param.Jinv[3 * i] * qw;
        F[b_w_exp_tilde_tmp + 5] += dt * b_param.Jinv[3 * i + 1] * qw;
        F[b_w_exp_tilde_tmp + 6] += dt * b_param.Jinv[3 * i + 2] * qw;
      }
      qy += b_x[i1 + 1] * b_param.g[i1];
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
    qw = b_x[0] * 0.0;
    c_a[0] = qw;
    c_a[3] = b_x[0] * -b_param.g[2];
    c_a[6] = b_x[0] * b_param.g[1];
    c_a[1] = b_x[0] * b_param.g[2];
    c_a[4] = qw;
    c_a[7] = b_x[0] * -b_param.g[0];
    c_a[2] = b_x[0] * -b_param.g[1];
    c_a[5] = b_x[0] * b_param.g[0];
    c_a[8] = qw;
    for (i = 0; i < 3; i++) {
      F[i + 7] = dt * (2.0 * (b_x[0] * b_param.g[i] - dv1[i]));
      qw = b_x[i + 1];
      c_w_exp_tilde_tmp = 11 * (i + 1);
      F[c_w_exp_tilde_tmp + 7] =
          dt * (2.0 * (((qy * b_b[3 * i] + b_x[1] * b_param.g[i]) -
                        b_param.g[0] * qw) +
                       c_a[3 * i]));
      b_w_exp_tilde_tmp = 3 * i + 1;
      F[c_w_exp_tilde_tmp + 8] =
          dt * (2.0 * (((qy * b_b[b_w_exp_tilde_tmp] + b_x[2] * b_param.g[i]) -
                        b_param.g[1] * qw) +
                       c_a[b_w_exp_tilde_tmp]));
      b_w_exp_tilde_tmp = 3 * i + 2;
      F[c_w_exp_tilde_tmp + 9] =
          dt * (2.0 * (((qy * b_b[b_w_exp_tilde_tmp] + b_x[3] * b_param.g[i]) -
                        b_param.g[2] * qw) +
                       c_a[b_w_exp_tilde_tmp]));
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
    expl_temp = 0.5 * (dt * dt);
    memset(&c_a[0], 0, 9U * sizeof(double));
    memset(&dv3[0], 0, 9U * sizeof(double));
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
    b_expl_temp = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      int d_w_exp_tilde_tmp;
      qw = c_a[3 * i1];
      qz = dv3[3 * i1];
      c_w_exp_tilde_tmp = 3 * i1 + 1;
      d_w_exp_tilde_tmp = 3 * i1 + 2;
      for (i = 0; i < 3; i++) {
        b_w_exp_tilde_tmp = i + 3 * i1;
        qy = w_exp_tilde_tmp[b_w_exp_tilde_tmp];
        n_idx_2 = ST[b_w_exp_tilde_tmp];
        qw += ST[3 * i] * qy;
        qz += 2.0 * w_exp_tilde_tmp[3 * i] * n_idx_2;
        b_w_exp_tilde_tmp = 3 * i + 1;
        c_a[c_w_exp_tilde_tmp] += ST[b_w_exp_tilde_tmp] * qy;
        dv3[c_w_exp_tilde_tmp] +=
            2.0 * w_exp_tilde_tmp[b_w_exp_tilde_tmp] * n_idx_2;
        b_w_exp_tilde_tmp = 3 * i + 2;
        c_a[d_w_exp_tilde_tmp] += ST[b_w_exp_tilde_tmp] * qy;
        dv3[d_w_exp_tilde_tmp] +=
            2.0 * w_exp_tilde_tmp[b_w_exp_tilde_tmp] * n_idx_2;
      }
      dv3[3 * i1] = qz;
      c_a[3 * i1] = qw;
      c_w_exp_tilde_tmp = 11 * (i1 + 4);
      F[c_w_exp_tilde_tmp + 7] = dt * ST[3 * i1] + expl_temp * (qw - qz);
      d_w_exp_tilde_tmp = 11 * (i1 + 7);
      F[d_w_exp_tilde_tmp + 7] = w_exp_tilde[3 * i1];
      b_w_exp_tilde_tmp = 3 * i1 + 1;
      F[c_w_exp_tilde_tmp + 8] =
          dt * ST[b_w_exp_tilde_tmp] +
          expl_temp * (c_a[b_w_exp_tilde_tmp] - dv3[b_w_exp_tilde_tmp]);
      F[d_w_exp_tilde_tmp + 8] = w_exp_tilde[b_w_exp_tilde_tmp];
      b_w_exp_tilde_tmp = 3 * i1 + 2;
      F[c_w_exp_tilde_tmp + 9] =
          dt * ST[b_w_exp_tilde_tmp] +
          expl_temp * (c_a[b_w_exp_tilde_tmp] - dv3[b_w_exp_tilde_tmp]);
      F[d_w_exp_tilde_tmp + 9] = w_exp_tilde[b_w_exp_tilde_tmp];
      qw = -b_x[i1 + 1];
      q[i1 + 1] = qw;
      b_expl_temp += qw * b_x[i1 + 7];
    }
    /*  Sola */
    /* %% for hardcoding: */
    /*  [2*qw*vx - 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz, 2*qx*vy -
     * 2*qw*vz - 2*qy*vx, 2*qw*vy + 2*qx*vz - 2*qz*vx] */
    /*  [2*qw*vy + 2*qx*vz - 2*qz*vx, 2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qx*vx +
     * 2*qy*vy + 2*qz*vz, 2*qy*vz - 2*qw*vx - 2*qz*vy] */
    /*  [2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qz*vx - 2*qx*vz - 2*qw*vy, 2*qw*vx -
     * 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz] */
    dn[0] = q[2] * b_x[9] - q[3] * b_x[8];
    dn[1] = q[3] * b_x[7] - q[1] * b_x[9];
    dn[2] = q[1] * b_x[8] - q[2] * b_x[7];
    for (i = 0; i < 3; i++) {
      qw = b_x[i + 7];
      c_dt[i] = dt * (2.0 * (q[0] * qw - dn[i]));
      qz = q[i + 1];
      c_w_exp_tilde_tmp = 3 * (i + 1);
      c_dt[c_w_exp_tilde_tmp] =
          dt * (2.0 * (((b_expl_temp * b_b[3 * i] + q[1] * qw) - b_x[7] * qz) +
                       q[0] * ST[3 * i]));
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_dt[c_w_exp_tilde_tmp + 1] =
          dt * (2.0 * (((b_expl_temp * b_b[b_w_exp_tilde_tmp] + q[2] * qw) -
                        b_x[8] * qz) +
                       q[0] * ST[b_w_exp_tilde_tmp]));
      b_w_exp_tilde_tmp = 3 * i + 2;
      c_dt[c_w_exp_tilde_tmp + 2] =
          dt * (2.0 * (((b_expl_temp * b_b[b_w_exp_tilde_tmp] + q[3] * qw) -
                        b_x[9] * qz) +
                       q[0] * ST[b_w_exp_tilde_tmp]));
    }
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
    qw = q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + q[3] * q[3]);
    qz = 2.0 * q[0];
    /*  Stevens */
    /* %% for hardcoding:  */
    /*  [qw^2 + qx^2 - qy^2 - qz^2,         2*qw*qz + 2*qx*qy,         2*qx*qz -
     * 2*qw*qy] */
    /*  [        2*qx*qy - 2*qw*qz, qw^2 - qx^2 + qy^2 - qz^2,         2*qw*qx +
     * 2*qy*qz] */
    /*  [        2*qw*qy + 2*qx*qz,         2*qy*qz - 2*qw*qx, qw^2 - qx^2 -
     * qy^2 + qz^2] */
    qy = qz * 0.0;
    c_a[0] = qy;
    c_a[3] = qz * -q[3];
    c_a[6] = qz * q[2];
    c_a[1] = qz * q[3];
    c_a[4] = qy;
    c_a[7] = qz * -q[1];
    c_a[2] = qz * -q[2];
    c_a[5] = qz * q[1];
    c_a[8] = qy;
    for (i = 0; i < 3; i++) {
      F[11 * (i + 7) + 10] =
          dt * ((qw * b_b[3 * i] + 2.0 * q[1] * q[i + 1]) - c_a[3 * i]);
    }
    /*  r_r = eye(3); */
    /*  column q */
    /*  column v */
    F[120] = 1.0;
    /*  column alt */
    /* %% discrete covariance */
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      for (i1 = 0; i1 < 11; i1++) {
        qw = P[i1 + 11 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_F[b_w_exp_tilde_tmp] += F[i2 + 11 * i1] * qw;
        }
      }
    }
    for (i1 = 0; i1 < 11; i1++) {
      for (i2 = 0; i2 < 11; i2++) {
        qw = 0.0;
        for (i = 0; i < 11; i++) {
          qw += b_F[i1 + 11 * i] * F[i2 + 11 * i];
        }
        b_w_exp_tilde_tmp = i1 + 11 * i2;
        P_pred[b_w_exp_tilde_tmp] = qw + Q[b_w_exp_tilde_tmp];
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
    for (i = 0; i < 3; i++) {
      c_w_exp_tilde_tmp = 11 * (i + 4);
      c_a[3 * i] = P_pred[c_w_exp_tilde_tmp + 4] + R[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_a[b_w_exp_tilde_tmp] =
          P_pred[c_w_exp_tilde_tmp + 5] + R[b_w_exp_tilde_tmp];
      b_w_exp_tilde_tmp = 3 * i + 2;
      c_a[b_w_exp_tilde_tmp] =
          P_pred[c_w_exp_tilde_tmp + 6] + R[b_w_exp_tilde_tmp];
    }
    inv(c_a, ST);
    memset(&K[0], 0, 33U * sizeof(double));
    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        qw = ST[i1 + 3 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          K[b_w_exp_tilde_tmp] += P_pred[i2 + 11 * (i1 + 4)] * qw;
        }
      }
    }
    /*  hardcoded P*H */
    d_eye(F);
    memcpy(&P[0], &F[0], 44U * sizeof(double));
    for (i = 0; i < 33; i++) {
      P[i + 44] = F[i + 44] - K[i];
    }
    memcpy(&P[77], &F[77], 44U * sizeof(double));
    /*  hardcoded K*H */
    /* %% correct covariance estimate */
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      for (i1 = 0; i1 < 11; i1++) {
        qw = P_pred[i1 + 11 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_F[b_w_exp_tilde_tmp] += P[i2 + 11 * i1] * qw;
        }
      }
    }
    memset(&b_K[0], 0, 33U * sizeof(double));
    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        qw = R[i1 + 3 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_K[b_w_exp_tilde_tmp] += K[i2 + 11 * i1] * qw;
        }
      }
    }
    memset(&F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      for (i1 = 0; i1 < 11; i1++) {
        qw = P[i + 11 * i1];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          F[b_w_exp_tilde_tmp] += b_F[i2 + 11 * i1] * qw;
        }
      }
    }
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        qw = K[i + 11 * i1];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_F[b_w_exp_tilde_tmp] += b_K[i2 + 11 * i1] * qw;
        }
      }
    }
    for (i = 0; i < 121; i++) {
      P[i] = F[i] + b_F[i];
    }
    /*  joseph stabilized */
    /* %% correct state estimate */
    qy = ((d1 / d4 * (board_gyro->meas[0] - c_b.board_gyro[0]) +
           d2 / d4 * (mti_gyro->meas[0] - c_b.mti_gyro[0])) +
          a[0] * (ad_gyro->meas[0] - c_b.ad_gyro[0])) -
         xhat[4];
    qz = d1 / d3;
    qw = d2 / d3;
    n_idx_2 = ((qz * (board_gyro->meas[1] - c_b.board_gyro[1]) +
                qw * (mti_gyro->meas[1] - c_b.mti_gyro[1])) +
               d * (ad_gyro->meas[1] - c_b.ad_gyro[1])) -
              xhat[5];
    qz = ((qz * (board_gyro->meas[2] - c_b.board_gyro[2]) +
           qw * (mti_gyro->meas[2] - c_b.mti_gyro[2])) +
          d * (ad_gyro->meas[2] - c_b.ad_gyro[2])) -
         xhat[6];
    for (i = 0; i < 11; i++) {
      b_x[i] = xhat[i] + ((K[i] * qy + K[i + 11] * n_idx_2) + K[i + 22] * qz);
    }
    /*  norms quaternion */
    /* %% inverse quaternion  */
    qz = c_norm(&b_x[0]);
    b_x[0] /= qz;
    b_x[1] /= qz;
    b_x[2] /= qz;
    b_x[3] /= qz;
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
  qz = b_norm(&b_x[7]);
  /*  return values */
  airdata->mach = qz / airdata->sonic_speed;
  /*  [-], Mach number */
  airdata->dynamic_pressure = 0.5 * airdata->density * (qz * qz);
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
