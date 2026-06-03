/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: navigation_codegen_entry.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 23:24:33
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

/* Function Definitions */
/*
 * Calls the pad and flight filters.
 *
 * Arguments    : double dt
 *                bool flight_phase
 *                const double x[11]
 *                const double P[121]
 *                struct1_T *bias
 *                const struct2_T *sens_filt
 *                const struct3_T *sens_input
 *                double x_ret[11]
 *                double P_ret[121]
 *                struct1_T *bias_ret
 *                struct2_T *sens_filt_ret
 * Return Type  : void
 */
void navigation_codegen_entry(double dt, bool flight_phase, const double x[11],
                              const double P[121], struct1_T *bias,
                              const struct2_T *sens_filt,
                              const struct3_T *sens_input, double x_ret[11],
                              double P_ret[121], struct1_T *bias_ret,
                              struct2_T *sens_filt_ret)
{
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
  static const double b[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double E[121];
  double F[121];
  double b_F[121];
  double K[33];
  double b_K[33];
  double aBuffer[16];
  double b_W_dt[16];
  double b_q[16];
  double c_W_dt[16];
  double x_pred[11];
  double ST[9];
  double c_a[9];
  double dv3[9];
  double w_exp_tilde[9];
  double c_q[4];
  double q[4];
  double a[3];
  double ad_accel_f[3];
  double b_dt[3];
  double b_w_exp_tilde[3];
  double board_accel_f[3];
  double dv1[3];
  double board_baro_f;
  double board_gyro_f_idx_0;
  double board_gyro_f_idx_1;
  double board_gyro_f_idx_2;
  double expl_temp;
  double mti_accel_f_idx_0;
  double mti_accel_f_idx_1;
  double mti_accel_f_idx_2;
  double mti_baro_f;
  double qw;
  double t1_density;
  int i;
  int i1;
  int i2;
  if (!isInitialized_GNC_codegen) {
    GNC_codegen_initialize();
  }
  memcpy(&x_ret[0], &x[0], 11U * sizeof(double));
  memcpy(&P_ret[0], &P[0], 121U * sizeof(double));
  *bias_ret = *bias;
  *sens_filt_ret = *sens_filt;
  /*     %% Pad filter iteration */
  if (!flight_phase) {
    /*  only before ignition (or if not run before) */
    /*  Computes on pad: inital state for flight filter, and bias values for the
     * sensors */
    /*  Outputs: initial state x, sensor biases bias */
    /*     %% Settings */
    /*  [s], low pass time constant */
    /*     %% parameters */
    board_accel_f[0] = sens_filt->board_accel_f[0];
    board_gyro_f_idx_0 = sens_filt->board_gyro_f[0];
    mti_accel_f_idx_0 = sens_filt->mti_accel_f[0];
    ad_accel_f[0] = sens_filt->ad_accel_f[0];
    board_accel_f[1] = sens_filt->board_accel_f[1];
    board_gyro_f_idx_1 = sens_filt->board_gyro_f[1];
    mti_accel_f_idx_1 = sens_filt->mti_accel_f[1];
    ad_accel_f[1] = sens_filt->ad_accel_f[1];
    board_accel_f[2] = sens_filt->board_accel_f[2];
    board_gyro_f_idx_2 = sens_filt->board_gyro_f[2];
    mti_accel_f_idx_2 = sens_filt->mti_accel_f[2];
    ad_accel_f[2] = sens_filt->ad_accel_f[2];
    board_baro_f = sens_filt->board_baro_f;
    mti_baro_f = sens_filt->mti_baro_f;
    /*     %% lowpass filter */
    /*  filtered = filtered + alpha*(measured-filtered); */
    /* %% lowpass filter function used in pad filter */
    if (sens_input->board_accel.status) {
      board_accel_f[0] = 0.0005 * sens_input->board_accel.meas[0] +
                         0.9995 * sens_filt->board_accel_f[0];
      board_accel_f[1] = 0.0005 * sens_input->board_accel.meas[1] +
                         0.9995 * sens_filt->board_accel_f[1];
      board_accel_f[2] = 0.0005 * sens_input->board_accel.meas[2] +
                         0.9995 * sens_filt->board_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->board_gyro.status) {
      board_gyro_f_idx_0 = 0.0005 * sens_input->board_gyro.meas[0] +
                           0.9995 * sens_filt->board_gyro_f[0];
      board_gyro_f_idx_1 = 0.0005 * sens_input->board_gyro.meas[1] +
                           0.9995 * sens_filt->board_gyro_f[1];
      board_gyro_f_idx_2 = 0.0005 * sens_input->board_gyro.meas[2] +
                           0.9995 * sens_filt->board_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->mti_accel.status) {
      mti_accel_f_idx_0 = 0.0005 * sens_input->mti_accel.meas[0] +
                          0.9995 * sens_filt->mti_accel_f[0];
      mti_accel_f_idx_1 = 0.0005 * sens_input->mti_accel.meas[1] +
                          0.9995 * sens_filt->mti_accel_f[1];
      mti_accel_f_idx_2 = 0.0005 * sens_input->mti_accel.meas[2] +
                          0.9995 * sens_filt->mti_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->mti_gyro.status) {
      sens_filt_ret->mti_gyro_f[0] = 0.0005 * sens_input->mti_gyro.meas[0] +
                                     0.9995 * sens_filt->mti_gyro_f[0];
      sens_filt_ret->mti_gyro_f[1] = 0.0005 * sens_input->mti_gyro.meas[1] +
                                     0.9995 * sens_filt->mti_gyro_f[1];
      sens_filt_ret->mti_gyro_f[2] = 0.0005 * sens_input->mti_gyro.meas[2] +
                                     0.9995 * sens_filt->mti_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->ad_accel.status) {
      ad_accel_f[0] = 0.0005 * sens_input->ad_accel.meas[0] +
                      0.9995 * sens_filt->ad_accel_f[0];
      ad_accel_f[1] = 0.0005 * sens_input->ad_accel.meas[1] +
                      0.9995 * sens_filt->ad_accel_f[1];
      ad_accel_f[2] = 0.0005 * sens_input->ad_accel.meas[2] +
                      0.9995 * sens_filt->ad_accel_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->ad_gyro.status) {
      sens_filt_ret->ad_gyro_f[0] = 0.0005 * sens_input->ad_gyro.meas[0] +
                                    0.9995 * sens_filt->ad_gyro_f[0];
      sens_filt_ret->ad_gyro_f[1] = 0.0005 * sens_input->ad_gyro.meas[1] +
                                    0.9995 * sens_filt->ad_gyro_f[1];
      sens_filt_ret->ad_gyro_f[2] = 0.0005 * sens_input->ad_gyro.meas[2] +
                                    0.9995 * sens_filt->ad_gyro_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->board_baro.status) {
      board_baro_f = 0.0005 * sens_input->board_baro.meas +
                     0.9995 * sens_filt->board_baro_f;
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->board_mag.status) {
      sens_filt_ret->board_mag_f[0] = 0.0005 * sens_input->board_mag.meas[0] +
                                      0.9995 * sens_filt->board_mag_f[0];
      sens_filt_ret->board_mag_f[1] = 0.0005 * sens_input->board_mag.meas[1] +
                                      0.9995 * sens_filt->board_mag_f[1];
      sens_filt_ret->board_mag_f[2] = 0.0005 * sens_input->board_mag.meas[2] +
                                      0.9995 * sens_filt->board_mag_f[2];
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->mti_baro.status) {
      mti_baro_f =
          0.0005 * sens_input->mti_baro.meas + 0.9995 * sens_filt->mti_baro_f;
    }
    /* %% lowpass filter function used in pad filter */
    if (sens_input->mti_mag.status) {
      sens_filt_ret->mti_mag_f[0] = 0.0005 * sens_input->mti_mag.meas[0] +
                                    0.9995 * sens_filt->mti_mag_f[0];
      sens_filt_ret->mti_mag_f[1] = 0.0005 * sens_input->mti_mag.meas[1] +
                                    0.9995 * sens_filt->mti_mag_f[1];
      sens_filt_ret->mti_mag_f[2] = 0.0005 * sens_input->mti_mag.meas[2] +
                                    0.9995 * sens_filt->mti_mag_f[2];
    }
    sens_filt_ret->board_baro_f = board_baro_f;
    sens_filt_ret->mti_baro_f = mti_baro_f;
    /*     %% Initial state determination     */
    /* %% Orientation */
    sens_filt_ret->board_accel_f[0] = board_accel_f[0];
    sens_filt_ret->board_gyro_f[0] = board_gyro_f_idx_0;
    sens_filt_ret->mti_accel_f[0] = mti_accel_f_idx_0;
    sens_filt_ret->ad_accel_f[0] = ad_accel_f[0];
    a[0] = 0.0;
    sens_filt_ret->board_accel_f[1] = board_accel_f[1];
    sens_filt_ret->board_gyro_f[1] = board_gyro_f_idx_1;
    sens_filt_ret->mti_accel_f[1] = mti_accel_f_idx_1;
    sens_filt_ret->ad_accel_f[1] = ad_accel_f[1];
    a[1] = 0.0;
    sens_filt_ret->board_accel_f[2] = board_accel_f[2];
    sens_filt_ret->board_gyro_f[2] = board_gyro_f_idx_2;
    sens_filt_ret->mti_accel_f[2] = mti_accel_f_idx_2;
    sens_filt_ret->ad_accel_f[2] = ad_accel_f[2];
    a[2] = 0.0;
    /*  acceleration  */
    if (sens_input->board_accel.status) {
      /*  only add alive IMUs to average */
      a[0] = board_accel_f[0];
      a[1] = board_accel_f[1];
      a[2] = board_accel_f[2];
    }
    if (sens_input->mti_accel.status) {
      a[0] += mti_accel_f_idx_0;
      a[1] += mti_accel_f_idx_1;
      a[2] += mti_accel_f_idx_2;
    }
    if (sens_input->ad_accel.status) {
      a[0] += ad_accel_f[0];
      a[1] += ad_accel_f[1];
      a[2] += ad_accel_f[2];
    }
    /* %% computes initial orientation of stationary body from gravity
     * acceleration  */
    /* %% Input: 3D acceleration vector */
    /* %% Output: Orientation quaternion */
    /* %% normed gravity vector in body-fixed frame */
    mti_accel_f_idx_2 = b_norm(a) + 1.0E-6;
    a[0] /= mti_accel_f_idx_2;
    a[1] /= mti_accel_f_idx_2;
    a[2] /= mti_accel_f_idx_2;
    /*  unit vector of gravity direction */
    /* %% determine initial orientation quaternion */
    qw = sqrt(0.5 * a[0] + 0.5);
    if (qw == 0.0) {
      /*  exact upside down case */
      mti_accel_f_idx_0 = 1.0;
      /*  either qy = 1 or qz = 1, this is arbitrary  */
      mti_accel_f_idx_2 = 0.0;
    } else {
      mti_accel_f_idx_0 = 0.5 * a[2] / qw;
      mti_accel_f_idx_2 = -0.5 * a[1] / qw;
    }
    q[0] = qw;
    q[1] = 0.0;
    q[2] = mti_accel_f_idx_0;
    q[3] = mti_accel_f_idx_2;
    mti_accel_f_idx_1 = c_norm(q);
    /*  a gets normed inside function */
    /* %% launch altitude */
    /* %% set constant initials */
    /*  stationary on rail */
    /* %% conconct state vector */
    expl_temp = qw / mti_accel_f_idx_1;
    q[0] = expl_temp;
    x_ret[0] = expl_temp;
    expl_temp = 0.0 / mti_accel_f_idx_1;
    q[1] = expl_temp;
    x_ret[1] = expl_temp;
    expl_temp = mti_accel_f_idx_0 / mti_accel_f_idx_1;
    q[2] = expl_temp;
    x_ret[2] = expl_temp;
    expl_temp = mti_accel_f_idx_2 / mti_accel_f_idx_1;
    q[3] = expl_temp;
    x_ret[3] = expl_temp;
    x_ret[10] = 420.0;
    /*     %% Bias determination */
    /* %% gyroscope */
    /* %% earth magnetic field */
    /*  computes rotation matrix from quaternion */
    /* %% norm quaternions */
    /*  q = quaternion_norm(q); */
    /* %% quaternion definition */
    /* %% skew symetric quaternion matrix */
    /* %% rotation matrix */
    x_ret[4] = 0.0;
    x_ret[7] = 0.0;
    bias->board_gyro[0] = board_gyro_f_idx_0;
    bias->mti_gyro[0] = sens_filt_ret->mti_gyro_f[0];
    bias->ad_gyro[0] = sens_filt_ret->ad_gyro_f[0];
    x_ret[5] = 0.0;
    x_ret[8] = 0.0;
    bias->board_gyro[1] = board_gyro_f_idx_1;
    bias->mti_gyro[1] = sens_filt_ret->mti_gyro_f[1];
    bias->ad_gyro[1] = sens_filt_ret->ad_gyro_f[1];
    x_ret[6] = 0.0;
    x_ret[9] = 0.0;
    bias->board_gyro[2] = board_gyro_f_idx_2;
    bias->mti_gyro[2] = sens_filt_ret->mti_gyro_f[2];
    bias->ad_gyro[2] = sens_filt_ret->ad_gyro_f[2];
    mti_accel_f_idx_2 =
        q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + expl_temp * expl_temp);
    mti_accel_f_idx_1 = 2.0 * q[0];
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
      ST[3 * i] = mti_accel_f_idx_2 * b[i] + qw * q[1];
      ST[3 * i + 1] = mti_accel_f_idx_2 * b[i + 3] + qw * q[2];
      ST[3 * i + 2] = mti_accel_f_idx_2 * b[i + 6] + qw * expl_temp;
    }
    mti_accel_f_idx_2 = mti_accel_f_idx_1 * 0.0;
    c_a[0] = mti_accel_f_idx_2;
    c_a[1] = mti_accel_f_idx_1 * -expl_temp;
    c_a[2] = mti_accel_f_idx_1 * q[2];
    c_a[3] = mti_accel_f_idx_1 * expl_temp;
    c_a[4] = mti_accel_f_idx_2;
    c_a[5] = mti_accel_f_idx_1 * -q[1];
    c_a[6] = mti_accel_f_idx_1 * -q[2];
    c_a[7] = mti_accel_f_idx_1 * q[1];
    c_a[8] = mti_accel_f_idx_2;
    for (i = 0; i < 9; i++) {
      ST[i] -= c_a[i];
    }
    /*  launch attitude */
    bias->board_mag_earth[0] = 0.0;
    bias->board_mag_earth[1] = 0.0;
    bias->board_mag_earth[2] = 0.0;
    for (i = 0; i < 3; i++) {
      mti_accel_f_idx_2 = sens_filt_ret->board_mag_f[i];
      bias->board_mag_earth[0] += ST[3 * i] * mti_accel_f_idx_2;
      bias->board_mag_earth[1] += ST[3 * i + 1] * mti_accel_f_idx_2;
      bias->board_mag_earth[2] += ST[3 * i + 2] * mti_accel_f_idx_2;
      bias->mti_mag_earth[i] = 0.0;
    }
    mti_accel_f_idx_2 = bias->mti_mag_earth[0];
    qw = bias->mti_mag_earth[1];
    mti_accel_f_idx_1 = bias->mti_mag_earth[2];
    for (i = 0; i < 3; i++) {
      mti_accel_f_idx_0 = sens_filt_ret->mti_mag_f[i];
      mti_accel_f_idx_2 += ST[3 * i] * mti_accel_f_idx_0;
      qw += ST[3 * i + 1] * mti_accel_f_idx_0;
      mti_accel_f_idx_1 += ST[3 * i + 2] * mti_accel_f_idx_0;
    }
    bias->mti_mag_earth[2] = mti_accel_f_idx_1;
    bias->mti_mag_earth[1] = qw;
    bias->mti_mag_earth[0] = mti_accel_f_idx_2;
    /* %% barometer */
    mti_accel_f_idx_2 =
        airdata_atmos(420.0, &mti_accel_f_idx_2, &t1_density, &qw,
                      &mti_accel_f_idx_0, &mti_accel_f_idx_1);
    /*  what the pressure should be at launch elevation */
    bias->board_baro = board_baro_f - mti_accel_f_idx_2;
    bias->mti_baro = mti_baro_f - mti_accel_f_idx_2;
    *bias_ret = *bias;
  }
  /*     %% Flight filter iteration */
  if (flight_phase) {
    double P_pred[121];
    double W_dt[16];
    double dv2[16];
    double c_dt[12];
    double w_exp_tilde_tmp[9];
    double dv[4];
    double b_a;
    double b_board_accel_f_tmp_tmp;
    double b_x;
    double board_accel_f_tmp_tmp;
    double c_board_accel_f_tmp_tmp;
    double d;
    double d1;
    double d2;
    double d3;
    double d4;
    double d5;
    double d6;
    double d7;
    double d8;
    double dphi_tmp;
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
    d = 9.9999999999999981E+9 * (double)sens_input->ad_gyro.status;
    board_accel_f_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_input->board_accel.status;
    b_board_accel_f_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_input->mti_accel.status;
    c_board_accel_f_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_input->ad_accel.status;
    mti_accel_f_idx_1 = (board_accel_f_tmp_tmp + b_board_accel_f_tmp_tmp) +
                        c_board_accel_f_tmp_tmp;
    board_accel_f[0] = mti_accel_f_idx_1;
    d1 = 9.9999999999999981E+9 * (double)sens_input->board_gyro.status;
    d2 = 9.9999999999999981E+9 * (double)sens_input->mti_gyro.status;
    d3 = d1 + d2;
    d4 = d3 + d;
    d /= d4;
    a[0] = d;
    board_accel_f[1] = mti_accel_f_idx_1;
    d = 0.0 / d3;
    board_accel_f[2] = mti_accel_f_idx_1;
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
    mti_accel_f_idx_1 = c_norm(&x[0]);
    q[0] = x[0] / mti_accel_f_idx_1;
    q[1] = x[1] / mti_accel_f_idx_1;
    q[2] = x[2] / mti_accel_f_idx_1;
    q[3] = x[3] / mti_accel_f_idx_1;
    /* %% incremental quaternion */
    dphi_tmp = b_norm(&x[4]);
    mti_accel_f_idx_2 = dphi_tmp * dt;
    board_baro_f = mti_accel_f_idx_2 / 2.0;
    if (dphi_tmp == 0.0) {
      ad_accel_f[0] = 0.0;
      ad_accel_f[1] = 0.0;
      ad_accel_f[2] = 0.0;
      board_gyro_f_idx_0 = 0.0;
      board_gyro_f_idx_1 = 0.0;
      board_gyro_f_idx_2 = 0.0;
    } else {
      board_gyro_f_idx_0 = x[4] / dphi_tmp;
      ad_accel_f[0] = board_gyro_f_idx_0;
      board_gyro_f_idx_1 = x[5] / dphi_tmp;
      ad_accel_f[1] = board_gyro_f_idx_1;
      board_gyro_f_idx_2 = x[6] / dphi_tmp;
      ad_accel_f[2] = board_gyro_f_idx_2;
    }
    mti_baro_f = sin(board_baro_f);
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
    ST[3] = -board_gyro_f_idx_2;
    ST[6] = board_gyro_f_idx_1;
    ST[1] = board_gyro_f_idx_2;
    ST[4] = 0.0;
    ST[7] = -board_gyro_f_idx_0;
    ST[2] = -board_gyro_f_idx_1;
    ST[5] = board_gyro_f_idx_0;
    ST[8] = 0.0;
    /* %% Rodrigues rotation formula */
    b_a = sin(mti_accel_f_idx_2);
    b_x = cos(mti_accel_f_idx_2);
    b_eye(w_exp_tilde_tmp);
    memset(&w_exp_tilde[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      mti_accel_f_idx_2 = w_exp_tilde[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        qw = ST[i1 + 3 * i];
        mti_accel_f_idx_2 += ST[3 * i1] * qw;
        w_exp_tilde[b_w_exp_tilde_tmp] += ST[3 * i1 + 1] * qw;
        w_exp_tilde[c_w_exp_tilde_tmp] += ST[3 * i1 + 2] * qw;
      }
      w_exp_tilde[3 * i] = mti_accel_f_idx_2;
    }
    for (i = 0; i < 9; i++) {
      w_exp_tilde[i] =
          (w_exp_tilde_tmp[i] - b_a * ST[i]) + (1.0 - b_x) * w_exp_tilde[i];
    }
    /*  aerodynamics model */
    /* %% air data   */
    /*  appends airadata struct with dynamic air data from velocity vector */
    /*  dynamic air data: dynamic pressure, mach number */
    mti_accel_f_idx_2 = b_norm(&x[7]);
    /*  return values */
    /*  [-], Mach number */
    /*  [Pa] */
    airdata_atmos(x[10], &qw, &t1_density, &mti_accel_f_idx_1,
                  &mti_accel_f_idx_0, &expl_temp);
    mti_accel_f_idx_0 =
        0.5 * t1_density * (mti_accel_f_idx_2 * mti_accel_f_idx_2);
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
    mti_accel_f_idx_2 =
        x[0] * x[0] - ((x[1] * x[1] + x[2] * x[2]) + x[3] * x[3]);
    mti_accel_f_idx_1 = 2.0 * x[0];
    for (i = 0; i < 3; i++) {
      qw = x[i + 1];
      ST[3 * i] = mti_accel_f_idx_2 * b[3 * i] + 2.0 * x[1] * qw;
      b_w_exp_tilde_tmp = 3 * i + 1;
      ST[b_w_exp_tilde_tmp] =
          mti_accel_f_idx_2 * b[b_w_exp_tilde_tmp] + 2.0 * x[2] * qw;
      b_w_exp_tilde_tmp = 3 * i + 2;
      ST[b_w_exp_tilde_tmp] =
          mti_accel_f_idx_2 * b[b_w_exp_tilde_tmp] + 2.0 * x[3] * qw;
    }
    mti_accel_f_idx_2 = mti_accel_f_idx_1 * 0.0;
    c_a[0] = mti_accel_f_idx_2;
    c_a[3] = mti_accel_f_idx_1 * -x[3];
    c_a[6] = mti_accel_f_idx_1 * x[2];
    c_a[1] = mti_accel_f_idx_1 * x[3];
    c_a[4] = mti_accel_f_idx_2;
    c_a[7] = mti_accel_f_idx_1 * -x[1];
    c_a[2] = mti_accel_f_idx_1 * -x[2];
    c_a[5] = mti_accel_f_idx_1 * x[1];
    c_a[8] = mti_accel_f_idx_2;
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
    dv[0] = cos(board_baro_f);
    memset(&c_a[0], 0, 9U * sizeof(double));
    memset(&b_w_exp_tilde[0], 0, 3U * sizeof(double));
    for (i = 0; i < 3; i++) {
      dv[i + 1] = ad_accel_f[i] * mti_baro_f;
      qw = c_a[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        mti_accel_f_idx_2 = param.J[i1 + 3 * i];
        qw += w_exp_tilde[3 * i1] * mti_accel_f_idx_2;
        c_a[b_w_exp_tilde_tmp] += w_exp_tilde[3 * i1 + 1] * mti_accel_f_idx_2;
        c_a[c_w_exp_tilde_tmp] += w_exp_tilde[3 * i1 + 2] * mti_accel_f_idx_2;
      }
      c_a[3 * i] = qw;
      mti_accel_f_idx_2 = x[i + 4];
      b_w_exp_tilde[0] += qw * mti_accel_f_idx_2;
      b_w_exp_tilde[1] += c_a[3 * i + 1] * mti_accel_f_idx_2;
      b_w_exp_tilde[2] += c_a[3 * i + 2] * mti_accel_f_idx_2;
    }
    ad_accel_f[0] = mti_accel_f_idx_0 * -0.0;
    ad_accel_f[1] =
        mti_accel_f_idx_0 * (-0.16182736457722724 * sin(b_atan2(x[9], x[7])));
    ad_accel_f[2] =
        mti_accel_f_idx_0 * (-0.16182736457722724 * -sin(b_atan2(x[8], x[7])));
    memset(&dv1[0], 0, 3U * sizeof(double));
    memset(&b_dt[0], 0, 3U * sizeof(double));
    mti_baro_f = dv1[0];
    t1_density = dv1[1];
    d5 = dv1[2];
    d6 = b_dt[0];
    d7 = b_dt[1];
    d8 = b_dt[2];
    mti_accel_f_idx_0 = x[7];
    board_gyro_f_idx_0 = x[8];
    expl_temp = x[9];
    for (i = 0; i < 3; i++) {
      mti_accel_f_idx_2 = param.Jinv[3 * i];
      qw = b_w_exp_tilde[i];
      mti_baro_f += mti_accel_f_idx_2 * qw;
      mti_accel_f_idx_1 = ad_accel_f[i];
      d6 += dt * mti_accel_f_idx_2 * mti_accel_f_idx_1;
      mti_accel_f_idx_2 = param.Jinv[3 * i + 1];
      t1_density += mti_accel_f_idx_2 * qw;
      d7 += dt * mti_accel_f_idx_2 * mti_accel_f_idx_1;
      mti_accel_f_idx_2 = param.Jinv[3 * i + 2];
      d5 += mti_accel_f_idx_2 * qw;
      d8 += dt * mti_accel_f_idx_2 * mti_accel_f_idx_1;
      mti_accel_f_idx_2 = board_accel_f[i];
      b_w_exp_tilde[i] = ((w_exp_tilde[i] * mti_accel_f_idx_0 +
                           w_exp_tilde[i + 3] * board_gyro_f_idx_0) +
                          w_exp_tilde[i + 6] * expl_temp) +
                         dt * ((board_accel_f_tmp_tmp / mti_accel_f_idx_2 *
                                    sens_input->board_accel.meas[i] +
                                b_board_accel_f_tmp_tmp / mti_accel_f_idx_2 *
                                    sens_input->mti_accel.meas[i]) +
                               c_board_accel_f_tmp_tmp / mti_accel_f_idx_2 *
                                   sens_input->ad_accel.meas[i]);
    }
    memset(&ad_accel_f[0], 0, 3U * sizeof(double));
    board_gyro_f_idx_1 = ad_accel_f[0];
    board_gyro_f_idx_2 = ad_accel_f[1];
    board_baro_f = ad_accel_f[2];
    mti_accel_f_idx_2 = x[7];
    qw = x[8];
    mti_accel_f_idx_1 = x[9];
    for (i = 0; i < 3; i++) {
      mti_accel_f_idx_0 = ST[3 * i];
      board_gyro_f_idx_0 = param.g[i];
      board_gyro_f_idx_1 += dt * mti_accel_f_idx_0 * board_gyro_f_idx_0;
      expl_temp = mti_accel_f_idx_0 * mti_accel_f_idx_2;
      mti_accel_f_idx_0 = ST[3 * i + 1];
      board_gyro_f_idx_2 += dt * mti_accel_f_idx_0 * board_gyro_f_idx_0;
      expl_temp += mti_accel_f_idx_0 * qw;
      mti_accel_f_idx_0 = ST[3 * i + 2];
      board_baro_f += dt * mti_accel_f_idx_0 * board_gyro_f_idx_0;
      expl_temp += mti_accel_f_idx_0 * mti_accel_f_idx_1;
      board_accel_f[i] = expl_temp;
    }
    memset(&c_q[0], 0, sizeof(double) << 2);
    mti_accel_f_idx_2 = c_q[0];
    qw = c_q[1];
    mti_accel_f_idx_1 = c_q[2];
    mti_accel_f_idx_0 = c_q[3];
    for (i = 0; i < 4; i++) {
      b_w_exp_tilde_tmp = i << 2;
      board_gyro_f_idx_0 = dv[i];
      mti_accel_f_idx_2 += b_q[b_w_exp_tilde_tmp] * board_gyro_f_idx_0;
      qw += b_q[b_w_exp_tilde_tmp + 1] * board_gyro_f_idx_0;
      mti_accel_f_idx_1 += b_q[b_w_exp_tilde_tmp + 2] * board_gyro_f_idx_0;
      mti_accel_f_idx_0 += b_q[b_w_exp_tilde_tmp + 3] * board_gyro_f_idx_0;
    }
    x_pred[0] = mti_accel_f_idx_2;
    x_pred[1] = qw;
    x_pred[2] = mti_accel_f_idx_1;
    x_pred[3] = mti_accel_f_idx_0;
    x_pred[4] = mti_baro_f + d6;
    x_pred[7] = b_w_exp_tilde[0] + board_gyro_f_idx_1;
    x_pred[5] = t1_density + d7;
    x_pred[8] = b_w_exp_tilde[1] + board_gyro_f_idx_2;
    x_pred[6] = d5 + d8;
    x_pred[9] = b_w_exp_tilde[2] + board_baro_f;
    x_pred[10] = x[10] + dt * board_accel_f[0];
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
    board_gyro_f_idx_1 = 0.5 * dt;
    qw = board_gyro_f_idx_1 * 0.0;
    W_dt[0] = qw;
    mti_accel_f_idx_1 = board_gyro_f_idx_1 * -x[4];
    W_dt[4] = mti_accel_f_idx_1;
    mti_accel_f_idx_2 = board_gyro_f_idx_1 * -x[5];
    W_dt[8] = mti_accel_f_idx_2;
    expl_temp = board_gyro_f_idx_1 * -x[6];
    W_dt[12] = expl_temp;
    mti_accel_f_idx_0 = board_gyro_f_idx_1 * x[4];
    W_dt[1] = mti_accel_f_idx_0;
    W_dt[5] = qw;
    board_gyro_f_idx_0 = board_gyro_f_idx_1 * x[6];
    W_dt[9] = board_gyro_f_idx_0;
    W_dt[13] = mti_accel_f_idx_2;
    mti_accel_f_idx_2 = board_gyro_f_idx_1 * x[5];
    W_dt[2] = mti_accel_f_idx_2;
    W_dt[6] = expl_temp;
    W_dt[10] = qw;
    W_dt[14] = mti_accel_f_idx_0;
    W_dt[3] = board_gyro_f_idx_0;
    W_dt[7] = mti_accel_f_idx_2;
    W_dt[11] = mti_accel_f_idx_1;
    W_dt[15] = qw;
    memset(&aBuffer[0], 0, sizeof(double) << 4);
    c_eye(dv2);
    memset(&b_W_dt[0], 0, sizeof(double) << 4);
    memset(&b_q[0], 0, sizeof(double) << 4);
    memset(&c_W_dt[0], 0, sizeof(double) << 4);
    for (i1 = 0; i1 < 4; i1++) {
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        qw = W_dt[i + b_w_exp_tilde_tmp];
        c_w_exp_tilde_tmp = i << 2;
        mti_accel_f_idx_2 = W_dt[c_w_exp_tilde_tmp] * qw;
        aBuffer[b_w_exp_tilde_tmp] += mti_accel_f_idx_2;
        b_W_dt[b_w_exp_tilde_tmp] += mti_accel_f_idx_2;
        b_q[b_w_exp_tilde_tmp] += mti_accel_f_idx_2;
        mti_accel_f_idx_2 = W_dt[c_w_exp_tilde_tmp + 1] * qw;
        aBuffer[b_w_exp_tilde_tmp + 1] += mti_accel_f_idx_2;
        b_W_dt[b_w_exp_tilde_tmp + 1] += mti_accel_f_idx_2;
        b_q[b_w_exp_tilde_tmp + 1] += mti_accel_f_idx_2;
        mti_accel_f_idx_2 = W_dt[c_w_exp_tilde_tmp + 2] * qw;
        aBuffer[b_w_exp_tilde_tmp + 2] += mti_accel_f_idx_2;
        b_W_dt[b_w_exp_tilde_tmp + 2] += mti_accel_f_idx_2;
        b_q[b_w_exp_tilde_tmp + 2] += mti_accel_f_idx_2;
        mti_accel_f_idx_2 = W_dt[c_w_exp_tilde_tmp + 3] * qw;
        aBuffer[b_w_exp_tilde_tmp + 3] += mti_accel_f_idx_2;
        b_W_dt[b_w_exp_tilde_tmp + 3] += mti_accel_f_idx_2;
        b_q[b_w_exp_tilde_tmp + 3] += mti_accel_f_idx_2;
      }
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        mti_accel_f_idx_2 = b_q[i + b_w_exp_tilde_tmp];
        c_w_exp_tilde_tmp = i << 2;
        c_W_dt[b_w_exp_tilde_tmp] +=
            W_dt[c_w_exp_tilde_tmp] * mti_accel_f_idx_2;
        c_W_dt[b_w_exp_tilde_tmp + 1] +=
            W_dt[c_w_exp_tilde_tmp + 1] * mti_accel_f_idx_2;
        c_W_dt[b_w_exp_tilde_tmp + 2] +=
            W_dt[c_w_exp_tilde_tmp + 2] * mti_accel_f_idx_2;
        c_W_dt[b_w_exp_tilde_tmp + 3] +=
            W_dt[c_w_exp_tilde_tmp + 3] * mti_accel_f_idx_2;
      }
    }
    memset(&b_q[0], 0, sizeof(double) << 4);
    for (i1 = 0; i1 < 4; i1++) {
      b_w_exp_tilde_tmp = i1 << 2;
      for (i = 0; i < 4; i++) {
        mti_accel_f_idx_2 = aBuffer[i + b_w_exp_tilde_tmp];
        c_w_exp_tilde_tmp = i << 2;
        b_q[b_w_exp_tilde_tmp] +=
            aBuffer[c_w_exp_tilde_tmp] * mti_accel_f_idx_2;
        b_q[b_w_exp_tilde_tmp + 1] +=
            aBuffer[c_w_exp_tilde_tmp + 1] * mti_accel_f_idx_2;
        b_q[b_w_exp_tilde_tmp + 2] +=
            aBuffer[c_w_exp_tilde_tmp + 2] * mti_accel_f_idx_2;
        b_q[b_w_exp_tilde_tmp + 3] +=
            aBuffer[c_w_exp_tilde_tmp + 3] * mti_accel_f_idx_2;
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
    qw = board_gyro_f_idx_1 * q[0];
    b_q[0] = qw;
    mti_accel_f_idx_1 = board_gyro_f_idx_1 * -q[1];
    b_q[4] = mti_accel_f_idx_1;
    mti_accel_f_idx_0 = board_gyro_f_idx_1 * -q[2];
    b_q[8] = mti_accel_f_idx_0;
    mti_accel_f_idx_2 = board_gyro_f_idx_1 * -q[3];
    b_q[12] = mti_accel_f_idx_2;
    expl_temp = board_gyro_f_idx_1 * q[1];
    b_q[1] = expl_temp;
    b_q[5] = qw;
    b_q[9] = mti_accel_f_idx_2;
    mti_accel_f_idx_2 = board_gyro_f_idx_1 * q[2];
    b_q[13] = mti_accel_f_idx_2;
    b_q[2] = mti_accel_f_idx_2;
    mti_accel_f_idx_2 = board_gyro_f_idx_1 * q[3];
    b_q[6] = mti_accel_f_idx_2;
    b_q[10] = qw;
    b_q[14] = mti_accel_f_idx_1;
    b_q[3] = mti_accel_f_idx_2;
    b_q[7] = mti_accel_f_idx_0;
    b_q[11] = expl_temp;
    b_q[15] = qw;
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
    airdata_atmos(x[10], &mti_accel_f_idx_2, &t1_density, &qw,
                  &mti_accel_f_idx_0, &mti_accel_f_idx_1);
    /*  computes matrix exponential of rotation */
    /* %% incremental angle */
    /* %% normed skew-symmetric matrix */
    if (dphi_tmp == 0.0) {
      board_gyro_f_idx_0 = 0.0;
      board_gyro_f_idx_1 = 0.0;
      board_gyro_f_idx_2 = 0.0;
    } else {
      board_gyro_f_idx_0 = x[4] / dphi_tmp;
      board_gyro_f_idx_1 = x[5] / dphi_tmp;
      board_gyro_f_idx_2 = x[6] / dphi_tmp;
    }
    /*  skew symmetric matrix / cross-product jacobian */
    ST[0] = 0.0;
    ST[3] = -board_gyro_f_idx_2;
    ST[6] = board_gyro_f_idx_1;
    ST[1] = board_gyro_f_idx_2;
    ST[4] = 0.0;
    ST[7] = -board_gyro_f_idx_0;
    ST[2] = -board_gyro_f_idx_1;
    ST[5] = board_gyro_f_idx_0;
    ST[8] = 0.0;
    /* %% Rodrigues rotation formula */
    memset(&w_exp_tilde[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      mti_accel_f_idx_2 = w_exp_tilde[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        qw = ST[i1 + 3 * i];
        mti_accel_f_idx_2 += ST[3 * i1] * qw;
        w_exp_tilde[b_w_exp_tilde_tmp] += ST[3 * i1 + 1] * qw;
        w_exp_tilde[c_w_exp_tilde_tmp] += ST[3 * i1 + 2] * qw;
      }
      w_exp_tilde[3 * i] = mti_accel_f_idx_2;
    }
    for (i = 0; i < 9; i++) {
      w_exp_tilde[i] =
          (w_exp_tilde_tmp[i] - b_a * ST[i]) + (1.0 - b_x) * w_exp_tilde[i];
    }
    memset(&dv3[0], 0, 9U * sizeof(double));
    for (i = 0; i < 3; i++) {
      qw = dv3[3 * i];
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_w_exp_tilde_tmp = 3 * i + 2;
      for (i1 = 0; i1 < 3; i1++) {
        mti_accel_f_idx_1 = w_exp_tilde[i1 + 3 * i];
        qw += b_param.Jinv[3 * i1] * mti_accel_f_idx_1;
        dv3[b_w_exp_tilde_tmp] += b_param.Jinv[3 * i1 + 1] * mti_accel_f_idx_1;
        dv3[c_w_exp_tilde_tmp] += b_param.Jinv[3 * i1 + 2] * mti_accel_f_idx_1;
        F[(i1 + 11 * (i + 4)) + 4] = 0.0;
      }
      dv3[3 * i] = qw;
    }
    c_a[1] = t1_density * (-0.08091368228861362 * x[9]);
    mti_accel_f_idx_2 = t1_density * -0.0;
    c_a[4] = mti_accel_f_idx_2;
    c_a[7] = t1_density * (-0.08091368228861362 * x[7]);
    c_a[2] = t1_density * (-0.08091368228861362 * -x[8]);
    c_a[5] = t1_density * (-0.08091368228861362 * -x[7]);
    c_a[8] = mti_accel_f_idx_2;
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
    mti_accel_f_idx_0 = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      c_a[3 * i1] = mti_accel_f_idx_2;
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
      mti_accel_f_idx_0 += x[i1 + 1] * b_param.g[i1];
    }
    /*  Sola */
    /* %% for hardcoding: */
    /*  [2*qw*vx - 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz, 2*qx*vy -
     * 2*qw*vz - 2*qy*vx, 2*qw*vy + 2*qx*vz - 2*qz*vx] */
    /*  [2*qw*vy + 2*qx*vz - 2*qz*vx, 2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qx*vx +
     * 2*qy*vy + 2*qz*vz, 2*qy*vz - 2*qw*vx - 2*qz*vy] */
    /*  [2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qz*vx - 2*qx*vz - 2*qw*vy, 2*qw*vx -
     * 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz] */
    ad_accel_f[0] = x[2] * b_param.g[2] - b_param.g[1] * x[3];
    ad_accel_f[1] = b_param.g[0] * x[3] - x[1] * b_param.g[2];
    ad_accel_f[2] = x[1] * b_param.g[1] - b_param.g[0] * x[2];
    mti_accel_f_idx_2 = x[0] * 0.0;
    c_a[0] = mti_accel_f_idx_2;
    c_a[3] = x[0] * -b_param.g[2];
    c_a[6] = x[0] * b_param.g[1];
    c_a[1] = x[0] * b_param.g[2];
    c_a[4] = mti_accel_f_idx_2;
    c_a[7] = x[0] * -b_param.g[0];
    c_a[2] = x[0] * -b_param.g[1];
    c_a[5] = x[0] * b_param.g[0];
    c_a[8] = mti_accel_f_idx_2;
    for (i = 0; i < 3; i++) {
      F[i + 7] = dt * (2.0 * (x[0] * b_param.g[i] - ad_accel_f[i]));
      mti_accel_f_idx_2 = x[i + 1];
      c_w_exp_tilde_tmp = 11 * (i + 1);
      F[c_w_exp_tilde_tmp + 7] =
          dt * (2.0 * (((mti_accel_f_idx_0 * b[3 * i] + x[1] * b_param.g[i]) -
                        b_param.g[0] * mti_accel_f_idx_2) +
                       c_a[3 * i]));
      b_w_exp_tilde_tmp = 3 * i + 1;
      F[c_w_exp_tilde_tmp + 8] =
          dt *
          (2.0 *
           (((mti_accel_f_idx_0 * b[b_w_exp_tilde_tmp] + x[2] * b_param.g[i]) -
             b_param.g[1] * mti_accel_f_idx_2) +
            c_a[b_w_exp_tilde_tmp]));
      b_w_exp_tilde_tmp = 3 * i + 2;
      F[c_w_exp_tilde_tmp + 9] =
          dt *
          (2.0 *
           (((mti_accel_f_idx_0 * b[b_w_exp_tilde_tmp] + x[3] * b_param.g[i]) -
             b_param.g[2] * mti_accel_f_idx_2) +
            c_a[b_w_exp_tilde_tmp]));
    }
    /*  jacobian of {exp(-tilde(w)*dt)*v} wrt w    */
    /*  skew symmetric matrix / cross-product jacobian */
    ST[0] = 0.0;
    ST[3] = -x[9];
    ST[6] = x[8];
    ST[1] = x[9];
    ST[4] = 0.0;
    ST[7] = -x[7];
    ST[2] = -x[8];
    ST[5] = x[7];
    ST[8] = 0.0;
    /*  skew symmetric matrix / cross-product jacobian */
    w_exp_tilde_tmp[0] = 0.0;
    w_exp_tilde_tmp[3] = -x[6];
    w_exp_tilde_tmp[6] = x[5];
    w_exp_tilde_tmp[1] = x[6];
    w_exp_tilde_tmp[4] = 0.0;
    w_exp_tilde_tmp[7] = -x[4];
    w_exp_tilde_tmp[2] = -x[5];
    w_exp_tilde_tmp[5] = x[4];
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
    q[0] = x[0];
    /*  computes rotation matrix from quaternion */
    /* %% norm quaternions */
    /*  q = quaternion_norm(q); */
    /* %% quaternion definition */
    /* %% rotation matrix */
    /*  S = (qw^2-qv'*qv)*eye(3) + 2*qv*qv' - 2*qw*q_tilde; % Stevens */
    /* %% Jacobian of rotation wrt quaternion */
    board_gyro_f_idx_0 = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      int d_w_exp_tilde_tmp;
      mti_accel_f_idx_2 = c_a[3 * i1];
      qw = dv3[3 * i1];
      c_w_exp_tilde_tmp = 3 * i1 + 1;
      d_w_exp_tilde_tmp = 3 * i1 + 2;
      for (i = 0; i < 3; i++) {
        b_w_exp_tilde_tmp = i + 3 * i1;
        mti_accel_f_idx_0 = w_exp_tilde_tmp[b_w_exp_tilde_tmp];
        mti_accel_f_idx_1 = ST[b_w_exp_tilde_tmp];
        mti_accel_f_idx_2 += ST[3 * i] * mti_accel_f_idx_0;
        qw += 2.0 * w_exp_tilde_tmp[3 * i] * mti_accel_f_idx_1;
        b_w_exp_tilde_tmp = 3 * i + 1;
        c_a[c_w_exp_tilde_tmp] += ST[b_w_exp_tilde_tmp] * mti_accel_f_idx_0;
        dv3[c_w_exp_tilde_tmp] +=
            2.0 * w_exp_tilde_tmp[b_w_exp_tilde_tmp] * mti_accel_f_idx_1;
        b_w_exp_tilde_tmp = 3 * i + 2;
        c_a[d_w_exp_tilde_tmp] += ST[b_w_exp_tilde_tmp] * mti_accel_f_idx_0;
        dv3[d_w_exp_tilde_tmp] +=
            2.0 * w_exp_tilde_tmp[b_w_exp_tilde_tmp] * mti_accel_f_idx_1;
      }
      dv3[3 * i1] = qw;
      c_a[3 * i1] = mti_accel_f_idx_2;
      c_w_exp_tilde_tmp = 11 * (i1 + 4);
      F[c_w_exp_tilde_tmp + 7] =
          dt * ST[3 * i1] + expl_temp * (mti_accel_f_idx_2 - qw);
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
      qw = -x[i1 + 1];
      q[i1 + 1] = qw;
      board_gyro_f_idx_0 += qw * x[i1 + 7];
    }
    /*  Sola */
    /* %% for hardcoding: */
    /*  [2*qw*vx - 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz, 2*qx*vy -
     * 2*qw*vz - 2*qy*vx, 2*qw*vy + 2*qx*vz - 2*qz*vx] */
    /*  [2*qw*vy + 2*qx*vz - 2*qz*vx, 2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qx*vx +
     * 2*qy*vy + 2*qz*vz, 2*qy*vz - 2*qw*vx - 2*qz*vy] */
    /*  [2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qz*vx - 2*qx*vz - 2*qw*vy, 2*qw*vx -
     * 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz] */
    ad_accel_f[0] = q[2] * x[9] - q[3] * x[8];
    ad_accel_f[1] = q[3] * x[7] - q[1] * x[9];
    ad_accel_f[2] = q[1] * x[8] - q[2] * x[7];
    for (i = 0; i < 3; i++) {
      qw = x[i + 7];
      c_dt[i] = dt * (2.0 * (q[0] * qw - ad_accel_f[i]));
      mti_accel_f_idx_1 = q[i + 1];
      c_w_exp_tilde_tmp = 3 * (i + 1);
      c_dt[c_w_exp_tilde_tmp] =
          dt * (2.0 * (((board_gyro_f_idx_0 * b[3 * i] + q[1] * qw) -
                        x[7] * mti_accel_f_idx_1) +
                       q[0] * ST[3 * i]));
      b_w_exp_tilde_tmp = 3 * i + 1;
      c_dt[c_w_exp_tilde_tmp + 1] =
          dt *
          (2.0 * (((board_gyro_f_idx_0 * b[b_w_exp_tilde_tmp] + q[2] * qw) -
                   x[8] * mti_accel_f_idx_1) +
                  q[0] * ST[b_w_exp_tilde_tmp]));
      b_w_exp_tilde_tmp = 3 * i + 2;
      c_dt[c_w_exp_tilde_tmp + 2] =
          dt *
          (2.0 * (((board_gyro_f_idx_0 * b[b_w_exp_tilde_tmp] + q[3] * qw) -
                   x[9] * mti_accel_f_idx_1) +
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
    mti_accel_f_idx_1 = 2.0 * q[0];
    /*  Stevens */
    /* %% for hardcoding:  */
    /*  [qw^2 + qx^2 - qy^2 - qz^2,         2*qw*qz + 2*qx*qy,         2*qx*qz -
     * 2*qw*qy] */
    /*  [        2*qx*qy - 2*qw*qz, qw^2 - qx^2 + qy^2 - qz^2,         2*qw*qx +
     * 2*qy*qz] */
    /*  [        2*qw*qy + 2*qx*qz,         2*qy*qz - 2*qw*qx, qw^2 - qx^2 -
     * qy^2 + qz^2] */
    mti_accel_f_idx_0 = mti_accel_f_idx_1 * 0.0;
    c_a[0] = mti_accel_f_idx_0;
    c_a[3] = mti_accel_f_idx_1 * -q[3];
    c_a[6] = mti_accel_f_idx_1 * q[2];
    c_a[1] = mti_accel_f_idx_1 * q[3];
    c_a[4] = mti_accel_f_idx_0;
    c_a[7] = mti_accel_f_idx_1 * -q[1];
    c_a[2] = mti_accel_f_idx_1 * -q[2];
    c_a[5] = mti_accel_f_idx_1 * q[1];
    c_a[8] = mti_accel_f_idx_0;
    for (i = 0; i < 3; i++) {
      F[11 * (i + 7) + 10] =
          dt * ((qw * b[3 * i] + 2.0 * q[1] * q[i + 1]) - c_a[3 * i]);
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
        mti_accel_f_idx_2 = P[i1 + 11 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_F[b_w_exp_tilde_tmp] += F[i2 + 11 * i1] * mti_accel_f_idx_2;
        }
      }
    }
    for (i1 = 0; i1 < 11; i1++) {
      for (i2 = 0; i2 < 11; i2++) {
        mti_accel_f_idx_2 = 0.0;
        for (i = 0; i < 11; i++) {
          mti_accel_f_idx_2 += b_F[i1 + 11 * i] * F[i2 + 11 * i];
        }
        b_w_exp_tilde_tmp = i1 + 11 * i2;
        P_pred[b_w_exp_tilde_tmp] = mti_accel_f_idx_2 + Q[b_w_exp_tilde_tmp];
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
        mti_accel_f_idx_2 = ST[i1 + 3 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          K[b_w_exp_tilde_tmp] +=
              P_pred[i2 + 11 * (i1 + 4)] * mti_accel_f_idx_2;
        }
      }
    }
    /*  hardcoded P*H */
    d_eye(F);
    memcpy(&E[0], &F[0], 44U * sizeof(double));
    for (i = 0; i < 33; i++) {
      E[i + 44] = F[i + 44] - K[i];
    }
    memcpy(&E[77], &F[77], 44U * sizeof(double));
    /*  hardcoded K*H */
    /* %% correct covariance estimate */
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      for (i1 = 0; i1 < 11; i1++) {
        mti_accel_f_idx_2 = P_pred[i1 + 11 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_F[b_w_exp_tilde_tmp] += E[i2 + 11 * i1] * mti_accel_f_idx_2;
        }
      }
    }
    memset(&b_K[0], 0, 33U * sizeof(double));
    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        mti_accel_f_idx_2 = R[i1 + 3 * i];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_K[b_w_exp_tilde_tmp] += K[i2 + 11 * i1] * mti_accel_f_idx_2;
        }
      }
    }
    memset(&P_ret[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      for (i1 = 0; i1 < 11; i1++) {
        mti_accel_f_idx_2 = E[i + 11 * i1];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          P_ret[b_w_exp_tilde_tmp] += b_F[i2 + 11 * i1] * mti_accel_f_idx_2;
        }
      }
    }
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i = 0; i < 11; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        mti_accel_f_idx_2 = K[i + 11 * i1];
        for (i2 = 0; i2 < 11; i2++) {
          b_w_exp_tilde_tmp = i2 + 11 * i;
          b_F[b_w_exp_tilde_tmp] += b_K[i2 + 11 * i1] * mti_accel_f_idx_2;
        }
      }
    }
    for (i = 0; i < 121; i++) {
      P_ret[i] += b_F[i];
    }
    /*  joseph stabilized */
    /* %% correct state estimate */
    mti_accel_f_idx_0 =
        ((d1 / d4 * (sens_input->board_gyro.meas[0] - bias->board_gyro[0]) +
          d2 / d4 * (sens_input->mti_gyro.meas[0] - bias->mti_gyro[0])) +
         a[0] * (sens_input->ad_gyro.meas[0] - bias->ad_gyro[0])) -
        x_pred[4];
    qw = d1 / d3;
    mti_accel_f_idx_2 = d2 / d3;
    mti_accel_f_idx_1 =
        ((qw * (sens_input->board_gyro.meas[1] - bias->board_gyro[1]) +
          mti_accel_f_idx_2 *
              (sens_input->mti_gyro.meas[1] - bias->mti_gyro[1])) +
         d * (sens_input->ad_gyro.meas[1] - bias->ad_gyro[1])) -
        x_pred[5];
    mti_accel_f_idx_2 =
        ((qw * (sens_input->board_gyro.meas[2] - bias->board_gyro[2]) +
          mti_accel_f_idx_2 *
              (sens_input->mti_gyro.meas[2] - bias->mti_gyro[2])) +
         d * (sens_input->ad_gyro.meas[2] - bias->ad_gyro[2])) -
        x_pred[6];
    for (i = 0; i < 11; i++) {
      x_ret[i] = x_pred[i] +
                 ((K[i] * mti_accel_f_idx_0 + K[i + 11] * mti_accel_f_idx_1) +
                  K[i + 22] * mti_accel_f_idx_2);
    }
    /*  norms quaternion */
    /* %% inverse quaternion  */
    mti_accel_f_idx_2 = c_norm(&x_ret[0]);
    x_ret[0] /= mti_accel_f_idx_2;
    x_ret[1] /= mti_accel_f_idx_2;
    x_ret[2] /= mti_accel_f_idx_2;
    x_ret[3] /= mti_accel_f_idx_2;
    /*  norm quaternions */
    /*     %% Correction steps, sequential for each additional sensor */
    /* %% Correction with barometer, magnetometer */
    /* %% R is a square matrix (size length of sensor vector), tuning for
     * expected measurement noise magnitude E(noise) */
    /* %% Barometer */
    if (sens_input->board_baro.status) {
      /*  only correct with alive IMUs */
      /* %% y = [ P(1) ] */
      memcpy(&x_pred[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&F[0], &P_ret[0], 121U * sizeof(double));
      b_ekf_correct(x_pred, F, sens_input->board_baro.meas, bias->board_baro,
                    x_ret, P_ret);
    }
    if (sens_input->mti_baro.status) {
      /*  only correct with alive IMUs */
      /* %% y = [ P(1) ] */
      memcpy(&x_pred[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&F[0], &P_ret[0], 121U * sizeof(double));
      b_ekf_correct(x_pred, F, sens_input->mti_baro.meas, bias->mti_baro, x_ret,
                    P_ret);
    }
    /* %% Magnetometer */
    if (sens_input->board_mag.status) {
      /* %% y = [  Mag(3) ] */
      memcpy(&x_pred[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&F[0], &P_ret[0], 121U * sizeof(double));
      ekf_correct(x_pred, F, sens_input->board_mag.meas, bias->board_mag_earth,
                  b, x_ret, P_ret);
    }
    if (sens_input->mti_mag.status) {
      /* %% y = [  Mag(3) ] */
      memcpy(&x_pred[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&F[0], &P_ret[0], 121U * sizeof(double));
      ekf_correct(x_pred, F, sens_input->mti_mag.meas, bias->mti_mag_earth, b,
                  x_ret, P_ret);
    }
  }
}

/*
 * File trailer for navigation_codegen_entry.c
 *
 * [EOF]
 */
