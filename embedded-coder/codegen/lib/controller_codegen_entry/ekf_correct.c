/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ekf_correct.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

/* Include Files */
#include "ekf_correct.h"
#include "airdata_atmos.h"
#include "eye.h"
#include "inv.h"
#include "mpower.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * Computes EKF correction step for other sensor data.
 *  Input function handles: measurement model, model jacobian;
 *  Input variables: old state x, old covariance P, measurement y, sensor bias
 * b; Input parameters: sensor noise weighting R; Outputs: new state x, new
 * covariance P
 *
 * Arguments    : const double x[11]
 *                const double b_P[121]
 *                double y
 *                double b
 *                double x_new[11]
 *                double P_new[121]
 * Return Type  : void
 */
void b_ekf_correct(const double x[11], const double b_P[121], double y,
                   double b, double x_new[11], double P_new[121])
{
  __m128d b_r;
  __m128d r1;
  double E[121];
  double b_E[121];
  double H[11];
  double b_H[11];
  double altitude_ratio;
  double b_b;
  double expl_temp;
  double layer_idx_1;
  double layer_idx_2;
  double layer_idx_3;
  int i;
  int i1;
  int i2;
  int i3;
  int layer_idx_0;
  /*     %% Correction */
  /*  computes a-posteriori state and covariance estimates */
  /*  Uses discrete-time measurement model and analytical Jacobian */
  /* %% compute expected measurement and difference to measured values */
  /*  Computes measurement prediction using current state and sensor biases */
  /*     %% decomp */
  /*  decompose state vector: [q(4); w(3); v(3); alt] */
  /*     %% atmosphere model */
  /*  [Pa], measured atmospheric pressure */
  /*     %% measurement prediction */
  layer_idx_2 = airdata_atmos(x[10], &layer_idx_2, &layer_idx_3,
                              &altitude_ratio, &expl_temp, &layer_idx_1);
  b_b = y - (layer_idx_2 + b);
  /* %% compute Jacobian: H = dh/dx */
  /*  H = jacobian(@model_measurement, 0, x, b);  */
  /*  Computes measurement prediction using current state and sensor biases */
  /*     %% decomp */
  /*  decompose state vector: [q(4); w(3); v(3); alt;] */
  /*     %% Initialize */
  /*  Jacobian is of size: length(measurement) rows, length(x) columns */
  /*  fill with zeros at first */
  memset(&H[0], 0, 11U * sizeof(double));
  /*     %% atmosphere model */
  /*  [P, ~, ~, ~] = model_airdata(alt); */
  /*  computes atmospheric data and gradient w.r.t. altitude from altitude,
   * according to US standard atmosphere  */
  /*  atmospheric data: static pressure, temperature, density, local speed of
   * sound */
  /*  calculations found in Stengel 2004, pp. 30 */
  /*  parameters */
  /*  adiabatic index */
  /*  specific gas constant for air */
  /*  troposphere */
  /*  tropopause */
  /*  stratosphere */
  /*  stratosphere 2 */
  /*  base height, P_base, T_base, lapse rate; */
  /*  mean earth radius */
  /*  zero height gravity */
  /*  geopotential altitude */
  altitude_ratio = 6.356766E+6 / (6.356766E+6 - x[10]);
  expl_temp = altitude_ratio * x[10];
  /*  select atmosphere behaviour from table */
  layer_idx_0 = 0;
  layer_idx_1 = 101325.0;
  layer_idx_2 = 288.15;
  layer_idx_3 = 0.0065;
  if (expl_temp > 11000.0) {
    if (expl_temp < 20000.0) {
      layer_idx_0 = 11000;
      layer_idx_1 = 22632.1;
      layer_idx_2 = 216.65;
      layer_idx_3 = 0.0;
    } else if (expl_temp < 32000.0) {
      layer_idx_0 = 20000;
      layer_idx_1 = 5474.9;
      layer_idx_2 = 216.65;
      layer_idx_3 = -0.001;
    } else {
      layer_idx_0 = 32000;
      layer_idx_1 = 868.02;
      layer_idx_2 = 228.65;
      layer_idx_3 = -0.0028;
    }
  }
  /*  base height */
  /*  base pressure */
  /*  base temperature */
  /*  temperature lapse rate */
  if (layer_idx_3 == 0.0) {
    layer_idx_2 *= 287.0579;
    layer_idx_2 = -layer_idx_1 * 9.81 / layer_idx_2 *
                  (altitude_ratio * altitude_ratio) *
                  exp(-9.81 * (expl_temp - (double)layer_idx_0) / layer_idx_2);
  } else {
    layer_idx_2 = -layer_idx_1 * 9.81 / (layer_idx_2 * 287.0579) *
                  (altitude_ratio * altitude_ratio) *
                  mpower(1.0 - layer_idx_3 / layer_idx_2 *
                                   (expl_temp - (double)layer_idx_0),
                         9.81 / (287.0579 * layer_idx_3) - 1.0);
  }
  /*  return values */
  /*  need to initialize struct field before calling airdata_dynamic */
  H[10] = layer_idx_2;
  /*     %% measurement prediction */
  /*  y = [W; M; P]; */
  /* %% compute Kalman gain (and helper matrices) */
  memset(&b_H[0], 0, 11U * sizeof(double));
  layer_idx_2 = 0.0;
  memset(&x_new[0], 0, 11U * sizeof(double));
  for (i = 0; i < 11; i++) {
    layer_idx_3 = b_H[i];
    for (i1 = 0; i1 < 11; i1++) {
      layer_idx_3 += H[i1] * b_P[i1 + 11 * i];
    }
    b_H[i] = layer_idx_3;
    layer_idx_2 += layer_idx_3 * H[i];
    for (i1 = 0; i1 <= 8; i1 += 2) {
      b_r = _mm_loadu_pd(&x_new[i1]);
      _mm_storeu_pd(&x_new[i1],
                    _mm_add_pd(b_r, _mm_mul_pd(_mm_loadu_pd(&b_P[i1 + 11 * i]),
                                               _mm_set1_pd(H[i]))));
    }
    x_new[10] += b_P[11 * i + 10] * H[i];
  }
  layer_idx_2 = 1.0 / (layer_idx_2 + 100.0);
  for (i1 = 0; i1 <= 8; i1 += 2) {
    b_r = _mm_loadu_pd(&x_new[i1]);
    _mm_storeu_pd(&x_new[i1], _mm_mul_pd(b_r, _mm_set1_pd(layer_idx_2)));
  }
  x_new[10] *= layer_idx_2;
  d_eye(E);
  for (i = 0; i < 11; i++) {
    for (i1 = 0; i1 <= 8; i1 += 2) {
      b_r = _mm_loadu_pd(&x_new[i1]);
      layer_idx_0 = i1 + 11 * i;
      r1 = _mm_loadu_pd(&E[layer_idx_0]);
      _mm_storeu_pd(
          &b_E[layer_idx_0],
          _mm_sub_pd(r1, _mm_add_pd(_mm_set1_pd(0.0),
                                    _mm_mul_pd(b_r, _mm_set1_pd(H[i])))));
    }
    layer_idx_0 = 11 * i + 10;
    b_E[layer_idx_0] = E[layer_idx_0] - x_new[10] * H[i];
  }
  /* %% correct covariance estimate */
  memset(&E[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    layer_idx_0 = 11 * i1 + 10;
    for (i = 0; i < 11; i++) {
      layer_idx_2 = b_P[i + 11 * i1];
      for (i2 = 0; i2 <= 8; i2 += 2) {
        b_r = _mm_loadu_pd(&b_E[i2 + 11 * i]);
        i3 = i2 + 11 * i1;
        r1 = _mm_loadu_pd(&E[i3]);
        _mm_storeu_pd(
            &E[i3], _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(layer_idx_2))));
      }
      E[layer_idx_0] += b_E[11 * i + 10] * layer_idx_2;
    }
  }
  memset(&P_new[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    layer_idx_0 = 11 * i1 + 10;
    for (i = 0; i < 11; i++) {
      layer_idx_2 = b_E[i1 + 11 * i];
      for (i2 = 0; i2 <= 8; i2 += 2) {
        b_r = _mm_loadu_pd(&E[i2 + 11 * i]);
        i3 = i2 + 11 * i1;
        r1 = _mm_loadu_pd(&P_new[i3]);
        _mm_storeu_pd(
            &P_new[i3],
            _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(layer_idx_2))));
      }
      P_new[layer_idx_0] += E[11 * i + 10] * layer_idx_2;
    }
  }
  for (i1 = 0; i1 < 11; i1++) {
    for (i = 0; i <= 8; i += 2) {
      b_r = _mm_loadu_pd(&x_new[i]);
      _mm_storeu_pd(&E[i + 11 * i1],
                    _mm_mul_pd(_mm_mul_pd(b_r, _mm_set1_pd(100.0)),
                               _mm_set1_pd(x_new[i1])));
    }
    E[11 * i1 + 10] = x_new[10] * 100.0 * x_new[i1];
  }
  for (i1 = 0; i1 <= 118; i1 += 2) {
    b_r = _mm_loadu_pd(&P_new[i1]);
    r1 = _mm_loadu_pd(&E[i1]);
    _mm_storeu_pd(&P_new[i1], _mm_add_pd(b_r, r1));
  }
  P_new[120] += E[120];
  /*  joseph stabilized */
  /* %% correct state estimate */
  for (i1 = 0; i1 <= 8; i1 += 2) {
    b_r = _mm_loadu_pd(&x_new[i1]);
    _mm_storeu_pd(&x_new[i1], _mm_add_pd(_mm_loadu_pd(&x[i1]),
                                         _mm_mul_pd(b_r, _mm_set1_pd(b_b))));
  }
  x_new[10] = x[10] + x_new[10] * b_b;
  /*  norms quaternion */
  /* %% inverse quaternion  */
  b_r = _mm_loadu_pd(&x_new[0]);
  r1 = _mm_set1_pd(c_norm(&x_new[0]));
  _mm_storeu_pd(&x_new[0], _mm_div_pd(b_r, r1));
  b_r = _mm_loadu_pd(&x_new[2]);
  _mm_storeu_pd(&x_new[2], _mm_div_pd(b_r, r1));
  /*  norm quaternions */
}

/*
 * Computes EKF correction step for other sensor data.
 *  Input function handles: measurement model, model jacobian;
 *  Input variables: old state x, old covariance P, measurement y, sensor bias
 * b; Input parameters: sensor noise weighting R; Outputs: new state x, new
 * covariance P
 *
 * Arguments    : const double x[11]
 *                const double b_P[121]
 *                const double y[3]
 *                const double b[3]
 *                const double R[9]
 *                double x_new[11]
 *                double P_new[121]
 * Return Type  : void
 */
void ekf_correct(const double x[11], const double b_P[121], const double y[3],
                 const double b[3], const double R[9], double x_new[11],
                 double P_new[121])
{
  static const signed char c_b[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  __m128d b_r;
  __m128d r1;
  double E[121];
  double b_E[121];
  double H[33];
  double K[33];
  double c_P[33];
  double y_tmp[33];
  double b_b[9];
  double c_a[9];
  double c_x[3];
  double H_tmp;
  double a;
  double b_a;
  double b_x;
  double x_tmp;
  int b_H_tmp;
  int c_H_tmp;
  int i;
  int i1;
  int i2;
  /*     %% Correction */
  /*  computes a-posteriori state and covariance estimates */
  /*  Uses discrete-time measurement model and analytical Jacobian */
  /* %% compute expected measurement and difference to measured values */
  /*  Computes measurement prediction using current state and sensor biases */
  /*     %% decomp */
  /*  decompose state vector: [q(4); w(3); v(3); alt; Cl; delta] */
  /*     %% magnetic field model */
  /*  computes rotation matrix from quaternion */
  /* %% norm quaternions */
  /*  q = quaternion_norm(q); */
  /* %% quaternion definition */
  /* %% skew symetric quaternion matrix */
  /* %% rotation matrix */
  a = x[0] * x[0] - ((x[1] * x[1] + x[2] * x[2]) + x[3] * x[3]);
  b_a = 2.0 * x[0];
  /*  Stevens */
  /* %% for hardcoding:  */
  /*  [qw^2 + qx^2 - qy^2 - qz^2,         2*qw*qz + 2*qx*qy,         2*qx*qz -
   * 2*qw*qy] */
  /*  [        2*qx*qy - 2*qw*qz, qw^2 - qx^2 + qy^2 - qz^2,         2*qw*qx +
   * 2*qy*qz] */
  /*  [        2*qw*qy + 2*qx*qz,         2*qy*qz - 2*qw*qx, qw^2 - qx^2 - qy^2
   * + qz^2] */
  /*  [Gauss], measured Earth magnetic field in body frame */
  /*  TODO: add iron corrections */
  /*     %% measurement prediction */
  /* %% compute Jacobian: H = dh/dx */
  /*  H = jacobian(@model_measurement, 0, x, b);  */
  /*  Computes measurement prediction using current state and sensor biases */
  /*     %% decomp */
  /*  decompose state vector: [q(4); w(3); v(3); alt] */
  /*     %% Initialize */
  /*  Jacobian is of size: length(measurement) rows, length(x) columns */
  /*  fill with zeros at first */
  memset(&H[0], 0, 33U * sizeof(double));
  /*     %% magnetic field model */
  /*  S = quaternion_rotmatrix(q); */
  /*  M = S * M_E; % Earth magnetic field in body frame */
  /*  TODO: add iron corrections */
  /*  computes rotation matrix from quaternion */
  /* %% norm quaternions */
  /*  q = quaternion_norm(q); */
  /* %% quaternion definition */
  /* %% rotation matrix */
  /*  S = (qw^2-qv'*qv)*eye(3) + 2*qv*qv' - 2*qw*q_tilde; % Stevens */
  /* %% Jacobian of rotation wrt quaternion */
  b_x = (b[0] * x[1] + b[1] * x[2]) + b[2] * x[3];
  c_x[0] = x[2] * b[2] - b[1] * x[3];
  c_x[1] = b[0] * x[3] - x[1] * b[2];
  c_x[2] = x[1] * b[1] - b[0] * x[2];
  x_tmp = x[0] * 0.0;
  c_a[0] = x_tmp;
  c_a[3] = x[0] * -b[2];
  c_a[6] = x[0] * b[1];
  c_a[1] = x[0] * b[2];
  c_a[4] = x_tmp;
  c_a[7] = x[0] * -b[0];
  c_a[2] = x[0] * -b[1];
  c_a[5] = x[0] * b[0];
  c_a[8] = x_tmp;
  for (i = 0; i < 3; i++) {
    H[i] = 2.0 * (x[0] * b[i] - c_x[i]);
    H_tmp = x[i + 1];
    b_H_tmp = 3 * (i + 1);
    H[b_H_tmp] =
        2.0 * (((b_x * (double)c_b[3 * i] + x[1] * b[i]) - b[0] * H_tmp) +
               c_a[3 * i]);
    c_H_tmp = 3 * i + 1;
    H[b_H_tmp + 1] =
        2.0 * (((b_x * (double)c_b[c_H_tmp] + x[2] * b[i]) - b[1] * H_tmp) +
               c_a[c_H_tmp]);
    c_H_tmp = 3 * i + 2;
    H[b_H_tmp + 2] =
        2.0 * (((b_x * (double)c_b[c_H_tmp] + x[3] * b[i]) - b[2] * H_tmp) +
               c_a[c_H_tmp]);
  }
  /*  Sola */
  /* %% for hardcoding: */
  /*  [2*qw*vx - 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz, 2*qx*vy -
   * 2*qw*vz - 2*qy*vx, 2*qw*vy + 2*qx*vz - 2*qz*vx] */
  /*  [2*qw*vy + 2*qx*vz - 2*qz*vx, 2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qx*vx +
   * 2*qy*vy + 2*qz*vz, 2*qy*vz - 2*qw*vx - 2*qz*vy] */
  /*  [2*qw*vz - 2*qx*vy + 2*qy*vx, 2*qz*vx - 2*qx*vz - 2*qw*vy, 2*qw*vx -
   * 2*qy*vz + 2*qz*vy, 2*qx*vx + 2*qy*vy + 2*qz*vz] */
  /*     %% measurement prediction */
  /* %% compute Kalman gain (and helper matrices) */
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 11; i1++) {
      y_tmp[i1 + 11 * i] = H[i + 3 * i1];
    }
  }
  memset(&K[0], 0, 33U * sizeof(double));
  for (i = 0; i < 11; i++) {
    c_H_tmp = 3 * i + 2;
    for (i1 = 0; i1 < 11; i1++) {
      x_tmp = b_P[i1 + 11 * i];
      b_r = _mm_loadu_pd(&H[3 * i1]);
      r1 = _mm_loadu_pd(&K[3 * i]);
      _mm_storeu_pd(&K[3 * i],
                    _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(x_tmp))));
      K[c_H_tmp] += H[3 * i1 + 2] * x_tmp;
    }
  }
  memset(&c_P[0], 0, 33U * sizeof(double));
  for (i2 = 0; i2 < 3; i2++) {
    for (i1 = 0; i1 < 3; i1++) {
      x_tmp = 0.0;
      for (i = 0; i < 11; i++) {
        x_tmp += K[i2 + 3 * i] * y_tmp[i + 11 * i1];
      }
      c_H_tmp = i2 + 3 * i1;
      b_b[c_H_tmp] = x_tmp + R[c_H_tmp];
    }
    c_H_tmp = 11 * i2 + 10;
    for (i1 = 0; i1 < 11; i1++) {
      x_tmp = y_tmp[i1 + 11 * i2];
      for (i = 0; i <= 8; i += 2) {
        b_H_tmp = i + 11 * i2;
        b_r = _mm_loadu_pd(&c_P[b_H_tmp]);
        _mm_storeu_pd(
            &c_P[b_H_tmp],
            _mm_add_pd(b_r, _mm_mul_pd(_mm_loadu_pd(&b_P[i + 11 * i1]),
                                       _mm_set1_pd(x_tmp))));
      }
      c_P[c_H_tmp] += b_P[11 * i1 + 10] * x_tmp;
    }
  }
  inv(b_b, c_a);
  memset(&K[0], 0, 33U * sizeof(double));
  for (i1 = 0; i1 < 3; i1++) {
    c_H_tmp = 11 * i1 + 10;
    for (i = 0; i < 3; i++) {
      x_tmp = c_a[i + 3 * i1];
      for (i2 = 0; i2 <= 8; i2 += 2) {
        b_r = _mm_loadu_pd(&c_P[i2 + 11 * i]);
        b_H_tmp = i2 + 11 * i1;
        r1 = _mm_loadu_pd(&K[b_H_tmp]);
        _mm_storeu_pd(&K[b_H_tmp],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(x_tmp))));
      }
      K[c_H_tmp] += c_P[11 * i + 10] * x_tmp;
    }
  }
  d_eye(E);
  for (i = 0; i < 11; i++) {
    x_tmp = K[i];
    H_tmp = K[i + 11];
    b_x = K[i + 22];
    for (i1 = 0; i1 < 11; i1++) {
      c_H_tmp = i + 11 * i1;
      b_E[c_H_tmp] = E[c_H_tmp] - ((x_tmp * H[3 * i1] + H_tmp * H[3 * i1 + 1]) +
                                   b_x * H[3 * i1 + 2]);
    }
  }
  /* %% correct covariance estimate */
  memset(&E[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    c_H_tmp = 11 * i1 + 10;
    for (i = 0; i < 11; i++) {
      x_tmp = b_P[i + 11 * i1];
      for (i2 = 0; i2 <= 8; i2 += 2) {
        b_r = _mm_loadu_pd(&b_E[i2 + 11 * i]);
        b_H_tmp = i2 + 11 * i1;
        r1 = _mm_loadu_pd(&E[b_H_tmp]);
        _mm_storeu_pd(&E[b_H_tmp],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(x_tmp))));
      }
      E[c_H_tmp] += b_E[11 * i + 10] * x_tmp;
    }
  }
  memset(&y_tmp[0], 0, 33U * sizeof(double));
  for (i1 = 0; i1 < 3; i1++) {
    c_H_tmp = 11 * i1 + 10;
    for (i = 0; i < 3; i++) {
      x_tmp = R[i + 3 * i1];
      for (i2 = 0; i2 <= 8; i2 += 2) {
        b_r = _mm_loadu_pd(&K[i2 + 11 * i]);
        b_H_tmp = i2 + 11 * i1;
        r1 = _mm_loadu_pd(&y_tmp[b_H_tmp]);
        _mm_storeu_pd(&y_tmp[b_H_tmp],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(x_tmp))));
      }
      y_tmp[c_H_tmp] += K[11 * i + 10] * x_tmp;
    }
  }
  memset(&P_new[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    c_H_tmp = 11 * i1 + 10;
    for (i = 0; i < 11; i++) {
      x_tmp = b_E[i1 + 11 * i];
      for (i2 = 0; i2 <= 8; i2 += 2) {
        b_r = _mm_loadu_pd(&E[i2 + 11 * i]);
        b_H_tmp = i2 + 11 * i1;
        r1 = _mm_loadu_pd(&P_new[b_H_tmp]);
        _mm_storeu_pd(&P_new[b_H_tmp],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(x_tmp))));
      }
      P_new[c_H_tmp] += E[11 * i + 10] * x_tmp;
    }
  }
  memset(&E[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    c_H_tmp = 11 * i1 + 10;
    for (i = 0; i < 3; i++) {
      x_tmp = K[i1 + 11 * i];
      for (i2 = 0; i2 <= 8; i2 += 2) {
        b_r = _mm_loadu_pd(&y_tmp[i2 + 11 * i]);
        b_H_tmp = i2 + 11 * i1;
        r1 = _mm_loadu_pd(&E[b_H_tmp]);
        _mm_storeu_pd(&E[b_H_tmp],
                      _mm_add_pd(r1, _mm_mul_pd(b_r, _mm_set1_pd(x_tmp))));
      }
      E[c_H_tmp] += y_tmp[11 * i + 10] * x_tmp;
    }
  }
  for (i = 0; i <= 118; i += 2) {
    b_r = _mm_loadu_pd(&P_new[i]);
    r1 = _mm_loadu_pd(&E[i]);
    _mm_storeu_pd(&P_new[i], _mm_add_pd(b_r, r1));
  }
  P_new[120] += E[120];
  /*  joseph stabilized */
  /* %% correct state estimate */
  for (i = 0; i < 3; i++) {
    x_tmp = x[i + 1];
    c_a[3 * i] = a * (double)c_b[3 * i] + 2.0 * x[1] * x_tmp;
    c_H_tmp = 3 * i + 1;
    c_a[c_H_tmp] = a * (double)c_b[c_H_tmp] + 2.0 * x[2] * x_tmp;
    c_H_tmp = 3 * i + 2;
    c_a[c_H_tmp] = a * (double)c_b[c_H_tmp] + 2.0 * x[3] * x_tmp;
  }
  x_tmp = b_a * 0.0;
  b_b[0] = x_tmp;
  b_b[3] = b_a * -x[3];
  b_b[6] = b_a * x[2];
  b_b[1] = b_a * x[3];
  b_b[4] = x_tmp;
  b_b[7] = b_a * -x[1];
  b_b[2] = b_a * -x[2];
  b_b[5] = b_a * x[1];
  b_b[8] = x_tmp;
  b_r = _mm_loadu_pd(&c_a[0]);
  r1 = _mm_loadu_pd(&b_b[0]);
  _mm_storeu_pd(&c_a[0], _mm_sub_pd(b_r, r1));
  b_r = _mm_loadu_pd(&c_a[2]);
  r1 = _mm_loadu_pd(&b_b[2]);
  _mm_storeu_pd(&c_a[2], _mm_sub_pd(b_r, r1));
  b_r = _mm_loadu_pd(&c_a[4]);
  r1 = _mm_loadu_pd(&b_b[4]);
  _mm_storeu_pd(&c_a[4], _mm_sub_pd(b_r, r1));
  b_r = _mm_loadu_pd(&c_a[6]);
  r1 = _mm_loadu_pd(&b_b[6]);
  _mm_storeu_pd(&c_a[6], _mm_sub_pd(b_r, r1));
  c_a[8] -= x_tmp;
  x_tmp = b[0];
  H_tmp = b[1];
  b_x = b[2];
  b_r = _mm_loadu_pd(&c_a[0]);
  r1 = _mm_mul_pd(b_r, _mm_set1_pd(x_tmp));
  b_r = _mm_loadu_pd(&c_a[3]);
  b_r = _mm_mul_pd(b_r, _mm_set1_pd(H_tmp));
  r1 = _mm_add_pd(r1, b_r);
  b_r = _mm_loadu_pd(&c_a[6]);
  b_r = _mm_mul_pd(b_r, _mm_set1_pd(b_x));
  b_r = _mm_add_pd(r1, b_r);
  r1 = _mm_loadu_pd(&y[0]);
  b_r = _mm_sub_pd(r1, b_r);
  _mm_storeu_pd(&c_x[0], b_r);
  c_x[2] = y[2] - ((c_a[2] * x_tmp + c_a[5] * H_tmp) + c_a[8] * b_x);
  x_tmp = c_x[0];
  H_tmp = c_x[1];
  b_x = c_x[2];
  for (i = 0; i <= 8; i += 2) {
    b_r = _mm_loadu_pd(&K[i]);
    r1 = _mm_mul_pd(b_r, _mm_set1_pd(x_tmp));
    b_r = _mm_loadu_pd(&K[i + 11]);
    b_r = _mm_mul_pd(b_r, _mm_set1_pd(H_tmp));
    r1 = _mm_add_pd(r1, b_r);
    b_r = _mm_loadu_pd(&K[i + 22]);
    b_r = _mm_mul_pd(b_r, _mm_set1_pd(b_x));
    b_r = _mm_add_pd(r1, b_r);
    r1 = _mm_loadu_pd(&x[i]);
    b_r = _mm_add_pd(r1, b_r);
    _mm_storeu_pd(&x_new[i], b_r);
  }
  x_new[10] = x[10] + ((K[10] * x_tmp + K[21] * H_tmp) + K[32] * b_x);
  /*  norms quaternion */
  /* %% inverse quaternion  */
  b_r = _mm_loadu_pd(&x_new[0]);
  r1 = _mm_set1_pd(c_norm(&x_new[0]));
  _mm_storeu_pd(&x_new[0], _mm_div_pd(b_r, r1));
  b_r = _mm_loadu_pd(&x_new[2]);
  _mm_storeu_pd(&x_new[2], _mm_div_pd(b_r, r1));
  /*  norm quaternions */
}

/*
 * File trailer for ekf_correct.c
 *
 * [EOF]
 */
