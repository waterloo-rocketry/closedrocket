/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: ekf_correct.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 23:24:33
 */

/* Include Files */
#include "ekf_correct.h"
#include "airdata_atmos.h"
#include "eye.h"
#include "inv.h"
#include "mpower.h"
#include "norm.h"
#include "rt_nonfinite.h"
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
 *                const double P[121]
 *                double y
 *                double b
 *                double x_new[11]
 *                double P_new[121]
 * Return Type  : void
 */
void b_ekf_correct(const double x[11], const double P[121], double y, double b,
                   double x_new[11], double P_new[121])
{
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
  for (i = 0; i < 11; i++) {
    layer_idx_3 = b_H[i];
    for (i1 = 0; i1 < 11; i1++) {
      layer_idx_3 += H[i1] * P[i1 + 11 * i];
    }
    b_H[i] = layer_idx_3;
    layer_idx_2 += layer_idx_3 * H[i];
  }
  layer_idx_2 = 1.0 / (layer_idx_2 + 100.0);
  memset(&x_new[0], 0, 11U * sizeof(double));
  for (i = 0; i < 11; i++) {
    for (i1 = 0; i1 < 11; i1++) {
      x_new[i1] += P[i1 + 11 * i] * H[i];
    }
  }
  for (i = 0; i < 11; i++) {
    x_new[i] *= layer_idx_2;
  }
  d_eye(E);
  for (i = 0; i < 11; i++) {
    for (i1 = 0; i1 < 11; i1++) {
      layer_idx_0 = i1 + 11 * i;
      b_E[layer_idx_0] = E[layer_idx_0] - x_new[i1] * H[i];
    }
  }
  /* %% correct covariance estimate */
  memset(&E[0], 0, 121U * sizeof(double));
  for (i = 0; i < 11; i++) {
    for (i1 = 0; i1 < 11; i1++) {
      layer_idx_2 = P[i1 + 11 * i];
      for (i2 = 0; i2 < 11; i2++) {
        layer_idx_0 = i2 + 11 * i;
        E[layer_idx_0] += b_E[i2 + 11 * i1] * layer_idx_2;
      }
    }
  }
  memset(&P_new[0], 0, 121U * sizeof(double));
  for (i = 0; i < 11; i++) {
    for (i1 = 0; i1 < 11; i1++) {
      layer_idx_2 = b_E[i + 11 * i1];
      for (i2 = 0; i2 < 11; i2++) {
        layer_idx_0 = i2 + 11 * i;
        P_new[layer_idx_0] += E[i2 + 11 * i1] * layer_idx_2;
      }
    }
  }
  for (i = 0; i < 11; i++) {
    for (i1 = 0; i1 < 11; i1++) {
      E[i1 + 11 * i] = x_new[i1] * 100.0 * x_new[i];
    }
  }
  for (i = 0; i < 121; i++) {
    P_new[i] += E[i];
  }
  /*  joseph stabilized */
  /* %% correct state estimate */
  for (i = 0; i < 11; i++) {
    x_new[i] = x[i] + x_new[i] * b_b;
  }
  /*  norms quaternion */
  /* %% inverse quaternion  */
  layer_idx_2 = c_norm(&x_new[0]);
  x_new[0] /= layer_idx_2;
  x_new[1] /= layer_idx_2;
  x_new[2] /= layer_idx_2;
  x_new[3] /= layer_idx_2;
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
 *                const double P[121]
 *                const double y[3]
 *                const double b[3]
 *                const double R[9]
 *                double x_new[11]
 *                double P_new[121]
 * Return Type  : void
 */
void ekf_correct(const double x[11], const double P[121], const double y[3],
                 const double b[3], const double R[9], double x_new[11],
                 double P_new[121])
{
  static const signed char c_b[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double E[121];
  double b_E[121];
  double H[33];
  double K[33];
  double b_P[33];
  double y_tmp[33];
  double b_b[9];
  double c_a[9];
  double c_x[3];
  double a;
  double b_a;
  double b_x;
  double d;
  double q_mag;
  int H_tmp;
  int b_H_tmp;
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
  q_mag = x[0] * 0.0;
  c_a[0] = q_mag;
  c_a[3] = x[0] * -b[2];
  c_a[6] = x[0] * b[1];
  c_a[1] = x[0] * b[2];
  c_a[4] = q_mag;
  c_a[7] = x[0] * -b[0];
  c_a[2] = x[0] * -b[1];
  c_a[5] = x[0] * b[0];
  c_a[8] = q_mag;
  for (i = 0; i < 3; i++) {
    H[i] = 2.0 * (x[0] * b[i] - c_x[i]);
    q_mag = x[i + 1];
    H_tmp = 3 * (i + 1);
    H[H_tmp] =
        2.0 * (((b_x * (double)c_b[3 * i] + x[1] * b[i]) - b[0] * q_mag) +
               c_a[3 * i]);
    b_H_tmp = 3 * i + 1;
    H[H_tmp + 1] =
        2.0 * (((b_x * (double)c_b[b_H_tmp] + x[2] * b[i]) - b[1] * q_mag) +
               c_a[b_H_tmp]);
    b_H_tmp = 3 * i + 2;
    H[H_tmp + 2] =
        2.0 * (((b_x * (double)c_b[b_H_tmp] + x[3] * b[i]) - b[2] * q_mag) +
               c_a[b_H_tmp]);
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
  for (i1 = 0; i1 < 11; i1++) {
    q_mag = K[3 * i1];
    b_H_tmp = 3 * i1 + 1;
    H_tmp = 3 * i1 + 2;
    for (i = 0; i < 11; i++) {
      d = P[i + 11 * i1];
      q_mag += H[3 * i] * d;
      K[b_H_tmp] += H[3 * i + 1] * d;
      K[H_tmp] += H[3 * i + 2] * d;
    }
    K[3 * i1] = q_mag;
  }
  memset(&b_P[0], 0, 33U * sizeof(double));
  for (i2 = 0; i2 < 3; i2++) {
    for (i1 = 0; i1 < 3; i1++) {
      q_mag = 0.0;
      for (i = 0; i < 11; i++) {
        q_mag += K[i2 + 3 * i] * y_tmp[i + 11 * i1];
      }
      b_H_tmp = i2 + 3 * i1;
      b_b[b_H_tmp] = q_mag + R[b_H_tmp];
    }
    for (i1 = 0; i1 < 11; i1++) {
      q_mag = y_tmp[i1 + 11 * i2];
      for (i = 0; i < 11; i++) {
        b_H_tmp = i + 11 * i2;
        b_P[b_H_tmp] += P[i + 11 * i1] * q_mag;
      }
    }
  }
  inv(b_b, c_a);
  memset(&K[0], 0, 33U * sizeof(double));
  for (i1 = 0; i1 < 3; i1++) {
    for (i = 0; i < 3; i++) {
      q_mag = c_a[i + 3 * i1];
      for (i2 = 0; i2 < 11; i2++) {
        b_H_tmp = i2 + 11 * i1;
        K[b_H_tmp] += b_P[i2 + 11 * i] * q_mag;
      }
    }
  }
  d_eye(E);
  for (i = 0; i < 11; i++) {
    q_mag = K[i];
    d = K[i + 11];
    b_x = K[i + 22];
    for (i1 = 0; i1 < 11; i1++) {
      b_H_tmp = i + 11 * i1;
      b_E[b_H_tmp] = E[b_H_tmp] - ((q_mag * H[3 * i1] + d * H[3 * i1 + 1]) +
                                   b_x * H[3 * i1 + 2]);
    }
  }
  /* %% correct covariance estimate */
  memset(&E[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    for (i = 0; i < 11; i++) {
      q_mag = P[i + 11 * i1];
      for (i2 = 0; i2 < 11; i2++) {
        b_H_tmp = i2 + 11 * i1;
        E[b_H_tmp] += b_E[i2 + 11 * i] * q_mag;
      }
    }
  }
  memset(&y_tmp[0], 0, 33U * sizeof(double));
  for (i1 = 0; i1 < 3; i1++) {
    for (i = 0; i < 3; i++) {
      q_mag = R[i + 3 * i1];
      for (i2 = 0; i2 < 11; i2++) {
        b_H_tmp = i2 + 11 * i1;
        y_tmp[b_H_tmp] += K[i2 + 11 * i] * q_mag;
      }
    }
  }
  memset(&P_new[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    for (i = 0; i < 11; i++) {
      q_mag = b_E[i1 + 11 * i];
      for (i2 = 0; i2 < 11; i2++) {
        b_H_tmp = i2 + 11 * i1;
        P_new[b_H_tmp] += E[i2 + 11 * i] * q_mag;
      }
    }
  }
  memset(&E[0], 0, 121U * sizeof(double));
  for (i1 = 0; i1 < 11; i1++) {
    for (i = 0; i < 3; i++) {
      q_mag = K[i1 + 11 * i];
      for (i2 = 0; i2 < 11; i2++) {
        b_H_tmp = i2 + 11 * i1;
        E[b_H_tmp] += y_tmp[i2 + 11 * i] * q_mag;
      }
    }
  }
  for (i = 0; i < 121; i++) {
    P_new[i] += E[i];
  }
  /*  joseph stabilized */
  /* %% correct state estimate */
  for (i = 0; i < 3; i++) {
    q_mag = x[i + 1];
    c_a[3 * i] = a * (double)c_b[3 * i] + 2.0 * x[1] * q_mag;
    b_H_tmp = 3 * i + 1;
    c_a[b_H_tmp] = a * (double)c_b[b_H_tmp] + 2.0 * x[2] * q_mag;
    b_H_tmp = 3 * i + 2;
    c_a[b_H_tmp] = a * (double)c_b[b_H_tmp] + 2.0 * x[3] * q_mag;
  }
  q_mag = b_a * 0.0;
  b_b[0] = q_mag;
  b_b[3] = b_a * -x[3];
  b_b[6] = b_a * x[2];
  b_b[1] = b_a * x[3];
  b_b[4] = q_mag;
  b_b[7] = b_a * -x[1];
  b_b[2] = b_a * -x[2];
  b_b[5] = b_a * x[1];
  b_b[8] = q_mag;
  for (i = 0; i < 9; i++) {
    c_a[i] -= b_b[i];
  }
  q_mag = b[0];
  d = b[1];
  b_x = b[2];
  for (i = 0; i < 3; i++) {
    c_x[i] = y[i] - ((c_a[i] * q_mag + c_a[i + 3] * d) + c_a[i + 6] * b_x);
  }
  q_mag = c_x[0];
  d = c_x[1];
  b_x = c_x[2];
  for (i = 0; i < 11; i++) {
    x_new[i] = x[i] + ((K[i] * q_mag + K[i + 11] * d) + K[i + 22] * b_x);
  }
  /*  norms quaternion */
  /* %% inverse quaternion  */
  q_mag = c_norm(&x_new[0]);
  x_new[0] /= q_mag;
  x_new[1] /= q_mag;
  x_new[2] /= q_mag;
  x_new[3] /= q_mag;
  /*  norm quaternions */
}

/*
 * File trailer for ekf_correct.c
 *
 * [EOF]
 */
