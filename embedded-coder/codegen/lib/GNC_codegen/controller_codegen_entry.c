/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

/* Include Files */
#include "controller_codegen_entry.h"
#include "GNC_codegen_data.h"
#include "GNC_codegen_initialize.h"
#include "diag.h"
#include "eye.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * u : control command, desired canard angle (rad)
 *  r : roll angle target (rad)
 *  time : current time stamp (s)
 *  xR : roll state [roll angle (rad), roll rate (rad/s)]
 *  pdyn : dynamic pressure (Pa)
 *
 * Arguments    : double b_time
 *                const double xR[2]
 *                double pdyn
 *                double delta
 *                double *u
 *                double *b_r
 *                double *C_l_delta
 * Return Type  : void
 */
void controller_codegen_entry(double b_time, const double xR[2], double pdyn,
                              double delta, double *u, double *b_r,
                              double *C_l_delta)
{
  static double P_minus[4];
  static double w_old;
  double Q[4];
  double dv[4];
  double dv1[4];
  double K[2];
  double L_delta;
  double b_delta;
  double b_tmp_tmp;
  double c_r;
  double pdyn_params;
  double r_idx_0;
  int i;
  int i1;
  if (!isInitialized_GNC_codegen) {
    GNC_codegen_initialize();
  }
  /*     %% Constants */
  /*  [rad], limit output to this angle */
  /*  [rad/s^2 (angular accelaration) / rad (canard angle)] limit roll control
   * derivative for low authority conditions */
  /*  [Pa] deactivate at low authority near apogee */
  /*     %% Reference signal */
  /* %% Generates reference signal r (rad) for roll program */
  /* %% includes multiple roll angle steps */
  if ((b_time >= 22.0) && (b_time < 27.0)) {
    *b_r = 0.5;
  } else if ((b_time >= 27.0) && (b_time < 32.0)) {
    *b_r = -0.5;
  } else if ((b_time >= 32.0) && (b_time < 39.0)) {
    *b_r = 0.5;
  } else {
    *b_r = 0.0;
  }
  /*  r = 0; % deactivate roll program */
  /*     %% controller algorithm */
  /* %% Computes control signal of the adaptive LQR controller. */
  /* %% Coefficient Estimation */
  pdyn_params = pdyn * 0.00061367415999999994;
  b_delta = delta;
  /*  estimates the canard aerodynamic coefficients from canard angle, roll
   * rates, air data */
  /*  coeffs : canard coefficients C_l_delta and C_l_0 */
  /*  time : current time stamp (s) */
  /*  w : angular rate measurement (rad/s) */
  /*  d : canard angle measurement or command (rad) */
  /*  pdyn_params : dynamic pressure * constant parameters (pressure * area *
   * arm / inertia) */
  /*  dw/dt = c * Cl * d + C0 * c = phi_k' * p */
  /*  C_l_delta = canard lift coeff */
  /*  C_l_0 = rocket induced angular acceleration / (rho * area * arm) */
  /*     %% tuning parameters */
  /* %% covariance */
  diag(Q);
  /* %% small angle cutoff */
  if (fabs(delta) < 0.005) {
    /*  prevents high noise density for small delta from affecting estimate */
    b_delta = 0.0;
    /*  probably should make more rigorous */
  }
  /* %% lowpass */
  /*  time constant */
  /*     %% initialize */
  if (!P_minus_not_empty) {
    P_minus[0] = Q[0];
    P_minus[1] = Q[1];
    P_minus[2] = Q[2];
    P_minus[3] = Q[3];
    P_minus_not_empty = true;
  }
  if (!w_old_not_empty) {
    w_old = xR[1];
    w_old_not_empty = true;
  }
  /*     %% lowpass command and measurement */
  b_delta = 0.75 * d_old + 0.25 * b_delta;
  w_dot_old = 0.75 * w_dot_old + 0.25 * (xR[1] - w_old) / (b_time - t);
  /*     %% Kalman filter */
  r_idx_0 = pdyn_params * b_delta;
  /*  regression  */
  Q[0] += P_minus[0];
  Q[1] += P_minus[1];
  Q[2] += P_minus[2];
  Q[3] += P_minus[3];
  /*  covariance prediction */
  memset(&K[0], 0, sizeof(double) << 1);
  L_delta = r_idx_0 * Q[0];
  b_tmp_tmp = pdyn_params * Q[3];
  c_r = ((K[0] + L_delta) + pdyn_params * Q[1]) * r_idx_0 +
        ((K[1] + r_idx_0 * Q[2]) + b_tmp_tmp) * pdyn_params;
  /*  correction gain. the stuff inside in brackets is just a scalar so you can
   * just divide */
  K[0] = (L_delta + Q[2] * pdyn_params) / (c_r + 1.0);
  K[1] = (Q[1] * r_idx_0 + b_tmp_tmp) / (c_r + 1.0);
  L_delta = w_dot_old - (r_idx_0 * c[0] + pdyn_params * c[1]);
  /*  coefficient correction */
  eye(P_minus);
  c[0] += K[0] * L_delta;
  dv[0] = P_minus[0] - K[0] * r_idx_0;
  dv[1] = P_minus[1] - K[1] * r_idx_0;
  c[1] += K[1] * L_delta;
  dv[2] = P_minus[2] - K[0] * pdyn_params;
  dv[3] = P_minus[3] - K[1] * pdyn_params;
  memset(&P_minus[0], 0, sizeof(double) << 2);
  L_delta = dv[0];
  c_r = dv[1];
  r_idx_0 = dv[2];
  b_tmp_tmp = dv[3];
  for (i = 0; i < 2; i++) {
    i1 = i << 1;
    t = Q[i1];
    w_old = P_minus[i1] + L_delta * t;
    d_old = P_minus[i1 + 1] + c_r * t;
    t = Q[i1 + 1];
    w_old += r_idx_0 * t;
    P_minus[i1] = w_old;
    d_old += b_tmp_tmp * t;
    P_minus[i1 + 1] = d_old;
  }
  memset(&dv1[0], 0, sizeof(double) << 2);
  L_delta = P_minus[0];
  c_r = P_minus[1];
  r_idx_0 = P_minus[2];
  b_tmp_tmp = P_minus[3];
  for (i = 0; i < 2; i++) {
    t = dv[i];
    i1 = i << 1;
    w_old = dv1[i1] + L_delta * t;
    d_old = dv1[i1 + 1] + c_r * t;
    Q[i1] = K[0] * K[i];
    t = dv[i + 2];
    w_old += r_idx_0 * t;
    dv1[i1] = w_old;
    d_old += b_tmp_tmp * t;
    dv1[i1 + 1] = d_old;
    Q[i1 + 1] = K[1] * K[i];
  }
  P_minus[0] = dv1[0] + Q[0];
  P_minus[1] = dv1[1] + Q[1];
  P_minus[2] = dv1[2] + Q[2];
  P_minus[3] = dv1[3] + Q[3];
  /*  covariance correction. Joseph form for numerical stability */
  /*     %% update for next cycle */
  t = b_time;
  w_old = xR[1];
  d_old = b_delta;
  *C_l_delta = c[0];
  L_delta = c[0] * pdyn_params;
  if (fabs(L_delta) < 10.0) {
    if (L_delta >= 0.0) {
      L_delta = 10.0;
    } else {
      L_delta = -10.0;
    }
  }
  /* %% Control Law, Feedback + Feedforward tracking */
  /*  computes the optimal control signal for a flight condition  */
  /*  u : control signal, desired canard angle (rad) */
  /*  xR : roll state [roll angle (rad); roll rate (rad/s)] */
  /*  r : reference signal, desired roll angle (rad) */
  /*  L_delta : roll acceleration control derivative (rad/s^2 / rad) */
  /*     %% tuning parameters */
  /* %% Q_phi : weight of angle error, Q_omega : weight of rate error */
  /* %% Mode 1: Tracking (low rates), Mode 2: Damping only (high rates) */
  /* %% thresholds */
  /*  w < w_low: fully mode 1 */
  /*  w > w_high: fully mode 2 */
  /*     %% Control mode switching: linear crossfade */
  r_idx_0 = fmax(0.0, fmin(1.0, (fabs(xR[1]) - 0.5) / 0.5));
  /*  clamp to [0,1] */
  /*     %% feedback gains */
  L_delta = -1.0 / L_delta;
  b_tmp_tmp = (1.0 - r_idx_0) * 5.0;
  c_r = sqrt(b_tmp_tmp);
  K[0] = L_delta * c_r;
  /*     %% feedforward gain */
  /*  simplifies to this */
  /*     %% control command */
  /* %% limit output to allowable angle */
  *u = fmin(
      fmax((K[0] * xR[0] +
            L_delta * sqrt(2.0 * c_r + (b_tmp_tmp + r_idx_0 * 20.0)) * xR[1]) +
               -K[0] * *b_r,
           -0.17453292519943295),
      0.17453292519943295);
  /*  upper bounds */
  if (pdyn < 500.0) {
    /*  disable during low control authority */
    *u = 0.0;
  }
}

/*
 * File trailer for controller_codegen_entry.c
 *
 * [EOF]
 */
