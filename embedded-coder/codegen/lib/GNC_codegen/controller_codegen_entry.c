/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 01-Jun-2026 00:25:13
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
 *                double dt_ctrl
 *                const double xR[2]
 *                double pdyn
 *                double delta
 *                double w_old
 *                const double coeffs[2]
 *                const double P_minus[4]
 *                double d_old
 *                double w_dot_old
 *                double *u
 *                double *b_r
 *                double coeffs_ret[2]
 *                double *w_old_ret
 *                double P_minus_ret[4]
 *                double *d_old_ret
 *                double *w_dot_old_ret
 * Return Type  : void
 */
void controller_codegen_entry(double b_time, double dt_ctrl, const double xR[2],
                              double pdyn, double delta, double w_old,
                              const double coeffs[2], const double P_minus[4],
                              double d_old, double w_dot_old, double *u,
                              double *b_r, double coeffs_ret[2],
                              double *w_old_ret, double P_minus_ret[4],
                              double *d_old_ret, double *w_dot_old_ret)
{
  double P[4];
  double dv[4];
  double dv1[4];
  double K[2];
  double L_delta;
  double b_tmp_tmp;
  double c_delta;
  double c_r;
  double d;
  double d1;
  double d2;
  double pdyn_params;
  double r_idx_0;
  double w_dot;
  int P_minus_ret_tmp;
  int i;
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
  c_delta = delta / 2.0;
  /* %% small angle cutoff */
  if (fabs(c_delta) < 0.005) {
    /*  prevents high noise density for small delta from affecting estimate */
    c_delta = 0.0;
    /*  probably should make more rigorous */
  }
  /* %% lowpass */
  /*  time constant */
  /*     %% lowpass command and measurement */
  c_delta = 0.75 * d_old + 0.25 * c_delta;
  w_dot = 0.75 * w_dot_old + 0.25 * (xR[1] - w_old) / dt_ctrl;
  /*     %% Kalman filter */
  r_idx_0 = pdyn_params * c_delta;
  /*  regression  */
  diag(P);
  P[0] += P_minus[0];
  P[1] += P_minus[1];
  P[2] += P_minus[2];
  P[3] += P_minus[3];
  /*  covariance prediction */
  memset(&K[0], 0, sizeof(double) << 1);
  L_delta = r_idx_0 * P[0];
  b_tmp_tmp = pdyn_params * P[3];
  c_r = ((K[0] + L_delta) + pdyn_params * P[1]) * r_idx_0 +
        ((K[1] + r_idx_0 * P[2]) + b_tmp_tmp) * pdyn_params;
  /*  correction gain. the stuff inside in brackets is just a scalar so you can
   * just divide */
  K[0] = (L_delta + P[2] * pdyn_params) / (c_r + 1.0);
  K[1] = (P[1] * r_idx_0 + b_tmp_tmp) / (c_r + 1.0);
  L_delta = w_dot - (r_idx_0 * coeffs[0] + pdyn_params * coeffs[1]);
  /*  coefficient correction */
  eye(dv);
  coeffs_ret[0] = coeffs[0] + K[0] * L_delta;
  dv1[0] = dv[0] - K[0] * r_idx_0;
  dv1[1] = dv[1] - K[1] * r_idx_0;
  coeffs_ret[1] = coeffs[1] + K[1] * L_delta;
  dv1[2] = dv[2] - K[0] * pdyn_params;
  dv1[3] = dv[3] - K[1] * pdyn_params;
  memset(&dv[0], 0, sizeof(double) << 2);
  L_delta = dv1[0];
  b_tmp_tmp = dv1[1];
  c_r = dv1[2];
  r_idx_0 = dv1[3];
  for (i = 0; i < 2; i++) {
    P_minus_ret_tmp = i << 1;
    d = P[P_minus_ret_tmp];
    d1 = dv[P_minus_ret_tmp] + L_delta * d;
    d2 = dv[P_minus_ret_tmp + 1] + b_tmp_tmp * d;
    d = P[P_minus_ret_tmp + 1];
    d1 += c_r * d;
    dv[P_minus_ret_tmp] = d1;
    d2 += r_idx_0 * d;
    dv[P_minus_ret_tmp + 1] = d2;
  }
  memset(&P_minus_ret[0], 0, sizeof(double) << 2);
  L_delta = dv[0];
  b_tmp_tmp = dv[1];
  c_r = dv[2];
  r_idx_0 = dv[3];
  for (i = 0; i < 2; i++) {
    d = dv1[i];
    P_minus_ret_tmp = i << 1;
    d1 = P_minus_ret[P_minus_ret_tmp] + L_delta * d;
    d2 = P_minus_ret[P_minus_ret_tmp + 1] + b_tmp_tmp * d;
    P[P_minus_ret_tmp] = K[0] * K[i];
    d = dv1[i + 2];
    d1 += c_r * d;
    P_minus_ret[P_minus_ret_tmp] = d1;
    d2 += r_idx_0 * d;
    P_minus_ret[P_minus_ret_tmp + 1] = d2;
    P[P_minus_ret_tmp + 1] = K[1] * K[i];
  }
  P_minus_ret[0] += P[0];
  P_minus_ret[1] += P[1];
  P_minus_ret[2] += P[2];
  P_minus_ret[3] += P[3];
  /*  covariance correction. Joseph form for numerical stability */
  /*     %% update for next cycle */
  *w_old_ret = xR[1];
  *d_old_ret = c_delta;
  *w_dot_old_ret = w_dot;
  L_delta = 2.0 * coeffs_ret[0] * pdyn_params;
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
           -0.3490658503988659),
      0.3490658503988659);
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
