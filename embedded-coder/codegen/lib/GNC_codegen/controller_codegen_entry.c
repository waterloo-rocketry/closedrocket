/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: controller_codegen_entry.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 23:24:33
 */

/* Include Files */
#include "controller_codegen_entry.h"
#include "GNC_codegen_data.h"
#include "GNC_codegen_initialize.h"
#include "GNC_codegen_types.h"
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
 *                const struct0_T *ctrl_mem_in
 *                double *u
 *                double *b_r
 *                struct0_T *ctrl_mem_out
 * Return Type  : void
 */
void controller_codegen_entry(double b_time, double dt_ctrl, const double xR[2],
                              double pdyn, double delta,
                              const struct0_T *ctrl_mem_in, double *u,
                              double *b_r, struct0_T *ctrl_mem_out)
{
  double P[4];
  double dv[4];
  double dv1[4];
  double dv2[4];
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
  c_delta = 0.75 * ctrl_mem_in->d_old + 0.25 * c_delta;
  w_dot = 0.75 * ctrl_mem_in->w_dot_old +
          0.25 * (xR[1] - ctrl_mem_in->w_old) / dt_ctrl;
  /*     %% Kalman filter */
  r_idx_0 = pdyn_params * c_delta;
  /*  regression  */
  diag(P);
  P[0] += ctrl_mem_in->P_minus[0];
  P[1] += ctrl_mem_in->P_minus[1];
  P[2] += ctrl_mem_in->P_minus[2];
  P[3] += ctrl_mem_in->P_minus[3];
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
  L_delta = w_dot - (r_idx_0 * ctrl_mem_in->coeffs[0] +
                     pdyn_params * ctrl_mem_in->coeffs[1]);
  /*  coefficient correction */
  eye(dv);
  ctrl_mem_out->coeffs[0] = ctrl_mem_in->coeffs[0] + K[0] * L_delta;
  dv1[0] = dv[0] - K[0] * r_idx_0;
  dv1[1] = dv[1] - K[1] * r_idx_0;
  ctrl_mem_out->coeffs[1] = ctrl_mem_in->coeffs[1] + K[1] * L_delta;
  dv1[2] = dv[2] - K[0] * pdyn_params;
  dv1[3] = dv[3] - K[1] * pdyn_params;
  memset(&dv[0], 0, sizeof(double) << 2);
  L_delta = dv1[0];
  b_tmp_tmp = dv1[1];
  c_r = dv1[2];
  r_idx_0 = dv1[3];
  for (i = 0; i < 2; i++) {
    i1 = i << 1;
    d = P[i1];
    d1 = dv[i1] + L_delta * d;
    d2 = dv[i1 + 1] + b_tmp_tmp * d;
    d = P[i1 + 1];
    d1 += c_r * d;
    dv[i1] = d1;
    d2 += r_idx_0 * d;
    dv[i1 + 1] = d2;
  }
  memset(&dv2[0], 0, sizeof(double) << 2);
  L_delta = dv[0];
  b_tmp_tmp = dv[1];
  c_r = dv[2];
  r_idx_0 = dv[3];
  for (i = 0; i < 2; i++) {
    d = dv1[i];
    i1 = i << 1;
    d1 = dv2[i1] + L_delta * d;
    d2 = dv2[i1 + 1] + b_tmp_tmp * d;
    P[i1] = K[0] * K[i];
    d = dv1[i + 2];
    d1 += c_r * d;
    dv2[i1] = d1;
    d2 += r_idx_0 * d;
    dv2[i1 + 1] = d2;
    P[i1 + 1] = K[1] * K[i];
  }
  ctrl_mem_out->P_minus[0] = dv2[0] + P[0];
  ctrl_mem_out->P_minus[1] = dv2[1] + P[1];
  ctrl_mem_out->P_minus[2] = dv2[2] + P[2];
  ctrl_mem_out->P_minus[3] = dv2[3] + P[3];
  /*  covariance correction. Joseph form for numerical stability */
  /*     %% update for next cycle */
  ctrl_mem_out->w_old = xR[1];
  ctrl_mem_out->d_old = c_delta;
  ctrl_mem_out->w_dot_old = w_dot;
  L_delta = 2.0 * ctrl_mem_out->coeffs[0] * pdyn_params;
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
