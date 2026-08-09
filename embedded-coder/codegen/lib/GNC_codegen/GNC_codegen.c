#include "GNC_codegen.h"
#include "GNC_codegen_types.h"
#include <math.h>
#include <string.h>

static const double dv[9] = {0.64, 0.0, 0.0, 0.0, 189.5, 0.0, 0.0, 0.0, 189.5};

static const double dv1[9] = {1.5625,
                              0.0,
                              0.0,
                              0.0,
                              0.0052770448548812663,
                              0.0,
                              0.0,
                              0.0,
                              0.0052770448548812663};

static double airdata_atmos(double altitude, double *airdata_temperature,
                            double *airdata_density,
                            double *airdata_sonic_speed, double *airdata_mach,
                            double *airdata_dynamic_pressure);

static double b_norm(const double x[3]);

static double c_norm(const double x[4]);

static void controller_codegen_entry_init(GNC_codegenStackData *SD);

static void dynamics_init(GNC_codegenStackData *SD);

static void dynamics_jacobian_init(GNC_codegenStackData *SD);

static void ekf_correct(const double x[11], const double P[121], double y,
                        double b, double x_new[11], double P_new[121]);

static void ekf_prefilter_imu_init(GNC_codegenStackData *SD);

static void mrdiv(const double A[33], const double B[9], double Y[33]);

static void pad_filter_init(GNC_codegenStackData *SD);

static double airdata_atmos(double altitude, double *airdata_temperature,
                            double *airdata_density,
                            double *airdata_sonic_speed, double *airdata_mach,
                            double *airdata_dynamic_pressure) {
  double airdata_pressure;
  double layer_idx_1;
  double layer_idx_2;
  double layer_idx_3;
  double pressure;
  double temperature;
  int layer_idx_0;
  altitude = (6.356766E+6 * altitude) / (6.356766E+6 - altitude);
  layer_idx_0 = 0;
  layer_idx_1 = 101325.0;
  layer_idx_2 = 288.15;
  layer_idx_3 = 0.0065;
  if (altitude > 11000.0) {
    if (altitude < 20000.0) {
      layer_idx_0 = 11000;
      layer_idx_2 = 216.65;
      layer_idx_3 = 0.0;
      pressure = 22632.1 * exp((-9.81 * (altitude - 11000.0)) / 62191.094035);
    } else {
      if (altitude < 32000.0) {
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
      pressure = layer_idx_1 * pow(1.0 - ((layer_idx_3 / layer_idx_2) *
                                          (altitude - ((double)layer_idx_0))),
                                   9.81 / (287.0579 * layer_idx_3));
    }
  } else {
    pressure = layer_idx_1 * pow(1.0 - ((layer_idx_3 / layer_idx_2) *
                                        (altitude - ((double)layer_idx_0))),
                                 9.81 / (287.0579 * layer_idx_3));
  }
  temperature =
      layer_idx_2 - (layer_idx_3 * (altitude - ((double)layer_idx_0)));
  *airdata_density = pressure / (287.0579 * temperature);
  *airdata_sonic_speed = sqrt(401.88106 * temperature);
  airdata_pressure = pressure;
  *airdata_temperature = temperature;
  *airdata_mach = 0.0;
  *airdata_dynamic_pressure = 0.0;
  return airdata_pressure;
}

static double b_norm(const double x[3]) {
  double absxk;
  double scale;
  double t;
  double y;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = ((y * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = ((y * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  return scale * sqrt(y);
}

static double c_norm(const double x[4]) {
  double absxk;
  double scale;
  double t;
  double y;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }
  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = ((y * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = ((y * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[3]);
  if (absxk > scale) {
    t = scale / absxk;
    y = ((y * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  return scale * sqrt(y);
}

static void controller_codegen_entry_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->param.Cn_alpha = 10.0;
  SD->pd->param.Cn_omega = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->param.J[i] = dv[i];
    SD->pd->param.Jinv[i] = dv1[i];
  }
  SD->pd->param.altitude_initial = 420.0;
  SD->pd->param.c_aero = -0.035602020206989993;
  SD->pd->param.c_aero_damping = 0.039162222227688996;
  SD->pd->param.c_canard = 0.00054680328244827419;
  SD->pd->param.d_ad[0] = 1.72;
  SD->pd->param.d_board[0] = 1.71;
  SD->pd->param.d_mti[0] = 1.72;
  SD->pd->param.g[0] = -9.81;
  SD->pd->param.d_ad[1] = -0.03;
  SD->pd->param.d_board[1] = -0.03;
  SD->pd->param.d_mti[1] = 0.0;
  SD->pd->param.g[1] = 0.0;
  SD->pd->param.d_ad[2] = 0.0;
  SD->pd->param.d_board[2] = -0.01;
  SD->pd->param.d_mti[2] = 0.0;
  SD->pd->param.g[2] = 0.0;
}

static void dynamics_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->d_param.Cn_alpha = 10.0;
  SD->pd->d_param.Cn_omega = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->d_param.J[i] = dv[i];
    SD->pd->d_param.Jinv[i] = dv1[i];
  }
  SD->pd->d_param.altitude_initial = 420.0;
  SD->pd->d_param.c_aero = -0.035602020206989993;
  SD->pd->d_param.c_aero_damping = 0.039162222227688996;
  SD->pd->d_param.c_canard = 0.00054680328244827419;
  SD->pd->d_param.d_ad[0] = 1.72;
  SD->pd->d_param.d_board[0] = 1.71;
  SD->pd->d_param.d_mti[0] = 1.72;
  SD->pd->d_param.g[0] = -9.81;
  SD->pd->d_param.d_ad[1] = -0.03;
  SD->pd->d_param.d_board[1] = -0.03;
  SD->pd->d_param.d_mti[1] = 0.0;
  SD->pd->d_param.g[1] = 0.0;
  SD->pd->d_param.d_ad[2] = 0.0;
  SD->pd->d_param.d_board[2] = -0.01;
  SD->pd->d_param.d_mti[2] = 0.0;
  SD->pd->d_param.g[2] = 0.0;
}

static void dynamics_jacobian_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->e_param.Cn_alpha = 10.0;
  SD->pd->e_param.Cn_omega = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->e_param.J[i] = dv[i];
    SD->pd->e_param.Jinv[i] = dv1[i];
  }
  SD->pd->e_param.altitude_initial = 420.0;
  SD->pd->e_param.c_aero = -0.035602020206989993;
  SD->pd->e_param.c_aero_damping = 0.039162222227688996;
  SD->pd->e_param.c_canard = 0.00054680328244827419;
  SD->pd->e_param.d_ad[0] = 1.72;
  SD->pd->e_param.d_board[0] = 1.71;
  SD->pd->e_param.d_mti[0] = 1.72;
  SD->pd->e_param.g[0] = -9.81;
  SD->pd->e_param.d_ad[1] = -0.03;
  SD->pd->e_param.d_board[1] = -0.03;
  SD->pd->e_param.d_mti[1] = 0.0;
  SD->pd->e_param.g[1] = 0.0;
  SD->pd->e_param.d_ad[2] = 0.0;
  SD->pd->e_param.d_board[2] = -0.01;
  SD->pd->e_param.d_mti[2] = 0.0;
  SD->pd->e_param.g[2] = 0.0;
}

static void ekf_correct(const double x[11], const double P[121], double y,
                        double b, double x_new[11], double P_new[121]) {
  double E[121];
  double b_E[121];
  double b_K[121];
  double c_E[121];
  double H[11];
  double K[11];
  double b_H[11];
  double absxk;
  double airdata_altitude_pressure;
  double altitude;
  double altitude_ratio;
  double b_b;
  double b_expl_temp;
  double c_H;
  double c_expl_temp;
  double d_expl_temp;
  double e_expl_temp;
  double expl_temp;
  double layer_idx_1;
  double layer_idx_2;
  double layer_idx_3;
  double q_mag;
  double scale;
  double t;
  double t0_pressure;
  int b_i;
  int b_k;
  int i;
  int i1;
  int i10;
  int i12;
  int i13;
  int i14;
  int i15;
  int i17;
  int i18;
  int i19;
  int i2;
  int i3;
  int i5;
  int i6;
  int i7;
  int i8;
  int i9;
  int k;
  int layer_idx_0;
  signed char b_I[121];
  t0_pressure = airdata_atmos(x[10], &expl_temp, &b_expl_temp, &c_expl_temp,
                              &d_expl_temp, &e_expl_temp);
  b_b = y - (t0_pressure + b);
  memset(&H[0], 0, 11U * (sizeof(double)));
  altitude_ratio = 6.356766E+6 / (6.356766E+6 - x[10]);
  altitude = altitude_ratio * x[10];
  layer_idx_0 = 0;
  layer_idx_1 = 101325.0;
  layer_idx_2 = 288.15;
  layer_idx_3 = 0.0065;
  if (altitude > 11000.0) {
    if (altitude < 20000.0) {
      airdata_altitude_pressure =
          (-3.5699790210323479 * (altitude_ratio * altitude_ratio)) *
          exp((-9.81 * (altitude - 11000.0)) / 62191.094035);
    } else {
      if (altitude < 32000.0) {
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
      airdata_altitude_pressure =
          ((((-layer_idx_1) * 9.81) / (layer_idx_2 * 287.0579)) *
           (altitude_ratio * altitude_ratio)) *
          pow(1.0 - ((layer_idx_3 / layer_idx_2) *
                     (altitude - ((double)layer_idx_0))),
              (9.81 / (287.0579 * layer_idx_3)) - 1.0);
    }
  } else {
    airdata_altitude_pressure =
        ((((-layer_idx_1) * 9.81) / (layer_idx_2 * 287.0579)) *
         (altitude_ratio * altitude_ratio)) *
        pow(1.0 - ((layer_idx_3 / layer_idx_2) *
                   (altitude - ((double)layer_idx_0))),
            (9.81 / (287.0579 * layer_idx_3)) - 1.0);
  }
  H[10] = airdata_altitude_pressure;
  memset(&b_H[0], 0, 11U * (sizeof(double)));
  c_H = 0.0;
  i = 0;
  for (i1 = 0; i1 < 11; i1++) {
    double d;
    d = b_H[i1];
    for (i3 = 0; i3 < 11; i3++) {
      d += H[i3] * P[i3 + i];
    }
    b_H[i1] = d;
    c_H += d * H[i1];
    i += 11;
  }
  for (i2 = 0; i2 < 11; i2++) {
    double d1;
    int i4;
    d1 = 0.0;
    i4 = 0;
    for (i5 = 0; i5 < 11; i5++) {
      d1 += P[i4 + i2] * H[i5];
      i4 += 11;
    }
    K[i2] = d1 / (c_H + 100.0);
  }
  memset(&b_I[0], 0, 121U * (sizeof(signed char)));
  k = 0;
  for (b_k = 0; b_k < 11; b_k++) {
    b_I[k] = (signed char)1;
    k += 12;
  }
  i6 = 0;
  for (i7 = 0; i7 < 11; i7++) {
    for (i8 = 0; i8 < 11; i8++) {
      int E_tmp;
      E_tmp = i8 + i6;
      E[E_tmp] = ((double)b_I[E_tmp]) - (K[i8] * H[i7]);
    }
    i6 += 11;
  }
  memset(&b_E[0], 0, 121U * (sizeof(double)));
  i9 = 0;
  for (i10 = 0; i10 < 11; i10++) {
    int i11;
    i11 = 0;
    for (i12 = 0; i12 < 11; i12++) {
      double d2;
      d2 = P[i12 + i9];
      for (i15 = 0; i15 < 11; i15++) {
        int b_E_tmp;
        b_E_tmp = i15 + i9;
        b_E[b_E_tmp] += E[i15 + i11] * d2;
      }
      i11 += 11;
    }
    i9 += 11;
  }
  memset(&c_E[0], 0, 121U * (sizeof(double)));
  i13 = 0;
  for (i14 = 0; i14 < 11; i14++) {
    int i16;
    i16 = 0;
    for (i18 = 0; i18 < 11; i18++) {
      double d3;
      d3 = E[i16 + i14];
      for (i19 = 0; i19 < 11; i19++) {
        int c_E_tmp;
        c_E_tmp = i19 + i13;
        c_E[c_E_tmp] += b_E[i19 + i16] * d3;
      }
      b_K[i18 + i13] = (K[i18] * 100.0) * K[i14];
      i16 += 11;
    }
    i13 += 11;
  }
  for (i17 = 0; i17 < 121; i17++) {
    P_new[i17] = c_E[i17] + b_K[i17];
  }
  for (b_i = 0; b_i < 11; b_i++) {
    x_new[b_i] = x[b_i] + (K[b_i] * b_b);
  }
  scale = 3.3121686421112381E-170;
  absxk = fabs(x_new[0]);
  if (absxk > 3.3121686421112381E-170) {
    q_mag = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    q_mag = t * t;
  }
  absxk = fabs(x_new[1]);
  if (absxk > scale) {
    t = scale / absxk;
    q_mag = ((q_mag * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    q_mag += t * t;
  }
  absxk = fabs(x_new[2]);
  if (absxk > scale) {
    t = scale / absxk;
    q_mag = ((q_mag * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    q_mag += t * t;
  }
  absxk = fabs(x_new[3]);
  if (absxk > scale) {
    t = scale / absxk;
    q_mag = ((q_mag * t) * t) + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    q_mag += t * t;
  }
  q_mag = scale * sqrt(q_mag);
  x_new[0] /= q_mag;
  x_new[1] /= q_mag;
  x_new[2] /= q_mag;
  x_new[3] /= q_mag;
}

static void ekf_prefilter_imu_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->c_param.Cn_alpha = 10.0;
  SD->pd->c_param.Cn_omega = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->c_param.J[i] = dv[i];
    SD->pd->c_param.Jinv[i] = dv1[i];
  }
  SD->pd->c_param.altitude_initial = 420.0;
  SD->pd->c_param.c_aero = -0.035602020206989993;
  SD->pd->c_param.c_aero_damping = 0.039162222227688996;
  SD->pd->c_param.c_canard = 0.00054680328244827419;
  SD->pd->c_param.d_ad[0] = 1.72;
  SD->pd->c_param.d_board[0] = 1.71;
  SD->pd->c_param.d_mti[0] = 1.72;
  SD->pd->c_param.g[0] = -9.81;
  SD->pd->c_param.d_ad[1] = -0.03;
  SD->pd->c_param.d_board[1] = -0.03;
  SD->pd->c_param.d_mti[1] = 0.0;
  SD->pd->c_param.g[1] = 0.0;
  SD->pd->c_param.d_ad[2] = 0.0;
  SD->pd->c_param.d_board[2] = -0.01;
  SD->pd->c_param.d_mti[2] = 0.0;
  SD->pd->c_param.g[2] = 0.0;
}

static void mrdiv(const double A[33], const double B[9], double Y[33]) {
  double b_A[9];
  double a21;
  double maxval;
  int k;
  int r1;
  int r2;
  int r3;
  memcpy(&b_A[0], &B[0], 9U * (sizeof(double)));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(B[0]);
  a21 = fabs(B[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }
  if (fabs(B[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }
  b_A[r2] = B[r2] / B[r1];
  b_A[r3] /= b_A[r1];
  b_A[r2 + 3] -= b_A[r2] * b_A[r1 + 3];
  b_A[r3 + 3] -= b_A[r3] * b_A[r1 + 3];
  b_A[r2 + 6] -= b_A[r2] * b_A[r1 + 6];
  b_A[r3 + 6] -= b_A[r3] * b_A[r1 + 6];
  if (fabs(b_A[r3 + 3]) > fabs(b_A[r2 + 3])) {
    int rtemp;
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }
  b_A[r3 + 3] /= b_A[r2 + 3];
  b_A[r3 + 6] -= b_A[r3 + 3] * b_A[r2 + 6];
  for (k = 0; k < 11; k++) {
    int Y_tmp;
    int b_Y_tmp;
    int c_Y_tmp;
    Y_tmp = k + (11 * r1);
    Y[Y_tmp] = A[k] / b_A[r1];
    b_Y_tmp = k + (11 * r2);
    Y[b_Y_tmp] = A[k + 11] - (Y[Y_tmp] * b_A[r1 + 3]);
    c_Y_tmp = k + (11 * r3);
    Y[c_Y_tmp] = A[k + 22] - (Y[Y_tmp] * b_A[r1 + 6]);
    Y[b_Y_tmp] /= b_A[r2 + 3];
    Y[c_Y_tmp] -= Y[b_Y_tmp] * b_A[r2 + 6];
    Y[c_Y_tmp] /= b_A[r3 + 6];
    Y[b_Y_tmp] -= Y[c_Y_tmp] * b_A[r3 + 3];
    Y[Y_tmp] -= Y[c_Y_tmp] * b_A[r3];
    Y[Y_tmp] -= Y[b_Y_tmp] * b_A[r2];
  }
}

static void pad_filter_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->b_param.Cn_alpha = 10.0;
  SD->pd->b_param.Cn_omega = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->b_param.J[i] = dv[i];
    SD->pd->b_param.Jinv[i] = dv1[i];
  }
  SD->pd->b_param.altitude_initial = 420.0;
  SD->pd->b_param.c_aero = -0.035602020206989993;
  SD->pd->b_param.c_aero_damping = 0.039162222227688996;
  SD->pd->b_param.c_canard = 0.00054680328244827419;
  SD->pd->b_param.d_ad[0] = 1.72;
  SD->pd->b_param.d_board[0] = 1.71;
  SD->pd->b_param.d_mti[0] = 1.72;
  SD->pd->b_param.g[0] = -9.81;
  SD->pd->b_param.d_ad[1] = -0.03;
  SD->pd->b_param.d_board[1] = -0.03;
  SD->pd->b_param.d_mti[1] = 0.0;
  SD->pd->b_param.g[1] = 0.0;
  SD->pd->b_param.d_ad[2] = 0.0;
  SD->pd->b_param.d_board[2] = -0.01;
  SD->pd->b_param.d_mti[2] = 0.0;
  SD->pd->b_param.g[2] = 0.0;
}

void GNC_codegen_initialize(GNC_codegenStackData *SD) {
  controller_codegen_entry_init(SD);
  pad_filter_init(SD);
  ekf_prefilter_imu_init(SD);
  dynamics_init(SD);
  dynamics_jacobian_init(SD);
}

void GNC_codegen_terminate(void) {}

void controller_codegen_entry(GNC_codegenStackData *SD, double b_time,
                              double dt_ctrl, const double where_it_is[2],
                              double pdyn, double delta_encoder,
                              struct0_T *ctrl_mem, double *u_motor,
                              double where_it_isnt[2], bool *w_status_ctrl) {
  static const double b_dv[5] = {-2.0943951023931953, 0.0, 0.0, 0.0, 0.0};
  static const double b_dv1[5] = {0.0, 0.0, 1.5, -1.5, 0.0};
  static const signed char b_iv[5] = {7, 13, 21, 31, 41};
  double P[4];
  double b_K[4];
  double dv2[4];
  double dv3[4];
  double dv4[4];
  double K[2];
  double b_r[2];
  double L_delta;
  double a;
  double b;
  double b_x;
  double blend;
  double c_r;
  double d1;
  double d10;
  double d11;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d8;
  double d9;
  double delta;
  double delta_lp;
  double deviation_idx_0;
  double deviation_idx_1;
  double pdyn_params;
  double r_idx_0;
  double r_phi;
  double r_w;
  double w_dot_lp;
  double x_tmp;
  int i1;
  int i2;
  int i3;
  int i4;
  int step_idx;
  *w_status_ctrl = true;
  r_phi = 0.0;
  r_w = 0.0;
  for (step_idx = 0; step_idx < 5; step_idx++) {
    int i;
    i = (int)b_iv[step_idx];
    if (b_time >= ((double)i)) {
      double d;
      double q;
      double r;
      double x;
      d = b_dv[step_idx];
      x = (b_dv1[step_idx] + ((b_time - ((double)i)) * d)) + 3.1415926535897931;
      q = fabs(x / 6.2831853071795862);
      if (fabs(q - floor(q + 0.5)) > (2.2204460492503131E-16 * q)) {
        r = fmod(x, 6.2831853071795862);
      } else {
        r = 0.0;
      }
      if (r == 0.0) {
        r = 0.0;
      } else if (r < 0.0) {
        r += 6.2831853071795862;
      }
      r_phi = r - 3.1415926535897931;
      r_w = d;
    }
  }
  where_it_isnt[0] = r_phi;
  where_it_isnt[1] = r_w;
  delta = delta_encoder / -2.0;
  pdyn_params = pdyn * SD->pd->param.c_canard;
  if (fabs(delta) < 0.005) {
    delta = 0.0;
  }
  delta_lp = (0.75 * ctrl_mem->delta_lp) + (0.25 * delta);
  w_dot_lp = (0.75 * ctrl_mem->w_dot_lp) +
             ((0.25 * (where_it_is[1] - ctrl_mem->w)) / dt_ctrl);
  r_idx_0 = pdyn_params * delta_lp;
  P[0] = ctrl_mem->P[0] + 5.0E-5;
  P[1] = ctrl_mem->P[1];
  P[2] = ctrl_mem->P[2];
  P[3] = ctrl_mem->P[3] + 5.0E-9;
  memset(&b_r[0], 0, (sizeof(double)) << 1);
  d1 = r_idx_0 * P[0];
  d2 = pdyn_params * P[3];
  c_r = (((b_r[0] + d1) + (pdyn_params * P[1])) * r_idx_0) +
        (((b_r[1] + (r_idx_0 * P[2])) + d2) * pdyn_params);
  K[0] = (d1 + (P[2] * pdyn_params)) / (c_r + 1.0);
  K[1] = ((P[1] * r_idx_0) + d2) / (c_r + 1.0);
  b = w_dot_lp -
      ((r_idx_0 * ctrl_mem->coeffs[0]) + (pdyn_params * ctrl_mem->coeffs[1]));
  ctrl_mem->coeffs[0] += K[0] * b;
  ctrl_mem->coeffs[1] += K[1] * b;
  dv2[0] = 1.0 - (K[0] * r_idx_0);
  dv2[1] = 0.0 - (K[1] * r_idx_0);
  dv2[2] = 0.0 - (K[0] * pdyn_params);
  dv2[3] = 1.0 - (K[1] * pdyn_params);
  memset(&dv3[0], 0, (sizeof(double)) << 2);
  i1 = 0;
  d3 = dv2[0];
  d4 = dv2[1];
  d5 = dv2[2];
  d6 = dv2[3];
  for (i2 = 0; i2 < 2; i2++) {
    double d7;
    d7 = P[i1];
    dv3[i1] += d3 * d7;
    dv3[i1 + 1] += d4 * d7;
    d7 = P[i1 + 1];
    dv3[i1] += d5 * d7;
    dv3[i1 + 1] += d6 * d7;
    i1 += 2;
  }
  memset(&dv4[0], 0, (sizeof(double)) << 2);
  i3 = 0;
  d8 = dv3[0];
  d9 = dv3[1];
  d10 = dv3[2];
  d11 = dv3[3];
  for (i4 = 0; i4 < 2; i4++) {
    double d12;
    d12 = dv2[i4];
    dv4[i3] += d8 * d12;
    dv4[i3 + 1] += d9 * d12;
    b_K[i3] = K[0] * K[i4];
    d12 = dv2[i4 + 2];
    dv4[i3] += d10 * d12;
    dv4[i3 + 1] += d11 * d12;
    b_K[i3 + 1] = K[1] * K[i4];
    i3 += 2;
  }
  ctrl_mem->P[0] = dv4[0] + b_K[0];
  ctrl_mem->P[1] = dv4[1] + b_K[1];
  ctrl_mem->P[2] = dv4[2] + b_K[2];
  ctrl_mem->P[3] = dv4[3] + b_K[3];
  ctrl_mem->w = where_it_is[1];
  ctrl_mem->delta_lp = delta_lp;
  ctrl_mem->w_dot_lp = w_dot_lp;
  L_delta = ctrl_mem->coeffs[0] * pdyn_params;
  if (fabs(L_delta) < 10.0) {
    if (L_delta >= 0.0) {
      L_delta = 10.0;
    } else {
      L_delta = -10.0;
    }
  }
  deviation_idx_0 = where_it_is[0] - r_phi;
  deviation_idx_1 = where_it_is[1] - r_w;
  blend = fmax(0.0, fmin(1.0, (fabs(deviation_idx_1) - 0.5) / 0.5));
  a = -1.0 / L_delta;
  x_tmp = (1.0 - blend) * 5.0;
  b_x = sqrt(x_tmp);
  if (deviation_idx_0 > 3.1415926535897931) {
    deviation_idx_0 -= 6.2831853071795862;
  } else if (deviation_idx_0 < -3.1415926535897931) {
    deviation_idx_0 += 6.2831853071795862;
  }
  *u_motor = fmin(fmax(((a * b_x) * deviation_idx_0) +
                           ((a * sqrt((2.0 * b_x) + (x_tmp + (blend * 20.0)))) *
                            deviation_idx_1),
                       -0.17453292519943295),
                  0.17453292519943295) *
             -2.0;
  if (pdyn < 500.0) {
    *u_motor = 0.0;
    *w_status_ctrl = false;
  }
  if (b_time < 7.0) {
    *u_motor = 0.0;
    *w_status_ctrl = false;
  }
}

void navigation_codegen_entry(GNC_codegenStackData *SD, double dt,
                              bool flight_phase, double x[11], double P[121],
                              struct1_T *bias, struct2_T *sens_filt,
                              const struct3_T *sens_in, double *cov_norm,
                              double roll_state[2], double *pdyn,
                              bool *w_status_nav) {
  static const double Q[121] = {
      1.0E-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.01,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.01,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.01,    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.001};
  static const double R[9] = {1.0E-7, 0.0, 0.0, 0.0,   1.0E-7,
                              0.0,    0.0, 0.0, 1.0E-7};
  static const signed char d_b[16] = {1, 0, 0,  0, 0, -1, 0, 0,
                                      0, 0, -1, 0, 0, 0,  0, -1};
  static const signed char b_b[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  static const signed char c_b[9] = {0, 0, 0, 0, 1, 0, 0, 0, 1};
  double F[121];
  double b_E[121];
  double b_F[121];
  double b_P[121];
  double c_E[121];
  double c_K[121];
  double c_P[121];
  double b_K[33];
  double e_dt[12];
  double b_x[11];
  double c_x[11];
  double b_dv[9];
  double b_n_tilde[9];
  double b_skewed_exp_w_tmp[9];
  double b_w_exp_tilde[9];
  double b_w_w_tmp[9];
  double l_a[9];
  double w_tilde_sq[9];
  double c_q[4];
  double q[4];
  double ST[3];
  double b_dt[3];
  double c_dt[3];
  double c_w_exp_tilde[3];
  double c_y[3];
  double dv2[3];
  double airdata_density;
  double b_a;
  double b_expl_temp;
  double c_expl_temp;
  double d_expl_temp;
  double e_expl_temp;
  double expl_temp;
  double f_expl_temp;
  double g_expl_temp;
  double h_expl_temp;
  double i_expl_temp;
  double j_expl_temp;
  double k_expl_temp;
  double l_expl_temp;
  double t1_density;
  int b_i;
  int e_k;
  int h_k;
  int i10;
  int i100;
  int i101;
  int i102;
  int i11;
  int i13;
  int i15;
  int i17;
  int i19;
  int i2;
  int i20;
  int i23;
  int i25;
  int i28;
  int i29;
  int i30;
  int i32;
  int i33;
  int i34;
  int i36;
  int i38;
  int i4;
  int i42;
  int i43;
  int i44;
  int i45;
  int i46;
  int i47;
  int i48;
  int i5;
  int i51;
  int i53;
  int i56;
  int i58;
  int i6;
  int i62;
  int i64;
  int i66;
  int i67;
  int i7;
  int i70;
  int i71;
  int i73;
  int i75;
  int i77;
  int i79;
  int i80;
  int i82;
  int i83;
  int i84;
  int i86;
  int i88;
  int i9;
  int i91;
  int i92;
  int i94;
  int i95;
  int i96;
  int i97;
  int i98;
  int i99;
  int j;
  int j_k;
  signed char g_I[121];
  if (!((int)flight_phase)) {
    double a[3];
    double a_norm;
    double b_filtered;
    double filtered;
    double t1_pressure;
    if (sens_in->board_accel.status) {
      sens_filt->board_accel[0] = (0.0001 * sens_in->board_accel.meas[0]) +
                                  (0.9999 * sens_filt->board_accel[0]);
      sens_filt->board_accel[1] = (0.0001 * sens_in->board_accel.meas[1]) +
                                  (0.9999 * sens_filt->board_accel[1]);
      sens_filt->board_accel[2] = (0.0001 * sens_in->board_accel.meas[2]) +
                                  (0.9999 * sens_filt->board_accel[2]);
    }
    if (sens_in->board_gyro.status) {
      sens_filt->board_gyro[0] = (0.0001 * sens_in->board_gyro.meas[0]) +
                                 (0.9999 * sens_filt->board_gyro[0]);
      sens_filt->board_gyro[1] = (0.0001 * sens_in->board_gyro.meas[1]) +
                                 (0.9999 * sens_filt->board_gyro[1]);
      sens_filt->board_gyro[2] = (0.0001 * sens_in->board_gyro.meas[2]) +
                                 (0.9999 * sens_filt->board_gyro[2]);
    }
    if (sens_in->mti_accel.status) {
      sens_filt->mti_accel[0] = (0.0001 * sens_in->mti_accel.meas[0]) +
                                (0.9999 * sens_filt->mti_accel[0]);
      sens_filt->mti_accel[1] = (0.0001 * sens_in->mti_accel.meas[1]) +
                                (0.9999 * sens_filt->mti_accel[1]);
      sens_filt->mti_accel[2] = (0.0001 * sens_in->mti_accel.meas[2]) +
                                (0.9999 * sens_filt->mti_accel[2]);
    }
    if (sens_in->mti_gyro.status) {
      sens_filt->mti_gyro[0] = (0.0001 * sens_in->mti_gyro.meas[0]) +
                               (0.9999 * sens_filt->mti_gyro[0]);
      sens_filt->mti_gyro[1] = (0.0001 * sens_in->mti_gyro.meas[1]) +
                               (0.9999 * sens_filt->mti_gyro[1]);
      sens_filt->mti_gyro[2] = (0.0001 * sens_in->mti_gyro.meas[2]) +
                               (0.9999 * sens_filt->mti_gyro[2]);
    }
    if (sens_in->ad_gyro.status) {
      sens_filt->ad_gyro[0] = (0.0001 * sens_in->ad_gyro.meas[0]) +
                              (0.9999 * sens_filt->ad_gyro[0]);
      sens_filt->ad_gyro[1] = (0.0001 * sens_in->ad_gyro.meas[1]) +
                              (0.9999 * sens_filt->ad_gyro[1]);
      sens_filt->ad_gyro[2] = (0.0001 * sens_in->ad_gyro.meas[2]) +
                              (0.9999 * sens_filt->ad_gyro[2]);
    }
    filtered = sens_filt->board_baro;
    if (sens_in->board_baro.status) {
      filtered = (0.0005 * sens_in->board_baro.meas) +
                 (0.9995 * sens_filt->board_baro);
    }
    sens_filt->board_baro = filtered;
    b_filtered = sens_filt->mti_baro;
    if (sens_in->mti_baro.status) {
      b_filtered =
          (0.0005 * sens_in->mti_baro.meas) + (0.9995 * sens_filt->mti_baro);
    }
    sens_filt->mti_baro = b_filtered;
    a[0] = (sens_filt->board_accel[0] * ((double)sens_in->board_accel.status)) +
           (sens_filt->mti_accel[0] * ((double)sens_in->mti_accel.status));
    a[1] = (sens_filt->board_accel[1] * ((double)sens_in->board_accel.status)) +
           (sens_filt->mti_accel[1] * ((double)sens_in->mti_accel.status));
    a[2] = (sens_filt->board_accel[2] * ((double)sens_in->board_accel.status)) +
           (sens_filt->mti_accel[2] * ((double)sens_in->mti_accel.status));
    *w_status_nav = true;
    a_norm = b_norm(a);
    if (a_norm < 1.0E-6) {
      q[0] = 1.0;
      q[1] = 0.0;
      q[2] = 0.0;
      q[3] = 0.0;
      *w_status_nav = false;
    } else {
      double d;
      double qw;
      double qy;
      double qz;
      qw = sqrt((0.5 * (a[0] / a_norm)) + 0.5);
      if (qw == 0.0) {
        qy = 1.0;
        qz = 0.0;
      } else {
        qy = (0.5 * (a[2] / a_norm)) / qw;
        qz = (-0.5 * (a[1] / a_norm)) / qw;
      }
      q[0] = qw;
      q[1] = 0.0;
      q[2] = qy;
      q[3] = qz;
      d = c_norm(q);
      q[0] = qw / d;
      q[1] = 0.0 / d;
      q[2] = qy / d;
      q[3] = qz / d;
    }
    x[0] = q[0];
    x[1] = q[1];
    x[2] = q[2];
    x[3] = q[3];
    x[10] = SD->pd->b_param.altitude_initial;
    memset(&ST[0], 0, 3U * (sizeof(double)));
    memset(&ST[0], 0, 3U * (sizeof(double)));
    x[4] = 0.0;
    x[7] = 0.0;
    bias->board_gyro[0] =
        sens_filt->board_gyro[0] * ((double)sens_in->board_gyro.status);
    bias->mti_gyro[0] =
        sens_filt->mti_gyro[0] * ((double)sens_in->mti_gyro.status);
    bias->ad_gyro[0] =
        sens_filt->ad_gyro[0] * ((double)sens_in->ad_gyro.status);
    bias->board_mag_earth[0] = 0.0;
    bias->mti_mag_earth[0] = 0.0;
    x[5] = 0.0;
    x[8] = 0.0;
    bias->board_gyro[1] =
        sens_filt->board_gyro[1] * ((double)sens_in->board_gyro.status);
    bias->mti_gyro[1] =
        sens_filt->mti_gyro[1] * ((double)sens_in->mti_gyro.status);
    bias->ad_gyro[1] =
        sens_filt->ad_gyro[1] * ((double)sens_in->ad_gyro.status);
    bias->board_mag_earth[1] = 0.0;
    bias->mti_mag_earth[1] = 0.0;
    x[6] = 0.0;
    x[9] = 0.0;
    bias->board_gyro[2] =
        sens_filt->board_gyro[2] * ((double)sens_in->board_gyro.status);
    bias->mti_gyro[2] =
        sens_filt->mti_gyro[2] * ((double)sens_in->mti_gyro.status);
    bias->ad_gyro[2] =
        sens_filt->ad_gyro[2] * ((double)sens_in->ad_gyro.status);
    bias->board_mag_earth[2] = 0.0;
    bias->mti_mag_earth[2] = 0.0;
    t1_pressure =
        airdata_atmos(SD->pd->b_param.altitude_initial, &e_expl_temp,
                      &t1_density, &f_expl_temp, &g_expl_temp, &h_expl_temp);
    bias->board_baro =
        (filtered - t1_pressure) * ((double)sens_in->board_baro.status);
    bias->mti_baro =
        (b_filtered - t1_pressure) * ((double)sens_in->mti_baro.status);
  } else {
    double C_total_a[3];
    double C_total_w[3];
    double a[3];
    double C_total_a_tmp;
    double C_total_a_tmp_tmp;
    double C_total_w_tmp;
    double C_total_w_tmp_tmp;
    double b_C_total_a_tmp_tmp;
    double b_C_total_w_tmp_tmp;
    double w_idx_0;
    double w_idx_1;
    double w_idx_2;
    int i;
    int k;
    bool exitg1;
    bool status_fast;
    bool y;
    status_fast = false;
    a[0] = 0.0;
    w_idx_0 = 0.0;
    i = 1000000 * ((int)sens_in->ad_gyro.status);
    C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_in->board_accel.status);
    b_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_in->mti_accel.status);
    C_total_a_tmp = C_total_a_tmp_tmp + b_C_total_a_tmp_tmp;
    C_total_a[0] = C_total_a_tmp;
    C_total_w_tmp_tmp =
        9.9999999999999981E+9 * ((double)sens_in->board_gyro.status);
    b_C_total_w_tmp_tmp =
        9.9999999999999981E+9 * ((double)sens_in->mti_gyro.status);
    C_total_w_tmp = C_total_w_tmp_tmp + b_C_total_w_tmp_tmp;
    C_total_w[0] = C_total_w_tmp + ((double)i);
    a[1] = 0.0;
    w_idx_1 = 0.0;
    C_total_a[1] = C_total_a_tmp;
    C_total_w[1] = C_total_w_tmp;
    a[2] = 0.0;
    w_idx_2 = 0.0;
    C_total_a[2] = C_total_a_tmp;
    C_total_w[2] = C_total_w_tmp;
    y = false;
    k = 0;
    exitg1 = false;
    while ((!((int)exitg1)) && (k < 3)) {
      if (C_total_a[k] == 0.0) {
        y = true;
        exitg1 = true;
      } else {
        k++;
      }
    }
    if (!((int)y)) {
      int b_k;
      bool b_y;
      b_y = false;
      b_k = 0;
      exitg1 = false;
      while ((!((int)exitg1)) && (b_k < 3)) {
        if (C_total_w[b_k] == 0.0) {
          b_y = true;
          exitg1 = true;
        } else {
          b_k++;
        }
      }
      if (!((int)b_y)) {
        double w_tilde[9];
        double b_w_idx_1_tmp;
        double d1;
        double d10;
        double d11;
        double d2;
        double d3;
        double d4;
        double d6;
        double d7;
        double d8;
        double d9;
        double w_idx_1_tmp;
        int i1;
        status_fast = true;
        w_idx_0 = (((C_total_w_tmp_tmp / C_total_w[0]) *
                    (sens_in->board_gyro.meas[0] - bias->board_gyro[0])) +
                   ((b_C_total_w_tmp_tmp / C_total_w[0]) *
                    (sens_in->mti_gyro.meas[0] - bias->mti_gyro[0]))) +
                  ((((double)i) / C_total_w[0]) *
                   (sens_in->ad_gyro.meas[0] - bias->ad_gyro[0]));
        d1 = 0.0 / C_total_w_tmp;
        w_idx_1_tmp = C_total_w_tmp_tmp / C_total_w_tmp;
        b_w_idx_1_tmp = b_C_total_w_tmp_tmp / C_total_w_tmp;
        w_idx_1 = ((w_idx_1_tmp *
                    (sens_in->board_gyro.meas[1] - bias->board_gyro[1])) +
                   (b_w_idx_1_tmp *
                    (sens_in->mti_gyro.meas[1] - bias->mti_gyro[1]))) +
                  (d1 * (sens_in->ad_gyro.meas[1] - bias->ad_gyro[1]));
        w_idx_2 = ((w_idx_1_tmp *
                    (sens_in->board_gyro.meas[2] - bias->board_gyro[2])) +
                   (b_w_idx_1_tmp *
                    (sens_in->mti_gyro.meas[2] - bias->mti_gyro[2]))) +
                  (d1 * (sens_in->ad_gyro.meas[2] - bias->ad_gyro[2]));
        w_tilde[0] = 0.0;
        w_tilde[3] = -w_idx_2;
        w_tilde[6] = w_idx_1;
        w_tilde[1] = w_idx_2;
        w_tilde[4] = 0.0;
        w_tilde[7] = -w_idx_0;
        w_tilde[2] = -w_idx_1;
        w_tilde[5] = w_idx_0;
        w_tilde[8] = 0.0;
        memset(&w_tilde_sq[0], 0, 9U * (sizeof(double)));
        i1 = 0;
        for (i2 = 0; i2 < 3; i2++) {
          int i3;
          i3 = 0;
          for (i4 = 0; i4 < 3; i4++) {
            double d5;
            d5 = w_tilde[i4 + i1];
            w_tilde_sq[i1] += w_tilde[i3] * d5;
            w_tilde_sq[i1 + 1] += w_tilde[i3 + 1] * d5;
            w_tilde_sq[i1 + 2] += w_tilde[i3 + 2] * d5;
            i3 += 3;
          }
          i1 += 3;
        }
        d2 = SD->pd->c_param.d_board[0];
        d3 = SD->pd->c_param.d_mti[0];
        d4 = SD->pd->c_param.d_ad[0];
        d6 = SD->pd->c_param.d_board[1];
        d7 = SD->pd->c_param.d_mti[1];
        d8 = SD->pd->c_param.d_ad[1];
        d9 = SD->pd->c_param.d_board[2];
        d10 = SD->pd->c_param.d_mti[2];
        d11 = SD->pd->c_param.d_ad[2];
        for (i6 = 0; i6 < 3; i6++) {
          double d12;
          double d13;
          double d14;
          double d15;
          double d17;
          d12 = w_tilde_sq[i6];
          d13 = d12 * d2;
          d14 = d12 * d3;
          d15 = d12 * d4;
          d12 = w_tilde_sq[i6 + 3];
          d13 += d12 * d6;
          d14 += d12 * d7;
          d15 += d12 * d8;
          d12 = w_tilde_sq[i6 + 6];
          d13 += d12 * d9;
          d14 += d12 * d10;
          d15 += d12 * d11;
          d17 = C_total_a[i6];
          a[i6] = (((C_total_a_tmp_tmp / d17) *
                    (sens_in->board_accel.meas[i6] - d13)) +
                   ((b_C_total_a_tmp_tmp / d17) *
                    (sens_in->mti_accel.meas[i6] - d14))) +
                  ((0.0 / d17) * (sens_in->ad_accel.meas[i6] - d15));
        }
      }
    }
    *w_status_nav = status_fast;
    if (status_fast) {
      double E[121];
      double P_pred[121];
      double K[33];
      double b_q[16];
      double d_I[16];
      double dq[16];
      double d_dt[12];
      double x_pred[11];
      double S[9];
      double b_P_pred[9];
      double dv3[9];
      double f_I[9];
      double h_a[9];
      double n_tilde[9];
      double skewed_exp_w_tmp[9];
      double w_exp_tilde[9];
      double w_w_tmp[9];
      double f_x[4];
      double r_q_tmp[4];
      double b_S[3];
      double b_dv1[3];
      double c_r_q_tmp[3];
      double dn[3];
      double i_x[3];
      double a_tmp;
      double b;
      double b_a_tmp;
      double b_dphi_tmp;
      double b_r_q_tmp;
      double c_a;
      double d20;
      double d21;
      double d22;
      double d23;
      double d24;
      double d25;
      double d26;
      double d27;
      double d28;
      double d29;
      double d30;
      double d31;
      double d34;
      double d36;
      double d37;
      double d39;
      double d43;
      double d44;
      double d45;
      double d_a;
      double d_x;
      double dphi;
      double dphi_tmp;
      double dq_idx_0;
      double dq_idx_1;
      double dq_idx_2;
      double dq_idx_3;
      double e_a;
      double e_x;
      double f_a;
      double g_a;
      double g_x;
      double h_x;
      double i_a;
      double j_a;
      double k_a;
      double m_a;
      double n_a;
      double n_idx_0;
      double n_idx_1;
      double n_idx_2;
      double q_mag;
      int c_k;
      int d_k;
      int f_k;
      int g_k;
      int i12;
      int i16;
      int i18;
      int i21;
      int i22;
      int i26;
      int i27;
      int i31;
      int i37;
      int i39;
      int i41;
      int i49;
      int i52;
      int i54;
      int i55;
      int i60;
      int i61;
      int i63;
      int i69;
      int i72;
      int i76;
      int i78;
      int i81;
      int i85;
      int i89;
      int i90;
      int i_k;
      signed char c_I[16];
      signed char b_I[9];
      signed char e_I[9];
      q_mag = c_norm(&x[0]);
      q[0] = x[0] / q_mag;
      q[1] = x[1] / q_mag;
      q[2] = x[2] / q_mag;
      q[3] = x[3] / q_mag;
      dphi_tmp = b_norm(&x[4]);
      b_dphi_tmp = dphi_tmp * dt;
      dphi = b_dphi_tmp / 2.0;
      if (dphi_tmp == 0.0) {
        dn[0] = 0.0;
        dn[1] = 0.0;
        dn[2] = 0.0;
        n_idx_0 = 0.0;
        n_idx_1 = 0.0;
        n_idx_2 = 0.0;
      } else {
        dn[0] = x[4] / dphi_tmp;
        dn[1] = x[5] / dphi_tmp;
        dn[2] = x[6] / dphi_tmp;
        n_idx_0 = x[4] / dphi_tmp;
        n_idx_1 = x[5] / dphi_tmp;
        n_idx_2 = x[6] / dphi_tmp;
      }
      b = sin(dphi);
      d_x = cos(dphi);
      n_tilde[0] = 0.0;
      n_tilde[3] = -n_idx_2;
      n_tilde[6] = n_idx_1;
      n_tilde[1] = n_idx_2;
      n_tilde[4] = 0.0;
      n_tilde[7] = -n_idx_0;
      n_tilde[2] = -n_idx_1;
      n_tilde[5] = n_idx_0;
      n_tilde[8] = 0.0;
      c_a = sin(b_dphi_tmp);
      e_x = cos(b_dphi_tmp);
      for (i5 = 0; i5 < 9; i5++) {
        b_I[i5] = (signed char)0;
      }
      memset(&b_n_tilde[0], 0, 9U * (sizeof(double)));
      c_k = 0;
      d_k = 0;
      for (e_k = 0; e_k < 3; e_k++) {
        int i8;
        b_I[c_k] = (signed char)1;
        i8 = 0;
        for (i9 = 0; i9 < 3; i9++) {
          double d16;
          d16 = n_tilde[i9 + d_k];
          b_n_tilde[d_k] += n_tilde[i8] * d16;
          b_n_tilde[d_k + 1] += n_tilde[i8 + 1] * d16;
          b_n_tilde[d_k + 2] += n_tilde[i8 + 2] * d16;
          i8 += 3;
        }
        c_k += 4;
        d_k += 3;
      }
      for (i7 = 0; i7 < 9; i7++) {
        w_exp_tilde[i7] = (((double)b_I[i7]) - (c_a * n_tilde[i7])) +
                          ((1.0 - e_x) * b_n_tilde[i7]);
      }
      airdata_atmos(x[10], &i_expl_temp, &airdata_density, &j_expl_temp,
                    &k_expl_temp, &l_expl_temp);
      a_tmp = 0.5 * airdata_density;
      d_a = (a_tmp * SD->pd->d_param.c_aero) * SD->pd->d_param.Cn_alpha;
      b_a_tmp = (-0.5 * airdata_density) * b_norm(&x[7]);
      e_a =
          (b_a_tmp * SD->pd->d_param.c_aero_damping) * SD->pd->d_param.Cn_omega;
      f_a = (x[0] * x[0]) - (((x[1] * x[1]) + (x[2] * x[2])) + (x[3] * x[3]));
      g_a = 2.0 * x[0];
      for (i10 = 0; i10 < 3; i10++) {
        double c_a_tmp;
        int d_a_tmp;
        int e_a_tmp;
        c_a_tmp = x[i10 + 1];
        h_a[3 * i10] =
            (f_a * ((double)b_b[3 * i10])) + ((2.0 * x[1]) * c_a_tmp);
        d_a_tmp = (3 * i10) + 1;
        h_a[d_a_tmp] =
            (f_a * ((double)b_b[d_a_tmp])) + ((2.0 * x[2]) * c_a_tmp);
        e_a_tmp = (3 * i10) + 2;
        h_a[e_a_tmp] =
            (f_a * ((double)b_b[e_a_tmp])) + ((2.0 * x[3]) * c_a_tmp);
      }
      b_dv[0] = 0.0;
      b_dv[3] = g_a * (-x[3]);
      b_dv[6] = g_a * x[2];
      b_dv[1] = g_a * x[3];
      b_dv[4] = 0.0;
      b_dv[7] = g_a * (-x[1]);
      b_dv[2] = g_a * (-x[2]);
      b_dv[5] = g_a * x[1];
      b_dv[8] = 0.0;
      for (i11 = 0; i11 < 9; i11++) {
        S[i11] = h_a[i11] - b_dv[i11];
      }
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
      f_x[0] = d_x;
      memset(&b_w_exp_tilde[0], 0, 9U * (sizeof(double)));
      memset(&c_w_exp_tilde[0], 0, 3U * (sizeof(double)));
      b_dv1[0] = 0.0;
      b_dv1[1] = d_a * (x[7] * x[9]);
      b_dv1[2] = d_a * ((-x[7]) * x[8]);
      dv2[0] = 0.0;
      dv2[1] = e_a * x[5];
      dv2[2] = e_a * x[6];
      i12 = 0;
      for (i13 = 0; i13 < 3; i13++) {
        int i14;
        f_x[i13 + 1] = dn[i13] * b;
        i14 = 0;
        for (i15 = 0; i15 < 3; i15++) {
          double d19;
          d19 = SD->pd->d_param.J[i15 + i12];
          b_w_exp_tilde[i12] += w_exp_tilde[i14] * d19;
          b_w_exp_tilde[i12 + 1] += w_exp_tilde[i14 + 1] * d19;
          b_w_exp_tilde[i12 + 2] += w_exp_tilde[i14 + 2] * d19;
          i14 += 3;
        }
        double d18;
        d18 = x[i13 + 4];
        c_w_exp_tilde[0] += b_w_exp_tilde[i12] * d18;
        c_w_exp_tilde[1] += b_w_exp_tilde[i12 + 1] * d18;
        c_w_exp_tilde[2] += b_w_exp_tilde[i12 + 2] * d18;
        b_dv1[i13] += dv2[i13];
        i12 += 3;
      }
      memset(&dv2[0], 0, 3U * (sizeof(double)));
      memset(&b_dt[0], 0, 3U * (sizeof(double)));
      memset(&c_dt[0], 0, 3U * (sizeof(double)));
      i16 = 0;
      d20 = dv2[0];
      d21 = dv2[1];
      d22 = dv2[2];
      d23 = b_dt[0];
      d24 = b_dt[1];
      d25 = b_dt[2];
      d26 = x[7];
      d27 = x[8];
      d28 = x[9];
      d29 = c_dt[0];
      d30 = c_dt[1];
      d31 = c_dt[2];
      for (i17 = 0; i17 < 3; i17++) {
        double d32;
        double d33;
        double d35;
        double d38;
        double d40;
        d32 = SD->pd->d_param.Jinv[i16];
        d33 = c_w_exp_tilde[i17];
        d20 += d32 * d33;
        d35 = b_dv1[i17];
        d23 += (dt * d32) * d35;
        d38 = S[i16];
        d29 += (dt * d38) * SD->pd->d_param.g[i17];
        d40 = d38 * d26;
        d32 = SD->pd->d_param.Jinv[i16 + 1];
        d21 += d32 * d33;
        d24 += (dt * d32) * d35;
        d38 = S[i16 + 1];
        d30 += (dt * d38) * SD->pd->d_param.g[i17];
        d40 += d38 * d27;
        d32 = SD->pd->d_param.Jinv[i16 + 2];
        d22 += d32 * d33;
        d25 += (dt * d32) * d35;
        d38 = S[i16 + 2];
        d31 += (dt * d38) * SD->pd->d_param.g[i17];
        d40 += d38 * d28;
        c_w_exp_tilde[i17] =
            (((w_exp_tilde[i17] * d26) + (w_exp_tilde[i17 + 3] * d27)) +
             (w_exp_tilde[i17 + 6] * d28)) +
            (dt * a[i17]);
        b_S[i17] = d40;
        i16 += 3;
      }
      memset(&c_q[0], 0, (sizeof(double)) << 2);
      i18 = 0;
      d34 = c_q[0];
      d36 = c_q[1];
      d37 = c_q[2];
      d39 = c_q[3];
      for (i19 = 0; i19 < 4; i19++) {
        double d41;
        d41 = f_x[i19];
        d34 += b_q[i18] * d41;
        d36 += b_q[i18 + 1] * d41;
        d37 += b_q[i18 + 2] * d41;
        d39 += b_q[i18 + 3] * d41;
        i18 += 4;
      }
      x_pred[0] = d34;
      x_pred[1] = d36;
      x_pred[2] = d37;
      x_pred[3] = d39;
      x_pred[4] = d20 + d23;
      x_pred[7] = c_w_exp_tilde[0] + d29;
      x_pred[5] = d21 + d24;
      x_pred[8] = c_w_exp_tilde[1] + d30;
      x_pred[6] = d22 + d25;
      x_pred[9] = c_w_exp_tilde[2] + d31;
      x_pred[10] = x[10] + (dt * b_S[0]);
      memset(&F[0], 0, 121U * (sizeof(double)));
      if (dphi_tmp > 0.0) {
        dq_idx_0 = d_x;
        dq_idx_1 = (x[4] / dphi_tmp) * b;
        dq_idx_2 = (x[5] / dphi_tmp) * b;
        dq_idx_3 = (x[6] / dphi_tmp) * b;
      } else {
        dq_idx_0 = 1.0;
        dq_idx_1 = 0.0;
        dq_idx_2 = 0.0;
        dq_idx_3 = 0.0;
      }
      for (i20 = 0; i20 < 16; i20++) {
        c_I[i20] = (signed char)0;
      }
      c_I[0] = (signed char)1;
      c_I[5] = (signed char)1;
      c_I[10] = (signed char)1;
      c_I[15] = (signed char)1;
      dq[0] = dq_idx_0;
      dq[4] = -dq_idx_1;
      dq[8] = -dq_idx_2;
      dq[12] = -dq_idx_3;
      dq[1] = dq_idx_1;
      dq[5] = dq_idx_0;
      dq[9] = dq_idx_3;
      dq[13] = -dq_idx_2;
      dq[2] = dq_idx_2;
      dq[6] = -dq_idx_3;
      dq[10] = dq_idx_0;
      dq[14] = dq_idx_1;
      dq[3] = dq_idx_3;
      dq[7] = dq_idx_2;
      dq[11] = -dq_idx_1;
      dq[15] = dq_idx_0;
      i21 = 0;
      i22 = 0;
      for (i23 = 0; i23 < 4; i23++) {
        int i24;
        d_I[i21] = (((double)c_I[i21]) - (q[0] * q[i23])) / q_mag;
        F[i22] = 0.0;
        d_I[i21 + 1] = (((double)c_I[i21 + 1]) - (q[1] * q[i23])) / q_mag;
        F[i22 + 1] = 0.0;
        d_I[i21 + 2] = (((double)c_I[i21 + 2]) - (q[2] * q[i23])) / q_mag;
        F[i22 + 2] = 0.0;
        d_I[i21 + 3] = (((double)c_I[i21 + 3]) - (q[3] * q[i23])) / q_mag;
        F[i22 + 3] = 0.0;
        i24 = 0;
        for (i25 = 0; i25 < 4; i25++) {
          double d42;
          d42 = d_I[i25 + i21];
          F[i22] += dq[i24] * d42;
          F[i22 + 1] += dq[i24 + 1] * d42;
          F[i22 + 2] += dq[i24 + 2] * d42;
          F[i22 + 3] += dq[i24 + 3] * d42;
          i24 += 4;
        }
        i21 += 4;
        i22 += 11;
      }
      i_a = 0.5 * dt;
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
      i26 = 0;
      i27 = 0;
      for (i28 = 0; i28 < 3; i28++) {
        F[i26 + 44] = i_a * b_q[i27 + 4];
        F[i26 + 45] = i_a * b_q[i27 + 5];
        F[i26 + 46] = i_a * b_q[i27 + 6];
        F[i26 + 47] = i_a * b_q[i27 + 7];
        i26 += 11;
        i27 += 4;
      }
      j_a = (a_tmp * SD->pd->e_param.c_aero) * SD->pd->e_param.Cn_alpha;
      k_a =
          (b_a_tmp * SD->pd->e_param.c_aero_damping) * SD->pd->e_param.Cn_omega;
      if (dphi_tmp == 0.0) {
        n_idx_0 = 0.0;
        n_idx_1 = 0.0;
        n_idx_2 = 0.0;
      } else {
        n_idx_0 = x[4] / dphi_tmp;
        n_idx_1 = x[5] / dphi_tmp;
        n_idx_2 = x[6] / dphi_tmp;
      }
      n_tilde[0] = 0.0;
      n_tilde[3] = -n_idx_2;
      n_tilde[6] = n_idx_1;
      n_tilde[1] = n_idx_2;
      n_tilde[4] = 0.0;
      n_tilde[7] = -n_idx_0;
      n_tilde[2] = -n_idx_1;
      n_tilde[5] = n_idx_0;
      n_tilde[8] = 0.0;
      for (i29 = 0; i29 < 9; i29++) {
        b_I[i29] = (signed char)0;
      }
      b_I[0] = (signed char)1;
      b_I[4] = (signed char)1;
      b_I[8] = (signed char)1;
      for (i30 = 0; i30 < 9; i30++) {
        w_w_tmp[i30] = dt * SD->pd->e_param.Jinv[i30];
      }
      memset(&c_y[0], 0, 3U * (sizeof(double)));
      i31 = 0;
      d43 = c_y[0];
      d44 = c_y[1];
      d45 = c_y[2];
      for (i32 = 0; i32 < 3; i32++) {
        double d46;
        d46 = x[i32 + 4];
        d43 += SD->pd->e_param.J[i31] * d46;
        d44 += SD->pd->e_param.J[i31 + 1] * d46;
        d45 += SD->pd->e_param.J[i31 + 2] * d46;
        i31 += 3;
      }
      h_a[0] = 0.0;
      h_a[3] = -x[6];
      h_a[6] = x[5];
      h_a[1] = x[6];
      h_a[4] = 0.0;
      h_a[7] = -x[4];
      h_a[2] = -x[5];
      h_a[5] = x[4];
      h_a[8] = 0.0;
      for (i33 = 0; i33 < 9; i33++) {
        e_I[i33] = (signed char)0;
      }
      b_dv[0] = 0.0;
      b_dv[3] = -d45;
      b_dv[6] = d44;
      b_dv[1] = d45;
      b_dv[4] = 0.0;
      b_dv[7] = -d43;
      b_dv[2] = -d44;
      b_dv[5] = d43;
      b_dv[8] = 0.0;
      memset(&l_a[0], 0, 9U * (sizeof(double)));
      f_k = 0;
      g_k = 0;
      for (h_k = 0; h_k < 3; h_k++) {
        int i35;
        e_I[f_k] = (signed char)1;
        i35 = 0;
        for (i36 = 0; i36 < 3; i36++) {
          double d47;
          d47 = SD->pd->e_param.J[i36 + g_k];
          l_a[g_k] += h_a[i35] * d47;
          l_a[g_k + 1] += h_a[i35 + 1] * d47;
          l_a[g_k + 2] += h_a[i35 + 2] * d47;
          i35 += 3;
        }
        f_k += 4;
        g_k += 3;
      }
      for (i34 = 0; i34 < 9; i34++) {
        b_dv[i34] -= l_a[i34];
      }
      memset(&b_w_w_tmp[0], 0, 9U * (sizeof(double)));
      i37 = 0;
      for (i38 = 0; i38 < 3; i38++) {
        int i40;
        i40 = 0;
        for (i42 = 0; i42 < 3; i42++) {
          double d48;
          int I_tmp;
          d48 = k_a * ((double)c_b[i42 + i37]);
          b_w_w_tmp[i37] += w_w_tmp[i40] * d48;
          b_w_w_tmp[i37 + 1] += w_w_tmp[i40 + 1] * d48;
          b_w_w_tmp[i37 + 2] += w_w_tmp[i40 + 2] * d48;
          I_tmp = i40 + i38;
          f_I[I_tmp] =
              ((double)e_I[I_tmp]) + (((w_w_tmp[i38] * b_dv[i40]) +
                                       (w_w_tmp[i38 + 3] * b_dv[i40 + 1])) +
                                      (w_w_tmp[i38 + 6] * b_dv[i40 + 2]));
          i40 += 3;
        }
        i37 += 3;
      }
      i39 = 0;
      i41 = 0;
      for (i43 = 0; i43 < 3; i43++) {
        F[i39 + 48] = f_I[i41] + b_w_w_tmp[i41];
        F[i39 + 49] = f_I[i41 + 1] + b_w_w_tmp[i41 + 1];
        F[i39 + 50] = f_I[i41 + 2] + b_w_w_tmp[i41 + 2];
        i39 += 11;
        i41 += 3;
      }
      b_dv[1] = j_a * x[9];
      b_dv[4] = 0.0;
      b_dv[7] = j_a * x[7];
      b_dv[2] = j_a * (-x[8]);
      b_dv[5] = j_a * (-x[7]);
      b_dv[8] = 0.0;
      g_x = 0.0;
      for (i44 = 0; i44 < 3; i44++) {
        double d49;
        double d50;
        double d51;
        int F_tmp;
        b_dv[3 * i44] = 0.0;
        F_tmp = 11 * (i44 + 7);
        d49 = 0.0;
        d50 = 0.0;
        d51 = 0.0;
        for (i45 = 0; i45 < 3; i45++) {
          double d52;
          d52 = b_dv[i45 + (3 * i44)];
          d49 += w_w_tmp[3 * i45] * d52;
          d50 += w_w_tmp[(3 * i45) + 1] * d52;
          d51 += w_w_tmp[(3 * i45) + 2] * d52;
        }
        F[F_tmp + 6] = d51;
        F[F_tmp + 5] = d50;
        F[F_tmp + 4] = d49;
        g_x += x[i44 + 1] * SD->pd->e_param.g[i44];
      }
      h_x = x[0];
      i_x[0] = (x[2] * SD->pd->e_param.g[2]) - (SD->pd->e_param.g[1] * x[3]);
      i_x[1] = (SD->pd->e_param.g[0] * x[3]) - (x[1] * SD->pd->e_param.g[2]);
      i_x[2] = (x[1] * SD->pd->e_param.g[1]) - (SD->pd->e_param.g[0] * x[2]);
      dv3[0] = 0.0;
      dv3[3] = x[0] * (-SD->pd->e_param.g[2]);
      dv3[6] = x[0] * SD->pd->e_param.g[1];
      dv3[1] = x[0] * SD->pd->e_param.g[2];
      dv3[4] = 0.0;
      dv3[7] = x[0] * (-SD->pd->e_param.g[0]);
      dv3[2] = x[0] * (-SD->pd->e_param.g[1]);
      dv3[5] = x[0] * SD->pd->e_param.g[0];
      dv3[8] = 0.0;
      for (i46 = 0; i46 < 3; i46++) {
        double b_F_tmp;
        int c_F_tmp;
        int d_F_tmp;
        int e_F_tmp;
        F[i46 + 7] = dt * (2.0 * ((h_x * SD->pd->e_param.g[i46]) - i_x[i46]));
        b_F_tmp = x[i46 + 1];
        c_F_tmp = 11 * (i46 + 1);
        F[c_F_tmp + 7] = dt * (2.0 * ((((g_x * ((double)b_b[3 * i46])) +
                                        (x[1] * SD->pd->e_param.g[i46])) -
                                       (SD->pd->e_param.g[0] * b_F_tmp)) +
                                      dv3[3 * i46]));
        d_F_tmp = (3 * i46) + 1;
        F[c_F_tmp + 8] = dt * (2.0 * ((((g_x * ((double)b_b[d_F_tmp])) +
                                        (x[2] * SD->pd->e_param.g[i46])) -
                                       (SD->pd->e_param.g[1] * b_F_tmp)) +
                                      dv3[d_F_tmp]));
        e_F_tmp = (3 * i46) + 2;
        F[c_F_tmp + 9] = dt * (2.0 * ((((g_x * ((double)b_b[e_F_tmp])) +
                                        (x[3] * SD->pd->e_param.g[i46])) -
                                       (SD->pd->e_param.g[2] * b_F_tmp)) +
                                      dv3[e_F_tmp]));
      }
      skewed_exp_w_tmp[0] = 0.0;
      skewed_exp_w_tmp[3] = -x[9];
      skewed_exp_w_tmp[6] = x[8];
      skewed_exp_w_tmp[1] = x[9];
      skewed_exp_w_tmp[4] = 0.0;
      skewed_exp_w_tmp[7] = -x[7];
      skewed_exp_w_tmp[2] = -x[8];
      skewed_exp_w_tmp[5] = x[7];
      skewed_exp_w_tmp[8] = 0.0;
      m_a = 0.5 * (dt * dt);
      memset(&b_skewed_exp_w_tmp[0], 0, 9U * (sizeof(double)));
      memset(&b_dv[0], 0, 9U * (sizeof(double)));
      memset(&b_n_tilde[0], 0, 9U * (sizeof(double)));
      r_q_tmp[0] = x[0];
      b_r_q_tmp = 0.0;
      for (i47 = 0; i47 < 3; i47++) {
        double d53;
        double d54;
        double d56;
        int b_n_tilde_tmp;
        int b_skewed_exp_w_tmp_tmp;
        int n_tilde_tmp;
        int skewed_exp_w_tmp_tmp;
        d53 = b_skewed_exp_w_tmp[3 * i47];
        d54 = b_dv[3 * i47];
        skewed_exp_w_tmp_tmp = (3 * i47) + 1;
        b_skewed_exp_w_tmp_tmp = (3 * i47) + 2;
        for (i48 = 0; i48 < 3; i48++) {
          double d55;
          double d57;
          int c_skewed_exp_w_tmp_tmp;
          int d_skewed_exp_w_tmp_tmp;
          int i50;
          i50 = i48 + (3 * i47);
          d55 = h_a[i50];
          d57 = skewed_exp_w_tmp[i50];
          d53 += skewed_exp_w_tmp[3 * i48] * d55;
          d54 += (2.0 * h_a[3 * i48]) * d57;
          c_skewed_exp_w_tmp_tmp = (3 * i48) + 1;
          b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] +=
              skewed_exp_w_tmp[c_skewed_exp_w_tmp_tmp] * d55;
          b_dv[skewed_exp_w_tmp_tmp] +=
              (2.0 * h_a[c_skewed_exp_w_tmp_tmp]) * d57;
          d_skewed_exp_w_tmp_tmp = (3 * i48) + 2;
          b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] +=
              skewed_exp_w_tmp[d_skewed_exp_w_tmp_tmp] * d55;
          b_dv[b_skewed_exp_w_tmp_tmp] +=
              (2.0 * h_a[d_skewed_exp_w_tmp_tmp]) * d57;
        }
        b_dv[3 * i47] = d54;
        b_skewed_exp_w_tmp[3 * i47] = d53;
        d56 = b_n_tilde[3 * i47];
        n_tilde_tmp = (3 * i47) + 1;
        b_n_tilde_tmp = (3 * i47) + 2;
        for (i53 = 0; i53 < 3; i53++) {
          double d58;
          int f_F_tmp;
          f_F_tmp = i53 + (3 * i47);
          F[(i53 + (11 * (i47 + 4))) + 7] =
              (dt * skewed_exp_w_tmp[f_F_tmp]) +
              (m_a * (b_skewed_exp_w_tmp[f_F_tmp] - b_dv[f_F_tmp]));
          d58 = n_tilde[f_F_tmp];
          d56 += n_tilde[3 * i53] * d58;
          b_n_tilde[n_tilde_tmp] += n_tilde[(3 * i53) + 1] * d58;
          b_n_tilde[b_n_tilde_tmp] += n_tilde[(3 * i53) + 2] * d58;
        }
        double d59;
        int g_F_tmp;
        int h_F_tmp;
        int i_F_tmp;
        b_n_tilde[3 * i47] = d56;
        g_F_tmp = 11 * (i47 + 7);
        F[g_F_tmp + 7] = (((double)b_I[3 * i47]) - (c_a * n_tilde[3 * i47])) +
                         ((1.0 - e_x) * d56);
        h_F_tmp = (3 * i47) + 1;
        F[g_F_tmp + 8] = (((double)b_I[h_F_tmp]) - (c_a * n_tilde[h_F_tmp])) +
                         ((1.0 - e_x) * b_n_tilde[h_F_tmp]);
        i_F_tmp = (3 * i47) + 2;
        F[g_F_tmp + 9] = (((double)b_I[i_F_tmp]) - (c_a * n_tilde[i_F_tmp])) +
                         ((1.0 - e_x) * b_n_tilde[i_F_tmp]);
        d59 = -x[i47 + 1];
        r_q_tmp[i47 + 1] = d59;
        b_r_q_tmp += d59 * x[i47 + 7];
      }
      c_r_q_tmp[0] = (r_q_tmp[2] * x[9]) - (r_q_tmp[3] * x[8]);
      c_r_q_tmp[1] = (r_q_tmp[3] * x[7]) - (r_q_tmp[1] * x[9]);
      c_r_q_tmp[2] = (r_q_tmp[1] * x[8]) - (r_q_tmp[2] * x[7]);
      i49 = 0;
      for (i51 = 0; i51 < 3; i51++) {
        double b_dt_tmp;
        double dt_tmp;
        dt_tmp = x[i51 + 7];
        d_dt[i51] = dt * (2.0 * ((r_q_tmp[0] * dt_tmp) - c_r_q_tmp[i51]));
        b_dt_tmp = r_q_tmp[i51 + 1];
        d_dt[i49 + 3] =
            dt * (2.0 *
                  ((((b_r_q_tmp * ((double)b_b[i49])) + (r_q_tmp[1] * dt_tmp)) -
                    (x[7] * b_dt_tmp)) +
                   (r_q_tmp[0] * skewed_exp_w_tmp[i49])));
        d_dt[i49 + 4] = dt * (2.0 * ((((b_r_q_tmp * ((double)b_b[i49 + 1])) +
                                       (r_q_tmp[2] * dt_tmp)) -
                                      (x[8] * b_dt_tmp)) +
                                     (r_q_tmp[0] * skewed_exp_w_tmp[i49 + 1])));
        d_dt[i49 + 5] = dt * (2.0 * ((((b_r_q_tmp * ((double)b_b[i49 + 2])) +
                                       (r_q_tmp[3] * dt_tmp)) -
                                      (x[9] * b_dt_tmp)) +
                                     (r_q_tmp[0] * skewed_exp_w_tmp[i49 + 2])));
        i49 += 3;
      }
      memset(&e_dt[0], 0, 12U * (sizeof(double)));
      i52 = 0;
      i54 = 0;
      i55 = 0;
      for (i56 = 0; i56 < 4; i56++) {
        int i57;
        i57 = 0;
        for (i58 = 0; i58 < 4; i58++) {
          int i59;
          i59 = (int)d_b[i58 + i55];
          e_dt[i54] += d_dt[i57] * ((double)i59);
          e_dt[i54 + 1] += d_dt[i57 + 1] * ((double)i59);
          e_dt[i54 + 2] += d_dt[i57 + 2] * ((double)i59);
          i57 += 3;
        }
        F[i52 + 10] = e_dt[i54];
        i52 += 11;
        i54 += 3;
        i55 += 4;
      }
      double o_a;
      n_a = (r_q_tmp[0] * r_q_tmp[0]) -
            (((r_q_tmp[1] * r_q_tmp[1]) + (r_q_tmp[2] * r_q_tmp[2])) +
             (r_q_tmp[3] * r_q_tmp[3]));
      o_a = 2.0 * r_q_tmp[0];
      b_dv[0] = 0.0;
      b_dv[3] = o_a * (-r_q_tmp[3]);
      b_dv[6] = o_a * r_q_tmp[2];
      b_dv[1] = o_a * r_q_tmp[3];
      b_dv[4] = 0.0;
      b_dv[7] = o_a * (-r_q_tmp[1]);
      b_dv[2] = o_a * (-r_q_tmp[2]);
      b_dv[5] = o_a * r_q_tmp[1];
      b_dv[8] = 0.0;
      i60 = 0;
      i61 = 0;
      for (i62 = 0; i62 < 3; i62++) {
        double f_a_tmp;
        f_a_tmp = r_q_tmp[i62 + 1];
        h_a[i60] = (n_a * ((double)b_b[i60])) + ((2.0 * r_q_tmp[1]) * f_a_tmp);
        h_a[i60 + 1] =
            (n_a * ((double)b_b[i60 + 1])) + ((2.0 * r_q_tmp[2]) * f_a_tmp);
        h_a[i60 + 2] =
            (n_a * ((double)b_b[i60 + 2])) + ((2.0 * r_q_tmp[3]) * f_a_tmp);
        F[i61 + 87] = dt * (h_a[i60] - b_dv[i60]);
        i60 += 3;
        i61 += 11;
      }
      F[120] = 1.0;
      memset(&b_F[0], 0, 121U * (sizeof(double)));
      i63 = 0;
      for (i64 = 0; i64 < 11; i64++) {
        int i65;
        i65 = 0;
        for (i67 = 0; i67 < 11; i67++) {
          double d60;
          d60 = P[i67 + i63];
          for (i71 = 0; i71 < 11; i71++) {
            int j_F_tmp;
            j_F_tmp = i71 + i63;
            b_F[j_F_tmp] += F[i71 + i65] * d60;
          }
          i65 += 11;
        }
        i63 += 11;
      }
      for (i66 = 0; i66 < 11; i66++) {
        int i68;
        i68 = 0;
        for (i70 = 0; i70 < 11; i70++) {
          double d61;
          int i74;
          d61 = 0.0;
          i74 = 0;
          for (i75 = 0; i75 < 11; i75++) {
            d61 += b_F[i74 + i66] * F[i74 + i70];
            i74 += 11;
          }
          int P_pred_tmp;
          P_pred_tmp = i68 + i66;
          P_pred[P_pred_tmp] = d61 + Q[P_pred_tmp];
          i68 += 11;
        }
      }
      i69 = 0;
      i72 = 0;
      for (i73 = 0; i73 < 3; i73++) {
        b_P_pred[i69] = P_pred[i72 + 48] + R[i69];
        b_P_pred[i69 + 1] = P_pred[i72 + 49] + R[i69 + 1];
        b_P_pred[i69 + 2] = P_pred[i72 + 50] + R[i69 + 2];
        i69 += 3;
        i72 += 11;
      }
      mrdiv(&P_pred[44], b_P_pred, K);
      memset(&g_I[0], 0, 121U * (sizeof(signed char)));
      i_k = 0;
      for (j_k = 0; j_k < 11; j_k++) {
        g_I[i_k] = (signed char)1;
        i_k += 12;
      }
      i76 = 0;
      for (i77 = 0; i77 < 4; i77++) {
        for (i79 = 0; i79 < 11; i79++) {
          int E_tmp;
          E_tmp = i79 + i76;
          E[E_tmp] = (double)g_I[E_tmp];
        }
        i76 += 11;
      }
      i78 = 0;
      for (i80 = 0; i80 < 3; i80++) {
        for (i82 = 0; i82 < 11; i82++) {
          int b_E_tmp;
          b_E_tmp = i82 + i78;
          E[b_E_tmp + 44] = ((double)g_I[b_E_tmp + 44]) - K[b_E_tmp];
        }
        i78 += 11;
      }
      i81 = 0;
      for (i83 = 0; i83 < 4; i83++) {
        for (i84 = 0; i84 < 11; i84++) {
          int c_E_tmp;
          c_E_tmp = (i84 + i81) + 77;
          E[c_E_tmp] = (double)g_I[c_E_tmp];
        }
        i81 += 11;
      }
      memset(&b_E[0], 0, 121U * (sizeof(double)));
      i85 = 0;
      for (i86 = 0; i86 < 11; i86++) {
        int i87;
        i87 = 0;
        for (i88 = 0; i88 < 11; i88++) {
          double d62;
          d62 = P_pred[i88 + i85];
          for (i92 = 0; i92 < 11; i92++) {
            int d_E_tmp;
            d_E_tmp = i92 + i85;
            b_E[d_E_tmp] += E[i92 + i87] * d62;
          }
          i87 += 11;
        }
        i85 += 11;
      }
      memset(&b_K[0], 0, 33U * (sizeof(double)));
      i89 = 0;
      i90 = 0;
      for (i91 = 0; i91 < 3; i91++) {
        int i93;
        i93 = 0;
        for (i94 = 0; i94 < 3; i94++) {
          double d63;
          d63 = R[i94 + i89];
          for (i96 = 0; i96 < 11; i96++) {
            int K_tmp;
            K_tmp = i96 + i90;
            b_K[K_tmp] += K[i96 + i93] * d63;
          }
          i93 += 11;
        }
        i89 += 3;
        i90 += 11;
      }
      memset(&c_E[0], 0, 121U * (sizeof(double)));
      memset(&c_K[0], 0, 121U * (sizeof(double)));
      for (i95 = 0; i95 < 11; i95++) {
        for (i98 = 0; i98 < 11; i98++) {
          double d64;
          d64 = E[i95 + (11 * i98)];
          for (i100 = 0; i100 < 11; i100++) {
            int e_E_tmp;
            e_E_tmp = i100 + (11 * i95);
            c_E[e_E_tmp] += b_E[i100 + (11 * i98)] * d64;
          }
        }
        for (i99 = 0; i99 < 3; i99++) {
          double d65;
          d65 = K[i95 + (11 * i99)];
          for (i102 = 0; i102 < 11; i102++) {
            int b_K_tmp;
            b_K_tmp = i102 + (11 * i95);
            c_K[b_K_tmp] += b_K[i102 + (11 * i99)] * d65;
          }
        }
      }
      for (i97 = 0; i97 < 121; i97++) {
        P[i97] = c_E[i97] + c_K[i97];
      }
      w_idx_0 -= x_pred[4];
      w_idx_1 -= x_pred[5];
      w_idx_2 -= x_pred[6];
      for (i101 = 0; i101 < 11; i101++) {
        x[i101] =
            x_pred[i101] + (((K[i101] * w_idx_0) + (K[i101 + 11] * w_idx_1)) +
                            (K[i101 + 22] * w_idx_2));
      }
      double b_q_mag;
      b_q_mag = c_norm(&x[0]);
      x[0] /= b_q_mag;
      x[1] /= b_q_mag;
      x[2] /= b_q_mag;
      x[3] /= b_q_mag;
    }
    if (sens_in->board_baro.status) {
      memcpy(&b_x[0], &x[0], 11U * (sizeof(double)));
      memcpy(&b_P[0], &P[0], 121U * (sizeof(double)));
      ekf_correct(b_x, b_P, sens_in->board_baro.meas, bias->board_baro, x, P);
    }
    if (sens_in->mti_baro.status) {
      memcpy(&c_x[0], &x[0], 11U * (sizeof(double)));
      memcpy(&c_P[0], &P[0], 121U * (sizeof(double)));
      ekf_correct(c_x, c_P, sens_in->mti_baro.meas, bias->mti_baro, x, P);
    }
  }
  b_a = b_norm(&x[7]);
  airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp,
                &d_expl_temp);
  *pdyn = (0.5 * t1_density) * (b_a * b_a);
  *cov_norm = 0.0;
  for (b_i = 0; b_i < 11; b_i++) {
    double s;
    s = 0.0;
    for (j = 0; j < 11; j++) {
      s += fabs(P[b_i + (11 * j)]);
    }
    if (s > (*cov_norm)) {
      *cov_norm = s;
    }
  }
  roll_state[0] =
      atan2(2.0 * ((x[2] * x[3]) + (x[0] * x[1])),
            (((x[0] * x[0]) - (x[1] * x[1])) - (x[2] * x[2])) + (x[3] * x[3]));
  roll_state[1] = x[4];
}
