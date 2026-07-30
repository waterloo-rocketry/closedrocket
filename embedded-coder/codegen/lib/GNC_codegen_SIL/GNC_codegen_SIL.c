#include "GNC_codegen_SIL.h"
#include "GNC_codegen_SIL_types.h"
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

static void b_ekf_correct(const double x[11], const double P[121], double y,
                          double b, double x_new[11], double P_new[121]);

static double b_norm(const double x[3]);

static double c_norm(const double x[4]);

static void controller_codegen_entry_init(GNC_codegen_SILStackData *b_SD);

static void dynamics_init(GNC_codegen_SILStackData *b_SD);

static void dynamics_jacobian_init(GNC_codegen_SILStackData *b_SD);

static void ekf_correct(const double x[11], const double P[121],
                        const double y[3], const double b[3], const double R[9],
                        double x_new[11], double P_new[121]);

static void ekf_prefilter_imu_init(GNC_codegen_SILStackData *b_SD);

static void mrdiv(const double A[33], const double B[9], double Y[33]);

static void pad_filter_init(GNC_codegen_SILStackData *b_SD);

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
  altitude = 6.356766E+6 * altitude / (6.356766E+6 - altitude);
  layer_idx_0 = 0;
  layer_idx_1 = 101325.0;
  layer_idx_2 = 288.15;
  layer_idx_3 = 0.0065;
  if (altitude > 11000.0) {
    if (altitude < 20000.0) {
      layer_idx_0 = 11000;
      layer_idx_2 = 216.65;
      layer_idx_3 = 0.0;
      pressure = 22632.1 * exp(-9.81 * (altitude - 11000.0) / 62191.094035);
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
      pressure = layer_idx_1 * pow(1.0 - layer_idx_3 / layer_idx_2 *
                                             (altitude - (double)layer_idx_0),
                                   9.81 / (287.0579 * layer_idx_3));
    }
  } else {
    pressure = layer_idx_1 * pow(1.0 - layer_idx_3 / layer_idx_2 *
                                           (altitude - (double)layer_idx_0),
                                 9.81 / (287.0579 * layer_idx_3));
  }
  temperature = layer_idx_2 - layer_idx_3 * (altitude - (double)layer_idx_0);
  *airdata_density = pressure / (287.0579 * temperature);
  *airdata_sonic_speed = sqrt(401.88106 * temperature);
  airdata_pressure = pressure;
  *airdata_temperature = temperature;
  *airdata_mach = 0.0;
  *airdata_dynamic_pressure = 0.0;
  return airdata_pressure;
}

static void b_ekf_correct(const double x[11], const double P[121], double y,
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
  int i;
  int i1;
  int i10;
  int i11;
  int i12;
  int i2;
  int i3;
  int i4;
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
  memset(&H[0], 0, 11U * sizeof(double));
  altitude_ratio = 6.356766E+6 / (6.356766E+6 - x[10]);
  altitude = altitude_ratio * x[10];
  layer_idx_0 = 0;
  layer_idx_1 = 101325.0;
  layer_idx_2 = 288.15;
  layer_idx_3 = 0.0065;
  if (altitude > 11000.0) {
    if (altitude < 20000.0) {
      airdata_altitude_pressure =
          -3.5699790210323479 * (altitude_ratio * altitude_ratio) *
          exp(-9.81 * (altitude - 11000.0) / 62191.094035);
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
          -layer_idx_1 * 9.81 / (layer_idx_2 * 287.0579) *
          (altitude_ratio * altitude_ratio) *
          pow(1.0 -
                  layer_idx_3 / layer_idx_2 * (altitude - (double)layer_idx_0),
              9.81 / (287.0579 * layer_idx_3) - 1.0);
    }
  } else {
    airdata_altitude_pressure =
        -layer_idx_1 * 9.81 / (layer_idx_2 * 287.0579) *
        (altitude_ratio * altitude_ratio) *
        pow(1.0 - layer_idx_3 / layer_idx_2 * (altitude - (double)layer_idx_0),
            9.81 / (287.0579 * layer_idx_3) - 1.0);
  }
  H[10] = airdata_altitude_pressure;
  memset(&b_H[0], 0, 11U * sizeof(double));
  c_H = 0.0;
  for (i = 0; i < 11; i++) {
    double d;
    d = b_H[i];
    for (i2 = 0; i2 < 11; i2++) {
      d += H[i2] * P[i2 + 11 * i];
    }
    b_H[i] = d;
    c_H += d * H[i];
  }
  for (i1 = 0; i1 < 11; i1++) {
    double d1;
    d1 = 0.0;
    for (i3 = 0; i3 < 11; i3++) {
      d1 += P[i1 + 11 * i3] * H[i3];
    }
    K[i1] = d1 / (c_H + 100.0);
  }
  memset(&b_I[0], 0, 121U * sizeof(signed char));
  for (k = 0; k < 11; k++) {
    b_I[k + 11 * k] = 1;
  }
  for (i4 = 0; i4 < 11; i4++) {
    for (i5 = 0; i5 < 11; i5++) {
      int E_tmp;
      E_tmp = i5 + 11 * i4;
      E[E_tmp] = (double)b_I[E_tmp] - K[i5] * H[i4];
    }
  }
  memset(&b_E[0], 0, 121U * sizeof(double));
  for (i6 = 0; i6 < 11; i6++) {
    for (i7 = 0; i7 < 11; i7++) {
      double d2;
      d2 = P[i7 + 11 * i6];
      for (i9 = 0; i9 < 11; i9++) {
        int b_E_tmp;
        b_E_tmp = i9 + 11 * i6;
        b_E[b_E_tmp] += E[i9 + 11 * i7] * d2;
      }
    }
  }
  memset(&c_E[0], 0, 121U * sizeof(double));
  for (i8 = 0; i8 < 11; i8++) {
    for (i10 = 0; i10 < 11; i10++) {
      double d3;
      d3 = E[i8 + 11 * i10];
      for (i12 = 0; i12 < 11; i12++) {
        int c_E_tmp;
        c_E_tmp = i12 + 11 * i8;
        c_E[c_E_tmp] += b_E[i12 + 11 * i10] * d3;
      }
      b_K[i10 + 11 * i8] = K[i10] * 100.0 * K[i8];
    }
  }
  for (i11 = 0; i11 < 121; i11++) {
    P_new[i11] = c_E[i11] + b_K[i11];
  }
  for (b_i = 0; b_i < 11; b_i++) {
    x_new[b_i] = x[b_i] + K[b_i] * b_b;
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
    q_mag = q_mag * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    q_mag += t * t;
  }
  absxk = fabs(x_new[2]);
  if (absxk > scale) {
    t = scale / absxk;
    q_mag = q_mag * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    q_mag += t * t;
  }
  absxk = fabs(x_new[3]);
  if (absxk > scale) {
    t = scale / absxk;
    q_mag = q_mag * t * t + 1.0;
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
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
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
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  absxk = fabs(x[3]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }
  return scale * sqrt(y);
}

static void controller_codegen_entry_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->param.J[i] = dv[i];
    b_SD->pd->param.Jinv[i] = dv1[i];
  }
  b_SD->pd->param.altitude_initial = 420.0;
  b_SD->pd->param.c_aero = -0.035602020206989993;
  b_SD->pd->param.c_canard = 0.00054680328244827419;
  b_SD->pd->param.d_ad[0] = 1.72;
  b_SD->pd->param.d_board[0] = 1.71;
  b_SD->pd->param.d_mti[0] = 1.72;
  b_SD->pd->param.g[0] = -9.81;
  b_SD->pd->param.d_ad[1] = -0.03;
  b_SD->pd->param.d_board[1] = -0.03;
  b_SD->pd->param.d_mti[1] = 0.0;
  b_SD->pd->param.g[1] = 0.0;
  b_SD->pd->param.d_ad[2] = 0.0;
  b_SD->pd->param.d_board[2] = -0.01;
  b_SD->pd->param.d_mti[2] = 0.0;
  b_SD->pd->param.g[2] = 0.0;
}

static void dynamics_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->d_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->d_param.J[i] = dv[i];
    b_SD->pd->d_param.Jinv[i] = dv1[i];
  }
  b_SD->pd->d_param.altitude_initial = 420.0;
  b_SD->pd->d_param.c_aero = -0.035602020206989993;
  b_SD->pd->d_param.c_canard = 0.00054680328244827419;
  b_SD->pd->d_param.d_ad[0] = 1.72;
  b_SD->pd->d_param.d_board[0] = 1.71;
  b_SD->pd->d_param.d_mti[0] = 1.72;
  b_SD->pd->d_param.g[0] = -9.81;
  b_SD->pd->d_param.d_ad[1] = -0.03;
  b_SD->pd->d_param.d_board[1] = -0.03;
  b_SD->pd->d_param.d_mti[1] = 0.0;
  b_SD->pd->d_param.g[1] = 0.0;
  b_SD->pd->d_param.d_ad[2] = 0.0;
  b_SD->pd->d_param.d_board[2] = -0.01;
  b_SD->pd->d_param.d_mti[2] = 0.0;
  b_SD->pd->d_param.g[2] = 0.0;
}

static void dynamics_jacobian_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->e_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->e_param.J[i] = dv[i];
    b_SD->pd->e_param.Jinv[i] = dv1[i];
  }
  b_SD->pd->e_param.altitude_initial = 420.0;
  b_SD->pd->e_param.c_aero = -0.035602020206989993;
  b_SD->pd->e_param.c_canard = 0.00054680328244827419;
  b_SD->pd->e_param.d_ad[0] = 1.72;
  b_SD->pd->e_param.d_board[0] = 1.71;
  b_SD->pd->e_param.d_mti[0] = 1.72;
  b_SD->pd->e_param.g[0] = -9.81;
  b_SD->pd->e_param.d_ad[1] = -0.03;
  b_SD->pd->e_param.d_board[1] = -0.03;
  b_SD->pd->e_param.d_mti[1] = 0.0;
  b_SD->pd->e_param.g[1] = 0.0;
  b_SD->pd->e_param.d_ad[2] = 0.0;
  b_SD->pd->e_param.d_board[2] = -0.01;
  b_SD->pd->e_param.d_mti[2] = 0.0;
  b_SD->pd->e_param.g[2] = 0.0;
}

static void ekf_correct(const double x[11], const double P[121],
                        const double y[3], const double b[3], const double R[9],
                        double x_new[11], double P_new[121]) {
  static const signed char b_b[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double E[121];
  double b_E[121];
  double c_E[121];
  double c_K[121];
  double H[33];
  double K[33];
  double b_H[33];
  double b_K[33];
  double b_P[33];
  double y_tmp[33];
  double b_dv[9];
  double c_H[9];
  double c_a[9];
  double b_y[3];
  double c_x[3];
  double a;
  double absxk;
  double b_a;
  double b_x;
  double d11;
  double d12;
  double d13;
  double d14;
  double d15;
  double d16;
  double q_mag;
  double scale;
  double t;
  int i;
  int i1;
  int i10;
  int i11;
  int i12;
  int i13;
  int i14;
  int i15;
  int i16;
  int i17;
  int i18;
  int i19;
  int i2;
  int i20;
  int i21;
  int i22;
  int i23;
  int i24;
  int i25;
  int i26;
  int i27;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int i9;
  int k;
  signed char b_I[121];
  a = x[0] * x[0] - ((x[1] * x[1] + x[2] * x[2]) + x[3] * x[3]);
  b_a = 2.0 * x[0];
  memset(&H[0], 0, 33U * sizeof(double));
  b_x = (b[0] * x[1] + b[1] * x[2]) + b[2] * x[3];
  c_x[0] = x[2] * b[2] - b[1] * x[3];
  c_x[1] = b[0] * x[3] - x[1] * b[2];
  c_x[2] = x[1] * b[1] - b[0] * x[2];
  b_dv[0] = 0.0;
  b_dv[3] = x[0] * -b[2];
  b_dv[6] = x[0] * b[1];
  b_dv[1] = x[0] * b[2];
  b_dv[4] = 0.0;
  b_dv[7] = x[0] * -b[0];
  b_dv[2] = x[0] * -b[1];
  b_dv[5] = x[0] * b[0];
  b_dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    double H_tmp;
    int b_H_tmp;
    int c_H_tmp;
    int d_H_tmp;
    H[i] = 2.0 * (x[0] * b[i] - c_x[i]);
    H_tmp = x[i + 1];
    b_H_tmp = 3 * (i + 1);
    H[b_H_tmp] =
        2.0 * (((b_x * (double)b_b[3 * i] + x[1] * b[i]) - b[0] * H_tmp) +
               b_dv[3 * i]);
    c_H_tmp = 3 * i + 1;
    H[b_H_tmp + 1] =
        2.0 * (((b_x * (double)b_b[c_H_tmp] + x[2] * b[i]) - b[1] * H_tmp) +
               b_dv[c_H_tmp]);
    d_H_tmp = 3 * i + 2;
    H[b_H_tmp + 2] =
        2.0 * (((b_x * (double)b_b[d_H_tmp] + x[3] * b[i]) - b[2] * H_tmp) +
               b_dv[d_H_tmp]);
  }
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 11; i2++) {
      y_tmp[i2 + 11 * i1] = H[i1 + 3 * i2];
    }
  }
  memset(&b_H[0], 0, 33U * sizeof(double));
  for (i3 = 0; i3 < 11; i3++) {
    double d;
    int e_H_tmp;
    int f_H_tmp;
    d = b_H[3 * i3];
    e_H_tmp = 3 * i3 + 1;
    f_H_tmp = 3 * i3 + 2;
    for (i5 = 0; i5 < 11; i5++) {
      double d1;
      d1 = P[i5 + 11 * i3];
      d += H[3 * i5] * d1;
      b_H[e_H_tmp] += H[3 * i5 + 1] * d1;
      b_H[f_H_tmp] += H[3 * i5 + 2] * d1;
    }
    b_H[3 * i3] = d;
  }
  memset(&b_P[0], 0, 33U * sizeof(double));
  for (i4 = 0; i4 < 3; i4++) {
    for (i6 = 0; i6 < 11; i6++) {
      double d2;
      d2 = y_tmp[i6 + 11 * i4];
      for (i8 = 0; i8 < 11; i8++) {
        int P_tmp;
        P_tmp = i8 + 11 * i4;
        b_P[P_tmp] += P[i8 + 11 * i6] * d2;
      }
    }
    for (i7 = 0; i7 < 3; i7++) {
      double d3;
      d3 = 0.0;
      for (i9 = 0; i9 < 11; i9++) {
        d3 += b_H[i4 + 3 * i9] * y_tmp[i9 + 11 * i7];
      }
      int g_H_tmp;
      g_H_tmp = i4 + 3 * i7;
      c_H[g_H_tmp] = d3 + R[g_H_tmp];
    }
  }
  mrdiv(b_P, c_H, K);
  memset(&b_I[0], 0, 121U * sizeof(signed char));
  for (k = 0; k < 11; k++) {
    b_I[k + 11 * k] = 1;
  }
  for (i10 = 0; i10 < 11; i10++) {
    double d4;
    double d5;
    double d6;
    d4 = K[i10];
    d5 = K[i10 + 11];
    d6 = K[i10 + 22];
    for (i12 = 0; i12 < 11; i12++) {
      int E_tmp;
      E_tmp = i10 + 11 * i12;
      E[E_tmp] = (double)b_I[E_tmp] - ((d4 * H[3 * i12] + d5 * H[3 * i12 + 1]) +
                                       d6 * H[3 * i12 + 2]);
    }
  }
  memset(&b_E[0], 0, 121U * sizeof(double));
  for (i11 = 0; i11 < 11; i11++) {
    for (i13 = 0; i13 < 11; i13++) {
      double d7;
      d7 = P[i13 + 11 * i11];
      for (i15 = 0; i15 < 11; i15++) {
        int b_E_tmp;
        b_E_tmp = i15 + 11 * i11;
        b_E[b_E_tmp] += E[i15 + 11 * i13] * d7;
      }
    }
  }
  memset(&b_K[0], 0, 33U * sizeof(double));
  for (i14 = 0; i14 < 3; i14++) {
    for (i16 = 0; i16 < 3; i16++) {
      double d8;
      d8 = R[i16 + 3 * i14];
      for (i17 = 0; i17 < 11; i17++) {
        int K_tmp;
        K_tmp = i17 + 11 * i14;
        b_K[K_tmp] += K[i17 + 11 * i16] * d8;
      }
    }
  }
  memset(&c_E[0], 0, 121U * sizeof(double));
  memset(&c_K[0], 0, 121U * sizeof(double));
  for (i18 = 0; i18 < 11; i18++) {
    for (i19 = 0; i19 < 11; i19++) {
      double d9;
      d9 = E[i18 + 11 * i19];
      for (i22 = 0; i22 < 11; i22++) {
        int c_E_tmp;
        c_E_tmp = i22 + 11 * i18;
        c_E[c_E_tmp] += b_E[i22 + 11 * i19] * d9;
      }
    }
    for (i21 = 0; i21 < 3; i21++) {
      double d10;
      d10 = K[i18 + 11 * i21];
      for (i24 = 0; i24 < 11; i24++) {
        int b_K_tmp;
        b_K_tmp = i24 + 11 * i18;
        c_K[b_K_tmp] += b_K[i24 + 11 * i21] * d10;
      }
    }
  }
  for (i20 = 0; i20 < 121; i20++) {
    P_new[i20] = c_E[i20] + c_K[i20];
  }
  for (i23 = 0; i23 < 3; i23++) {
    double a_tmp;
    int b_a_tmp;
    int c_a_tmp;
    a_tmp = x[i23 + 1];
    c_a[3 * i23] = a * (double)b_b[3 * i23] + 2.0 * x[1] * a_tmp;
    b_a_tmp = 3 * i23 + 1;
    c_a[b_a_tmp] = a * (double)b_b[b_a_tmp] + 2.0 * x[2] * a_tmp;
    c_a_tmp = 3 * i23 + 2;
    c_a[c_a_tmp] = a * (double)b_b[c_a_tmp] + 2.0 * x[3] * a_tmp;
  }
  b_dv[0] = 0.0;
  b_dv[3] = b_a * -x[3];
  b_dv[6] = b_a * x[2];
  b_dv[1] = b_a * x[3];
  b_dv[4] = 0.0;
  b_dv[7] = b_a * -x[1];
  b_dv[2] = b_a * -x[2];
  b_dv[5] = b_a * x[1];
  b_dv[8] = 0.0;
  for (i25 = 0; i25 < 9; i25++) {
    c_a[i25] -= b_dv[i25];
  }
  d11 = b[0];
  d12 = b[1];
  d13 = b[2];
  for (i26 = 0; i26 < 3; i26++) {
    b_y[i26] =
        y[i26] - ((c_a[i26] * d11 + c_a[i26 + 3] * d12) + c_a[i26 + 6] * d13);
  }
  d14 = b_y[0];
  d15 = b_y[1];
  d16 = b_y[2];
  for (i27 = 0; i27 < 11; i27++) {
    x_new[i27] =
        x[i27] + ((K[i27] * d14 + K[i27 + 11] * d15) + K[i27 + 22] * d16);
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
    q_mag = q_mag * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    q_mag += t * t;
  }
  absxk = fabs(x_new[2]);
  if (absxk > scale) {
    t = scale / absxk;
    q_mag = q_mag * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    q_mag += t * t;
  }
  absxk = fabs(x_new[3]);
  if (absxk > scale) {
    t = scale / absxk;
    q_mag = q_mag * t * t + 1.0;
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

static void ekf_prefilter_imu_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->c_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->c_param.J[i] = dv[i];
    b_SD->pd->c_param.Jinv[i] = dv1[i];
  }
  b_SD->pd->c_param.altitude_initial = 420.0;
  b_SD->pd->c_param.c_aero = -0.035602020206989993;
  b_SD->pd->c_param.c_canard = 0.00054680328244827419;
  b_SD->pd->c_param.d_ad[0] = 1.72;
  b_SD->pd->c_param.d_board[0] = 1.71;
  b_SD->pd->c_param.d_mti[0] = 1.72;
  b_SD->pd->c_param.g[0] = -9.81;
  b_SD->pd->c_param.d_ad[1] = -0.03;
  b_SD->pd->c_param.d_board[1] = -0.03;
  b_SD->pd->c_param.d_mti[1] = 0.0;
  b_SD->pd->c_param.g[1] = 0.0;
  b_SD->pd->c_param.d_ad[2] = 0.0;
  b_SD->pd->c_param.d_board[2] = -0.01;
  b_SD->pd->c_param.d_mti[2] = 0.0;
  b_SD->pd->c_param.g[2] = 0.0;
}

static void mrdiv(const double A[33], const double B[9], double Y[33]) {
  double b_A[9];
  double a21;
  double maxval;
  int k;
  int r1;
  int r2;
  int r3;
  memcpy(&b_A[0], &B[0], 9U * sizeof(double));
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
    Y_tmp = k + 11 * r1;
    Y[Y_tmp] = A[k] / b_A[r1];
    b_Y_tmp = k + 11 * r2;
    Y[b_Y_tmp] = A[k + 11] - Y[Y_tmp] * b_A[r1 + 3];
    c_Y_tmp = k + 11 * r3;
    Y[c_Y_tmp] = A[k + 22] - Y[Y_tmp] * b_A[r1 + 6];
    Y[b_Y_tmp] /= b_A[r2 + 3];
    Y[c_Y_tmp] -= Y[b_Y_tmp] * b_A[r2 + 6];
    Y[c_Y_tmp] /= b_A[r3 + 6];
    Y[b_Y_tmp] -= Y[c_Y_tmp] * b_A[r3 + 3];
    Y[Y_tmp] -= Y[c_Y_tmp] * b_A[r3];
    Y[Y_tmp] -= Y[b_Y_tmp] * b_A[r2];
  }
}

static void pad_filter_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->b_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->b_param.J[i] = dv[i];
    b_SD->pd->b_param.Jinv[i] = dv1[i];
  }
  b_SD->pd->b_param.altitude_initial = 420.0;
  b_SD->pd->b_param.c_aero = -0.035602020206989993;
  b_SD->pd->b_param.c_canard = 0.00054680328244827419;
  b_SD->pd->b_param.d_ad[0] = 1.72;
  b_SD->pd->b_param.d_board[0] = 1.71;
  b_SD->pd->b_param.d_mti[0] = 1.72;
  b_SD->pd->b_param.g[0] = -9.81;
  b_SD->pd->b_param.d_ad[1] = -0.03;
  b_SD->pd->b_param.d_board[1] = -0.03;
  b_SD->pd->b_param.d_mti[1] = 0.0;
  b_SD->pd->b_param.g[1] = 0.0;
  b_SD->pd->b_param.d_ad[2] = 0.0;
  b_SD->pd->b_param.d_board[2] = -0.01;
  b_SD->pd->b_param.d_mti[2] = 0.0;
  b_SD->pd->b_param.g[2] = 0.0;
}

void GNC_codegen_SIL_initialize(GNC_codegen_SILStackData *b_SD) {
  controller_codegen_entry_init(b_SD);
  pad_filter_init(b_SD);
  ekf_prefilter_imu_init(b_SD);
  dynamics_init(b_SD);
  dynamics_jacobian_init(b_SD);
}

void GNC_codegen_SIL_terminate(void) {}

void controller_codegen_entry(GNC_codegen_SILStackData *b_SD, double b_time,
                              double dt_ctrl, const double where_it_is[2],
                              double pdyn, double delta_encoder,
                              struct0_T *ctrl_mem, double *u_motor,
                              double where_it_isnt[2],
                              boolean_T *w_status_ctrl) {
  static const double b_dv[5] = {-1.5707963267948966, 0.0, 0.0, 0.0, 0.0};
  static const signed char b_iv[5] = {7, 15, 25, 35, 45};
  static const signed char b_iv1[5] = {0, 0, 1, -1, 0};
  double P[4];
  double b_K[4];
  double b_dv1[4];
  double dv2[4];
  double dv3[4];
  double K[2];
  double b_r[2];
  double L_delta;
  double a;
  double b;
  double b_x;
  double blend;
  double c_r;
  double d1;
  double d11;
  double d12;
  double d13;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
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
  int i3;
  int step_idx;
  *w_status_ctrl = true;
  r_phi = 0.0;
  r_w = 0.0;
  for (step_idx = 0; step_idx < 5; step_idx++) {
    int i;
    i = b_iv[step_idx];
    if (b_time >= i) {
      double d;
      double q;
      double r;
      double x;
      d = b_dv[step_idx];
      x = ((double)b_iv1[step_idx] + (b_time - (double)i) * d) +
          3.1415926535897931;
      q = fabs(x / 6.2831853071795862);
      if (fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q) {
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
  pdyn_params = pdyn * b_SD->pd->param.c_canard;
  if (fabs(delta) < 0.005) {
    delta = 0.0;
  }
  delta_lp = 0.75 * ctrl_mem->delta_lp + 0.25 * delta;
  w_dot_lp = 0.75 * ctrl_mem->w_dot_lp +
             0.25 * (where_it_is[1] - ctrl_mem->w) / dt_ctrl;
  r_idx_0 = pdyn_params * delta_lp;
  P[0] = ctrl_mem->P[0] + 5.0E-5;
  P[1] = ctrl_mem->P[1];
  P[2] = ctrl_mem->P[2];
  P[3] = ctrl_mem->P[3] + 5.0E-9;
  memset(&b_r[0], 0, sizeof(double) << 1);
  d1 = r_idx_0 * P[0];
  d2 = pdyn_params * P[3];
  c_r = ((b_r[0] + d1) + pdyn_params * P[1]) * r_idx_0 +
        ((b_r[1] + r_idx_0 * P[2]) + d2) * pdyn_params;
  K[0] = (d1 + P[2] * pdyn_params) / (c_r + 1.0);
  K[1] = (P[1] * r_idx_0 + d2) / (c_r + 1.0);
  b = w_dot_lp -
      (r_idx_0 * ctrl_mem->coeffs[0] + pdyn_params * ctrl_mem->coeffs[1]);
  ctrl_mem->coeffs[0] += K[0] * b;
  ctrl_mem->coeffs[1] += K[1] * b;
  b_dv1[0] = 1.0 - K[0] * r_idx_0;
  b_dv1[1] = 0.0 - K[1] * r_idx_0;
  b_dv1[2] = 0.0 - K[0] * pdyn_params;
  b_dv1[3] = 1.0 - K[1] * pdyn_params;
  memset(&dv2[0], 0, sizeof(double) << 2);
  d3 = b_dv1[0];
  d4 = b_dv1[1];
  d5 = b_dv1[2];
  d6 = b_dv1[3];
  for (i1 = 0; i1 < 2; i1++) {
    double d10;
    double d7;
    double d8;
    int i2;
    i2 = i1 << 1;
    d7 = P[i2];
    d8 = dv2[i2] + d3 * d7;
    d10 = dv2[i2 + 1] + d4 * d7;
    d7 = P[i2 + 1];
    d8 += d5 * d7;
    dv2[i2] = d8;
    d10 += d6 * d7;
    dv2[i2 + 1] = d10;
  }
  memset(&dv3[0], 0, sizeof(double) << 2);
  d9 = dv2[0];
  d11 = dv2[1];
  d12 = dv2[2];
  d13 = dv2[3];
  for (i3 = 0; i3 < 2; i3++) {
    double d14;
    double d15;
    double d16;
    int i4;
    d14 = b_dv1[i3];
    i4 = i3 << 1;
    d15 = dv3[i4] + d9 * d14;
    d16 = dv3[i4 + 1] + d11 * d14;
    b_K[i4] = K[0] * K[i3];
    d14 = b_dv1[i3 + 2];
    d15 += d12 * d14;
    dv3[i4] = d15;
    d16 += d13 * d14;
    dv3[i4 + 1] = d16;
    b_K[i4 + 1] = K[1] * K[i3];
  }
  ctrl_mem->P[0] = dv3[0] + b_K[0];
  ctrl_mem->P[1] = dv3[1] + b_K[1];
  ctrl_mem->P[2] = dv3[2] + b_K[2];
  ctrl_mem->P[3] = dv3[3] + b_K[3];
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
  *u_motor = fmin(fmax(a * b_x * deviation_idx_0 +
                           a * sqrt(2.0 * b_x + (x_tmp + blend * 20.0)) *
                               deviation_idx_1,
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

void navigation_codegen_entry(GNC_codegen_SILStackData *b_SD, double dt,
                              boolean_T flight_phase, double x[11],
                              double P[121], struct1_T *bias,
                              struct2_T *sens_filt, struct3_T *sens_in,
                              double *cov_norm, double roll_state[2],
                              double *pdyn, boolean_T *w_status_nav) {
  static const double Q[121] = {
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      1.0E-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.001,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.001,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.001,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0001,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.001};
  static const double R[9] = {1.0E-9, 0.0, 0.0, 0.0,   1.0E-9,
                              0.0,    0.0, 0.0, 1.0E-9};
  static const double b[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double F[121];
  double b_E[121];
  double b_F[121];
  double b_P[121];
  double c_E[121];
  double c_K[121];
  double c_P[121];
  double d_P[121];
  double e_P[121];
  double b_K[33];
  double b_W_dt[16];
  double c_b[16];
  double d_b[16];
  double d_x[11];
  double e_x[11];
  double f_x[11];
  double g_x[11];
  double b_dv[9];
  double b_n_tilde[9];
  double b_w_exp_tilde[9];
  double c_skewed_exp_w_tmp[9];
  double w_tilde_sq[9];
  double c_q[4];
  double q[4];
  double b_ST[3];
  double b_dt[3];
  double c_dt[3];
  double c_w_exp_tilde[3];
  double dv3[3];
  double b_a;
  double b_expl_temp;
  double b_x;
  double c_expl_temp;
  double c_x;
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
  double m_expl_temp;
  double n_expl_temp;
  double o_expl_temp;
  double p_expl_temp;
  double t1_density;
  int b_i;
  int c_k;
  int d_k;
  int e_k;
  int f_k;
  int i1;
  int i10;
  int i11;
  int i12;
  int i13;
  int i14;
  int i15;
  int i16;
  int i19;
  int i2;
  int i21;
  int i22;
  int i23;
  int i24;
  int i25;
  int i26;
  int i27;
  int i28;
  int i29;
  int i3;
  int i32;
  int i33;
  int i34;
  int i35;
  int i36;
  int i37;
  int i38;
  int i4;
  int i40;
  int i41;
  int i42;
  int i43;
  int i44;
  int i45;
  int i46;
  int i47;
  int i48;
  int i49;
  int i5;
  int i50;
  int i51;
  int i52;
  int i53;
  int i54;
  int i55;
  int i56;
  int i57;
  int i58;
  int i59;
  int i6;
  int i60;
  int i61;
  int i62;
  int i63;
  int i7;
  int i8;
  int i9;
  int j;
  signed char d_I[121];
  b_x = b_norm(sens_in->board_mag.meas);
  if (fabs(b_x) >= 1.0E-6) {
    sens_in->board_mag.meas[0] /= b_x;
    sens_in->board_mag.meas[1] /= b_x;
    sens_in->board_mag.meas[2] /= b_x;
  }
  c_x = b_norm(sens_in->mti_mag.meas);
  if (fabs(c_x) >= 1.0E-6) {
    sens_in->mti_mag.meas[0] /= c_x;
    sens_in->mti_mag.meas[1] /= c_x;
    sens_in->mti_mag.meas[2] /= c_x;
  }
  if (!flight_phase) {
    double ST[9];
    double e_a[9];
    double a[3];
    double a_norm;
    double b_filtered;
    double c_a;
    double d19;
    double d20;
    double d21;
    double d24;
    double d25;
    double d26;
    double d_a;
    double filtered;
    if (sens_in->board_accel.status) {
      sens_filt->board_accel[0] = 0.0001 * sens_in->board_accel.meas[0] +
                                  0.9999 * sens_filt->board_accel[0];
      sens_filt->board_accel[1] = 0.0001 * sens_in->board_accel.meas[1] +
                                  0.9999 * sens_filt->board_accel[1];
      sens_filt->board_accel[2] = 0.0001 * sens_in->board_accel.meas[2] +
                                  0.9999 * sens_filt->board_accel[2];
    }
    if (sens_in->board_gyro.status) {
      sens_filt->board_gyro[0] = 0.0001 * sens_in->board_gyro.meas[0] +
                                 0.9999 * sens_filt->board_gyro[0];
      sens_filt->board_gyro[1] = 0.0001 * sens_in->board_gyro.meas[1] +
                                 0.9999 * sens_filt->board_gyro[1];
      sens_filt->board_gyro[2] = 0.0001 * sens_in->board_gyro.meas[2] +
                                 0.9999 * sens_filt->board_gyro[2];
    }
    if (sens_in->mti_accel.status) {
      sens_filt->mti_accel[0] = 0.0001 * sens_in->mti_accel.meas[0] +
                                0.9999 * sens_filt->mti_accel[0];
      sens_filt->mti_accel[1] = 0.0001 * sens_in->mti_accel.meas[1] +
                                0.9999 * sens_filt->mti_accel[1];
      sens_filt->mti_accel[2] = 0.0001 * sens_in->mti_accel.meas[2] +
                                0.9999 * sens_filt->mti_accel[2];
    }
    if (sens_in->mti_gyro.status) {
      sens_filt->mti_gyro[0] =
          0.0001 * sens_in->mti_gyro.meas[0] + 0.9999 * sens_filt->mti_gyro[0];
      sens_filt->mti_gyro[1] =
          0.0001 * sens_in->mti_gyro.meas[1] + 0.9999 * sens_filt->mti_gyro[1];
      sens_filt->mti_gyro[2] =
          0.0001 * sens_in->mti_gyro.meas[2] + 0.9999 * sens_filt->mti_gyro[2];
    }
    if (sens_in->ad_accel.status) {
      sens_filt->ad_accel[0] =
          0.0001 * sens_in->ad_accel.meas[0] + 0.9999 * sens_filt->ad_accel[0];
      sens_filt->ad_accel[1] =
          0.0001 * sens_in->ad_accel.meas[1] + 0.9999 * sens_filt->ad_accel[1];
      sens_filt->ad_accel[2] =
          0.0001 * sens_in->ad_accel.meas[2] + 0.9999 * sens_filt->ad_accel[2];
    }
    if (sens_in->ad_gyro.status) {
      sens_filt->ad_gyro[0] =
          0.0001 * sens_in->ad_gyro.meas[0] + 0.9999 * sens_filt->ad_gyro[0];
      sens_filt->ad_gyro[1] =
          0.0001 * sens_in->ad_gyro.meas[1] + 0.9999 * sens_filt->ad_gyro[1];
      sens_filt->ad_gyro[2] =
          0.0001 * sens_in->ad_gyro.meas[2] + 0.9999 * sens_filt->ad_gyro[2];
    }
    filtered = sens_filt->board_baro;
    if (sens_in->board_baro.status) {
      filtered =
          0.0005 * sens_in->board_baro.meas + 0.9995 * sens_filt->board_baro;
    }
    sens_filt->board_baro = filtered;
    if (sens_in->board_mag.status) {
      sens_filt->board_mag[0] = 0.0005 * sens_in->board_mag.meas[0] +
                                0.9995 * sens_filt->board_mag[0];
      sens_filt->board_mag[1] = 0.0005 * sens_in->board_mag.meas[1] +
                                0.9995 * sens_filt->board_mag[1];
      sens_filt->board_mag[2] = 0.0005 * sens_in->board_mag.meas[2] +
                                0.9995 * sens_filt->board_mag[2];
    }
    b_filtered = sens_filt->mti_baro;
    if (sens_in->mti_baro.status) {
      b_filtered =
          0.0005 * sens_in->mti_baro.meas + 0.9995 * sens_filt->mti_baro;
    }
    sens_filt->mti_baro = b_filtered;
    if (sens_in->mti_mag.status) {
      sens_filt->mti_mag[0] =
          0.0005 * sens_in->mti_mag.meas[0] + 0.9995 * sens_filt->mti_mag[0];
      sens_filt->mti_mag[1] =
          0.0005 * sens_in->mti_mag.meas[1] + 0.9995 * sens_filt->mti_mag[1];
      sens_filt->mti_mag[2] =
          0.0005 * sens_in->mti_mag.meas[2] + 0.9995 * sens_filt->mti_mag[2];
    }
    a[0] = (sens_filt->board_accel[0] * (double)sens_in->board_accel.status +
            sens_filt->mti_accel[0] * (double)sens_in->mti_accel.status) +
           sens_filt->ad_accel[0] * (double)sens_in->ad_accel.status;
    a[1] = (sens_filt->board_accel[1] * (double)sens_in->board_accel.status +
            sens_filt->mti_accel[1] * (double)sens_in->mti_accel.status) +
           sens_filt->ad_accel[1] * (double)sens_in->ad_accel.status;
    a[2] = (sens_filt->board_accel[2] * (double)sens_in->board_accel.status +
            sens_filt->mti_accel[2] * (double)sens_in->mti_accel.status) +
           sens_filt->ad_accel[2] * (double)sens_in->ad_accel.status;
    *w_status_nav = true;
    a_norm = b_norm(a);
    if (a_norm < 1.0E-6) {
      q[0] = 1.0;
      q[1] = 0.0;
      q[2] = 0.0;
      q[3] = 0.0;
      *w_status_nav = false;
    } else {
      double d1;
      double qw;
      double qy;
      double qz;
      qw = sqrt(0.5 * (a[0] / a_norm) + 0.5);
      if (qw == 0.0) {
        qy = 1.0;
        qz = 0.0;
      } else {
        qy = 0.5 * (a[2] / a_norm) / qw;
        qz = -0.5 * (a[1] / a_norm) / qw;
      }
      q[0] = qw;
      q[1] = 0.0;
      q[2] = qy;
      q[3] = qz;
      d1 = c_norm(q);
      q[0] = qw / d1;
      q[1] = 0.0 / d1;
      q[2] = qy / d1;
      q[3] = qz / d1;
    }
    x[0] = q[0];
    x[1] = q[1];
    x[2] = q[2];
    x[3] = q[3];
    x[10] = b_SD->pd->b_param.altitude_initial;
    x[4] = 0.0;
    x[7] = 0.0;
    bias->board_gyro[0] =
        sens_filt->board_gyro[0] * (double)sens_in->board_gyro.status;
    bias->mti_gyro[0] =
        sens_filt->mti_gyro[0] * (double)sens_in->mti_gyro.status;
    bias->ad_gyro[0] = sens_filt->ad_gyro[0] * (double)sens_in->ad_gyro.status;
    x[5] = 0.0;
    x[8] = 0.0;
    bias->board_gyro[1] =
        sens_filt->board_gyro[1] * (double)sens_in->board_gyro.status;
    bias->mti_gyro[1] =
        sens_filt->mti_gyro[1] * (double)sens_in->mti_gyro.status;
    bias->ad_gyro[1] = sens_filt->ad_gyro[1] * (double)sens_in->ad_gyro.status;
    x[6] = 0.0;
    x[9] = 0.0;
    bias->board_gyro[2] =
        sens_filt->board_gyro[2] * (double)sens_in->board_gyro.status;
    bias->mti_gyro[2] =
        sens_filt->mti_gyro[2] * (double)sens_in->mti_gyro.status;
    bias->ad_gyro[2] = sens_filt->ad_gyro[2] * (double)sens_in->ad_gyro.status;
    c_a = q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + q[3] * q[3]);
    d_a = 2.0 * q[0];
    for (i3 = 0; i3 < 3; i3++) {
      double a_tmp;
      a_tmp = 2.0 * q[i3 + 1];
      e_a[3 * i3] = c_a * b[i3] + a_tmp * q[1];
      e_a[3 * i3 + 1] = c_a * b[i3 + 3] + a_tmp * q[2];
      e_a[3 * i3 + 2] = c_a * b[i3 + 6] + a_tmp * q[3];
    }
    b_dv[0] = 0.0;
    b_dv[1] = d_a * -q[3];
    b_dv[2] = d_a * q[2];
    b_dv[3] = d_a * q[3];
    b_dv[4] = 0.0;
    b_dv[5] = d_a * -q[1];
    b_dv[6] = d_a * -q[2];
    b_dv[7] = d_a * q[1];
    b_dv[8] = 0.0;
    for (i7 = 0; i7 < 9; i7++) {
      ST[i7] = e_a[i7] - b_dv[i7];
    }
    memset(&b_ST[0], 0, 3U * sizeof(double));
    d19 = b_ST[0];
    d20 = b_ST[1];
    d21 = b_ST[2];
    for (i9 = 0; i9 < 3; i9++) {
      double d23;
      d23 = sens_filt->board_mag[i9];
      d19 += ST[3 * i9] * d23;
      d20 += ST[3 * i9 + 1] * d23;
      d21 += ST[3 * i9 + 2] * d23;
    }
    bias->board_mag_earth[0] = d19 * (double)sens_in->board_mag.status;
    bias->board_mag_earth[1] = d20 * (double)sens_in->board_mag.status;
    bias->board_mag_earth[2] = d21 * (double)sens_in->board_mag.status;
    memset(&b_ST[0], 0, 3U * sizeof(double));
    d24 = b_ST[0];
    d25 = b_ST[1];
    d26 = b_ST[2];
    for (i11 = 0; i11 < 3; i11++) {
      double d27;
      d27 = sens_filt->mti_mag[i11];
      d24 += ST[3 * i11] * d27;
      d25 += ST[3 * i11 + 1] * d27;
      d26 += ST[3 * i11 + 2] * d27;
    }
    double t1_pressure;
    bias->mti_mag_earth[0] = d24 * (double)sens_in->mti_mag.status;
    bias->mti_mag_earth[1] = d25 * (double)sens_in->mti_mag.status;
    bias->mti_mag_earth[2] = d26 * (double)sens_in->mti_mag.status;
    t1_pressure =
        airdata_atmos(b_SD->pd->b_param.altitude_initial, &i_expl_temp,
                      &t1_density, &j_expl_temp, &k_expl_temp, &l_expl_temp);
    bias->board_baro =
        (filtered - t1_pressure) * (double)sens_in->board_baro.status;
    bias->mti_baro =
        (b_filtered - t1_pressure) * (double)sens_in->mti_baro.status;
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
    double c_C_total_a_tmp_tmp;
    double w_idx_0;
    double w_idx_1;
    double w_idx_2;
    int i;
    int k;
    boolean_T exitg1;
    boolean_T status_fast;
    boolean_T y;
    status_fast = false;
    a[0] = 0.0;
    w_idx_0 = 0.0;
    i = 1000000 * sens_in->ad_gyro.status;
    C_total_a_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_in->board_accel.status;
    b_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_in->mti_accel.status;
    c_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_in->ad_accel.status;
    C_total_a_tmp =
        (C_total_a_tmp_tmp + b_C_total_a_tmp_tmp) + c_C_total_a_tmp_tmp;
    C_total_a[0] = C_total_a_tmp;
    C_total_w_tmp_tmp =
        9.9999999999999981E+9 * (double)sens_in->board_gyro.status;
    b_C_total_w_tmp_tmp =
        9.9999999999999981E+9 * (double)sens_in->mti_gyro.status;
    C_total_w_tmp = C_total_w_tmp_tmp + b_C_total_w_tmp_tmp;
    C_total_w[0] = C_total_w_tmp + (double)i;
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
    while ((!exitg1) && (k < 3)) {
      if (C_total_a[k] == 0.0) {
        y = true;
        exitg1 = true;
      } else {
        k++;
      }
    }
    if (!y) {
      int b_k;
      boolean_T b_y;
      b_y = false;
      b_k = 0;
      exitg1 = false;
      while ((!exitg1) && (b_k < 3)) {
        if (C_total_w[b_k] == 0.0) {
          b_y = true;
          exitg1 = true;
        } else {
          b_k++;
        }
      }
      if (!b_y) {
        double w_tilde[9];
        double b_w_idx_1_tmp;
        double d;
        double d10;
        double d11;
        double d12;
        double d3;
        double d4;
        double d5;
        double d6;
        double d7;
        double d9;
        double w_idx_1_tmp;
        status_fast = true;
        w_idx_0 = (C_total_w_tmp_tmp / C_total_w[0] *
                       (sens_in->board_gyro.meas[0] - bias->board_gyro[0]) +
                   b_C_total_w_tmp_tmp / C_total_w[0] *
                       (sens_in->mti_gyro.meas[0] - bias->mti_gyro[0])) +
                  (double)i / C_total_w[0] *
                      (sens_in->ad_gyro.meas[0] - bias->ad_gyro[0]);
        d = 0.0 / C_total_w_tmp;
        w_idx_1_tmp = C_total_w_tmp_tmp / C_total_w_tmp;
        b_w_idx_1_tmp = b_C_total_w_tmp_tmp / C_total_w_tmp;
        w_idx_1 =
            (w_idx_1_tmp * (sens_in->board_gyro.meas[1] - bias->board_gyro[1]) +
             b_w_idx_1_tmp * (sens_in->mti_gyro.meas[1] - bias->mti_gyro[1])) +
            d * (sens_in->ad_gyro.meas[1] - bias->ad_gyro[1]);
        w_idx_2 =
            (w_idx_1_tmp * (sens_in->board_gyro.meas[2] - bias->board_gyro[2]) +
             b_w_idx_1_tmp * (sens_in->mti_gyro.meas[2] - bias->mti_gyro[2])) +
            d * (sens_in->ad_gyro.meas[2] - bias->ad_gyro[2]);
        w_tilde[0] = 0.0;
        w_tilde[3] = -w_idx_2;
        w_tilde[6] = w_idx_1;
        w_tilde[1] = w_idx_2;
        w_tilde[4] = 0.0;
        w_tilde[7] = -w_idx_0;
        w_tilde[2] = -w_idx_1;
        w_tilde[5] = w_idx_0;
        w_tilde[8] = 0.0;
        memset(&w_tilde_sq[0], 0, 9U * sizeof(double));
        for (i1 = 0; i1 < 3; i1++) {
          double d2;
          int b_w_tilde_sq_tmp;
          int w_tilde_sq_tmp;
          d2 = w_tilde_sq[3 * i1];
          w_tilde_sq_tmp = 3 * i1 + 1;
          b_w_tilde_sq_tmp = 3 * i1 + 2;
          for (i2 = 0; i2 < 3; i2++) {
            double d8;
            d8 = w_tilde[i2 + 3 * i1];
            d2 += w_tilde[3 * i2] * d8;
            w_tilde_sq[w_tilde_sq_tmp] += w_tilde[3 * i2 + 1] * d8;
            w_tilde_sq[b_w_tilde_sq_tmp] += w_tilde[3 * i2 + 2] * d8;
          }
          w_tilde_sq[3 * i1] = d2;
        }
        d3 = b_SD->pd->c_param.d_board[0];
        d4 = b_SD->pd->c_param.d_mti[0];
        d5 = b_SD->pd->c_param.d_ad[0];
        d6 = b_SD->pd->c_param.d_board[1];
        d7 = b_SD->pd->c_param.d_mti[1];
        d9 = b_SD->pd->c_param.d_ad[1];
        d10 = b_SD->pd->c_param.d_board[2];
        d11 = b_SD->pd->c_param.d_mti[2];
        d12 = b_SD->pd->c_param.d_ad[2];
        for (i5 = 0; i5 < 3; i5++) {
          double d13;
          double d14;
          double d16;
          double d17;
          double d22;
          d13 = w_tilde_sq[i5];
          d14 = d13 * d3;
          d16 = d13 * d4;
          d17 = d13 * d5;
          d13 = w_tilde_sq[i5 + 3];
          d14 += d13 * d6;
          d16 += d13 * d7;
          d17 += d13 * d9;
          d13 = w_tilde_sq[i5 + 6];
          d14 += d13 * d10;
          d16 += d13 * d11;
          d17 += d13 * d12;
          d22 = C_total_a[i5];
          a[i5] =
              (C_total_a_tmp_tmp / d22 * (sens_in->board_accel.meas[i5] - d14) +
               b_C_total_a_tmp_tmp / d22 *
                   (sens_in->mti_accel.meas[i5] - d16)) +
              c_C_total_a_tmp_tmp / d22 * (sens_in->ad_accel.meas[i5] - d17);
        }
      }
    }
    *w_status_nav = status_fast;
    if (status_fast) {
      double E[121];
      double P_pred[121];
      double K[33];
      double W_dt[16];
      double b_q[16];
      double m_a[16];
      double d_dt[12];
      double x_pred[11];
      double S[9];
      double b_P_pred[9];
      double b_skewed_exp_w_tmp[9];
      double dv4[9];
      double e_a[9];
      double n_tilde[9];
      double skewed_exp_w_tmp[9];
      double w_exp_tilde[9];
      double b_dv1[4];
      double r_q_tmp[4];
      double b_S[3];
      double c_r_q_tmp[3];
      double dn[3];
      double dv2[3];
      double k_x[3];
      double b_b;
      double b_dphi_tmp;
      double b_r_q_tmp;
      double d31;
      double d32;
      double d33;
      double d34;
      double d35;
      double d36;
      double d37;
      double d38;
      double d39;
      double d40;
      double d41;
      double d42;
      double d45;
      double d46;
      double d48;
      double d49;
      double dphi;
      double dphi_tmp;
      double f_a;
      double h_a;
      double h_x;
      double i_a;
      double i_x;
      double j_a;
      double j_x;
      double k_a;
      double l_a;
      double n_a;
      double n_idx_0;
      double n_idx_1;
      double n_idx_2;
      double o_a;
      double p_a;
      double q_mag;
      signed char c_I[16];
      signed char b_I[9];
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
      b_b = sin(dphi);
      n_tilde[0] = 0.0;
      n_tilde[3] = -n_idx_2;
      n_tilde[6] = n_idx_1;
      n_tilde[1] = n_idx_2;
      n_tilde[4] = 0.0;
      n_tilde[7] = -n_idx_0;
      n_tilde[2] = -n_idx_1;
      n_tilde[5] = n_idx_0;
      n_tilde[8] = 0.0;
      f_a = sin(b_dphi_tmp);
      h_x = cos(b_dphi_tmp);
      for (i4 = 0; i4 < 9; i4++) {
        b_I[i4] = 0;
      }
      memset(&b_n_tilde[0], 0, 9U * sizeof(double));
      for (c_k = 0; c_k < 3; c_k++) {
        double d15;
        int b_n_tilde_tmp;
        int n_tilde_tmp;
        b_I[c_k + 3 * c_k] = 1;
        d15 = b_n_tilde[3 * c_k];
        n_tilde_tmp = 3 * c_k + 1;
        b_n_tilde_tmp = 3 * c_k + 2;
        for (i8 = 0; i8 < 3; i8++) {
          double d18;
          d18 = n_tilde[i8 + 3 * c_k];
          d15 += n_tilde[3 * i8] * d18;
          b_n_tilde[n_tilde_tmp] += n_tilde[3 * i8 + 1] * d18;
          b_n_tilde[b_n_tilde_tmp] += n_tilde[3 * i8 + 2] * d18;
        }
        b_n_tilde[3 * c_k] = d15;
      }
      for (i6 = 0; i6 < 9; i6++) {
        w_exp_tilde[i6] =
            ((double)b_I[i6] - f_a * n_tilde[i6]) + (1.0 - h_x) * b_n_tilde[i6];
      }
      double g_a;
      g_a = b_norm(&x[7]);
      airdata_atmos(x[10], &e_expl_temp, &t1_density, &f_expl_temp,
                    &g_expl_temp, &h_expl_temp);
      h_a = 0.5 * t1_density * (g_a * g_a);
      i_a = b_SD->pd->d_param.c_aero * b_SD->pd->d_param.Cn_alpha;
      j_a = x[0] * x[0] - ((x[1] * x[1] + x[2] * x[2]) + x[3] * x[3]);
      k_a = 2.0 * x[0];
      for (i10 = 0; i10 < 3; i10++) {
        double b_a_tmp;
        int c_a_tmp;
        int d_a_tmp;
        b_a_tmp = x[i10 + 1];
        e_a[3 * i10] = j_a * b[3 * i10] + 2.0 * x[1] * b_a_tmp;
        c_a_tmp = 3 * i10 + 1;
        e_a[c_a_tmp] = j_a * b[c_a_tmp] + 2.0 * x[2] * b_a_tmp;
        d_a_tmp = 3 * i10 + 2;
        e_a[d_a_tmp] = j_a * b[d_a_tmp] + 2.0 * x[3] * b_a_tmp;
      }
      b_dv[0] = 0.0;
      b_dv[3] = k_a * -x[3];
      b_dv[6] = k_a * x[2];
      b_dv[1] = k_a * x[3];
      b_dv[4] = 0.0;
      b_dv[7] = k_a * -x[1];
      b_dv[2] = k_a * -x[2];
      b_dv[5] = k_a * x[1];
      b_dv[8] = 0.0;
      for (i12 = 0; i12 < 9; i12++) {
        S[i12] = e_a[i12] - b_dv[i12];
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
      b_dv1[0] = cos(dphi);
      memset(&b_w_exp_tilde[0], 0, 9U * sizeof(double));
      memset(&c_w_exp_tilde[0], 0, 3U * sizeof(double));
      for (i13 = 0; i13 < 3; i13++) {
        double d28;
        int b_w_exp_tilde_tmp;
        int w_exp_tilde_tmp;
        b_dv1[i13 + 1] = dn[i13] * b_b;
        d28 = b_w_exp_tilde[3 * i13];
        w_exp_tilde_tmp = 3 * i13 + 1;
        b_w_exp_tilde_tmp = 3 * i13 + 2;
        for (i14 = 0; i14 < 3; i14++) {
          double d29;
          d29 = b_SD->pd->d_param.J[i14 + 3 * i13];
          d28 += w_exp_tilde[3 * i14] * d29;
          b_w_exp_tilde[w_exp_tilde_tmp] += w_exp_tilde[3 * i14 + 1] * d29;
          b_w_exp_tilde[b_w_exp_tilde_tmp] += w_exp_tilde[3 * i14 + 2] * d29;
        }
        double d30;
        b_w_exp_tilde[3 * i13] = d28;
        d30 = x[i13 + 4];
        c_w_exp_tilde[0] += d28 * d30;
        c_w_exp_tilde[1] += b_w_exp_tilde[3 * i13 + 1] * d30;
        c_w_exp_tilde[2] += b_w_exp_tilde[3 * i13 + 2] * d30;
      }
      dv2[0] = 0.0;
      dv2[1] = h_a * (i_a * sin(atan2(x[9], x[7])));
      dv2[2] = h_a * (i_a * -sin(atan2(x[8], x[7])));
      memset(&dv3[0], 0, 3U * sizeof(double));
      memset(&b_dt[0], 0, 3U * sizeof(double));
      memset(&c_dt[0], 0, 3U * sizeof(double));
      d31 = dv3[0];
      d32 = dv3[1];
      d33 = dv3[2];
      d34 = b_dt[0];
      d35 = b_dt[1];
      d36 = b_dt[2];
      d37 = x[7];
      d38 = x[8];
      d39 = x[9];
      d40 = c_dt[0];
      d41 = c_dt[1];
      d42 = c_dt[2];
      for (i15 = 0; i15 < 3; i15++) {
        double d43;
        double d44;
        double d47;
        double d50;
        double d51;
        double d52;
        int i17;
        int i18;
        d43 = b_SD->pd->d_param.Jinv[3 * i15];
        d44 = c_w_exp_tilde[i15];
        d31 += d43 * d44;
        d47 = dv2[i15];
        d34 += dt * d43 * d47;
        d50 = S[3 * i15];
        d51 = b_SD->pd->d_param.g[i15];
        d40 += dt * d50 * d51;
        d52 = d50 * d37;
        i17 = 3 * i15 + 1;
        d43 = b_SD->pd->d_param.Jinv[i17];
        d32 += d43 * d44;
        d35 += dt * d43 * d47;
        d50 = S[i17];
        d41 += dt * d50 * d51;
        d52 += d50 * d38;
        i18 = 3 * i15 + 2;
        d43 = b_SD->pd->d_param.Jinv[i18];
        d33 += d43 * d44;
        d36 += dt * d43 * d47;
        d50 = S[i18];
        d42 += dt * d50 * d51;
        d52 += d50 * d39;
        c_w_exp_tilde[i15] =
            ((w_exp_tilde[i15] * d37 + w_exp_tilde[i15 + 3] * d38) +
             w_exp_tilde[i15 + 6] * d39) +
            dt * a[i15];
        b_S[i15] = d52;
      }
      memset(&c_q[0], 0, sizeof(double) << 2);
      d45 = c_q[0];
      d46 = c_q[1];
      d48 = c_q[2];
      d49 = c_q[3];
      for (i16 = 0; i16 < 4; i16++) {
        double d53;
        int q_tmp;
        q_tmp = i16 << 2;
        d53 = b_dv1[i16];
        d45 += b_q[q_tmp] * d53;
        d46 += b_q[q_tmp + 1] * d53;
        d48 += b_q[q_tmp + 2] * d53;
        d49 += b_q[q_tmp + 3] * d53;
      }
      x_pred[0] = d45;
      x_pred[1] = d46;
      x_pred[2] = d48;
      x_pred[3] = d49;
      x_pred[4] = d31 + d34;
      x_pred[7] = c_w_exp_tilde[0] + d40;
      x_pred[5] = d32 + d35;
      x_pred[8] = c_w_exp_tilde[1] + d41;
      x_pred[6] = d33 + d36;
      x_pred[9] = c_w_exp_tilde[2] + d42;
      x_pred[10] = x[10] + dt * b_S[0];
      memset(&F[0], 0, 121U * sizeof(double));
      l_a = 0.5 * dt;
      W_dt[0] = 0.0;
      W_dt[4] = l_a * -x[4];
      W_dt[8] = l_a * -x[5];
      W_dt[12] = l_a * -x[6];
      W_dt[1] = l_a * x[4];
      W_dt[5] = 0.0;
      W_dt[9] = l_a * x[6];
      W_dt[13] = l_a * -x[5];
      W_dt[2] = l_a * x[5];
      W_dt[6] = l_a * -x[6];
      W_dt[10] = 0.0;
      W_dt[14] = l_a * x[4];
      W_dt[3] = l_a * x[6];
      W_dt[7] = l_a * x[5];
      W_dt[11] = l_a * -x[4];
      W_dt[15] = 0.0;
      memset(&c_b[0], 0, sizeof(double) << 4);
      for (i19 = 0; i19 < 4; i19++) {
        int i20;
        i20 = i19 << 2;
        for (i22 = 0; i22 < 4; i22++) {
          double d54;
          int b_tmp;
          d54 = W_dt[i22 + i20];
          b_tmp = i22 << 2;
          c_b[i20] += W_dt[b_tmp] * d54;
          c_b[i20 + 1] += W_dt[b_tmp + 1] * d54;
          c_b[i20 + 2] += W_dt[b_tmp + 2] * d54;
          c_b[i20 + 3] += W_dt[b_tmp + 3] * d54;
        }
      }
      for (i21 = 0; i21 < 16; i21++) {
        c_I[i21] = 0;
      }
      memset(&b_W_dt[0], 0, sizeof(double) << 4);
      memset(&d_b[0], 0, sizeof(double) << 4);
      for (d_k = 0; d_k < 4; d_k++) {
        double d55;
        double d56;
        double d57;
        double d58;
        double d59;
        double d60;
        double d61;
        double d62;
        int I_tmp;
        I_tmp = d_k << 2;
        c_I[d_k + I_tmp] = 1;
        d55 = b_W_dt[I_tmp];
        d56 = b_W_dt[I_tmp + 1];
        d57 = b_W_dt[I_tmp + 2];
        d58 = b_W_dt[I_tmp + 3];
        d59 = d_b[I_tmp];
        d60 = d_b[I_tmp + 1];
        d61 = d_b[I_tmp + 2];
        d62 = d_b[I_tmp + 3];
        for (i23 = 0; i23 < 4; i23++) {
          double d63;
          int W_dt_tmp;
          d63 = c_b[i23 + I_tmp];
          W_dt_tmp = i23 << 2;
          d55 += W_dt[W_dt_tmp] * d63;
          d59 += c_b[W_dt_tmp] * d63;
          d56 += W_dt[W_dt_tmp + 1] * d63;
          d60 += c_b[W_dt_tmp + 1] * d63;
          d57 += W_dt[W_dt_tmp + 2] * d63;
          d61 += c_b[W_dt_tmp + 2] * d63;
          d58 += W_dt[W_dt_tmp + 3] * d63;
          d62 += c_b[W_dt_tmp + 3] * d63;
        }
        d_b[I_tmp + 3] = d62;
        d_b[I_tmp + 2] = d61;
        d_b[I_tmp + 1] = d60;
        d_b[I_tmp] = d59;
        b_W_dt[I_tmp + 3] = d58;
        b_W_dt[I_tmp + 2] = d57;
        b_W_dt[I_tmp + 1] = d56;
        b_W_dt[I_tmp] = d55;
        F[11 * d_k] = ((((double)c_I[I_tmp] + W_dt[I_tmp]) + 0.5 * c_b[I_tmp]) +
                       0.16666666666666666 * d55) +
                      0.041666666666666664 * d59;
        F[11 * d_k + 1] = ((((double)c_I[I_tmp + 1] + W_dt[I_tmp + 1]) +
                            0.5 * c_b[I_tmp + 1]) +
                           0.16666666666666666 * d56) +
                          0.041666666666666664 * d60;
        F[11 * d_k + 2] = ((((double)c_I[I_tmp + 2] + W_dt[I_tmp + 2]) +
                            0.5 * c_b[I_tmp + 2]) +
                           0.16666666666666666 * d57) +
                          0.041666666666666664 * d61;
        F[11 * d_k + 3] = ((((double)c_I[I_tmp + 3] + W_dt[I_tmp + 3]) +
                            0.5 * c_b[I_tmp + 3]) +
                           0.16666666666666666 * d58) +
                          0.041666666666666664 * d62;
      }
      double e_a_tmp;
      double f_a_tmp;
      double g_a_tmp;
      double h_a_tmp;
      double i_a_tmp;
      double j_a_tmp;
      double k_a_tmp;
      e_a_tmp = l_a * q[0];
      m_a[0] = e_a_tmp;
      f_a_tmp = l_a * -q[1];
      m_a[4] = f_a_tmp;
      g_a_tmp = l_a * -q[2];
      m_a[8] = g_a_tmp;
      h_a_tmp = l_a * -q[3];
      m_a[12] = h_a_tmp;
      i_a_tmp = l_a * q[1];
      m_a[1] = i_a_tmp;
      m_a[5] = e_a_tmp;
      m_a[9] = h_a_tmp;
      j_a_tmp = l_a * q[2];
      m_a[13] = j_a_tmp;
      m_a[2] = j_a_tmp;
      k_a_tmp = l_a * q[3];
      m_a[6] = k_a_tmp;
      m_a[10] = e_a_tmp;
      m_a[14] = f_a_tmp;
      m_a[3] = k_a_tmp;
      m_a[7] = g_a_tmp;
      m_a[11] = i_a_tmp;
      m_a[15] = e_a_tmp;
      for (i24 = 0; i24 < 3; i24++) {
        int F_tmp;
        int b_F_tmp;
        F_tmp = (i24 + 1) << 2;
        b_F_tmp = 11 * (i24 + 4);
        F[b_F_tmp] = m_a[F_tmp];
        F[b_F_tmp + 1] = m_a[F_tmp + 1];
        F[b_F_tmp + 2] = m_a[F_tmp + 2];
        F[b_F_tmp + 3] = m_a[F_tmp + 3];
      }
      n_a = 0.5 * b_SD->pd->e_param.c_aero * b_SD->pd->e_param.Cn_alpha;
      airdata_atmos(x[10], &m_expl_temp, &t1_density, &n_expl_temp,
                    &o_expl_temp, &p_expl_temp);
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
      for (i25 = 0; i25 < 9; i25++) {
        b_I[i25] = 0;
      }
      memset(&b_n_tilde[0], 0, 9U * sizeof(double));
      for (e_k = 0; e_k < 3; e_k++) {
        double d64;
        int c_n_tilde_tmp;
        int d_n_tilde_tmp;
        b_I[e_k + 3 * e_k] = 1;
        d64 = b_n_tilde[3 * e_k];
        c_n_tilde_tmp = 3 * e_k + 1;
        d_n_tilde_tmp = 3 * e_k + 2;
        for (i27 = 0; i27 < 3; i27++) {
          double d65;
          d65 = n_tilde[i27 + 3 * e_k];
          d64 += n_tilde[3 * i27] * d65;
          b_n_tilde[c_n_tilde_tmp] += n_tilde[3 * i27 + 1] * d65;
          b_n_tilde[d_n_tilde_tmp] += n_tilde[3 * i27 + 2] * d65;
        }
        b_n_tilde[3 * e_k] = d64;
      }
      for (i26 = 0; i26 < 9; i26++) {
        w_exp_tilde[i26] = ((double)b_I[i26] - f_a * n_tilde[i26]) +
                           (1.0 - h_x) * b_n_tilde[i26];
      }
      memset(&b_dv[0], 0, 9U * sizeof(double));
      for (i28 = 0; i28 < 3; i28++) {
        double d66;
        int i30;
        int i31;
        d66 = b_dv[3 * i28];
        i30 = 3 * i28 + 1;
        i31 = 3 * i28 + 2;
        for (i33 = 0; i33 < 3; i33++) {
          double d68;
          d68 = w_exp_tilde[i33 + 3 * i28];
          d66 += b_SD->pd->e_param.Jinv[3 * i33] * d68;
          b_dv[i30] += b_SD->pd->e_param.Jinv[3 * i33 + 1] * d68;
          b_dv[i31] += b_SD->pd->e_param.Jinv[3 * i33 + 2] * d68;
          F[(i33 + 11 * (i28 + 4)) + 4] = 0.0;
        }
        b_dv[3 * i28] = d66;
      }
      for (i29 = 0; i29 < 3; i29++) {
        int F_tmp_tmp;
        F_tmp_tmp = 11 * (i29 + 4);
        for (i32 = 0; i32 < 3; i32++) {
          double d67;
          d67 = b_SD->pd->e_param.J[i32 + 3 * i29];
          F[F_tmp_tmp + 4] += b_dv[3 * i32] * d67;
          F[F_tmp_tmp + 5] += b_dv[3 * i32 + 1] * d67;
          F[F_tmp_tmp + 6] += b_dv[3 * i32 + 2] * d67;
        }
      }
      b_dv[1] = t1_density * (n_a * x[9]);
      b_dv[4] = 0.0;
      b_dv[7] = t1_density * (n_a * x[7]);
      b_dv[2] = t1_density * (n_a * -x[8]);
      b_dv[5] = t1_density * (n_a * -x[7]);
      b_dv[8] = 0.0;
      i_x = 0.0;
      for (i34 = 0; i34 < 3; i34++) {
        double d69;
        double d70;
        double d71;
        int c_F_tmp;
        b_dv[3 * i34] = 0.0;
        c_F_tmp = 11 * (i34 + 7);
        d69 = 0.0;
        d70 = 0.0;
        d71 = 0.0;
        for (i35 = 0; i35 < 3; i35++) {
          double d72;
          d72 = b_dv[i35 + 3 * i34];
          d69 += dt * b_SD->pd->e_param.Jinv[3 * i35] * d72;
          d70 += dt * b_SD->pd->e_param.Jinv[3 * i35 + 1] * d72;
          d71 += dt * b_SD->pd->e_param.Jinv[3 * i35 + 2] * d72;
        }
        F[c_F_tmp + 6] = d71;
        F[c_F_tmp + 5] = d70;
        F[c_F_tmp + 4] = d69;
        i_x += x[i34 + 1] * b_SD->pd->e_param.g[i34];
      }
      j_x = x[0];
      k_x[0] = x[2] * b_SD->pd->e_param.g[2] - b_SD->pd->e_param.g[1] * x[3];
      k_x[1] = b_SD->pd->e_param.g[0] * x[3] - x[1] * b_SD->pd->e_param.g[2];
      k_x[2] = x[1] * b_SD->pd->e_param.g[1] - b_SD->pd->e_param.g[0] * x[2];
      dv4[0] = 0.0;
      dv4[3] = x[0] * -b_SD->pd->e_param.g[2];
      dv4[6] = x[0] * b_SD->pd->e_param.g[1];
      dv4[1] = x[0] * b_SD->pd->e_param.g[2];
      dv4[4] = 0.0;
      dv4[7] = x[0] * -b_SD->pd->e_param.g[0];
      dv4[2] = x[0] * -b_SD->pd->e_param.g[1];
      dv4[5] = x[0] * b_SD->pd->e_param.g[0];
      dv4[8] = 0.0;
      skewed_exp_w_tmp[0] = 0.0;
      skewed_exp_w_tmp[3] = -x[9];
      skewed_exp_w_tmp[6] = x[8];
      skewed_exp_w_tmp[1] = x[9];
      skewed_exp_w_tmp[4] = 0.0;
      skewed_exp_w_tmp[7] = -x[7];
      skewed_exp_w_tmp[2] = -x[8];
      skewed_exp_w_tmp[5] = x[7];
      skewed_exp_w_tmp[8] = 0.0;
      b_skewed_exp_w_tmp[0] = 0.0;
      b_skewed_exp_w_tmp[3] = -x[6];
      b_skewed_exp_w_tmp[6] = x[5];
      b_skewed_exp_w_tmp[1] = x[6];
      b_skewed_exp_w_tmp[4] = 0.0;
      b_skewed_exp_w_tmp[7] = -x[4];
      b_skewed_exp_w_tmp[2] = -x[5];
      b_skewed_exp_w_tmp[5] = x[4];
      b_skewed_exp_w_tmp[8] = 0.0;
      o_a = 0.5 * (dt * dt);
      memset(&c_skewed_exp_w_tmp[0], 0, 9U * sizeof(double));
      memset(&b_dv[0], 0, 9U * sizeof(double));
      r_q_tmp[0] = x[0];
      b_r_q_tmp = 0.0;
      for (i36 = 0; i36 < 3; i36++) {
        double d73;
        double d74;
        double d_F_tmp;
        int e_F_tmp;
        int f_F_tmp;
        int g_F_tmp;
        F[i36 + 7] = dt * (2.0 * (j_x * b_SD->pd->e_param.g[i36] - k_x[i36]));
        d_F_tmp = x[i36 + 1];
        e_F_tmp = 11 * (i36 + 1);
        F[e_F_tmp + 7] =
            dt * (2.0 * (((i_x * b[3 * i36] + x[1] * b_SD->pd->e_param.g[i36]) -
                          b_SD->pd->e_param.g[0] * d_F_tmp) +
                         dv4[3 * i36]));
        f_F_tmp = 3 * i36 + 1;
        F[e_F_tmp + 8] =
            dt * (2.0 * (((i_x * b[f_F_tmp] + x[2] * b_SD->pd->e_param.g[i36]) -
                          b_SD->pd->e_param.g[1] * d_F_tmp) +
                         dv4[f_F_tmp]));
        g_F_tmp = 3 * i36 + 2;
        F[e_F_tmp + 9] =
            dt * (2.0 * (((i_x * b[g_F_tmp] + x[3] * b_SD->pd->e_param.g[i36]) -
                          b_SD->pd->e_param.g[2] * d_F_tmp) +
                         dv4[g_F_tmp]));
        d73 = c_skewed_exp_w_tmp[3 * i36];
        d74 = b_dv[3 * i36];
        for (i38 = 0; i38 < 3; i38++) {
          double d75;
          double d76;
          int b_skewed_exp_w_tmp_tmp;
          int i39;
          int skewed_exp_w_tmp_tmp;
          i39 = i38 + 3 * i36;
          d75 = b_skewed_exp_w_tmp[i39];
          d76 = skewed_exp_w_tmp[i39];
          d73 += skewed_exp_w_tmp[3 * i38] * d75;
          d74 += 2.0 * b_skewed_exp_w_tmp[3 * i38] * d76;
          skewed_exp_w_tmp_tmp = 3 * i38 + 1;
          c_skewed_exp_w_tmp[f_F_tmp] +=
              skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d75;
          b_dv[f_F_tmp] += 2.0 * b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d76;
          b_skewed_exp_w_tmp_tmp = 3 * i38 + 2;
          c_skewed_exp_w_tmp[g_F_tmp] +=
              skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d75;
          b_dv[g_F_tmp] +=
              2.0 * b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d76;
        }
        int h_F_tmp;
        int i_F_tmp;
        b_dv[3 * i36] = d74;
        c_skewed_exp_w_tmp[3 * i36] = d73;
        h_F_tmp = 11 * (i36 + 4);
        F[h_F_tmp + 7] = dt * skewed_exp_w_tmp[3 * i36] + o_a * (d73 - d74);
        i_F_tmp = 11 * (i36 + 7);
        F[i_F_tmp + 7] = w_exp_tilde[3 * i36];
        F[h_F_tmp + 8] = dt * skewed_exp_w_tmp[f_F_tmp] +
                         o_a * (c_skewed_exp_w_tmp[f_F_tmp] - b_dv[f_F_tmp]);
        F[i_F_tmp + 8] = w_exp_tilde[f_F_tmp];
        F[h_F_tmp + 9] = dt * skewed_exp_w_tmp[g_F_tmp] +
                         o_a * (c_skewed_exp_w_tmp[g_F_tmp] - b_dv[g_F_tmp]);
        F[i_F_tmp + 9] = w_exp_tilde[g_F_tmp];
        r_q_tmp[i36 + 1] = -d_F_tmp;
        b_r_q_tmp += -d_F_tmp * x[i36 + 7];
      }
      c_r_q_tmp[0] = r_q_tmp[2] * x[9] - r_q_tmp[3] * x[8];
      c_r_q_tmp[1] = r_q_tmp[3] * x[7] - r_q_tmp[1] * x[9];
      c_r_q_tmp[2] = r_q_tmp[1] * x[8] - r_q_tmp[2] * x[7];
      for (i37 = 0; i37 < 3; i37++) {
        double b_dt_tmp;
        double dt_tmp;
        int c_dt_tmp;
        int d_dt_tmp;
        int e_dt_tmp;
        dt_tmp = x[i37 + 7];
        d_dt[i37] = dt * (2.0 * (r_q_tmp[0] * dt_tmp - c_r_q_tmp[i37]));
        b_dt_tmp = r_q_tmp[i37 + 1];
        c_dt_tmp = 3 * (i37 + 1);
        d_dt[c_dt_tmp] =
            dt * (2.0 * (((b_r_q_tmp * b[3 * i37] + r_q_tmp[1] * dt_tmp) -
                          x[7] * b_dt_tmp) +
                         r_q_tmp[0] * skewed_exp_w_tmp[3 * i37]));
        d_dt_tmp = 3 * i37 + 1;
        d_dt[c_dt_tmp + 1] =
            dt * (2.0 * (((b_r_q_tmp * b[d_dt_tmp] + r_q_tmp[2] * dt_tmp) -
                          x[8] * b_dt_tmp) +
                         r_q_tmp[0] * skewed_exp_w_tmp[d_dt_tmp]));
        e_dt_tmp = 3 * i37 + 2;
        d_dt[c_dt_tmp + 2] =
            dt * (2.0 * (((b_r_q_tmp * b[e_dt_tmp] + r_q_tmp[3] * dt_tmp) -
                          x[9] * b_dt_tmp) +
                         r_q_tmp[0] * skewed_exp_w_tmp[e_dt_tmp]));
      }
      double q_a;
      F[10] = d_dt[0];
      F[21] = d_dt[3];
      F[32] = d_dt[6];
      F[43] = d_dt[9];
      p_a = r_q_tmp[0] * r_q_tmp[0] -
            ((r_q_tmp[1] * r_q_tmp[1] + r_q_tmp[2] * r_q_tmp[2]) +
             r_q_tmp[3] * r_q_tmp[3]);
      q_a = 2.0 * r_q_tmp[0];
      b_dv[0] = 0.0;
      b_dv[3] = q_a * -r_q_tmp[3];
      b_dv[6] = q_a * r_q_tmp[2];
      b_dv[1] = q_a * r_q_tmp[3];
      b_dv[4] = 0.0;
      b_dv[7] = q_a * -r_q_tmp[1];
      b_dv[2] = q_a * -r_q_tmp[2];
      b_dv[5] = q_a * r_q_tmp[1];
      b_dv[8] = 0.0;
      for (i40 = 0; i40 < 3; i40++) {
        F[11 * (i40 + 7) + 10] =
            dt * ((p_a * b[3 * i40] + 2.0 * r_q_tmp[1] * r_q_tmp[i40 + 1]) -
                  b_dv[3 * i40]);
      }
      F[120] = 1.0;
      memset(&b_F[0], 0, 121U * sizeof(double));
      for (i41 = 0; i41 < 11; i41++) {
        for (i42 = 0; i42 < 11; i42++) {
          double d77;
          d77 = P[i42 + 11 * i41];
          for (i45 = 0; i45 < 11; i45++) {
            int j_F_tmp;
            j_F_tmp = i45 + 11 * i41;
            b_F[j_F_tmp] += F[i45 + 11 * i42] * d77;
          }
        }
      }
      for (i43 = 0; i43 < 11; i43++) {
        for (i44 = 0; i44 < 11; i44++) {
          double d78;
          d78 = 0.0;
          for (i47 = 0; i47 < 11; i47++) {
            d78 += b_F[i43 + 11 * i47] * F[i44 + 11 * i47];
          }
          int c_P_pred_tmp;
          c_P_pred_tmp = i43 + 11 * i44;
          P_pred[c_P_pred_tmp] = d78 + Q[c_P_pred_tmp];
        }
      }
      for (i46 = 0; i46 < 3; i46++) {
        int P_pred_tmp;
        int b_P_pred_tmp;
        int d_P_pred_tmp;
        P_pred_tmp = 11 * (i46 + 4);
        b_P_pred[3 * i46] = P_pred[P_pred_tmp + 4] + R[3 * i46];
        b_P_pred_tmp = 3 * i46 + 1;
        b_P_pred[b_P_pred_tmp] = P_pred[P_pred_tmp + 5] + R[b_P_pred_tmp];
        d_P_pred_tmp = 3 * i46 + 2;
        b_P_pred[d_P_pred_tmp] = P_pred[P_pred_tmp + 6] + R[d_P_pred_tmp];
      }
      mrdiv(&P_pred[44], b_P_pred, K);
      memset(&d_I[0], 0, 121U * sizeof(signed char));
      for (f_k = 0; f_k < 11; f_k++) {
        d_I[f_k + 11 * f_k] = 1;
      }
      for (i48 = 0; i48 < 44; i48++) {
        E[i48] = d_I[i48];
      }
      for (i49 = 0; i49 < 33; i49++) {
        E[i49 + 44] = (double)d_I[i49 + 44] - K[i49];
      }
      for (i50 = 0; i50 < 44; i50++) {
        E[i50 + 77] = d_I[i50 + 77];
      }
      memset(&b_E[0], 0, 121U * sizeof(double));
      for (i51 = 0; i51 < 11; i51++) {
        for (i52 = 0; i52 < 11; i52++) {
          double d79;
          d79 = P_pred[i52 + 11 * i51];
          for (i54 = 0; i54 < 11; i54++) {
            int E_tmp;
            E_tmp = i54 + 11 * i51;
            b_E[E_tmp] += E[i54 + 11 * i52] * d79;
          }
        }
      }
      memset(&b_K[0], 0, 33U * sizeof(double));
      for (i53 = 0; i53 < 3; i53++) {
        for (i55 = 0; i55 < 3; i55++) {
          double d80;
          d80 = R[i55 + 3 * i53];
          for (i56 = 0; i56 < 11; i56++) {
            int K_tmp;
            K_tmp = i56 + 11 * i53;
            b_K[K_tmp] += K[i56 + 11 * i55] * d80;
          }
        }
      }
      memset(&c_E[0], 0, 121U * sizeof(double));
      memset(&c_K[0], 0, 121U * sizeof(double));
      for (i57 = 0; i57 < 11; i57++) {
        for (i59 = 0; i59 < 11; i59++) {
          double d81;
          d81 = E[i57 + 11 * i59];
          for (i61 = 0; i61 < 11; i61++) {
            int b_E_tmp;
            b_E_tmp = i61 + 11 * i57;
            c_E[b_E_tmp] += b_E[i61 + 11 * i59] * d81;
          }
        }
        for (i60 = 0; i60 < 3; i60++) {
          double d82;
          d82 = K[i57 + 11 * i60];
          for (i63 = 0; i63 < 11; i63++) {
            int b_K_tmp;
            b_K_tmp = i63 + 11 * i57;
            c_K[b_K_tmp] += b_K[i63 + 11 * i60] * d82;
          }
        }
      }
      for (i58 = 0; i58 < 121; i58++) {
        P[i58] = c_E[i58] + c_K[i58];
      }
      w_idx_0 -= x_pred[4];
      w_idx_1 -= x_pred[5];
      w_idx_2 -= x_pred[6];
      for (i62 = 0; i62 < 11; i62++) {
        x[i62] = x_pred[i62] + ((K[i62] * w_idx_0 + K[i62 + 11] * w_idx_1) +
                                K[i62 + 22] * w_idx_2);
      }
      double b_q_mag;
      b_q_mag = c_norm(&x[0]);
      x[0] /= b_q_mag;
      x[1] /= b_q_mag;
      x[2] /= b_q_mag;
      x[3] /= b_q_mag;
    }
    if (sens_in->board_baro.status) {
      memcpy(&d_x[0], &x[0], 11U * sizeof(double));
      memcpy(&b_P[0], &P[0], 121U * sizeof(double));
      b_ekf_correct(d_x, b_P, sens_in->board_baro.meas, bias->board_baro, x, P);
    }
    if (sens_in->mti_baro.status) {
      memcpy(&e_x[0], &x[0], 11U * sizeof(double));
      memcpy(&c_P[0], &P[0], 121U * sizeof(double));
      b_ekf_correct(e_x, c_P, sens_in->mti_baro.meas, bias->mti_baro, x, P);
    }
    if (sens_in->board_mag.status) {
      memcpy(&f_x[0], &x[0], 11U * sizeof(double));
      memcpy(&d_P[0], &P[0], 121U * sizeof(double));
      ekf_correct(f_x, d_P, sens_in->board_mag.meas, bias->board_mag_earth, b,
                  x, P);
    }
    if (sens_in->mti_mag.status) {
      memcpy(&g_x[0], &x[0], 11U * sizeof(double));
      memcpy(&e_P[0], &P[0], 121U * sizeof(double));
      ekf_correct(g_x, e_P, sens_in->mti_mag.meas, bias->mti_mag_earth, b, x,
                  P);
    }
  }
  b_a = b_norm(&x[7]);
  airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp,
                &d_expl_temp);
  *pdyn = 0.5 * t1_density * (b_a * b_a);
  *cov_norm = 0.0;
  for (b_i = 0; b_i < 11; b_i++) {
    double s;
    s = 0.0;
    for (j = 0; j < 11; j++) {
      s += fabs(P[b_i + 11 * j]);
    }
    if (s > *cov_norm) {
      *cov_norm = s;
    }
  }
  roll_state[0] =
      atan2(2.0 * (x[2] * x[3] + x[0] * x[1]),
            ((x[0] * x[0] - x[1] * x[1]) - x[2] * x[2]) + x[3] * x[3]);
  roll_state[1] = x[4];
}
