#include "GNC_codegen_SIL.h"
#include "GNC_codegen_SIL_types.h"
#include <math.h>
#include <string.h>

static const double dv[9] = {0.46, 0.0, 0.0, 0.0, 49.5, 0.0, 0.0, 0.0, 49.5};

static const double dv1[9] = {2.1739130434782608,   0.0, 0.0, 0.0,
                              0.020202020202020204, 0.0, 0.0, 0.0,
                              0.020202020202020204};

static double airdata_atmos(double altitude, double *airdata_temperature,
                            double *airdata_density,
                            double *airdata_sonic_speed, double *airdata_mach,
                            double *airdata_dynamic_pressure);

static void b_ekf_correct(const double x[11], const double P[121], double y,
                          double b, double x_new[11], double P_new[121]);

static double b_norm(const double x[3]);

static void controller_codegen_entry_init(GNC_codegen_SILStackData *b_SD);

static void dynamics_init(GNC_codegen_SILStackData *b_SD);

static void dynamics_jacobian_init(GNC_codegen_SILStackData *b_SD);

static void ekf_correct(const double x[11], const double P[121],
                        const double y[3], const double b[3], const double R[9],
                        double x_new[11], double P_new[121]);

static void inv(const double x[9], double y[9]);

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
  double b_P[11];
  double absxk;
  double airdata_altitude_pressure;
  double altitude;
  double altitude_ratio;
  double b_b;
  double b_expl_temp;
  double c_H;
  double c_b;
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
  int i13;
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
    for (i1 = 0; i1 < 11; i1++) {
      d += H[i1] * P[i1 + 11 * i];
    }
    b_H[i] = d;
    c_H += d * H[i];
  }
  c_b = 1.0 / (c_H + 100.0);
  memset(&b_P[0], 0, 11U * sizeof(double));
  for (i2 = 0; i2 < 11; i2++) {
    for (i3 = 0; i3 < 11; i3++) {
      b_P[i3] += P[i3 + 11 * i2] * H[i2];
    }
  }
  for (i4 = 0; i4 < 11; i4++) {
    K[i4] = b_P[i4] * c_b;
  }
  memset(&b_I[0], 0, 121U * sizeof(signed char));
  for (k = 0; k < 11; k++) {
    b_I[k + 11 * k] = 1;
  }
  for (i5 = 0; i5 < 11; i5++) {
    for (i6 = 0; i6 < 11; i6++) {
      int E_tmp;
      E_tmp = i6 + 11 * i5;
      E[E_tmp] = (double)b_I[E_tmp] - K[i6] * H[i5];
    }
  }
  memset(&b_E[0], 0, 121U * sizeof(double));
  for (i7 = 0; i7 < 11; i7++) {
    for (i8 = 0; i8 < 11; i8++) {
      double d1;
      d1 = P[i8 + 11 * i7];
      for (i10 = 0; i10 < 11; i10++) {
        int b_E_tmp;
        b_E_tmp = i10 + 11 * i7;
        b_E[b_E_tmp] += E[i10 + 11 * i8] * d1;
      }
    }
  }
  memset(&c_E[0], 0, 121U * sizeof(double));
  for (i9 = 0; i9 < 11; i9++) {
    for (i11 = 0; i11 < 11; i11++) {
      double d2;
      d2 = E[i9 + 11 * i11];
      for (i13 = 0; i13 < 11; i13++) {
        int c_E_tmp;
        c_E_tmp = i13 + 11 * i9;
        c_E[c_E_tmp] += b_E[i13 + 11 * i11] * d2;
      }
      b_K[i11 + 11 * i9] = K[i11] * 100.0 * K[i9];
    }
  }
  for (i12 = 0; i12 < 121; i12++) {
    P_new[i12] = c_E[i12] + b_K[i12];
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

static void controller_codegen_entry_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->param.J[i] = dv[i];
    b_SD->pd->param.Jinv[i] = dv1[i];
  }
  b_SD->pd->param.c_aero = -0.016182736457722724;
  b_SD->pd->param.c_canard = 0.00061367415999999994;
  b_SD->pd->param.elevation = 420.0;
  b_SD->pd->param.g[0] = -9.81;
  b_SD->pd->param.g[1] = 0.0;
  b_SD->pd->param.g[2] = 0.0;
}

static void dynamics_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->c_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->c_param.J[i] = dv[i];
    b_SD->pd->c_param.Jinv[i] = dv1[i];
  }
  b_SD->pd->c_param.c_aero = -0.016182736457722724;
  b_SD->pd->c_param.c_canard = 0.00061367415999999994;
  b_SD->pd->c_param.elevation = 420.0;
  b_SD->pd->c_param.g[0] = -9.81;
  b_SD->pd->c_param.g[1] = 0.0;
  b_SD->pd->c_param.g[2] = 0.0;
}

static void dynamics_jacobian_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->d_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->d_param.J[i] = dv[i];
    b_SD->pd->d_param.Jinv[i] = dv1[i];
  }
  b_SD->pd->d_param.c_aero = -0.016182736457722724;
  b_SD->pd->d_param.c_canard = 0.00061367415999999994;
  b_SD->pd->d_param.elevation = 420.0;
  b_SD->pd->d_param.g[0] = -9.81;
  b_SD->pd->d_param.g[1] = 0.0;
  b_SD->pd->d_param.g[2] = 0.0;
}

static void ekf_correct(const double x[11], const double P[121],
                        const double y[3], const double b[3], const double R[9],
                        double x_new[11], double P_new[121]) {
  static const signed char iv[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
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
  double d12;
  double d13;
  double d14;
  double d15;
  double d16;
  double d17;
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
  int i28;
  int i29;
  int i3;
  int i30;
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
        2.0 * (((b_x * (double)iv[3 * i] + x[1] * b[i]) - b[0] * H_tmp) +
               b_dv[3 * i]);
    c_H_tmp = 3 * i + 1;
    H[b_H_tmp + 1] =
        2.0 * (((b_x * (double)iv[c_H_tmp] + x[2] * b[i]) - b[1] * H_tmp) +
               b_dv[c_H_tmp]);
    d_H_tmp = 3 * i + 2;
    H[b_H_tmp + 2] =
        2.0 * (((b_x * (double)iv[d_H_tmp] + x[3] * b[i]) - b[2] * H_tmp) +
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
    for (i6 = 0; i6 < 3; i6++) {
      double d2;
      d2 = 0.0;
      for (i8 = 0; i8 < 11; i8++) {
        d2 += b_H[i4 + 3 * i8] * y_tmp[i8 + 11 * i6];
      }
      int g_H_tmp;
      g_H_tmp = i4 + 3 * i6;
      c_H[g_H_tmp] = d2 + R[g_H_tmp];
    }
    for (i7 = 0; i7 < 11; i7++) {
      double d3;
      d3 = y_tmp[i7 + 11 * i4];
      for (i10 = 0; i10 < 11; i10++) {
        int P_tmp;
        P_tmp = i10 + 11 * i4;
        b_P[P_tmp] += P[i10 + 11 * i7] * d3;
      }
    }
  }
  inv(c_H, b_dv);
  memset(&K[0], 0, 33U * sizeof(double));
  for (i9 = 0; i9 < 3; i9++) {
    for (i11 = 0; i11 < 3; i11++) {
      double d4;
      d4 = b_dv[i11 + 3 * i9];
      for (i12 = 0; i12 < 11; i12++) {
        int K_tmp;
        K_tmp = i12 + 11 * i9;
        K[K_tmp] += b_P[i12 + 11 * i11] * d4;
      }
    }
  }
  memset(&b_I[0], 0, 121U * sizeof(signed char));
  for (k = 0; k < 11; k++) {
    b_I[k + 11 * k] = 1;
  }
  for (i13 = 0; i13 < 11; i13++) {
    double d5;
    double d6;
    double d7;
    d5 = K[i13];
    d6 = K[i13 + 11];
    d7 = K[i13 + 22];
    for (i15 = 0; i15 < 11; i15++) {
      int E_tmp;
      E_tmp = i13 + 11 * i15;
      E[E_tmp] = (double)b_I[E_tmp] - ((d5 * H[3 * i15] + d6 * H[3 * i15 + 1]) +
                                       d7 * H[3 * i15 + 2]);
    }
  }
  memset(&b_E[0], 0, 121U * sizeof(double));
  for (i14 = 0; i14 < 11; i14++) {
    for (i16 = 0; i16 < 11; i16++) {
      double d8;
      d8 = P[i16 + 11 * i14];
      for (i18 = 0; i18 < 11; i18++) {
        int b_E_tmp;
        b_E_tmp = i18 + 11 * i14;
        b_E[b_E_tmp] += E[i18 + 11 * i16] * d8;
      }
    }
  }
  memset(&b_K[0], 0, 33U * sizeof(double));
  for (i17 = 0; i17 < 3; i17++) {
    for (i19 = 0; i19 < 3; i19++) {
      double d9;
      d9 = R[i19 + 3 * i17];
      for (i20 = 0; i20 < 11; i20++) {
        int b_K_tmp;
        b_K_tmp = i20 + 11 * i17;
        b_K[b_K_tmp] += K[i20 + 11 * i19] * d9;
      }
    }
  }
  memset(&c_E[0], 0, 121U * sizeof(double));
  memset(&c_K[0], 0, 121U * sizeof(double));
  for (i21 = 0; i21 < 11; i21++) {
    for (i22 = 0; i22 < 11; i22++) {
      double d10;
      d10 = E[i21 + 11 * i22];
      for (i25 = 0; i25 < 11; i25++) {
        int c_E_tmp;
        c_E_tmp = i25 + 11 * i21;
        c_E[c_E_tmp] += b_E[i25 + 11 * i22] * d10;
      }
    }
    for (i24 = 0; i24 < 3; i24++) {
      double d11;
      d11 = K[i21 + 11 * i24];
      for (i27 = 0; i27 < 11; i27++) {
        int c_K_tmp;
        c_K_tmp = i27 + 11 * i21;
        c_K[c_K_tmp] += b_K[i27 + 11 * i24] * d11;
      }
    }
  }
  for (i23 = 0; i23 < 121; i23++) {
    P_new[i23] = c_E[i23] + c_K[i23];
  }
  for (i26 = 0; i26 < 3; i26++) {
    double a_tmp;
    int b_a_tmp;
    int c_a_tmp;
    a_tmp = x[i26 + 1];
    c_a[3 * i26] = a * (double)iv[3 * i26] + 2.0 * x[1] * a_tmp;
    b_a_tmp = 3 * i26 + 1;
    c_a[b_a_tmp] = a * (double)iv[b_a_tmp] + 2.0 * x[2] * a_tmp;
    c_a_tmp = 3 * i26 + 2;
    c_a[c_a_tmp] = a * (double)iv[c_a_tmp] + 2.0 * x[3] * a_tmp;
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
  for (i28 = 0; i28 < 9; i28++) {
    c_a[i28] -= b_dv[i28];
  }
  d12 = b[0];
  d13 = b[1];
  d14 = b[2];
  for (i29 = 0; i29 < 3; i29++) {
    b_y[i29] =
        y[i29] - ((c_a[i29] * d12 + c_a[i29 + 3] * d13) + c_a[i29 + 6] * d14);
  }
  d15 = b_y[0];
  d16 = b_y[1];
  d17 = b_y[2];
  for (i30 = 0; i30 < 11; i30++) {
    x_new[i30] =
        x[i30] + ((K[i30] * d15 + K[i30 + 11] * d16) + K[i30 + 22] * d17);
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

static void inv(const double x[9], double y[9]) {
  double b_x[9];
  double absx11;
  double absx21;
  double absx31;
  double t2;
  double t3;
  int p1;
  int p2;
  int p3;
  memcpy(&b_x[0], &x[0], 9U * sizeof(double));
  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = fabs(x[0]);
  absx21 = fabs(x[1]);
  absx31 = fabs(x[2]);
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    b_x[0] = x[1];
    b_x[1] = x[0];
    b_x[3] = x[4];
    b_x[4] = x[3];
    b_x[6] = x[7];
    b_x[7] = x[6];
  } else if (absx31 > absx11) {
    p1 = 6;
    p3 = 0;
    b_x[0] = x[2];
    b_x[2] = x[0];
    b_x[3] = x[5];
    b_x[5] = x[3];
    b_x[6] = x[8];
    b_x[8] = x[6];
  }
  b_x[1] /= b_x[0];
  b_x[2] /= b_x[0];
  b_x[4] -= b_x[1] * b_x[3];
  b_x[5] -= b_x[2] * b_x[3];
  b_x[7] -= b_x[1] * b_x[6];
  b_x[8] -= b_x[2] * b_x[6];
  if (fabs(b_x[5]) > fabs(b_x[4])) {
    double t1;
    int itmp;
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    t1 = b_x[1];
    b_x[1] = b_x[2];
    b_x[2] = t1;
    t1 = b_x[4];
    b_x[4] = b_x[5];
    b_x[5] = t1;
    t1 = b_x[7];
    b_x[7] = b_x[8];
    b_x[8] = t1;
  }
  b_x[5] /= b_x[4];
  b_x[8] -= b_x[5] * b_x[7];
  t3 = (b_x[1] * b_x[5] - b_x[2]) / b_x[8];
  t2 = -(b_x[1] + b_x[7] * t3) / b_x[4];
  y[p1] = ((1.0 - b_x[3] * t2) - b_x[6] * t3) / b_x[0];
  y[p1 + 1] = t2;
  y[p1 + 2] = t3;
  t3 = -b_x[5] / b_x[8];
  t2 = (1.0 - b_x[7] * t3) / b_x[4];
  y[p2] = -(b_x[3] * t2 + b_x[6] * t3) / b_x[0];
  y[p2 + 1] = t2;
  y[p2 + 2] = t3;
  t3 = 1.0 / b_x[8];
  t2 = -b_x[7] * t3 / b_x[4];
  y[p3] = -(b_x[3] * t2 + b_x[6] * t3) / b_x[0];
  y[p3 + 1] = t2;
  y[p3 + 2] = t3;
}

static void pad_filter_init(GNC_codegen_SILStackData *b_SD) {
  int i;
  b_SD->pd->b_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    b_SD->pd->b_param.J[i] = dv[i];
    b_SD->pd->b_param.Jinv[i] = dv1[i];
  }
  b_SD->pd->b_param.c_aero = -0.016182736457722724;
  b_SD->pd->b_param.c_canard = 0.00061367415999999994;
  b_SD->pd->b_param.elevation = 420.0;
  b_SD->pd->b_param.g[0] = -9.81;
  b_SD->pd->b_param.g[1] = 0.0;
  b_SD->pd->b_param.g[2] = 0.0;
}

void GNC_codegen_SIL_initialize(GNC_codegen_SILStackData *b_SD) {
  controller_codegen_entry_init(b_SD);
  pad_filter_init(b_SD);
  dynamics_init(b_SD);
  dynamics_jacobian_init(b_SD);
}

void GNC_codegen_SIL_terminate(void) {}

void controller_codegen_entry(GNC_codegen_SILStackData *b_SD, double b_time,
                              double dt_ctrl, const double xR[2], double pdyn,
                              double delta, const struct0_T *ctrl_mem_in,
                              double *u, double *r, struct0_T *ctrl_mem_out) {
  double P[4];
  double b_K[4];
  double b_dv[4];
  double b_dv1[4];
  double dv2[4];
  double K[2];
  double b_r[2];
  double L_delta;
  double a;
  double b;
  double blend;
  double c_delta;
  double c_r;
  double d;
  double d1;
  double d10;
  double d11;
  double d12;
  double d2;
  double d3;
  double d4;
  double d5;
  double d8;
  double pdyn_params;
  double r_idx_0;
  double w_dot;
  double x;
  double x_tmp;
  int i;
  int i2;
  if ((b_time >= 22.0) && (b_time < 27.0)) {
    *r = 0.5;
  } else if ((b_time >= 27.0) && (b_time < 32.0)) {
    *r = -0.5;
  } else if ((b_time >= 32.0) && (b_time < 39.0)) {
    *r = 0.5;
  } else {
    *r = 0.0;
  }
  pdyn_params = pdyn * b_SD->pd->param.c_canard;
  c_delta = delta / 2.0;
  if (fabs(c_delta) < 0.005) {
    c_delta = 0.0;
  }
  c_delta = 0.75 * ctrl_mem_in->d_old + 0.25 * c_delta;
  w_dot = 0.75 * ctrl_mem_in->w_dot_old +
          0.25 * (xR[1] - ctrl_mem_in->w_old) / dt_ctrl;
  r_idx_0 = pdyn_params * c_delta;
  P[0] = ctrl_mem_in->P_minus[0] + 1.0E-5;
  P[1] = ctrl_mem_in->P_minus[1];
  P[2] = ctrl_mem_in->P_minus[2];
  P[3] = ctrl_mem_in->P_minus[3] + 1.0E-9;
  memset(&b_r[0], 0, sizeof(double) << 1);
  d = r_idx_0 * (ctrl_mem_in->P_minus[0] + 1.0E-5);
  d1 = pdyn_params * (ctrl_mem_in->P_minus[3] + 1.0E-9);
  c_r = ((b_r[0] + d) + pdyn_params * ctrl_mem_in->P_minus[1]) * r_idx_0 +
        ((b_r[1] + r_idx_0 * ctrl_mem_in->P_minus[2]) + d1) * pdyn_params;
  K[0] = (d + ctrl_mem_in->P_minus[2] * pdyn_params) / (c_r + 1.0);
  K[1] = (ctrl_mem_in->P_minus[1] * r_idx_0 + d1) / (c_r + 1.0);
  b = w_dot -
      (r_idx_0 * ctrl_mem_in->coeffs[0] + pdyn_params * ctrl_mem_in->coeffs[1]);
  ctrl_mem_out->coeffs[0] = ctrl_mem_in->coeffs[0] + K[0] * b;
  ctrl_mem_out->coeffs[1] = ctrl_mem_in->coeffs[1] + K[1] * b;
  b_dv[0] = 1.0 - K[0] * r_idx_0;
  b_dv[1] = 0.0 - K[1] * r_idx_0;
  b_dv[2] = 0.0 - K[0] * pdyn_params;
  b_dv[3] = 1.0 - K[1] * pdyn_params;
  memset(&b_dv1[0], 0, sizeof(double) << 2);
  d2 = b_dv[0];
  d3 = b_dv[1];
  d4 = b_dv[2];
  d5 = b_dv[3];
  for (i = 0; i < 2; i++) {
    double d6;
    double d7;
    double d9;
    int i1;
    i1 = i << 1;
    d6 = P[i1];
    d7 = b_dv1[i1] + d2 * d6;
    d9 = b_dv1[i1 + 1] + d3 * d6;
    d6 = P[i1 + 1];
    d7 += d4 * d6;
    b_dv1[i1] = d7;
    d9 += d5 * d6;
    b_dv1[i1 + 1] = d9;
  }
  memset(&dv2[0], 0, sizeof(double) << 2);
  d8 = b_dv1[0];
  d10 = b_dv1[1];
  d11 = b_dv1[2];
  d12 = b_dv1[3];
  for (i2 = 0; i2 < 2; i2++) {
    double d13;
    double d14;
    double d15;
    int i3;
    d13 = b_dv[i2];
    i3 = i2 << 1;
    d14 = dv2[i3] + d8 * d13;
    d15 = dv2[i3 + 1] + d10 * d13;
    b_K[i3] = K[0] * K[i2];
    d13 = b_dv[i2 + 2];
    d14 += d11 * d13;
    dv2[i3] = d14;
    d15 += d12 * d13;
    dv2[i3 + 1] = d15;
    b_K[i3 + 1] = K[1] * K[i2];
  }
  ctrl_mem_out->P_minus[0] = dv2[0] + b_K[0];
  ctrl_mem_out->P_minus[1] = dv2[1] + b_K[1];
  ctrl_mem_out->P_minus[2] = dv2[2] + b_K[2];
  ctrl_mem_out->P_minus[3] = dv2[3] + b_K[3];
  ctrl_mem_out->w_old = xR[1];
  ctrl_mem_out->d_old = c_delta;
  ctrl_mem_out->w_dot_old = w_dot;
  L_delta = ctrl_mem_out->coeffs[0] * pdyn_params / 2.0;
  if (fabs(L_delta) < 10.0) {
    if (L_delta >= 0.0) {
      L_delta = 10.0;
    } else {
      L_delta = -10.0;
    }
  }
  blend = fmax(0.0, fmin(1.0, (fabs(xR[1]) - 0.5) / 0.5));
  a = -1.0 / L_delta;
  x_tmp = (1.0 - blend) * 5.0;
  x = sqrt(x_tmp);
  K[0] = a * x;
  *u = fmin(
      fmax((K[0] * xR[0] + a * sqrt(2.0 * x + (x_tmp + blend * 20.0)) * xR[1]) +
               -K[0] * *r,
           -0.3490658503988659),
      0.3490658503988659);
  if (pdyn < 500.0) {
    *u = 0.0;
  }
}

void navigation_codegen_entry(GNC_codegen_SILStackData *b_SD, double dt,
                              boolean_T flight_phase, const double x[11],
                              const double P[121], const struct1_T *bias,
                              const struct2_T *sens_filt,
                              const struct3_T *sens_input, double x_ret[11],
                              double P_ret[121], struct1_T *bias_ret,
                              struct2_T *sens_filt_ret, struct6_T *airdata,
                              double roll_state[2]) {
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
  double F[121];
  double b_E[121];
  double b_F[121];
  double b_P_ret[121];
  double c_E[121];
  double c_K[121];
  double c_P_ret[121];
  double d_P_ret[121];
  double e_P_ret[121];
  double f_P_ret[121];
  double g_P_ret[121];
  double h_P_ret[121];
  double i_P_ret[121];
  double K[33];
  double b_K[33];
  double b_W_dt[16];
  double c_b[16];
  double d_b[16];
  double b_x_ret[11];
  double c_x_ret[11];
  double d_x_ret[11];
  double e_x_ret[11];
  double f_x_ret[11];
  double g_x_ret[11];
  double h_x_ret[11];
  double i_x_ret[11];
  double b_dv[9];
  double b_n_tilde[9];
  double b_w_exp_tilde[9];
  double c_skewed_exp_w_tmp[9];
  double c_q[4];
  double q[4];
  double b_dt[3];
  double c_dt[3];
  double c_w_exp_tilde[3];
  double dv3[3];
  double airspeed;
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
  int b_k;
  int c_k;
  int i;
  int i1;
  int i10;
  int i11;
  int i12;
  int i15;
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
  int i28;
  int i29;
  int i3;
  int i30;
  int i31;
  int i32;
  int i33;
  int i34;
  int i36;
  int i37;
  int i38;
  int i39;
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
  int i7;
  int i8;
  int i9;
  int k;
  signed char c_I[121];
  memcpy(&P_ret[0], &P[0], 121U * sizeof(double));
  *bias_ret = *bias;
  *sens_filt_ret = *sens_filt;
  if (!flight_phase) {
    double ST[9];
    double h_a[9];
    double a[3];
    double b_absxk;
    double b_scale;
    double b_t;
    double board_baro_f;
    double d11;
    double d12;
    double d5;
    double d6;
    double d9;
    double d_a;
    double f_a;
    double mti_baro_f;
    double qw;
    double qy;
    double qz;
    double y;
    board_baro_f = sens_filt->board_baro_f;
    mti_baro_f = sens_filt->mti_baro_f;
    if (sens_input->board_accel.status) {
      sens_filt_ret->board_accel_f[0] =
          0.0005 * sens_input->board_accel.meas[0] +
          0.9995 * sens_filt->board_accel_f[0];
      sens_filt_ret->board_accel_f[1] =
          0.0005 * sens_input->board_accel.meas[1] +
          0.9995 * sens_filt->board_accel_f[1];
      sens_filt_ret->board_accel_f[2] =
          0.0005 * sens_input->board_accel.meas[2] +
          0.9995 * sens_filt->board_accel_f[2];
    }
    if (sens_input->board_gyro.status) {
      sens_filt_ret->board_gyro_f[0] = 0.0005 * sens_input->board_gyro.meas[0] +
                                       0.9995 * sens_filt->board_gyro_f[0];
      sens_filt_ret->board_gyro_f[1] = 0.0005 * sens_input->board_gyro.meas[1] +
                                       0.9995 * sens_filt->board_gyro_f[1];
      sens_filt_ret->board_gyro_f[2] = 0.0005 * sens_input->board_gyro.meas[2] +
                                       0.9995 * sens_filt->board_gyro_f[2];
    }
    if (sens_input->mti_accel.status) {
      sens_filt_ret->mti_accel_f[0] = 0.0005 * sens_input->mti_accel.meas[0] +
                                      0.9995 * sens_filt->mti_accel_f[0];
      sens_filt_ret->mti_accel_f[1] = 0.0005 * sens_input->mti_accel.meas[1] +
                                      0.9995 * sens_filt->mti_accel_f[1];
      sens_filt_ret->mti_accel_f[2] = 0.0005 * sens_input->mti_accel.meas[2] +
                                      0.9995 * sens_filt->mti_accel_f[2];
    }
    if (sens_input->mti_gyro.status) {
      sens_filt_ret->mti_gyro_f[0] = 0.0005 * sens_input->mti_gyro.meas[0] +
                                     0.9995 * sens_filt->mti_gyro_f[0];
      sens_filt_ret->mti_gyro_f[1] = 0.0005 * sens_input->mti_gyro.meas[1] +
                                     0.9995 * sens_filt->mti_gyro_f[1];
      sens_filt_ret->mti_gyro_f[2] = 0.0005 * sens_input->mti_gyro.meas[2] +
                                     0.9995 * sens_filt->mti_gyro_f[2];
    }
    if (sens_input->ad_accel.status) {
      sens_filt_ret->ad_accel_f[0] = 0.0005 * sens_input->ad_accel.meas[0] +
                                     0.9995 * sens_filt->ad_accel_f[0];
      sens_filt_ret->ad_accel_f[1] = 0.0005 * sens_input->ad_accel.meas[1] +
                                     0.9995 * sens_filt->ad_accel_f[1];
      sens_filt_ret->ad_accel_f[2] = 0.0005 * sens_input->ad_accel.meas[2] +
                                     0.9995 * sens_filt->ad_accel_f[2];
    }
    if (sens_input->ad_gyro.status) {
      sens_filt_ret->ad_gyro_f[0] = 0.0005 * sens_input->ad_gyro.meas[0] +
                                    0.9995 * sens_filt->ad_gyro_f[0];
      sens_filt_ret->ad_gyro_f[1] = 0.0005 * sens_input->ad_gyro.meas[1] +
                                    0.9995 * sens_filt->ad_gyro_f[1];
      sens_filt_ret->ad_gyro_f[2] = 0.0005 * sens_input->ad_gyro.meas[2] +
                                    0.9995 * sens_filt->ad_gyro_f[2];
    }
    if (sens_input->board_baro.status) {
      board_baro_f = 0.0005 * sens_input->board_baro.meas +
                     0.9995 * sens_filt->board_baro_f;
    }
    if (sens_input->board_mag.status) {
      sens_filt_ret->board_mag_f[0] = 0.0005 * sens_input->board_mag.meas[0] +
                                      0.9995 * sens_filt->board_mag_f[0];
      sens_filt_ret->board_mag_f[1] = 0.0005 * sens_input->board_mag.meas[1] +
                                      0.9995 * sens_filt->board_mag_f[1];
      sens_filt_ret->board_mag_f[2] = 0.0005 * sens_input->board_mag.meas[2] +
                                      0.9995 * sens_filt->board_mag_f[2];
    }
    if (sens_input->mti_baro.status) {
      mti_baro_f =
          0.0005 * sens_input->mti_baro.meas + 0.9995 * sens_filt->mti_baro_f;
    }
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
    a[0] = 0.0;
    a[1] = 0.0;
    a[2] = 0.0;
    if (sens_input->board_accel.status) {
      a[0] = sens_filt_ret->board_accel_f[0];
      a[1] = sens_filt_ret->board_accel_f[1];
      a[2] = sens_filt_ret->board_accel_f[2];
    }
    if (sens_input->mti_accel.status) {
      a[0] += sens_filt_ret->mti_accel_f[0];
      a[1] += sens_filt_ret->mti_accel_f[1];
      a[2] += sens_filt_ret->mti_accel_f[2];
    }
    if (sens_input->ad_accel.status) {
      a[0] += sens_filt_ret->ad_accel_f[0];
      a[1] += sens_filt_ret->ad_accel_f[1];
      a[2] += sens_filt_ret->ad_accel_f[2];
    }
    d5 = b_norm(a) + 1.0E-6;
    qw = sqrt(0.5 * (a[0] / d5) + 0.5);
    if (qw == 0.0) {
      qy = 1.0;
      qz = 0.0;
    } else {
      qy = 0.5 * (a[2] / d5) / qw;
      qz = -0.5 * (a[1] / d5) / qw;
    }
    b_scale = 3.3121686421112381E-170;
    if (qw > 3.3121686421112381E-170) {
      y = 1.0;
      b_scale = qw;
    } else {
      b_t = qw / 3.3121686421112381E-170;
      y = b_t * b_t;
    }
    b_absxk = fabs(qy);
    if (b_absxk > b_scale) {
      b_t = b_scale / b_absxk;
      y = y * b_t * b_t + 1.0;
      b_scale = b_absxk;
    } else {
      b_t = b_absxk / b_scale;
      y += b_t * b_t;
    }
    b_absxk = fabs(qz);
    if (b_absxk > b_scale) {
      b_t = b_scale / b_absxk;
      y = y * b_t * b_t + 1.0;
      b_scale = b_absxk;
    } else {
      b_t = b_absxk / b_scale;
      y += b_t * b_t;
    }
    y = b_scale * sqrt(y);
    d6 = qw / y;
    q[0] = d6;
    x_ret[0] = d6;
    d6 = 0.0 / y;
    q[1] = d6;
    x_ret[1] = d6;
    d6 = qy / y;
    q[2] = d6;
    x_ret[2] = d6;
    d6 = qz / y;
    q[3] = d6;
    x_ret[3] = d6;
    x_ret[10] = b_SD->pd->b_param.elevation;
    x_ret[4] = 0.0;
    x_ret[7] = 0.0;
    bias_ret->board_gyro[0] = sens_filt_ret->board_gyro_f[0];
    bias_ret->mti_gyro[0] = sens_filt_ret->mti_gyro_f[0];
    bias_ret->ad_gyro[0] = sens_filt_ret->ad_gyro_f[0];
    x_ret[5] = 0.0;
    x_ret[8] = 0.0;
    bias_ret->board_gyro[1] = sens_filt_ret->board_gyro_f[1];
    bias_ret->mti_gyro[1] = sens_filt_ret->mti_gyro_f[1];
    bias_ret->ad_gyro[1] = sens_filt_ret->ad_gyro_f[1];
    x_ret[6] = 0.0;
    x_ret[9] = 0.0;
    bias_ret->board_gyro[2] = sens_filt_ret->board_gyro_f[2];
    bias_ret->mti_gyro[2] = sens_filt_ret->mti_gyro_f[2];
    bias_ret->ad_gyro[2] = sens_filt_ret->ad_gyro_f[2];
    d_a = q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + d6 * d6);
    f_a = 2.0 * q[0];
    for (i3 = 0; i3 < 3; i3++) {
      double a_tmp;
      a_tmp = 2.0 * q[i3 + 1];
      h_a[3 * i3] = d_a * b_b[i3] + a_tmp * q[1];
      h_a[3 * i3 + 1] = d_a * b_b[i3 + 3] + a_tmp * q[2];
      h_a[3 * i3 + 2] = d_a * b_b[i3 + 6] + a_tmp * d6;
    }
    b_dv[0] = 0.0;
    b_dv[1] = f_a * -d6;
    b_dv[2] = f_a * q[2];
    b_dv[3] = f_a * d6;
    b_dv[4] = 0.0;
    b_dv[5] = f_a * -q[1];
    b_dv[6] = f_a * -q[2];
    b_dv[7] = f_a * q[1];
    b_dv[8] = 0.0;
    for (i5 = 0; i5 < 9; i5++) {
      ST[i5] = h_a[i5] - b_dv[i5];
    }
    bias_ret->board_mag_earth[0] = 0.0;
    bias_ret->board_mag_earth[1] = 0.0;
    bias_ret->board_mag_earth[2] = 0.0;
    for (i7 = 0; i7 < 3; i7++) {
      double d10;
      d10 = sens_filt_ret->board_mag_f[i7];
      bias_ret->board_mag_earth[0] += ST[3 * i7] * d10;
      bias_ret->board_mag_earth[1] += ST[3 * i7 + 1] * d10;
      bias_ret->board_mag_earth[2] += ST[3 * i7 + 2] * d10;
      bias_ret->mti_mag_earth[i7] = 0.0;
    }
    d9 = bias_ret->mti_mag_earth[0];
    d11 = bias_ret->mti_mag_earth[1];
    d12 = bias_ret->mti_mag_earth[2];
    for (i8 = 0; i8 < 3; i8++) {
      double d13;
      d13 = sens_filt_ret->mti_mag_f[i8];
      d9 += ST[3 * i8] * d13;
      d11 += ST[3 * i8 + 1] * d13;
      d12 += ST[3 * i8 + 2] * d13;
    }
    double t1_pressure;
    bias_ret->mti_mag_earth[2] = d12;
    bias_ret->mti_mag_earth[1] = d11;
    bias_ret->mti_mag_earth[0] = d9;
    t1_pressure =
        airdata_atmos(b_SD->pd->b_param.elevation, &e_expl_temp, &t1_density,
                      &f_expl_temp, &g_expl_temp, &h_expl_temp);
    bias_ret->board_baro = board_baro_f - t1_pressure;
    bias_ret->mti_baro = mti_baro_f - t1_pressure;
  } else {
    double E[121];
    double P_pred[121];
    double W_dt[16];
    double b_q[16];
    double l_a[16];
    double d_dt[12];
    double x_pred[11];
    double S[9];
    double b_P_pred[9];
    double b_skewed_exp_w_tmp[9];
    double dv4[9];
    double h_a[9];
    double n_tilde[9];
    double skewed_exp_w_tmp[9];
    double w_exp_tilde[9];
    double b_dv1[4];
    double r_q_tmp[4];
    double C_total_a[3];
    double b_S[3];
    double c_r_q_tmp[3];
    double d_x[3];
    double dn[3];
    double dv2[3];
    double C_ad_w_idx_0;
    double C_total_a_tmp;
    double C_total_a_tmp_tmp;
    double absxk;
    double b;
    double b_C_total_a_tmp_tmp;
    double b_a;
    double b_dphi_tmp;
    double b_q_mag;
    double b_r_q_tmp;
    double b_x;
    double c_C_total_a_tmp_tmp;
    double c_absxk;
    double c_scale;
    double c_t;
    double c_x;
    double d;
    double d1;
    double d16;
    double d17;
    double d18;
    double d19;
    double d2;
    double d20;
    double d21;
    double d22;
    double d23;
    double d24;
    double d25;
    double d26;
    double d27;
    double d3;
    double d30;
    double d31;
    double d33;
    double d34;
    double d4;
    double d72;
    double d73;
    double d74;
    double dphi;
    double dphi_tmp;
    double e_a;
    double g_a;
    double i_a;
    double j_a;
    double k_a;
    double m_a;
    double n_a;
    double n_idx_0;
    double n_idx_1;
    double n_idx_2;
    double o_a;
    double q_mag;
    double scale;
    double t;
    signed char b_I[16];
    signed char w_exp_tilde_tmp[9];
    d = 9.9999999999999981E+9 * (double)sens_input->ad_gyro.status;
    C_total_a_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_input->board_accel.status;
    b_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_input->mti_accel.status;
    c_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * (double)sens_input->ad_accel.status;
    C_total_a_tmp =
        (C_total_a_tmp_tmp + b_C_total_a_tmp_tmp) + c_C_total_a_tmp_tmp;
    C_total_a[0] = C_total_a_tmp;
    d1 = 9.9999999999999981E+9 * (double)sens_input->board_gyro.status;
    d2 = 9.9999999999999981E+9 * (double)sens_input->mti_gyro.status;
    d3 = d1 + d2;
    d4 = d3 + d;
    d /= d4;
    C_ad_w_idx_0 = d;
    C_total_a[1] = C_total_a_tmp;
    d = 0.0 / d3;
    C_total_a[2] = C_total_a_tmp;
    scale = 3.3121686421112381E-170;
    absxk = fabs(x[0]);
    if (absxk > 3.3121686421112381E-170) {
      q_mag = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      q_mag = t * t;
    }
    absxk = fabs(x[1]);
    if (absxk > scale) {
      t = scale / absxk;
      q_mag = q_mag * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      q_mag += t * t;
    }
    absxk = fabs(x[2]);
    if (absxk > scale) {
      t = scale / absxk;
      q_mag = q_mag * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      q_mag += t * t;
    }
    absxk = fabs(x[3]);
    if (absxk > scale) {
      t = scale / absxk;
      q_mag = q_mag * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      q_mag += t * t;
    }
    q_mag = scale * sqrt(q_mag);
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
      double b_dn_tmp;
      double c_dn_tmp;
      double dn_tmp;
      dn_tmp = x[4] / dphi_tmp;
      dn[0] = dn_tmp;
      b_dn_tmp = x[5] / dphi_tmp;
      dn[1] = b_dn_tmp;
      c_dn_tmp = x[6] / dphi_tmp;
      dn[2] = c_dn_tmp;
      n_idx_0 = dn_tmp;
      n_idx_1 = b_dn_tmp;
      n_idx_2 = c_dn_tmp;
    }
    b = sin(dphi);
    n_tilde[0] = 0.0;
    n_tilde[3] = -n_idx_2;
    n_tilde[6] = n_idx_1;
    n_tilde[1] = n_idx_2;
    n_tilde[4] = 0.0;
    n_tilde[7] = -n_idx_0;
    n_tilde[2] = -n_idx_1;
    n_tilde[5] = n_idx_0;
    n_tilde[8] = 0.0;
    b_a = sin(b_dphi_tmp);
    b_x = cos(b_dphi_tmp);
    for (i = 0; i < 9; i++) {
      w_exp_tilde_tmp[i] = 0;
    }
    memset(&b_n_tilde[0], 0, 9U * sizeof(double));
    for (k = 0; k < 3; k++) {
      double d7;
      int b_n_tilde_tmp;
      int n_tilde_tmp;
      w_exp_tilde_tmp[k + 3 * k] = 1;
      d7 = b_n_tilde[3 * k];
      n_tilde_tmp = 3 * k + 1;
      b_n_tilde_tmp = 3 * k + 2;
      for (i2 = 0; i2 < 3; i2++) {
        double d8;
        d8 = n_tilde[i2 + 3 * k];
        d7 += n_tilde[3 * i2] * d8;
        b_n_tilde[n_tilde_tmp] += n_tilde[3 * i2 + 1] * d8;
        b_n_tilde[b_n_tilde_tmp] += n_tilde[3 * i2 + 2] * d8;
      }
      b_n_tilde[3 * k] = d7;
    }
    for (i1 = 0; i1 < 9; i1++) {
      w_exp_tilde[i1] = ((double)w_exp_tilde_tmp[i1] - b_a * n_tilde[i1]) +
                        (1.0 - b_x) * b_n_tilde[i1];
    }
    double c_a;
    c_a = b_norm(&x[7]);
    airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp,
                  &d_expl_temp);
    e_a = 0.5 * t1_density * (c_a * c_a);
    g_a = b_SD->pd->c_param.c_aero * b_SD->pd->c_param.Cn_alpha;
    i_a = x[0] * x[0] - ((x[1] * x[1] + x[2] * x[2]) + x[3] * x[3]);
    j_a = 2.0 * x[0];
    for (i4 = 0; i4 < 3; i4++) {
      double b_a_tmp;
      int c_a_tmp;
      int d_a_tmp;
      b_a_tmp = x[i4 + 1];
      h_a[3 * i4] = i_a * b_b[3 * i4] + 2.0 * x[1] * b_a_tmp;
      c_a_tmp = 3 * i4 + 1;
      h_a[c_a_tmp] = i_a * b_b[c_a_tmp] + 2.0 * x[2] * b_a_tmp;
      d_a_tmp = 3 * i4 + 2;
      h_a[d_a_tmp] = i_a * b_b[d_a_tmp] + 2.0 * x[3] * b_a_tmp;
    }
    b_dv[0] = 0.0;
    b_dv[3] = j_a * -x[3];
    b_dv[6] = j_a * x[2];
    b_dv[1] = j_a * x[3];
    b_dv[4] = 0.0;
    b_dv[7] = j_a * -x[1];
    b_dv[2] = j_a * -x[2];
    b_dv[5] = j_a * x[1];
    b_dv[8] = 0.0;
    for (i6 = 0; i6 < 9; i6++) {
      S[i6] = h_a[i6] - b_dv[i6];
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
    for (i9 = 0; i9 < 3; i9++) {
      double d14;
      int b_w_exp_tilde_tmp;
      int c_w_exp_tilde_tmp;
      b_dv1[i9 + 1] = dn[i9] * b;
      d14 = b_w_exp_tilde[3 * i9];
      b_w_exp_tilde_tmp = 3 * i9 + 1;
      c_w_exp_tilde_tmp = 3 * i9 + 2;
      for (i10 = 0; i10 < 3; i10++) {
        double d15;
        d15 = b_SD->pd->c_param.J[i10 + 3 * i9];
        d14 += w_exp_tilde[3 * i10] * d15;
        b_w_exp_tilde[b_w_exp_tilde_tmp] += w_exp_tilde[3 * i10 + 1] * d15;
        b_w_exp_tilde[c_w_exp_tilde_tmp] += w_exp_tilde[3 * i10 + 2] * d15;
      }
      double d_w_exp_tilde_tmp;
      b_w_exp_tilde[3 * i9] = d14;
      d_w_exp_tilde_tmp = x[i9 + 4];
      c_w_exp_tilde[0] += d14 * d_w_exp_tilde_tmp;
      c_w_exp_tilde[1] += b_w_exp_tilde[3 * i9 + 1] * d_w_exp_tilde_tmp;
      c_w_exp_tilde[2] += b_w_exp_tilde[3 * i9 + 2] * d_w_exp_tilde_tmp;
    }
    dv2[0] = 0.0;
    dv2[1] = e_a * (g_a * sin(atan2(x[9], x[7])));
    dv2[2] = e_a * (g_a * -sin(atan2(x[8], x[7])));
    memset(&dv3[0], 0, 3U * sizeof(double));
    memset(&b_dt[0], 0, 3U * sizeof(double));
    memset(&c_dt[0], 0, 3U * sizeof(double));
    d16 = dv3[0];
    d17 = dv3[1];
    d18 = dv3[2];
    d19 = b_dt[0];
    d20 = b_dt[1];
    d21 = b_dt[2];
    d22 = c_dt[0];
    d23 = c_dt[1];
    d24 = c_dt[2];
    d25 = x[7];
    d26 = x[8];
    d27 = x[9];
    for (i11 = 0; i11 < 3; i11++) {
      double d28;
      double d29;
      double d32;
      double d35;
      double d36;
      double d37;
      double d39;
      int i13;
      int i14;
      d28 = b_SD->pd->c_param.Jinv[3 * i11];
      d29 = c_w_exp_tilde[i11];
      d16 += d28 * d29;
      d32 = dv2[i11];
      d19 += dt * d28 * d32;
      d35 = S[3 * i11];
      d36 = b_SD->pd->c_param.g[i11];
      d22 += dt * d35 * d36;
      d37 = d35 * d25;
      i13 = 3 * i11 + 1;
      d28 = b_SD->pd->c_param.Jinv[i13];
      d17 += d28 * d29;
      d20 += dt * d28 * d32;
      d35 = S[i13];
      d23 += dt * d35 * d36;
      d37 += d35 * d26;
      i14 = 3 * i11 + 2;
      d28 = b_SD->pd->c_param.Jinv[i14];
      d18 += d28 * d29;
      d21 += dt * d28 * d32;
      d35 = S[i14];
      d24 += dt * d35 * d36;
      d37 += d35 * d27;
      d39 = C_total_a[i11];
      c_w_exp_tilde[i11] =
          ((w_exp_tilde[i11] * d25 + w_exp_tilde[i11 + 3] * d26) +
           w_exp_tilde[i11 + 6] * d27) +
          dt * ((C_total_a_tmp_tmp / d39 * sens_input->board_accel.meas[i11] +
                 b_C_total_a_tmp_tmp / d39 * sens_input->mti_accel.meas[i11]) +
                c_C_total_a_tmp_tmp / d39 * sens_input->ad_accel.meas[i11]);
      b_S[i11] = d37;
    }
    memset(&c_q[0], 0, sizeof(double) << 2);
    d30 = c_q[0];
    d31 = c_q[1];
    d33 = c_q[2];
    d34 = c_q[3];
    for (i12 = 0; i12 < 4; i12++) {
      double d38;
      int q_tmp;
      q_tmp = i12 << 2;
      d38 = b_dv1[i12];
      d30 += b_q[q_tmp] * d38;
      d31 += b_q[q_tmp + 1] * d38;
      d33 += b_q[q_tmp + 2] * d38;
      d34 += b_q[q_tmp + 3] * d38;
    }
    double W_dt_tmp;
    double b_W_dt_tmp;
    double c_W_dt_tmp;
    double d_W_dt_tmp;
    double e_W_dt_tmp;
    double f_W_dt_tmp;
    x_pred[0] = d30;
    x_pred[1] = d31;
    x_pred[2] = d33;
    x_pred[3] = d34;
    x_pred[4] = d16 + d19;
    x_pred[7] = c_w_exp_tilde[0] + d22;
    x_pred[5] = d17 + d20;
    x_pred[8] = c_w_exp_tilde[1] + d23;
    x_pred[6] = d18 + d21;
    x_pred[9] = c_w_exp_tilde[2] + d24;
    x_pred[10] = x[10] + dt * b_S[0];
    memset(&F[0], 0, 121U * sizeof(double));
    k_a = 0.5 * dt;
    W_dt[0] = 0.0;
    W_dt_tmp = k_a * -x[4];
    W_dt[4] = W_dt_tmp;
    b_W_dt_tmp = k_a * -x[5];
    W_dt[8] = b_W_dt_tmp;
    c_W_dt_tmp = k_a * -x[6];
    W_dt[12] = c_W_dt_tmp;
    d_W_dt_tmp = k_a * x[4];
    W_dt[1] = d_W_dt_tmp;
    W_dt[5] = 0.0;
    e_W_dt_tmp = k_a * x[6];
    W_dt[9] = e_W_dt_tmp;
    W_dt[13] = b_W_dt_tmp;
    f_W_dt_tmp = k_a * x[5];
    W_dt[2] = f_W_dt_tmp;
    W_dt[6] = c_W_dt_tmp;
    W_dt[10] = 0.0;
    W_dt[14] = d_W_dt_tmp;
    W_dt[3] = e_W_dt_tmp;
    W_dt[7] = f_W_dt_tmp;
    W_dt[11] = W_dt_tmp;
    W_dt[15] = 0.0;
    memset(&c_b[0], 0, sizeof(double) << 4);
    for (i15 = 0; i15 < 4; i15++) {
      int i16;
      i16 = i15 << 2;
      for (i18 = 0; i18 < 4; i18++) {
        double d40;
        int b_tmp;
        d40 = W_dt[i18 + i16];
        b_tmp = i18 << 2;
        c_b[i16] += W_dt[b_tmp] * d40;
        c_b[i16 + 1] += W_dt[b_tmp + 1] * d40;
        c_b[i16 + 2] += W_dt[b_tmp + 2] * d40;
        c_b[i16 + 3] += W_dt[b_tmp + 3] * d40;
      }
    }
    for (i17 = 0; i17 < 16; i17++) {
      b_I[i17] = 0;
    }
    memset(&b_W_dt[0], 0, sizeof(double) << 4);
    memset(&d_b[0], 0, sizeof(double) << 4);
    for (b_k = 0; b_k < 4; b_k++) {
      double d41;
      double d42;
      double d43;
      double d44;
      double d45;
      double d46;
      double d47;
      double d48;
      int I_tmp;
      I_tmp = b_k << 2;
      b_I[b_k + I_tmp] = 1;
      d41 = b_W_dt[I_tmp];
      d42 = b_W_dt[I_tmp + 1];
      d43 = b_W_dt[I_tmp + 2];
      d44 = b_W_dt[I_tmp + 3];
      d45 = d_b[I_tmp];
      d46 = d_b[I_tmp + 1];
      d47 = d_b[I_tmp + 2];
      d48 = d_b[I_tmp + 3];
      for (i19 = 0; i19 < 4; i19++) {
        double d49;
        int g_W_dt_tmp;
        d49 = c_b[i19 + I_tmp];
        g_W_dt_tmp = i19 << 2;
        d41 += W_dt[g_W_dt_tmp] * d49;
        d45 += c_b[g_W_dt_tmp] * d49;
        d42 += W_dt[g_W_dt_tmp + 1] * d49;
        d46 += c_b[g_W_dt_tmp + 1] * d49;
        d43 += W_dt[g_W_dt_tmp + 2] * d49;
        d47 += c_b[g_W_dt_tmp + 2] * d49;
        d44 += W_dt[g_W_dt_tmp + 3] * d49;
        d48 += c_b[g_W_dt_tmp + 3] * d49;
      }
      d_b[I_tmp + 3] = d48;
      d_b[I_tmp + 2] = d47;
      d_b[I_tmp + 1] = d46;
      d_b[I_tmp] = d45;
      b_W_dt[I_tmp + 3] = d44;
      b_W_dt[I_tmp + 2] = d43;
      b_W_dt[I_tmp + 1] = d42;
      b_W_dt[I_tmp] = d41;
      F[11 * b_k] = ((((double)b_I[I_tmp] + W_dt[I_tmp]) + 0.5 * c_b[I_tmp]) +
                     0.16666666666666666 * d41) +
                    0.041666666666666664 * d45;
      F[11 * b_k + 1] =
          ((((double)b_I[I_tmp + 1] + W_dt[I_tmp + 1]) + 0.5 * c_b[I_tmp + 1]) +
           0.16666666666666666 * d42) +
          0.041666666666666664 * d46;
      F[11 * b_k + 2] =
          ((((double)b_I[I_tmp + 2] + W_dt[I_tmp + 2]) + 0.5 * c_b[I_tmp + 2]) +
           0.16666666666666666 * d43) +
          0.041666666666666664 * d47;
      F[11 * b_k + 3] =
          ((((double)b_I[I_tmp + 3] + W_dt[I_tmp + 3]) + 0.5 * c_b[I_tmp + 3]) +
           0.16666666666666666 * d44) +
          0.041666666666666664 * d48;
    }
    double e_a_tmp;
    double f_a_tmp;
    double g_a_tmp;
    double h_a_tmp;
    double i_a_tmp;
    double j_a_tmp;
    double k_a_tmp;
    e_a_tmp = k_a * q[0];
    l_a[0] = e_a_tmp;
    f_a_tmp = k_a * -q[1];
    l_a[4] = f_a_tmp;
    g_a_tmp = k_a * -q[2];
    l_a[8] = g_a_tmp;
    h_a_tmp = k_a * -q[3];
    l_a[12] = h_a_tmp;
    i_a_tmp = k_a * q[1];
    l_a[1] = i_a_tmp;
    l_a[5] = e_a_tmp;
    l_a[9] = h_a_tmp;
    j_a_tmp = k_a * q[2];
    l_a[13] = j_a_tmp;
    l_a[2] = j_a_tmp;
    k_a_tmp = k_a * q[3];
    l_a[6] = k_a_tmp;
    l_a[10] = e_a_tmp;
    l_a[14] = f_a_tmp;
    l_a[3] = k_a_tmp;
    l_a[7] = g_a_tmp;
    l_a[11] = i_a_tmp;
    l_a[15] = e_a_tmp;
    for (i20 = 0; i20 < 3; i20++) {
      int F_tmp;
      int b_F_tmp;
      F_tmp = (i20 + 1) << 2;
      b_F_tmp = 11 * (i20 + 4);
      F[b_F_tmp] = l_a[F_tmp];
      F[b_F_tmp + 1] = l_a[F_tmp + 1];
      F[b_F_tmp + 2] = l_a[F_tmp + 2];
      F[b_F_tmp + 3] = l_a[F_tmp + 3];
    }
    m_a = 0.5 * b_SD->pd->d_param.c_aero * b_SD->pd->d_param.Cn_alpha;
    airdata_atmos(x[10], &i_expl_temp, &t1_density, &j_expl_temp, &k_expl_temp,
                  &l_expl_temp);
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
    memset(&b_n_tilde[0], 0, 9U * sizeof(double));
    for (i21 = 0; i21 < 3; i21++) {
      double d50;
      int c_n_tilde_tmp;
      int d_n_tilde_tmp;
      d50 = b_n_tilde[3 * i21];
      c_n_tilde_tmp = 3 * i21 + 1;
      d_n_tilde_tmp = 3 * i21 + 2;
      for (i23 = 0; i23 < 3; i23++) {
        double d51;
        d51 = n_tilde[i23 + 3 * i21];
        d50 += n_tilde[3 * i23] * d51;
        b_n_tilde[c_n_tilde_tmp] += n_tilde[3 * i23 + 1] * d51;
        b_n_tilde[d_n_tilde_tmp] += n_tilde[3 * i23 + 2] * d51;
      }
      b_n_tilde[3 * i21] = d50;
    }
    for (i22 = 0; i22 < 9; i22++) {
      w_exp_tilde[i22] = ((double)w_exp_tilde_tmp[i22] - b_a * n_tilde[i22]) +
                         (1.0 - b_x) * b_n_tilde[i22];
    }
    memset(&b_dv[0], 0, 9U * sizeof(double));
    for (i24 = 0; i24 < 3; i24++) {
      double d52;
      int i26;
      int i27;
      d52 = b_dv[3 * i24];
      i26 = 3 * i24 + 1;
      i27 = 3 * i24 + 2;
      for (i29 = 0; i29 < 3; i29++) {
        double d54;
        d54 = w_exp_tilde[i29 + 3 * i24];
        d52 += b_SD->pd->d_param.Jinv[3 * i29] * d54;
        b_dv[i26] += b_SD->pd->d_param.Jinv[3 * i29 + 1] * d54;
        b_dv[i27] += b_SD->pd->d_param.Jinv[3 * i29 + 2] * d54;
        F[(i29 + 11 * (i24 + 4)) + 4] = 0.0;
      }
      b_dv[3 * i24] = d52;
    }
    for (i25 = 0; i25 < 3; i25++) {
      int F_tmp_tmp;
      F_tmp_tmp = 11 * (i25 + 4);
      for (i28 = 0; i28 < 3; i28++) {
        double d53;
        d53 = b_SD->pd->d_param.J[i28 + 3 * i25];
        F[F_tmp_tmp + 4] += b_dv[3 * i28] * d53;
        F[F_tmp_tmp + 5] += b_dv[3 * i28 + 1] * d53;
        F[F_tmp_tmp + 6] += b_dv[3 * i28 + 2] * d53;
      }
    }
    b_dv[1] = t1_density * (m_a * x[9]);
    b_dv[4] = 0.0;
    b_dv[7] = t1_density * (m_a * x[7]);
    b_dv[2] = t1_density * (m_a * -x[8]);
    b_dv[5] = t1_density * (m_a * -x[7]);
    b_dv[8] = 0.0;
    c_x = 0.0;
    for (i30 = 0; i30 < 3; i30++) {
      double d55;
      double d56;
      double d57;
      int c_F_tmp;
      b_dv[3 * i30] = 0.0;
      c_F_tmp = 11 * (i30 + 7);
      d55 = 0.0;
      d56 = 0.0;
      d57 = 0.0;
      for (i31 = 0; i31 < 3; i31++) {
        double d58;
        d58 = b_dv[i31 + 3 * i30];
        d55 += dt * b_SD->pd->d_param.Jinv[3 * i31] * d58;
        d56 += dt * b_SD->pd->d_param.Jinv[3 * i31 + 1] * d58;
        d57 += dt * b_SD->pd->d_param.Jinv[3 * i31 + 2] * d58;
      }
      F[c_F_tmp + 6] = d57;
      F[c_F_tmp + 5] = d56;
      F[c_F_tmp + 4] = d55;
      c_x += x[i30 + 1] * b_SD->pd->d_param.g[i30];
    }
    d_x[0] = x[2] * b_SD->pd->d_param.g[2] - b_SD->pd->d_param.g[1] * x[3];
    d_x[1] = b_SD->pd->d_param.g[0] * x[3] - x[1] * b_SD->pd->d_param.g[2];
    d_x[2] = x[1] * b_SD->pd->d_param.g[1] - b_SD->pd->d_param.g[0] * x[2];
    dv4[0] = 0.0;
    dv4[3] = x[0] * -b_SD->pd->d_param.g[2];
    dv4[6] = x[0] * b_SD->pd->d_param.g[1];
    dv4[1] = x[0] * b_SD->pd->d_param.g[2];
    dv4[4] = 0.0;
    dv4[7] = x[0] * -b_SD->pd->d_param.g[0];
    dv4[2] = x[0] * -b_SD->pd->d_param.g[1];
    dv4[5] = x[0] * b_SD->pd->d_param.g[0];
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
    n_a = 0.5 * (dt * dt);
    memset(&c_skewed_exp_w_tmp[0], 0, 9U * sizeof(double));
    memset(&b_dv[0], 0, 9U * sizeof(double));
    r_q_tmp[0] = x[0];
    b_r_q_tmp = 0.0;
    for (i32 = 0; i32 < 3; i32++) {
      double d59;
      double d60;
      double d_F_tmp;
      int e_F_tmp;
      int f_F_tmp;
      int g_F_tmp;
      F[i32 + 7] = dt * (2.0 * (x[0] * b_SD->pd->d_param.g[i32] - d_x[i32]));
      d_F_tmp = x[i32 + 1];
      e_F_tmp = 11 * (i32 + 1);
      F[e_F_tmp + 7] =
          dt * (2.0 * (((c_x * b_b[3 * i32] + x[1] * b_SD->pd->d_param.g[i32]) -
                        b_SD->pd->d_param.g[0] * d_F_tmp) +
                       dv4[3 * i32]));
      f_F_tmp = 3 * i32 + 1;
      F[e_F_tmp + 8] =
          dt * (2.0 * (((c_x * b_b[f_F_tmp] + x[2] * b_SD->pd->d_param.g[i32]) -
                        b_SD->pd->d_param.g[1] * d_F_tmp) +
                       dv4[f_F_tmp]));
      g_F_tmp = 3 * i32 + 2;
      F[e_F_tmp + 9] =
          dt * (2.0 * (((c_x * b_b[g_F_tmp] + x[3] * b_SD->pd->d_param.g[i32]) -
                        b_SD->pd->d_param.g[2] * d_F_tmp) +
                       dv4[g_F_tmp]));
      d59 = c_skewed_exp_w_tmp[3 * i32];
      d60 = b_dv[3 * i32];
      for (i34 = 0; i34 < 3; i34++) {
        double d61;
        double d62;
        int b_skewed_exp_w_tmp_tmp;
        int i35;
        int skewed_exp_w_tmp_tmp;
        i35 = i34 + 3 * i32;
        d61 = b_skewed_exp_w_tmp[i35];
        d62 = skewed_exp_w_tmp[i35];
        d59 += skewed_exp_w_tmp[3 * i34] * d61;
        d60 += 2.0 * b_skewed_exp_w_tmp[3 * i34] * d62;
        skewed_exp_w_tmp_tmp = 3 * i34 + 1;
        c_skewed_exp_w_tmp[f_F_tmp] +=
            skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d61;
        b_dv[f_F_tmp] += 2.0 * b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d62;
        b_skewed_exp_w_tmp_tmp = 3 * i34 + 2;
        c_skewed_exp_w_tmp[g_F_tmp] +=
            skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d61;
        b_dv[g_F_tmp] += 2.0 * b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d62;
      }
      int h_F_tmp;
      int i_F_tmp;
      b_dv[3 * i32] = d60;
      c_skewed_exp_w_tmp[3 * i32] = d59;
      h_F_tmp = 11 * (i32 + 4);
      F[h_F_tmp + 7] = dt * skewed_exp_w_tmp[3 * i32] + n_a * (d59 - d60);
      i_F_tmp = 11 * (i32 + 7);
      F[i_F_tmp + 7] = w_exp_tilde[3 * i32];
      F[h_F_tmp + 8] = dt * skewed_exp_w_tmp[f_F_tmp] +
                       n_a * (c_skewed_exp_w_tmp[f_F_tmp] - b_dv[f_F_tmp]);
      F[i_F_tmp + 8] = w_exp_tilde[f_F_tmp];
      F[h_F_tmp + 9] = dt * skewed_exp_w_tmp[g_F_tmp] +
                       n_a * (c_skewed_exp_w_tmp[g_F_tmp] - b_dv[g_F_tmp]);
      F[i_F_tmp + 9] = w_exp_tilde[g_F_tmp];
      r_q_tmp[i32 + 1] = -d_F_tmp;
      b_r_q_tmp += -d_F_tmp * x[i32 + 7];
    }
    c_r_q_tmp[0] = r_q_tmp[2] * x[9] - r_q_tmp[3] * x[8];
    c_r_q_tmp[1] = r_q_tmp[3] * x[7] - r_q_tmp[1] * x[9];
    c_r_q_tmp[2] = r_q_tmp[1] * x[8] - r_q_tmp[2] * x[7];
    for (i33 = 0; i33 < 3; i33++) {
      double b_dt_tmp;
      double dt_tmp;
      int c_dt_tmp;
      int d_dt_tmp;
      int e_dt_tmp;
      dt_tmp = x[i33 + 7];
      d_dt[i33] = dt * (2.0 * (r_q_tmp[0] * dt_tmp - c_r_q_tmp[i33]));
      b_dt_tmp = r_q_tmp[i33 + 1];
      c_dt_tmp = 3 * (i33 + 1);
      d_dt[c_dt_tmp] =
          dt * (2.0 * (((b_r_q_tmp * b_b[3 * i33] + r_q_tmp[1] * dt_tmp) -
                        x[7] * b_dt_tmp) +
                       r_q_tmp[0] * skewed_exp_w_tmp[3 * i33]));
      d_dt_tmp = 3 * i33 + 1;
      d_dt[c_dt_tmp + 1] =
          dt * (2.0 * (((b_r_q_tmp * b_b[d_dt_tmp] + r_q_tmp[2] * dt_tmp) -
                        x[8] * b_dt_tmp) +
                       r_q_tmp[0] * skewed_exp_w_tmp[d_dt_tmp]));
      e_dt_tmp = 3 * i33 + 2;
      d_dt[c_dt_tmp + 2] =
          dt * (2.0 * (((b_r_q_tmp * b_b[e_dt_tmp] + r_q_tmp[3] * dt_tmp) -
                        x[9] * b_dt_tmp) +
                       r_q_tmp[0] * skewed_exp_w_tmp[e_dt_tmp]));
    }
    double p_a;
    F[10] = d_dt[0];
    F[21] = d_dt[3];
    F[32] = d_dt[6];
    F[43] = d_dt[9];
    o_a = r_q_tmp[0] * r_q_tmp[0] -
          ((r_q_tmp[1] * r_q_tmp[1] + r_q_tmp[2] * r_q_tmp[2]) +
           r_q_tmp[3] * r_q_tmp[3]);
    p_a = 2.0 * r_q_tmp[0];
    b_dv[0] = 0.0;
    b_dv[3] = p_a * -r_q_tmp[3];
    b_dv[6] = p_a * r_q_tmp[2];
    b_dv[1] = p_a * r_q_tmp[3];
    b_dv[4] = 0.0;
    b_dv[7] = p_a * -r_q_tmp[1];
    b_dv[2] = p_a * -r_q_tmp[2];
    b_dv[5] = p_a * r_q_tmp[1];
    b_dv[8] = 0.0;
    for (i36 = 0; i36 < 3; i36++) {
      F[11 * (i36 + 7) + 10] =
          dt * ((o_a * b_b[3 * i36] + 2.0 * r_q_tmp[1] * r_q_tmp[i36 + 1]) -
                b_dv[3 * i36]);
    }
    F[120] = 1.0;
    memset(&b_F[0], 0, 121U * sizeof(double));
    for (i37 = 0; i37 < 11; i37++) {
      for (i38 = 0; i38 < 11; i38++) {
        double d63;
        d63 = P[i38 + 11 * i37];
        for (i41 = 0; i41 < 11; i41++) {
          int j_F_tmp;
          j_F_tmp = i41 + 11 * i37;
          b_F[j_F_tmp] += F[i41 + 11 * i38] * d63;
        }
      }
    }
    for (i39 = 0; i39 < 11; i39++) {
      for (i40 = 0; i40 < 11; i40++) {
        double d64;
        d64 = 0.0;
        for (i43 = 0; i43 < 11; i43++) {
          d64 += b_F[i39 + 11 * i43] * F[i40 + 11 * i43];
        }
        int c_P_pred_tmp;
        c_P_pred_tmp = i39 + 11 * i40;
        P_pred[c_P_pred_tmp] = d64 + Q[c_P_pred_tmp];
      }
    }
    for (i42 = 0; i42 < 3; i42++) {
      int P_pred_tmp;
      int b_P_pred_tmp;
      int d_P_pred_tmp;
      P_pred_tmp = 11 * (i42 + 4);
      b_P_pred[3 * i42] = P_pred[P_pred_tmp + 4] + R[3 * i42];
      b_P_pred_tmp = 3 * i42 + 1;
      b_P_pred[b_P_pred_tmp] = P_pred[P_pred_tmp + 5] + R[b_P_pred_tmp];
      d_P_pred_tmp = 3 * i42 + 2;
      b_P_pred[d_P_pred_tmp] = P_pred[P_pred_tmp + 6] + R[d_P_pred_tmp];
    }
    inv(b_P_pred, b_dv);
    memset(&K[0], 0, 33U * sizeof(double));
    for (i44 = 0; i44 < 3; i44++) {
      for (i45 = 0; i45 < 3; i45++) {
        double d65;
        d65 = b_dv[i45 + 3 * i44];
        for (i46 = 0; i46 < 11; i46++) {
          int K_tmp;
          K_tmp = i46 + 11 * i44;
          K[K_tmp] += P_pred[i46 + 11 * (i45 + 4)] * d65;
        }
      }
    }
    memset(&c_I[0], 0, 121U * sizeof(signed char));
    for (c_k = 0; c_k < 11; c_k++) {
      c_I[c_k + 11 * c_k] = 1;
    }
    for (i47 = 0; i47 < 44; i47++) {
      E[i47] = c_I[i47];
    }
    for (i48 = 0; i48 < 33; i48++) {
      E[i48 + 44] = (double)c_I[i48 + 44] - K[i48];
    }
    for (i49 = 0; i49 < 44; i49++) {
      E[i49 + 77] = c_I[i49 + 77];
    }
    memset(&b_E[0], 0, 121U * sizeof(double));
    for (i50 = 0; i50 < 11; i50++) {
      for (i51 = 0; i51 < 11; i51++) {
        double d66;
        d66 = P_pred[i51 + 11 * i50];
        for (i53 = 0; i53 < 11; i53++) {
          int E_tmp;
          E_tmp = i53 + 11 * i50;
          b_E[E_tmp] += E[i53 + 11 * i51] * d66;
        }
      }
    }
    memset(&b_K[0], 0, 33U * sizeof(double));
    for (i52 = 0; i52 < 3; i52++) {
      for (i54 = 0; i54 < 3; i54++) {
        double d67;
        d67 = R[i54 + 3 * i52];
        for (i55 = 0; i55 < 11; i55++) {
          int b_K_tmp;
          b_K_tmp = i55 + 11 * i52;
          b_K[b_K_tmp] += K[i55 + 11 * i54] * d67;
        }
      }
    }
    memset(&c_E[0], 0, 121U * sizeof(double));
    memset(&c_K[0], 0, 121U * sizeof(double));
    for (i56 = 0; i56 < 11; i56++) {
      for (i57 = 0; i57 < 11; i57++) {
        double d68;
        d68 = E[i56 + 11 * i57];
        for (i60 = 0; i60 < 11; i60++) {
          int b_E_tmp;
          b_E_tmp = i60 + 11 * i56;
          c_E[b_E_tmp] += b_E[i60 + 11 * i57] * d68;
        }
      }
      for (i59 = 0; i59 < 3; i59++) {
        double d70;
        d70 = K[i56 + 11 * i59];
        for (i61 = 0; i61 < 11; i61++) {
          int c_K_tmp;
          c_K_tmp = i61 + 11 * i56;
          c_K[c_K_tmp] += b_K[i61 + 11 * i59] * d70;
        }
      }
    }
    for (i58 = 0; i58 < 121; i58++) {
      P_ret[i58] = c_E[i58] + c_K[i58];
    }
    double d69;
    double d71;
    d69 = d1 / d3;
    d71 = d2 / d3;
    d72 = ((d1 / d4 * (sens_input->board_gyro.meas[0] - bias->board_gyro[0]) +
            d2 / d4 * (sens_input->mti_gyro.meas[0] - bias->mti_gyro[0])) +
           C_ad_w_idx_0 * (sens_input->ad_gyro.meas[0] - bias->ad_gyro[0])) -
          x_pred[4];
    d73 = ((d69 * (sens_input->board_gyro.meas[1] - bias->board_gyro[1]) +
            d71 * (sens_input->mti_gyro.meas[1] - bias->mti_gyro[1])) +
           d * (sens_input->ad_gyro.meas[1] - bias->ad_gyro[1])) -
          x_pred[5];
    d74 = ((d69 * (sens_input->board_gyro.meas[2] - bias->board_gyro[2]) +
            d71 * (sens_input->mti_gyro.meas[2] - bias->mti_gyro[2])) +
           d * (sens_input->ad_gyro.meas[2] - bias->ad_gyro[2])) -
          x_pred[6];
    for (i62 = 0; i62 < 11; i62++) {
      x_ret[i62] = x_pred[i62] +
                   ((K[i62] * d72 + K[i62 + 11] * d73) + K[i62 + 22] * d74);
    }
    c_scale = 3.3121686421112381E-170;
    c_absxk = fabs(x_ret[0]);
    if (c_absxk > 3.3121686421112381E-170) {
      b_q_mag = 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / 3.3121686421112381E-170;
      b_q_mag = c_t * c_t;
    }
    c_absxk = fabs(x_ret[1]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = b_q_mag * c_t * c_t + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    c_absxk = fabs(x_ret[2]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = b_q_mag * c_t * c_t + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    c_absxk = fabs(x_ret[3]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = b_q_mag * c_t * c_t + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    b_q_mag = c_scale * sqrt(b_q_mag);
    x_ret[0] /= b_q_mag;
    x_ret[1] /= b_q_mag;
    x_ret[2] /= b_q_mag;
    x_ret[3] /= b_q_mag;
    if (sens_input->board_baro.status) {
      memcpy(&b_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&b_P_ret[0], &P_ret[0], 121U * sizeof(double));
      memcpy(&e_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&d_P_ret[0], &P_ret[0], 121U * sizeof(double));
      b_ekf_correct(e_x_ret, d_P_ret, sens_input->board_baro.meas,
                    bias->board_baro, x_ret, P_ret);
    }
    if (sens_input->mti_baro.status) {
      memcpy(&c_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&c_P_ret[0], &P_ret[0], 121U * sizeof(double));
      memcpy(&g_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&f_P_ret[0], &P_ret[0], 121U * sizeof(double));
      b_ekf_correct(g_x_ret, f_P_ret, sens_input->mti_baro.meas, bias->mti_baro,
                    x_ret, P_ret);
    }
    if (sens_input->board_mag.status) {
      memcpy(&d_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&e_P_ret[0], &P_ret[0], 121U * sizeof(double));
      memcpy(&h_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&h_P_ret[0], &P_ret[0], 121U * sizeof(double));
      ekf_correct(h_x_ret, h_P_ret, sens_input->board_mag.meas,
                  bias->board_mag_earth, b_b, x_ret, P_ret);
    }
    if (sens_input->mti_mag.status) {
      memcpy(&f_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&g_P_ret[0], &P_ret[0], 121U * sizeof(double));
      memcpy(&i_x_ret[0], &x_ret[0], 11U * sizeof(double));
      memcpy(&i_P_ret[0], &P_ret[0], 121U * sizeof(double));
      ekf_correct(i_x_ret, i_P_ret, sens_input->mti_mag.meas,
                  bias->mti_mag_earth, b_b, x_ret, P_ret);
    }
  }
  airdata->pressure = airdata_atmos(x_ret[10], &airdata->temperature,
                                    &airdata->density, &airdata->sonic_speed,
                                    &airdata->mach, &airdata->dynamic_pressure);
  airspeed = b_norm(&x_ret[7]);
  airdata->mach = airspeed / airdata->sonic_speed;
  airdata->dynamic_pressure = 0.5 * airdata->density * (airspeed * airspeed);
  roll_state[0] = atan2(
      2.0 * (x_ret[2] * x_ret[3] + x_ret[0] * x_ret[1]),
      ((x_ret[0] * x_ret[0] - x_ret[1] * x_ret[1]) - x_ret[2] * x_ret[2]) +
          x_ret[3] * x_ret[3]);
  roll_state[1] = x_ret[4];
}
