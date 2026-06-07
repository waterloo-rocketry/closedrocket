#include "GNC_codegen.h"
#include "GNC_codegen_types.h"
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

static void controller_codegen_entry_init(GNC_codegenStackData *SD);

static void dynamics_init(GNC_codegenStackData *SD);

static void dynamics_jacobian_init(GNC_codegenStackData *SD);

static void ekf_correct(const double x[11], const double P[121],
                        const double y[3], const double b[3], const double R[9],
                        double x_new[11], double P_new[121]);

static void inv(const double x[9], double y[9]);

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
  int b_k;
  int i;
  int i1;
  int i10;
  int i11;
  int i13;
  int i14;
  int i15;
  int i16;
  int i18;
  int i19;
  int i2;
  int i20;
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
    for (i2 = 0; i2 < 11; i2++) {
      d += H[i2] * P[i2 + i];
    }
    b_H[i1] = d;
    c_H += d * H[i1];
    i += 11;
  }
  c_b = 1.0 / (c_H + 100.0);
  memset(&b_P[0], 0, 11U * (sizeof(double)));
  i3 = 0;
  for (i4 = 0; i4 < 11; i4++) {
    for (i5 = 0; i5 < 11; i5++) {
      b_P[i5] += P[i5 + i3] * H[i4];
    }
    i3 += 11;
  }
  for (i6 = 0; i6 < 11; i6++) {
    K[i6] = b_P[i6] * c_b;
  }
  memset(&b_I[0], 0, 121U * (sizeof(signed char)));
  k = 0;
  for (b_k = 0; b_k < 11; b_k++) {
    b_I[k] = (signed char)1;
    k += 12;
  }
  i7 = 0;
  for (i8 = 0; i8 < 11; i8++) {
    for (i9 = 0; i9 < 11; i9++) {
      int E_tmp;
      E_tmp = i9 + i7;
      E[E_tmp] = ((double)b_I[E_tmp]) - (K[i9] * H[i8]);
    }
    i7 += 11;
  }
  memset(&b_E[0], 0, 121U * (sizeof(double)));
  i10 = 0;
  for (i11 = 0; i11 < 11; i11++) {
    int i12;
    i12 = 0;
    for (i13 = 0; i13 < 11; i13++) {
      double d1;
      d1 = P[i13 + i10];
      for (i16 = 0; i16 < 11; i16++) {
        int b_E_tmp;
        b_E_tmp = i16 + i10;
        b_E[b_E_tmp] += E[i16 + i12] * d1;
      }
      i12 += 11;
    }
    i10 += 11;
  }
  memset(&c_E[0], 0, 121U * (sizeof(double)));
  i14 = 0;
  for (i15 = 0; i15 < 11; i15++) {
    int i17;
    i17 = 0;
    for (i19 = 0; i19 < 11; i19++) {
      double d2;
      d2 = E[i17 + i15];
      for (i20 = 0; i20 < 11; i20++) {
        int c_E_tmp;
        c_E_tmp = i20 + i14;
        c_E[c_E_tmp] += b_E[i20 + i17] * d2;
      }
      b_K[i19 + i14] = (K[i19] * 100.0) * K[i15];
      i17 += 11;
    }
    i14 += 11;
  }
  for (i18 = 0; i18 < 121; i18++) {
    P_new[i18] = c_E[i18] + b_K[i18];
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

static void controller_codegen_entry_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->param.J[i] = dv[i];
    SD->pd->param.Jinv[i] = dv1[i];
  }
  SD->pd->param.c_aero = -0.016182736457722724;
  SD->pd->param.c_canard = 0.00061367415999999994;
  SD->pd->param.elevation = 420.0;
  SD->pd->param.g[0] = -9.81;
  SD->pd->param.g[1] = 0.0;
  SD->pd->param.g[2] = 0.0;
}

static void dynamics_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->c_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->c_param.J[i] = dv[i];
    SD->pd->c_param.Jinv[i] = dv1[i];
  }
  SD->pd->c_param.c_aero = -0.016182736457722724;
  SD->pd->c_param.c_canard = 0.00061367415999999994;
  SD->pd->c_param.elevation = 420.0;
  SD->pd->c_param.g[0] = -9.81;
  SD->pd->c_param.g[1] = 0.0;
  SD->pd->c_param.g[2] = 0.0;
}

static void dynamics_jacobian_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->d_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->d_param.J[i] = dv[i];
    SD->pd->d_param.Jinv[i] = dv1[i];
  }
  SD->pd->d_param.c_aero = -0.016182736457722724;
  SD->pd->d_param.c_canard = 0.00061367415999999994;
  SD->pd->d_param.elevation = 420.0;
  SD->pd->d_param.g[0] = -9.81;
  SD->pd->d_param.g[1] = 0.0;
  SD->pd->d_param.g[2] = 0.0;
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
  double d11;
  double d12;
  double d13;
  double d14;
  double d15;
  double d16;
  double q_mag;
  double scale;
  double t;
  int b_k;
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
  int i19;
  int i2;
  int i20;
  int i21;
  int i24;
  int i25;
  int i27;
  int i28;
  int i29;
  int i30;
  int i31;
  int i32;
  int i34;
  int i35;
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
  int i5;
  int i6;
  int i7;
  int i9;
  int k;
  signed char b_I[121];
  a = (x[0] * x[0]) - (((x[1] * x[1]) + (x[2] * x[2])) + (x[3] * x[3]));
  b_a = 2.0 * x[0];
  memset(&H[0], 0, 33U * (sizeof(double)));
  b_x = ((b[0] * x[1]) + (b[1] * x[2])) + (b[2] * x[3]);
  c_x[0] = (x[2] * b[2]) - (b[1] * x[3]);
  c_x[1] = (b[0] * x[3]) - (x[1] * b[2]);
  c_x[2] = (x[1] * b[1]) - (b[0] * x[2]);
  b_dv[0] = 0.0;
  b_dv[3] = x[0] * (-b[2]);
  b_dv[6] = x[0] * b[1];
  b_dv[1] = x[0] * b[2];
  b_dv[4] = 0.0;
  b_dv[7] = x[0] * (-b[0]);
  b_dv[2] = x[0] * (-b[1]);
  b_dv[5] = x[0] * b[0];
  b_dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    double H_tmp;
    int b_H_tmp;
    int c_H_tmp;
    int d_H_tmp;
    H[i] = 2.0 * ((x[0] * b[i]) - c_x[i]);
    H_tmp = x[i + 1];
    b_H_tmp = 3 * (i + 1);
    H[b_H_tmp] =
        2.0 *
        ((((b_x * ((double)iv[3 * i])) + (x[1] * b[i])) - (b[0] * H_tmp)) +
         b_dv[3 * i]);
    c_H_tmp = (3 * i) + 1;
    H[b_H_tmp + 1] =
        2.0 *
        ((((b_x * ((double)iv[c_H_tmp])) + (x[2] * b[i])) - (b[1] * H_tmp)) +
         b_dv[c_H_tmp]);
    d_H_tmp = (3 * i) + 2;
    H[b_H_tmp + 2] =
        2.0 *
        ((((b_x * ((double)iv[d_H_tmp])) + (x[3] * b[i])) - (b[2] * H_tmp)) +
         b_dv[d_H_tmp]);
  }
  i1 = 0;
  for (i2 = 0; i2 < 3; i2++) {
    int i3;
    i3 = 0;
    for (i4 = 0; i4 < 11; i4++) {
      y_tmp[i4 + i1] = H[i3 + i2];
      i3 += 3;
    }
    i1 += 11;
  }
  memset(&b_H[0], 0, 33U * (sizeof(double)));
  i5 = 0;
  i6 = 0;
  for (i7 = 0; i7 < 11; i7++) {
    int i8;
    i8 = 0;
    for (i9 = 0; i9 < 11; i9++) {
      double d;
      d = P[i9 + i5];
      b_H[i6] += H[i8] * d;
      b_H[i6 + 1] += H[i8 + 1] * d;
      b_H[i6 + 2] += H[i8 + 2] * d;
      i8 += 3;
    }
    i5 += 11;
    i6 += 3;
  }
  memset(&b_P[0], 0, 33U * (sizeof(double)));
  for (i10 = 0; i10 < 3; i10++) {
    for (i11 = 0; i11 < 3; i11++) {
      double d1;
      d1 = 0.0;
      for (i13 = 0; i13 < 11; i13++) {
        d1 += b_H[i10 + (3 * i13)] * y_tmp[i13 + (11 * i11)];
      }
      int e_H_tmp;
      e_H_tmp = i10 + (3 * i11);
      c_H[e_H_tmp] = d1 + R[e_H_tmp];
    }
    for (i12 = 0; i12 < 11; i12++) {
      double d2;
      d2 = y_tmp[i12 + (11 * i10)];
      for (i16 = 0; i16 < 11; i16++) {
        int P_tmp;
        P_tmp = i16 + (11 * i10);
        b_P[P_tmp] += P[i16 + (11 * i12)] * d2;
      }
    }
  }
  inv(c_H, b_dv);
  memset(&K[0], 0, 33U * (sizeof(double)));
  i14 = 0;
  i15 = 0;
  for (i17 = 0; i17 < 3; i17++) {
    int i18;
    i18 = 0;
    for (i19 = 0; i19 < 3; i19++) {
      double d3;
      d3 = b_dv[i19 + i14];
      for (i20 = 0; i20 < 11; i20++) {
        int K_tmp;
        K_tmp = i20 + i15;
        K[K_tmp] += b_P[i20 + i18] * d3;
      }
      i18 += 11;
    }
    i14 += 3;
    i15 += 11;
  }
  memset(&b_I[0], 0, 121U * (sizeof(signed char)));
  k = 0;
  for (b_k = 0; b_k < 11; b_k++) {
    b_I[k] = (signed char)1;
    k += 12;
  }
  for (i21 = 0; i21 < 11; i21++) {
    double d4;
    double d5;
    double d6;
    int i22;
    int i23;
    i22 = 0;
    i23 = 0;
    d4 = K[i21];
    d5 = K[i21 + 11];
    d6 = K[i21 + 22];
    for (i27 = 0; i27 < 11; i27++) {
      int E_tmp;
      E_tmp = i23 + i21;
      E[E_tmp] = ((double)b_I[E_tmp]) -
                 (((d4 * H[i22]) + (d5 * H[i22 + 1])) + (d6 * H[i22 + 2]));
      i22 += 3;
      i23 += 11;
    }
  }
  memset(&b_E[0], 0, 121U * (sizeof(double)));
  i24 = 0;
  for (i25 = 0; i25 < 11; i25++) {
    int i26;
    i26 = 0;
    for (i28 = 0; i28 < 11; i28++) {
      double d7;
      d7 = P[i28 + i24];
      for (i32 = 0; i32 < 11; i32++) {
        int b_E_tmp;
        b_E_tmp = i32 + i24;
        b_E[b_E_tmp] += E[i32 + i26] * d7;
      }
      i26 += 11;
    }
    i24 += 11;
  }
  memset(&b_K[0], 0, 33U * (sizeof(double)));
  i29 = 0;
  i30 = 0;
  for (i31 = 0; i31 < 3; i31++) {
    int i33;
    i33 = 0;
    for (i34 = 0; i34 < 3; i34++) {
      double d8;
      d8 = R[i34 + i29];
      for (i36 = 0; i36 < 11; i36++) {
        int b_K_tmp;
        b_K_tmp = i36 + i30;
        b_K[b_K_tmp] += K[i36 + i33] * d8;
      }
      i33 += 11;
    }
    i29 += 3;
    i30 += 11;
  }
  memset(&c_E[0], 0, 121U * (sizeof(double)));
  memset(&c_K[0], 0, 121U * (sizeof(double)));
  for (i35 = 0; i35 < 11; i35++) {
    for (i37 = 0; i37 < 11; i37++) {
      double d9;
      d9 = E[i35 + (11 * i37)];
      for (i40 = 0; i40 < 11; i40++) {
        int c_E_tmp;
        c_E_tmp = i40 + (11 * i35);
        c_E[c_E_tmp] += b_E[i40 + (11 * i37)] * d9;
      }
    }
    for (i39 = 0; i39 < 3; i39++) {
      double d10;
      d10 = K[i35 + (11 * i39)];
      for (i42 = 0; i42 < 11; i42++) {
        int c_K_tmp;
        c_K_tmp = i42 + (11 * i35);
        c_K[c_K_tmp] += b_K[i42 + (11 * i39)] * d10;
      }
    }
  }
  for (i38 = 0; i38 < 121; i38++) {
    P_new[i38] = c_E[i38] + c_K[i38];
  }
  for (i41 = 0; i41 < 3; i41++) {
    double a_tmp;
    int b_a_tmp;
    int c_a_tmp;
    a_tmp = x[i41 + 1];
    c_a[3 * i41] = (a * ((double)iv[3 * i41])) + ((2.0 * x[1]) * a_tmp);
    b_a_tmp = (3 * i41) + 1;
    c_a[b_a_tmp] = (a * ((double)iv[b_a_tmp])) + ((2.0 * x[2]) * a_tmp);
    c_a_tmp = (3 * i41) + 2;
    c_a[c_a_tmp] = (a * ((double)iv[c_a_tmp])) + ((2.0 * x[3]) * a_tmp);
  }
  b_dv[0] = 0.0;
  b_dv[3] = b_a * (-x[3]);
  b_dv[6] = b_a * x[2];
  b_dv[1] = b_a * x[3];
  b_dv[4] = 0.0;
  b_dv[7] = b_a * (-x[1]);
  b_dv[2] = b_a * (-x[2]);
  b_dv[5] = b_a * x[1];
  b_dv[8] = 0.0;
  for (i43 = 0; i43 < 9; i43++) {
    c_a[i43] -= b_dv[i43];
  }
  d11 = b[0];
  d12 = b[1];
  d13 = b[2];
  for (i44 = 0; i44 < 3; i44++) {
    b_y[i44] = y[i44] - (((c_a[i44] * d11) + (c_a[i44 + 3] * d12)) +
                         (c_a[i44 + 6] * d13));
  }
  d14 = b_y[0];
  d15 = b_y[1];
  d16 = b_y[2];
  for (i45 = 0; i45 < 11; i45++) {
    x_new[i45] =
        x[i45] + (((K[i45] * d14) + (K[i45 + 11] * d15)) + (K[i45 + 22] * d16));
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
  memcpy(&b_x[0], &x[0], 9U * (sizeof(double)));
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
  t3 = ((b_x[1] * b_x[5]) - b_x[2]) / b_x[8];
  t2 = (-(b_x[1] + (b_x[7] * t3))) / b_x[4];
  y[p1] = ((1.0 - (b_x[3] * t2)) - (b_x[6] * t3)) / b_x[0];
  y[p1 + 1] = t2;
  y[p1 + 2] = t3;
  t3 = (-b_x[5]) / b_x[8];
  t2 = (1.0 - (b_x[7] * t3)) / b_x[4];
  y[p2] = (-((b_x[3] * t2) + (b_x[6] * t3))) / b_x[0];
  y[p2 + 1] = t2;
  y[p2 + 2] = t3;
  t3 = 1.0 / b_x[8];
  t2 = ((-b_x[7]) * t3) / b_x[4];
  y[p3] = (-((b_x[3] * t2) + (b_x[6] * t3))) / b_x[0];
  y[p3 + 1] = t2;
  y[p3 + 2] = t3;
}

static void pad_filter_init(GNC_codegenStackData *SD) {
  int i;
  SD->pd->b_param.Cn_alpha = 10.0;
  for (i = 0; i < 9; i++) {
    SD->pd->b_param.J[i] = dv[i];
    SD->pd->b_param.Jinv[i] = dv1[i];
  }
  SD->pd->b_param.c_aero = -0.016182736457722724;
  SD->pd->b_param.c_canard = 0.00061367415999999994;
  SD->pd->b_param.elevation = 420.0;
  SD->pd->b_param.g[0] = -9.81;
  SD->pd->b_param.g[1] = 0.0;
  SD->pd->b_param.g[2] = 0.0;
}

void GNC_codegen_initialize(GNC_codegenStackData *SD) {
  controller_codegen_entry_init(SD);
  pad_filter_init(SD);
  dynamics_init(SD);
  dynamics_jacobian_init(SD);
}

void GNC_codegen_terminate(void) {}

void controller_codegen_entry(GNC_codegenStackData *SD, double b_time,
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
  double d2;
  double d3;
  double d4;
  double d5;
  double d7;
  double d8;
  double d9;
  double pdyn_params;
  double r_idx_0;
  double w_dot;
  double x;
  double x_tmp;
  int i;
  int i1;
  int i2;
  int i3;
  if ((b_time >= 22.0) && (b_time < 27.0)) {
    *r = 0.5;
  } else if ((b_time >= 27.0) && (b_time < 32.0)) {
    *r = -0.5;
  } else if ((b_time >= 32.0) && (b_time < 39.0)) {
    *r = 0.5;
  } else {
    *r = 0.0;
  }
  pdyn_params = pdyn * SD->pd->param.c_canard;
  c_delta = delta / 2.0;
  if (fabs(c_delta) < 0.005) {
    c_delta = 0.0;
  }
  c_delta = (0.75 * ctrl_mem_in->d_old) + (0.25 * c_delta);
  w_dot = (0.75 * ctrl_mem_in->w_dot_old) +
          ((0.25 * (xR[1] - ctrl_mem_in->w_old)) / dt_ctrl);
  r_idx_0 = pdyn_params * c_delta;
  P[0] = ctrl_mem_in->P_minus[0] + 1.0E-5;
  P[1] = ctrl_mem_in->P_minus[1];
  P[2] = ctrl_mem_in->P_minus[2];
  P[3] = ctrl_mem_in->P_minus[3] + 1.0E-9;
  memset(&b_r[0], 0, (sizeof(double)) << 1);
  d = r_idx_0 * (ctrl_mem_in->P_minus[0] + 1.0E-5);
  d1 = pdyn_params * (ctrl_mem_in->P_minus[3] + 1.0E-9);
  c_r = (((b_r[0] + d) + (pdyn_params * ctrl_mem_in->P_minus[1])) * r_idx_0) +
        (((b_r[1] + (r_idx_0 * ctrl_mem_in->P_minus[2])) + d1) * pdyn_params);
  K[0] = (d + (ctrl_mem_in->P_minus[2] * pdyn_params)) / (c_r + 1.0);
  K[1] = ((ctrl_mem_in->P_minus[1] * r_idx_0) + d1) / (c_r + 1.0);
  b = w_dot - ((r_idx_0 * ctrl_mem_in->coeffs[0]) +
               (pdyn_params * ctrl_mem_in->coeffs[1]));
  ctrl_mem_out->coeffs[0] = ctrl_mem_in->coeffs[0] + (K[0] * b);
  ctrl_mem_out->coeffs[1] = ctrl_mem_in->coeffs[1] + (K[1] * b);
  b_dv[0] = 1.0 - (K[0] * r_idx_0);
  b_dv[1] = 0.0 - (K[1] * r_idx_0);
  b_dv[2] = 0.0 - (K[0] * pdyn_params);
  b_dv[3] = 1.0 - (K[1] * pdyn_params);
  memset(&b_dv1[0], 0, (sizeof(double)) << 2);
  i = 0;
  d2 = b_dv[0];
  d3 = b_dv[1];
  d4 = b_dv[2];
  d5 = b_dv[3];
  for (i1 = 0; i1 < 2; i1++) {
    double d6;
    d6 = P[i];
    b_dv1[i] += d2 * d6;
    b_dv1[i + 1] += d3 * d6;
    d6 = P[i + 1];
    b_dv1[i] += d4 * d6;
    b_dv1[i + 1] += d5 * d6;
    i += 2;
  }
  memset(&dv2[0], 0, (sizeof(double)) << 2);
  i2 = 0;
  d7 = b_dv1[0];
  d8 = b_dv1[1];
  d9 = b_dv1[2];
  d10 = b_dv1[3];
  for (i3 = 0; i3 < 2; i3++) {
    double d11;
    d11 = b_dv[i3];
    dv2[i2] += d7 * d11;
    dv2[i2 + 1] += d8 * d11;
    b_K[i2] = K[0] * K[i3];
    d11 = b_dv[i3 + 2];
    dv2[i2] += d9 * d11;
    dv2[i2 + 1] += d10 * d11;
    b_K[i2 + 1] = K[1] * K[i3];
    i2 += 2;
  }
  ctrl_mem_out->P_minus[0] = dv2[0] + b_K[0];
  ctrl_mem_out->P_minus[1] = dv2[1] + b_K[1];
  ctrl_mem_out->P_minus[2] = dv2[2] + b_K[2];
  ctrl_mem_out->P_minus[3] = dv2[3] + b_K[3];
  ctrl_mem_out->w_old = xR[1];
  ctrl_mem_out->d_old = c_delta;
  ctrl_mem_out->w_dot_old = w_dot;
  L_delta = (ctrl_mem_out->coeffs[0] * pdyn_params) / 2.0;
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
  *u = fmin(fmax(((K[0] * xR[0]) +
                  ((a * sqrt((2.0 * x) + (x_tmp + (blend * 20.0)))) * xR[1])) +
                     ((-K[0]) * (*r)),
                 -0.3490658503988659),
            0.3490658503988659);
  if (pdyn < 500.0) {
    *u = 0.0;
  }
}

void navigation_codegen_entry(GNC_codegenStackData *SD, double dt,
                              bool flight_phase, const double x[11],
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
  int c_k;
  int f_k;
  int h_k;
  int i;
  int i1;
  int i10;
  int i100;
  int i101;
  int i102;
  int i103;
  int i104;
  int i12;
  int i14;
  int i16;
  int i18;
  int i20;
  int i22;
  int i24;
  int i25;
  int i29;
  int i3;
  int i30;
  int i33;
  int i35;
  int i37;
  int i38;
  int i41;
  int i45;
  int i46;
  int i48;
  int i49;
  int i5;
  int i50;
  int i51;
  int i53;
  int i54;
  int i58;
  int i6;
  int i60;
  int i62;
  int i63;
  int i66;
  int i67;
  int i69;
  int i7;
  int i71;
  int i74;
  int i76;
  int i77;
  int i79;
  int i8;
  int i81;
  int i82;
  int i84;
  int i85;
  int i86;
  int i88;
  int i90;
  int i93;
  int i94;
  int i96;
  int i97;
  int i98;
  int i99;
  signed char c_I[121];
  memcpy(&P_ret[0], &P[0], 121U * (sizeof(double)));
  *bias_ret = *bias;
  *sens_filt_ret = *sens_filt;
  if (!((int)flight_phase)) {
    double ST[9];
    double f_a[9];
    double a[3];
    double b_absxk;
    double b_scale;
    double b_t;
    double board_baro_f;
    double d10;
    double d11;
    double d5;
    double d6;
    double d9;
    double d_a;
    double e_a;
    double mti_baro_f;
    double qw;
    double qy;
    double qz;
    double y;
    int i11;
    int i4;
    int i9;
    board_baro_f = sens_filt->board_baro_f;
    mti_baro_f = sens_filt->mti_baro_f;
    if (sens_input->board_accel.status) {
      sens_filt_ret->board_accel_f[0] =
          (0.0005 * sens_input->board_accel.meas[0]) +
          (0.9995 * sens_filt->board_accel_f[0]);
      sens_filt_ret->board_accel_f[1] =
          (0.0005 * sens_input->board_accel.meas[1]) +
          (0.9995 * sens_filt->board_accel_f[1]);
      sens_filt_ret->board_accel_f[2] =
          (0.0005 * sens_input->board_accel.meas[2]) +
          (0.9995 * sens_filt->board_accel_f[2]);
    }
    if (sens_input->board_gyro.status) {
      sens_filt_ret->board_gyro_f[0] =
          (0.0005 * sens_input->board_gyro.meas[0]) +
          (0.9995 * sens_filt->board_gyro_f[0]);
      sens_filt_ret->board_gyro_f[1] =
          (0.0005 * sens_input->board_gyro.meas[1]) +
          (0.9995 * sens_filt->board_gyro_f[1]);
      sens_filt_ret->board_gyro_f[2] =
          (0.0005 * sens_input->board_gyro.meas[2]) +
          (0.9995 * sens_filt->board_gyro_f[2]);
    }
    if (sens_input->mti_accel.status) {
      sens_filt_ret->mti_accel_f[0] = (0.0005 * sens_input->mti_accel.meas[0]) +
                                      (0.9995 * sens_filt->mti_accel_f[0]);
      sens_filt_ret->mti_accel_f[1] = (0.0005 * sens_input->mti_accel.meas[1]) +
                                      (0.9995 * sens_filt->mti_accel_f[1]);
      sens_filt_ret->mti_accel_f[2] = (0.0005 * sens_input->mti_accel.meas[2]) +
                                      (0.9995 * sens_filt->mti_accel_f[2]);
    }
    if (sens_input->mti_gyro.status) {
      sens_filt_ret->mti_gyro_f[0] = (0.0005 * sens_input->mti_gyro.meas[0]) +
                                     (0.9995 * sens_filt->mti_gyro_f[0]);
      sens_filt_ret->mti_gyro_f[1] = (0.0005 * sens_input->mti_gyro.meas[1]) +
                                     (0.9995 * sens_filt->mti_gyro_f[1]);
      sens_filt_ret->mti_gyro_f[2] = (0.0005 * sens_input->mti_gyro.meas[2]) +
                                     (0.9995 * sens_filt->mti_gyro_f[2]);
    }
    if (sens_input->ad_accel.status) {
      sens_filt_ret->ad_accel_f[0] = (0.0005 * sens_input->ad_accel.meas[0]) +
                                     (0.9995 * sens_filt->ad_accel_f[0]);
      sens_filt_ret->ad_accel_f[1] = (0.0005 * sens_input->ad_accel.meas[1]) +
                                     (0.9995 * sens_filt->ad_accel_f[1]);
      sens_filt_ret->ad_accel_f[2] = (0.0005 * sens_input->ad_accel.meas[2]) +
                                     (0.9995 * sens_filt->ad_accel_f[2]);
    }
    if (sens_input->ad_gyro.status) {
      sens_filt_ret->ad_gyro_f[0] = (0.0005 * sens_input->ad_gyro.meas[0]) +
                                    (0.9995 * sens_filt->ad_gyro_f[0]);
      sens_filt_ret->ad_gyro_f[1] = (0.0005 * sens_input->ad_gyro.meas[1]) +
                                    (0.9995 * sens_filt->ad_gyro_f[1]);
      sens_filt_ret->ad_gyro_f[2] = (0.0005 * sens_input->ad_gyro.meas[2]) +
                                    (0.9995 * sens_filt->ad_gyro_f[2]);
    }
    if (sens_input->board_baro.status) {
      board_baro_f = (0.0005 * sens_input->board_baro.meas) +
                     (0.9995 * sens_filt->board_baro_f);
    }
    if (sens_input->board_mag.status) {
      sens_filt_ret->board_mag_f[0] = (0.0005 * sens_input->board_mag.meas[0]) +
                                      (0.9995 * sens_filt->board_mag_f[0]);
      sens_filt_ret->board_mag_f[1] = (0.0005 * sens_input->board_mag.meas[1]) +
                                      (0.9995 * sens_filt->board_mag_f[1]);
      sens_filt_ret->board_mag_f[2] = (0.0005 * sens_input->board_mag.meas[2]) +
                                      (0.9995 * sens_filt->board_mag_f[2]);
    }
    if (sens_input->mti_baro.status) {
      mti_baro_f = (0.0005 * sens_input->mti_baro.meas) +
                   (0.9995 * sens_filt->mti_baro_f);
    }
    if (sens_input->mti_mag.status) {
      sens_filt_ret->mti_mag_f[0] = (0.0005 * sens_input->mti_mag.meas[0]) +
                                    (0.9995 * sens_filt->mti_mag_f[0]);
      sens_filt_ret->mti_mag_f[1] = (0.0005 * sens_input->mti_mag.meas[1]) +
                                    (0.9995 * sens_filt->mti_mag_f[1]);
      sens_filt_ret->mti_mag_f[2] = (0.0005 * sens_input->mti_mag.meas[2]) +
                                    (0.9995 * sens_filt->mti_mag_f[2]);
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
    qw = sqrt((0.5 * (a[0] / d5)) + 0.5);
    if (qw == 0.0) {
      qy = 1.0;
      qz = 0.0;
    } else {
      qy = (0.5 * (a[2] / d5)) / qw;
      qz = (-0.5 * (a[1] / d5)) / qw;
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
      y = ((y * b_t) * b_t) + 1.0;
      b_scale = b_absxk;
    } else {
      b_t = b_absxk / b_scale;
      y += b_t * b_t;
    }
    b_absxk = fabs(qz);
    if (b_absxk > b_scale) {
      b_t = b_scale / b_absxk;
      y = ((y * b_t) * b_t) + 1.0;
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
    x_ret[10] = SD->pd->b_param.elevation;
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
    d_a = (q[0] * q[0]) - (((q[1] * q[1]) + (q[2] * q[2])) + (d6 * d6));
    e_a = 2.0 * q[0];
    i4 = 0;
    for (i5 = 0; i5 < 3; i5++) {
      double a_tmp;
      a_tmp = 2.0 * q[i5 + 1];
      f_a[i4] = (d_a * b_b[i5]) + (a_tmp * q[1]);
      f_a[i4 + 1] = (d_a * b_b[i5 + 3]) + (a_tmp * q[2]);
      f_a[i4 + 2] = (d_a * b_b[i5 + 6]) + (a_tmp * d6);
      i4 += 3;
    }
    b_dv[0] = 0.0;
    b_dv[1] = e_a * (-d6);
    b_dv[2] = e_a * q[2];
    b_dv[3] = e_a * d6;
    b_dv[4] = 0.0;
    b_dv[5] = e_a * (-q[1]);
    b_dv[6] = e_a * (-q[2]);
    b_dv[7] = e_a * q[1];
    b_dv[8] = 0.0;
    for (i7 = 0; i7 < 9; i7++) {
      ST[i7] = f_a[i7] - b_dv[i7];
    }
    bias_ret->board_mag_earth[0] = 0.0;
    bias_ret->board_mag_earth[1] = 0.0;
    bias_ret->board_mag_earth[2] = 0.0;
    i9 = 0;
    for (i10 = 0; i10 < 3; i10++) {
      double d8;
      d8 = sens_filt_ret->board_mag_f[i10];
      bias_ret->board_mag_earth[0] += ST[i9] * d8;
      bias_ret->board_mag_earth[1] += ST[i9 + 1] * d8;
      bias_ret->board_mag_earth[2] += ST[i9 + 2] * d8;
      bias_ret->mti_mag_earth[i10] = 0.0;
      i9 += 3;
    }
    i11 = 0;
    d9 = bias_ret->mti_mag_earth[0];
    d10 = bias_ret->mti_mag_earth[1];
    d11 = bias_ret->mti_mag_earth[2];
    for (i12 = 0; i12 < 3; i12++) {
      double d12;
      d12 = sens_filt_ret->mti_mag_f[i12];
      d9 += ST[i11] * d12;
      d10 += ST[i11 + 1] * d12;
      d11 += ST[i11 + 2] * d12;
      i11 += 3;
    }
    double t1_pressure;
    bias_ret->mti_mag_earth[2] = d11;
    bias_ret->mti_mag_earth[1] = d10;
    bias_ret->mti_mag_earth[0] = d9;
    t1_pressure =
        airdata_atmos(SD->pd->b_param.elevation, &e_expl_temp, &t1_density,
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
    double f_a[9];
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
    double d14;
    double d15;
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
    double d28;
    double d3;
    double d30;
    double d31;
    double d33;
    double d4;
    double d59;
    double d60;
    double d61;
    double dphi;
    double dphi_tmp;
    double g_a;
    double h_a;
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
    int b_k;
    int d_k;
    int e_k;
    int g_k;
    int i13;
    int i17;
    int i19;
    int i21;
    int i26;
    int i27;
    int i31;
    int i32;
    int i34;
    int i39;
    int i40;
    int i42;
    int i44;
    int i52;
    int i56;
    int i57;
    int i59;
    int i65;
    int i68;
    int i72;
    int i73;
    int i78;
    int i80;
    int i83;
    int i87;
    int i91;
    int i92;
    int k;
    signed char b_I[16];
    signed char w_exp_tilde_tmp[9];
    d = 9.9999999999999981E+9 * ((double)sens_input->ad_gyro.status);
    C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_input->board_accel.status);
    b_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_input->mti_accel.status);
    c_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_input->ad_accel.status);
    C_total_a_tmp =
        (C_total_a_tmp_tmp + b_C_total_a_tmp_tmp) + c_C_total_a_tmp_tmp;
    C_total_a[0] = C_total_a_tmp;
    d1 = 9.9999999999999981E+9 * ((double)sens_input->board_gyro.status);
    d2 = 9.9999999999999981E+9 * ((double)sens_input->mti_gyro.status);
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
      q_mag = ((q_mag * t) * t) + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      q_mag += t * t;
    }
    absxk = fabs(x[2]);
    if (absxk > scale) {
      t = scale / absxk;
      q_mag = ((q_mag * t) * t) + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      q_mag += t * t;
    }
    absxk = fabs(x[3]);
    if (absxk > scale) {
      t = scale / absxk;
      q_mag = ((q_mag * t) * t) + 1.0;
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
      w_exp_tilde_tmp[i] = (signed char)0;
    }
    memset(&b_n_tilde[0], 0, 9U * (sizeof(double)));
    k = 0;
    b_k = 0;
    for (c_k = 0; c_k < 3; c_k++) {
      int i2;
      w_exp_tilde_tmp[k] = (signed char)1;
      i2 = 0;
      for (i3 = 0; i3 < 3; i3++) {
        double d7;
        d7 = n_tilde[i3 + b_k];
        b_n_tilde[b_k] += n_tilde[i2] * d7;
        b_n_tilde[b_k + 1] += n_tilde[i2 + 1] * d7;
        b_n_tilde[b_k + 2] += n_tilde[i2 + 2] * d7;
        i2 += 3;
      }
      k += 4;
      b_k += 3;
    }
    for (i1 = 0; i1 < 9; i1++) {
      w_exp_tilde[i1] = (((double)w_exp_tilde_tmp[i1]) - (b_a * n_tilde[i1])) +
                        ((1.0 - b_x) * b_n_tilde[i1]);
    }
    double c_a;
    c_a = b_norm(&x[7]);
    airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp,
                  &d_expl_temp);
    g_a = (0.5 * t1_density) * (c_a * c_a);
    h_a = SD->pd->c_param.c_aero * SD->pd->c_param.Cn_alpha;
    i_a = (x[0] * x[0]) - (((x[1] * x[1]) + (x[2] * x[2])) + (x[3] * x[3]));
    j_a = 2.0 * x[0];
    for (i6 = 0; i6 < 3; i6++) {
      double b_a_tmp;
      int c_a_tmp;
      int d_a_tmp;
      b_a_tmp = x[i6 + 1];
      f_a[3 * i6] = (i_a * b_b[3 * i6]) + ((2.0 * x[1]) * b_a_tmp);
      c_a_tmp = (3 * i6) + 1;
      f_a[c_a_tmp] = (i_a * b_b[c_a_tmp]) + ((2.0 * x[2]) * b_a_tmp);
      d_a_tmp = (3 * i6) + 2;
      f_a[d_a_tmp] = (i_a * b_b[d_a_tmp]) + ((2.0 * x[3]) * b_a_tmp);
    }
    b_dv[0] = 0.0;
    b_dv[3] = j_a * (-x[3]);
    b_dv[6] = j_a * x[2];
    b_dv[1] = j_a * x[3];
    b_dv[4] = 0.0;
    b_dv[7] = j_a * (-x[1]);
    b_dv[2] = j_a * (-x[2]);
    b_dv[5] = j_a * x[1];
    b_dv[8] = 0.0;
    for (i8 = 0; i8 < 9; i8++) {
      S[i8] = f_a[i8] - b_dv[i8];
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
    memset(&b_w_exp_tilde[0], 0, 9U * (sizeof(double)));
    memset(&c_w_exp_tilde[0], 0, 3U * (sizeof(double)));
    i13 = 0;
    for (i14 = 0; i14 < 3; i14++) {
      int i15;
      b_dv1[i14 + 1] = dn[i14] * b;
      i15 = 0;
      for (i16 = 0; i16 < 3; i16++) {
        double d13;
        d13 = SD->pd->c_param.J[i16 + i13];
        b_w_exp_tilde[i13] += w_exp_tilde[i15] * d13;
        b_w_exp_tilde[i13 + 1] += w_exp_tilde[i15 + 1] * d13;
        b_w_exp_tilde[i13 + 2] += w_exp_tilde[i15 + 2] * d13;
        i15 += 3;
      }
      double b_w_exp_tilde_tmp;
      b_w_exp_tilde_tmp = x[i14 + 4];
      c_w_exp_tilde[0] += b_w_exp_tilde[i13] * b_w_exp_tilde_tmp;
      c_w_exp_tilde[1] += b_w_exp_tilde[i13 + 1] * b_w_exp_tilde_tmp;
      c_w_exp_tilde[2] += b_w_exp_tilde[i13 + 2] * b_w_exp_tilde_tmp;
      i13 += 3;
    }
    dv2[0] = 0.0;
    dv2[1] = g_a * (h_a * sin(atan2(x[9], x[7])));
    dv2[2] = g_a * (h_a * (-sin(atan2(x[8], x[7]))));
    memset(&dv3[0], 0, 3U * (sizeof(double)));
    memset(&b_dt[0], 0, 3U * (sizeof(double)));
    memset(&c_dt[0], 0, 3U * (sizeof(double)));
    i17 = 0;
    d14 = dv3[0];
    d15 = dv3[1];
    d16 = dv3[2];
    d17 = b_dt[0];
    d18 = b_dt[1];
    d19 = b_dt[2];
    d20 = c_dt[0];
    d21 = c_dt[1];
    d22 = c_dt[2];
    d23 = x[7];
    d24 = x[8];
    d25 = x[9];
    for (i18 = 0; i18 < 3; i18++) {
      double d26;
      double d27;
      double d29;
      double d32;
      double d34;
      double d36;
      d26 = SD->pd->c_param.Jinv[i17];
      d27 = c_w_exp_tilde[i18];
      d14 += d26 * d27;
      d29 = dv2[i18];
      d17 += (dt * d26) * d29;
      d32 = S[i17];
      d20 += (dt * d32) * SD->pd->c_param.g[i18];
      d34 = d32 * d23;
      d26 = SD->pd->c_param.Jinv[i17 + 1];
      d15 += d26 * d27;
      d18 += (dt * d26) * d29;
      d32 = S[i17 + 1];
      d21 += (dt * d32) * SD->pd->c_param.g[i18];
      d34 += d32 * d24;
      d26 = SD->pd->c_param.Jinv[i17 + 2];
      d16 += d26 * d27;
      d19 += (dt * d26) * d29;
      d32 = S[i17 + 2];
      d22 += (dt * d32) * SD->pd->c_param.g[i18];
      d34 += d32 * d25;
      d36 = C_total_a[i18];
      c_w_exp_tilde[i18] =
          (((w_exp_tilde[i18] * d23) + (w_exp_tilde[i18 + 3] * d24)) +
           (w_exp_tilde[i18 + 6] * d25)) +
          (dt *
           ((((C_total_a_tmp_tmp / d36) * sens_input->board_accel.meas[i18]) +
             ((b_C_total_a_tmp_tmp / d36) * sens_input->mti_accel.meas[i18])) +
            ((c_C_total_a_tmp_tmp / d36) * sens_input->ad_accel.meas[i18])));
      b_S[i18] = d34;
      i17 += 3;
    }
    memset(&c_q[0], 0, (sizeof(double)) << 2);
    i19 = 0;
    d28 = c_q[0];
    d30 = c_q[1];
    d31 = c_q[2];
    d33 = c_q[3];
    for (i20 = 0; i20 < 4; i20++) {
      double d35;
      d35 = b_dv1[i20];
      d28 += b_q[i19] * d35;
      d30 += b_q[i19 + 1] * d35;
      d31 += b_q[i19 + 2] * d35;
      d33 += b_q[i19 + 3] * d35;
      i19 += 4;
    }
    double W_dt_tmp;
    double b_W_dt_tmp;
    double c_W_dt_tmp;
    double d_W_dt_tmp;
    double e_W_dt_tmp;
    double f_W_dt_tmp;
    x_pred[0] = d28;
    x_pred[1] = d30;
    x_pred[2] = d31;
    x_pred[3] = d33;
    x_pred[4] = d14 + d17;
    x_pred[7] = c_w_exp_tilde[0] + d20;
    x_pred[5] = d15 + d18;
    x_pred[8] = c_w_exp_tilde[1] + d21;
    x_pred[6] = d16 + d19;
    x_pred[9] = c_w_exp_tilde[2] + d22;
    x_pred[10] = x[10] + (dt * b_S[0]);
    memset(&F[0], 0, 121U * (sizeof(double)));
    k_a = 0.5 * dt;
    W_dt[0] = 0.0;
    W_dt_tmp = k_a * (-x[4]);
    W_dt[4] = W_dt_tmp;
    b_W_dt_tmp = k_a * (-x[5]);
    W_dt[8] = b_W_dt_tmp;
    c_W_dt_tmp = k_a * (-x[6]);
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
    memset(&c_b[0], 0, (sizeof(double)) << 4);
    i21 = 0;
    for (i22 = 0; i22 < 4; i22++) {
      int i23;
      i23 = 0;
      for (i25 = 0; i25 < 4; i25++) {
        double d37;
        d37 = W_dt[i25 + i21];
        c_b[i21] += W_dt[i23] * d37;
        c_b[i21 + 1] += W_dt[i23 + 1] * d37;
        c_b[i21 + 2] += W_dt[i23 + 2] * d37;
        c_b[i21 + 3] += W_dt[i23 + 3] * d37;
        i23 += 4;
      }
      i21 += 4;
    }
    for (i24 = 0; i24 < 16; i24++) {
      b_I[i24] = (signed char)0;
    }
    memset(&b_W_dt[0], 0, (sizeof(double)) << 4);
    memset(&d_b[0], 0, (sizeof(double)) << 4);
    d_k = 0;
    e_k = 0;
    for (f_k = 0; f_k < 4; f_k++) {
      int i28;
      b_I[d_k] = (signed char)1;
      i28 = 0;
      for (i30 = 0; i30 < 4; i30++) {
        double d38;
        d38 = c_b[i30 + e_k];
        b_W_dt[e_k] += W_dt[i28] * d38;
        d_b[e_k] += c_b[i28] * d38;
        b_W_dt[e_k + 1] += W_dt[i28 + 1] * d38;
        d_b[e_k + 1] += c_b[i28 + 1] * d38;
        b_W_dt[e_k + 2] += W_dt[i28 + 2] * d38;
        d_b[e_k + 2] += c_b[i28 + 2] * d38;
        b_W_dt[e_k + 3] += W_dt[i28 + 3] * d38;
        d_b[e_k + 3] += c_b[i28 + 3] * d38;
        i28 += 4;
      }
      d_k += 5;
      e_k += 4;
    }
    i26 = 0;
    i27 = 0;
    for (i29 = 0; i29 < 4; i29++) {
      F[i26] = (((((double)b_I[i27]) + W_dt[i27]) + (0.5 * c_b[i27])) +
                (0.16666666666666666 * b_W_dt[i27])) +
               (0.041666666666666664 * d_b[i27]);
      F[i26 + 1] =
          (((((double)b_I[i27 + 1]) + W_dt[i27 + 1]) + (0.5 * c_b[i27 + 1])) +
           (0.16666666666666666 * b_W_dt[i27 + 1])) +
          (0.041666666666666664 * d_b[i27 + 1]);
      F[i26 + 2] =
          (((((double)b_I[i27 + 2]) + W_dt[i27 + 2]) + (0.5 * c_b[i27 + 2])) +
           (0.16666666666666666 * b_W_dt[i27 + 2])) +
          (0.041666666666666664 * d_b[i27 + 2]);
      F[i26 + 3] =
          (((((double)b_I[i27 + 3]) + W_dt[i27 + 3]) + (0.5 * c_b[i27 + 3])) +
           (0.16666666666666666 * b_W_dt[i27 + 3])) +
          (0.041666666666666664 * d_b[i27 + 3]);
      i26 += 11;
      i27 += 4;
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
    f_a_tmp = k_a * (-q[1]);
    l_a[4] = f_a_tmp;
    g_a_tmp = k_a * (-q[2]);
    l_a[8] = g_a_tmp;
    h_a_tmp = k_a * (-q[3]);
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
    i31 = 0;
    i32 = 0;
    for (i33 = 0; i33 < 3; i33++) {
      F[i31 + 44] = l_a[i32 + 4];
      F[i31 + 45] = l_a[i32 + 5];
      F[i31 + 46] = l_a[i32 + 6];
      F[i31 + 47] = l_a[i32 + 7];
      i31 += 11;
      i32 += 4;
    }
    m_a = (0.5 * SD->pd->d_param.c_aero) * SD->pd->d_param.Cn_alpha;
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
    memset(&b_n_tilde[0], 0, 9U * (sizeof(double)));
    i34 = 0;
    for (i35 = 0; i35 < 3; i35++) {
      int i36;
      i36 = 0;
      for (i38 = 0; i38 < 3; i38++) {
        double d39;
        d39 = n_tilde[i38 + i34];
        b_n_tilde[i34] += n_tilde[i36] * d39;
        b_n_tilde[i34 + 1] += n_tilde[i36 + 1] * d39;
        b_n_tilde[i34 + 2] += n_tilde[i36 + 2] * d39;
        i36 += 3;
      }
      i34 += 3;
    }
    for (i37 = 0; i37 < 9; i37++) {
      w_exp_tilde[i37] =
          (((double)w_exp_tilde_tmp[i37]) - (b_a * n_tilde[i37])) +
          ((1.0 - b_x) * b_n_tilde[i37]);
    }
    memset(&b_dv[0], 0, 9U * (sizeof(double)));
    i39 = 0;
    i40 = 0;
    for (i41 = 0; i41 < 3; i41++) {
      int i43;
      i43 = 0;
      for (i45 = 0; i45 < 3; i45++) {
        double d40;
        d40 = w_exp_tilde[i45 + i39];
        b_dv[i39] += SD->pd->d_param.Jinv[i43] * d40;
        b_dv[i39 + 1] += SD->pd->d_param.Jinv[i43 + 1] * d40;
        b_dv[i39 + 2] += SD->pd->d_param.Jinv[i43 + 2] * d40;
        F[(i45 + i40) + 48] = 0.0;
        i43 += 3;
      }
      i39 += 3;
      i40 += 11;
    }
    i42 = 0;
    i44 = 0;
    for (i46 = 0; i46 < 3; i46++) {
      int i47;
      i47 = 0;
      for (i48 = 0; i48 < 3; i48++) {
        double d41;
        d41 = SD->pd->d_param.J[i48 + i42];
        F[i44 + 48] += b_dv[i47] * d41;
        F[i44 + 49] += b_dv[i47 + 1] * d41;
        F[i44 + 50] += b_dv[i47 + 2] * d41;
        i47 += 3;
      }
      i42 += 3;
      i44 += 11;
    }
    b_dv[1] = t1_density * (m_a * x[9]);
    b_dv[4] = 0.0;
    b_dv[7] = t1_density * (m_a * x[7]);
    b_dv[2] = t1_density * (m_a * (-x[8]));
    b_dv[5] = t1_density * (m_a * (-x[7]));
    b_dv[8] = 0.0;
    c_x = 0.0;
    for (i49 = 0; i49 < 3; i49++) {
      double d42;
      double d43;
      double d44;
      int F_tmp;
      b_dv[3 * i49] = 0.0;
      F_tmp = 11 * (i49 + 7);
      d42 = 0.0;
      d43 = 0.0;
      d44 = 0.0;
      for (i50 = 0; i50 < 3; i50++) {
        double d45;
        d45 = b_dv[i50 + (3 * i49)];
        d42 += (dt * SD->pd->d_param.Jinv[3 * i50]) * d45;
        d43 += (dt * SD->pd->d_param.Jinv[(3 * i50) + 1]) * d45;
        d44 += (dt * SD->pd->d_param.Jinv[(3 * i50) + 2]) * d45;
      }
      F[F_tmp + 6] = d44;
      F[F_tmp + 5] = d43;
      F[F_tmp + 4] = d42;
      c_x += x[i49 + 1] * SD->pd->d_param.g[i49];
    }
    d_x[0] = (x[2] * SD->pd->d_param.g[2]) - (SD->pd->d_param.g[1] * x[3]);
    d_x[1] = (SD->pd->d_param.g[0] * x[3]) - (x[1] * SD->pd->d_param.g[2]);
    d_x[2] = (x[1] * SD->pd->d_param.g[1]) - (SD->pd->d_param.g[0] * x[2]);
    dv4[0] = 0.0;
    dv4[3] = x[0] * (-SD->pd->d_param.g[2]);
    dv4[6] = x[0] * SD->pd->d_param.g[1];
    dv4[1] = x[0] * SD->pd->d_param.g[2];
    dv4[4] = 0.0;
    dv4[7] = x[0] * (-SD->pd->d_param.g[0]);
    dv4[2] = x[0] * (-SD->pd->d_param.g[1]);
    dv4[5] = x[0] * SD->pd->d_param.g[0];
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
    memset(&c_skewed_exp_w_tmp[0], 0, 9U * (sizeof(double)));
    memset(&b_dv[0], 0, 9U * (sizeof(double)));
    r_q_tmp[0] = x[0];
    b_r_q_tmp = 0.0;
    for (i51 = 0; i51 < 3; i51++) {
      double b_F_tmp;
      double d46;
      double d47;
      int c_F_tmp;
      int d_F_tmp;
      int e_F_tmp;
      F[i51 + 7] = dt * (2.0 * ((x[0] * SD->pd->d_param.g[i51]) - d_x[i51]));
      b_F_tmp = x[i51 + 1];
      c_F_tmp = 11 * (i51 + 1);
      F[c_F_tmp + 7] =
          dt *
          (2.0 * ((((c_x * b_b[3 * i51]) + (x[1] * SD->pd->d_param.g[i51])) -
                   (SD->pd->d_param.g[0] * b_F_tmp)) +
                  dv4[3 * i51]));
      d_F_tmp = (3 * i51) + 1;
      F[c_F_tmp + 8] =
          dt *
          (2.0 * ((((c_x * b_b[d_F_tmp]) + (x[2] * SD->pd->d_param.g[i51])) -
                   (SD->pd->d_param.g[1] * b_F_tmp)) +
                  dv4[d_F_tmp]));
      e_F_tmp = (3 * i51) + 2;
      F[c_F_tmp + 9] =
          dt *
          (2.0 * ((((c_x * b_b[e_F_tmp]) + (x[3] * SD->pd->d_param.g[i51])) -
                   (SD->pd->d_param.g[2] * b_F_tmp)) +
                  dv4[e_F_tmp]));
      d46 = c_skewed_exp_w_tmp[3 * i51];
      d47 = b_dv[3 * i51];
      for (i54 = 0; i54 < 3; i54++) {
        double d48;
        double d49;
        int b_skewed_exp_w_tmp_tmp;
        int i55;
        int skewed_exp_w_tmp_tmp;
        i55 = i54 + (3 * i51);
        d48 = b_skewed_exp_w_tmp[i55];
        d49 = skewed_exp_w_tmp[i55];
        d46 += skewed_exp_w_tmp[3 * i54] * d48;
        d47 += (2.0 * b_skewed_exp_w_tmp[3 * i54]) * d49;
        skewed_exp_w_tmp_tmp = (3 * i54) + 1;
        c_skewed_exp_w_tmp[d_F_tmp] +=
            skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d48;
        b_dv[d_F_tmp] += (2.0 * b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp]) * d49;
        b_skewed_exp_w_tmp_tmp = (3 * i54) + 2;
        c_skewed_exp_w_tmp[e_F_tmp] +=
            skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d48;
        b_dv[e_F_tmp] +=
            (2.0 * b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp]) * d49;
      }
      int f_F_tmp;
      int g_F_tmp;
      b_dv[3 * i51] = d47;
      c_skewed_exp_w_tmp[3 * i51] = d46;
      f_F_tmp = 11 * (i51 + 4);
      F[f_F_tmp + 7] = (dt * skewed_exp_w_tmp[3 * i51]) + (n_a * (d46 - d47));
      g_F_tmp = 11 * (i51 + 7);
      F[g_F_tmp + 7] = w_exp_tilde[3 * i51];
      F[f_F_tmp + 8] = (dt * skewed_exp_w_tmp[d_F_tmp]) +
                       (n_a * (c_skewed_exp_w_tmp[d_F_tmp] - b_dv[d_F_tmp]));
      F[g_F_tmp + 8] = w_exp_tilde[d_F_tmp];
      F[f_F_tmp + 9] = (dt * skewed_exp_w_tmp[e_F_tmp]) +
                       (n_a * (c_skewed_exp_w_tmp[e_F_tmp] - b_dv[e_F_tmp]));
      F[g_F_tmp + 9] = w_exp_tilde[e_F_tmp];
      r_q_tmp[i51 + 1] = -b_F_tmp;
      b_r_q_tmp += (-b_F_tmp) * x[i51 + 7];
    }
    c_r_q_tmp[0] = (r_q_tmp[2] * x[9]) - (r_q_tmp[3] * x[8]);
    c_r_q_tmp[1] = (r_q_tmp[3] * x[7]) - (r_q_tmp[1] * x[9]);
    c_r_q_tmp[2] = (r_q_tmp[1] * x[8]) - (r_q_tmp[2] * x[7]);
    i52 = 0;
    for (i53 = 0; i53 < 3; i53++) {
      double b_dt_tmp;
      double dt_tmp;
      dt_tmp = x[i53 + 7];
      d_dt[i53] = dt * (2.0 * ((r_q_tmp[0] * dt_tmp) - c_r_q_tmp[i53]));
      b_dt_tmp = r_q_tmp[i53 + 1];
      d_dt[i52 + 3] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[i52]) + (r_q_tmp[1] * dt_tmp)) -
                        (x[7] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[i52])));
      d_dt[i52 + 4] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[i52 + 1]) + (r_q_tmp[2] * dt_tmp)) -
                        (x[8] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[i52 + 1])));
      d_dt[i52 + 5] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[i52 + 2]) + (r_q_tmp[3] * dt_tmp)) -
                        (x[9] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[i52 + 2])));
      i52 += 3;
    }
    double p_a;
    F[10] = d_dt[0];
    F[21] = d_dt[3];
    F[32] = d_dt[6];
    F[43] = d_dt[9];
    o_a = (r_q_tmp[0] * r_q_tmp[0]) -
          (((r_q_tmp[1] * r_q_tmp[1]) + (r_q_tmp[2] * r_q_tmp[2])) +
           (r_q_tmp[3] * r_q_tmp[3]));
    p_a = 2.0 * r_q_tmp[0];
    b_dv[0] = 0.0;
    b_dv[3] = p_a * (-r_q_tmp[3]);
    b_dv[6] = p_a * r_q_tmp[2];
    b_dv[1] = p_a * r_q_tmp[3];
    b_dv[4] = 0.0;
    b_dv[7] = p_a * (-r_q_tmp[1]);
    b_dv[2] = p_a * (-r_q_tmp[2]);
    b_dv[5] = p_a * r_q_tmp[1];
    b_dv[8] = 0.0;
    i56 = 0;
    i57 = 0;
    for (i58 = 0; i58 < 3; i58++) {
      double l_a_tmp;
      l_a_tmp = r_q_tmp[i58 + 1];
      f_a[i56] = (o_a * b_b[i56]) + ((2.0 * r_q_tmp[1]) * l_a_tmp);
      f_a[i56 + 1] = (o_a * b_b[i56 + 1]) + ((2.0 * r_q_tmp[2]) * l_a_tmp);
      f_a[i56 + 2] = (o_a * b_b[i56 + 2]) + ((2.0 * r_q_tmp[3]) * l_a_tmp);
      F[i57 + 87] = dt * (f_a[i56] - b_dv[i56]);
      i56 += 3;
      i57 += 11;
    }
    F[120] = 1.0;
    memset(&b_F[0], 0, 121U * (sizeof(double)));
    i59 = 0;
    for (i60 = 0; i60 < 11; i60++) {
      int i61;
      i61 = 0;
      for (i63 = 0; i63 < 11; i63++) {
        double d50;
        d50 = P[i63 + i59];
        for (i67 = 0; i67 < 11; i67++) {
          int h_F_tmp;
          h_F_tmp = i67 + i59;
          b_F[h_F_tmp] += F[i67 + i61] * d50;
        }
        i61 += 11;
      }
      i59 += 11;
    }
    for (i62 = 0; i62 < 11; i62++) {
      int i64;
      i64 = 0;
      for (i66 = 0; i66 < 11; i66++) {
        double d51;
        int i70;
        d51 = 0.0;
        i70 = 0;
        for (i71 = 0; i71 < 11; i71++) {
          d51 += b_F[i70 + i62] * F[i70 + i66];
          i70 += 11;
        }
        int P_pred_tmp;
        P_pred_tmp = i64 + i62;
        P_pred[P_pred_tmp] = d51 + Q[P_pred_tmp];
        i64 += 11;
      }
    }
    i65 = 0;
    i68 = 0;
    for (i69 = 0; i69 < 3; i69++) {
      b_P_pred[i65] = P_pred[i68 + 48] + R[i65];
      b_P_pred[i65 + 1] = P_pred[i68 + 49] + R[i65 + 1];
      b_P_pred[i65 + 2] = P_pred[i68 + 50] + R[i65 + 2];
      i65 += 3;
      i68 += 11;
    }
    inv(b_P_pred, b_dv);
    memset(&K[0], 0, 33U * (sizeof(double)));
    i72 = 0;
    i73 = 0;
    for (i74 = 0; i74 < 3; i74++) {
      int i75;
      i75 = 0;
      for (i76 = 0; i76 < 3; i76++) {
        double d52;
        d52 = b_dv[i76 + i72];
        for (i77 = 0; i77 < 11; i77++) {
          int K_tmp;
          K_tmp = i77 + i73;
          K[K_tmp] += P_pred[(i77 + i75) + 44] * d52;
        }
        i75 += 11;
      }
      i72 += 3;
      i73 += 11;
    }
    memset(&c_I[0], 0, 121U * (sizeof(signed char)));
    g_k = 0;
    for (h_k = 0; h_k < 11; h_k++) {
      c_I[g_k] = (signed char)1;
      g_k += 12;
    }
    i78 = 0;
    for (i79 = 0; i79 < 4; i79++) {
      for (i81 = 0; i81 < 11; i81++) {
        int E_tmp;
        E_tmp = i81 + i78;
        E[E_tmp] = (double)c_I[E_tmp];
      }
      i78 += 11;
    }
    i80 = 0;
    for (i82 = 0; i82 < 3; i82++) {
      for (i84 = 0; i84 < 11; i84++) {
        int b_E_tmp;
        b_E_tmp = i84 + i80;
        E[b_E_tmp + 44] = ((double)c_I[b_E_tmp + 44]) - K[b_E_tmp];
      }
      i80 += 11;
    }
    i83 = 0;
    for (i85 = 0; i85 < 4; i85++) {
      for (i86 = 0; i86 < 11; i86++) {
        int c_E_tmp;
        c_E_tmp = (i86 + i83) + 77;
        E[c_E_tmp] = (double)c_I[c_E_tmp];
      }
      i83 += 11;
    }
    memset(&b_E[0], 0, 121U * (sizeof(double)));
    i87 = 0;
    for (i88 = 0; i88 < 11; i88++) {
      int i89;
      i89 = 0;
      for (i90 = 0; i90 < 11; i90++) {
        double d53;
        d53 = P_pred[i90 + i87];
        for (i94 = 0; i94 < 11; i94++) {
          int d_E_tmp;
          d_E_tmp = i94 + i87;
          b_E[d_E_tmp] += E[i94 + i89] * d53;
        }
        i89 += 11;
      }
      i87 += 11;
    }
    memset(&b_K[0], 0, 33U * (sizeof(double)));
    i91 = 0;
    i92 = 0;
    for (i93 = 0; i93 < 3; i93++) {
      int i95;
      i95 = 0;
      for (i96 = 0; i96 < 3; i96++) {
        double d54;
        d54 = R[i96 + i91];
        for (i98 = 0; i98 < 11; i98++) {
          int b_K_tmp;
          b_K_tmp = i98 + i92;
          b_K[b_K_tmp] += K[i98 + i95] * d54;
        }
        i95 += 11;
      }
      i91 += 3;
      i92 += 11;
    }
    memset(&c_E[0], 0, 121U * (sizeof(double)));
    memset(&c_K[0], 0, 121U * (sizeof(double)));
    for (i97 = 0; i97 < 11; i97++) {
      for (i99 = 0; i99 < 11; i99++) {
        double d55;
        d55 = E[i97 + (11 * i99)];
        for (i102 = 0; i102 < 11; i102++) {
          int e_E_tmp;
          e_E_tmp = i102 + (11 * i97);
          c_E[e_E_tmp] += b_E[i102 + (11 * i99)] * d55;
        }
      }
      for (i101 = 0; i101 < 3; i101++) {
        double d57;
        d57 = K[i97 + (11 * i101)];
        for (i103 = 0; i103 < 11; i103++) {
          int c_K_tmp;
          c_K_tmp = i103 + (11 * i97);
          c_K[c_K_tmp] += b_K[i103 + (11 * i101)] * d57;
        }
      }
    }
    for (i100 = 0; i100 < 121; i100++) {
      P_ret[i100] = c_E[i100] + c_K[i100];
    }
    double d56;
    double d58;
    d56 = d1 / d3;
    d58 = d2 / d3;
    d59 =
        ((((d1 / d4) * (sens_input->board_gyro.meas[0] - bias->board_gyro[0])) +
          ((d2 / d4) * (sens_input->mti_gyro.meas[0] - bias->mti_gyro[0]))) +
         (C_ad_w_idx_0 * (sens_input->ad_gyro.meas[0] - bias->ad_gyro[0]))) -
        x_pred[4];
    d60 = (((d56 * (sens_input->board_gyro.meas[1] - bias->board_gyro[1])) +
            (d58 * (sens_input->mti_gyro.meas[1] - bias->mti_gyro[1]))) +
           (d * (sens_input->ad_gyro.meas[1] - bias->ad_gyro[1]))) -
          x_pred[5];
    d61 = (((d56 * (sens_input->board_gyro.meas[2] - bias->board_gyro[2])) +
            (d58 * (sens_input->mti_gyro.meas[2] - bias->mti_gyro[2]))) +
           (d * (sens_input->ad_gyro.meas[2] - bias->ad_gyro[2]))) -
          x_pred[6];
    for (i104 = 0; i104 < 11; i104++) {
      x_ret[i104] = x_pred[i104] + (((K[i104] * d59) + (K[i104 + 11] * d60)) +
                                    (K[i104 + 22] * d61));
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
      b_q_mag = ((b_q_mag * c_t) * c_t) + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    c_absxk = fabs(x_ret[2]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = ((b_q_mag * c_t) * c_t) + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    c_absxk = fabs(x_ret[3]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = ((b_q_mag * c_t) * c_t) + 1.0;
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
      memcpy(&b_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&b_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      memcpy(&e_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&d_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      b_ekf_correct(e_x_ret, d_P_ret, sens_input->board_baro.meas,
                    bias->board_baro, x_ret, P_ret);
    }
    if (sens_input->mti_baro.status) {
      memcpy(&c_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&c_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      memcpy(&g_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&f_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      b_ekf_correct(g_x_ret, f_P_ret, sens_input->mti_baro.meas, bias->mti_baro,
                    x_ret, P_ret);
    }
    if (sens_input->board_mag.status) {
      memcpy(&d_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&e_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      memcpy(&h_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&h_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      ekf_correct(h_x_ret, h_P_ret, sens_input->board_mag.meas,
                  bias->board_mag_earth, b_b, x_ret, P_ret);
    }
    if (sens_input->mti_mag.status) {
      memcpy(&f_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&g_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      memcpy(&i_x_ret[0], &x_ret[0], 11U * (sizeof(double)));
      memcpy(&i_P_ret[0], &P_ret[0], 121U * (sizeof(double)));
      ekf_correct(i_x_ret, i_P_ret, sens_input->mti_mag.meas,
                  bias->mti_mag_earth, b_b, x_ret, P_ret);
    }
  }
  airdata->pressure = airdata_atmos(x_ret[10], &airdata->temperature,
                                    &airdata->density, &airdata->sonic_speed,
                                    &airdata->mach, &airdata->dynamic_pressure);
  airspeed = b_norm(&x_ret[7]);
  airdata->mach = airspeed / airdata->sonic_speed;
  airdata->dynamic_pressure = (0.5 * airdata->density) * (airspeed * airspeed);
  roll_state[0] = atan2(2.0 * ((x_ret[2] * x_ret[3]) + (x_ret[0] * x_ret[1])),
                        (((x_ret[0] * x_ret[0]) - (x_ret[1] * x_ret[1])) -
                         (x_ret[2] * x_ret[2])) +
                            (x_ret[3] * x_ret[3]));
  roll_state[1] = x_ret[4];
}
