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

static void b_ekf_correct(const double x[11], const double P[121], double y,
                          double b, double x_new[11], double P_new[121]);

static double b_norm(const double x[3]);

static void controller_codegen_entry_init(GNC_codegenStackData *SD);

static void dynamics_init(GNC_codegenStackData *SD);

static void dynamics_jacobian_init(GNC_codegenStackData *SD);

static void ekf_correct(const double x[11], const double P[121],
                        const double y[3], const double b[3], const double R[9],
                        double x_new[11], double P_new[121]);

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
  SD->pd->param.altitude_initial = 420.0;
  SD->pd->param.c_aero = -0.035602020206989993;
  SD->pd->param.c_canard = 0.000649003220184504;
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
  SD->pd->c_param.altitude_initial = 420.0;
  SD->pd->c_param.c_aero = -0.035602020206989993;
  SD->pd->c_param.c_canard = 0.000649003220184504;
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
  SD->pd->d_param.altitude_initial = 420.0;
  SD->pd->d_param.c_aero = -0.035602020206989993;
  SD->pd->d_param.c_canard = 0.000649003220184504;
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
  double d10;
  double d11;
  double d12;
  double d13;
  double d14;
  double d15;
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
  int i18;
  int i19;
  int i2;
  int i21;
  int i22;
  int i23;
  int i24;
  int i25;
  int i26;
  int i28;
  int i29;
  int i30;
  int i31;
  int i32;
  int i33;
  int i34;
  int i35;
  int i36;
  int i37;
  int i38;
  int i39;
  int i4;
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
    for (i11 = 0; i11 < 11; i11++) {
      double d1;
      d1 = y_tmp[i11 + (11 * i10)];
      for (i13 = 0; i13 < 11; i13++) {
        int P_tmp;
        P_tmp = i13 + (11 * i10);
        b_P[P_tmp] += P[i13 + (11 * i11)] * d1;
      }
    }
    for (i12 = 0; i12 < 3; i12++) {
      double d2;
      d2 = 0.0;
      for (i14 = 0; i14 < 11; i14++) {
        d2 += b_H[i10 + (3 * i14)] * y_tmp[i14 + (11 * i12)];
      }
      int e_H_tmp;
      e_H_tmp = i10 + (3 * i12);
      c_H[e_H_tmp] = d2 + R[e_H_tmp];
    }
  }
  mrdiv(b_P, c_H, K);
  memset(&b_I[0], 0, 121U * (sizeof(signed char)));
  k = 0;
  for (b_k = 0; b_k < 11; b_k++) {
    b_I[k] = (signed char)1;
    k += 12;
  }
  for (i15 = 0; i15 < 11; i15++) {
    double d3;
    double d4;
    double d5;
    int i16;
    int i17;
    i16 = 0;
    i17 = 0;
    d3 = K[i15];
    d4 = K[i15 + 11];
    d5 = K[i15 + 22];
    for (i21 = 0; i21 < 11; i21++) {
      int E_tmp;
      E_tmp = i17 + i15;
      E[E_tmp] = ((double)b_I[E_tmp]) -
                 (((d3 * H[i16]) + (d4 * H[i16 + 1])) + (d5 * H[i16 + 2]));
      i16 += 3;
      i17 += 11;
    }
  }
  memset(&b_E[0], 0, 121U * (sizeof(double)));
  i18 = 0;
  for (i19 = 0; i19 < 11; i19++) {
    int i20;
    i20 = 0;
    for (i22 = 0; i22 < 11; i22++) {
      double d6;
      d6 = P[i22 + i18];
      for (i26 = 0; i26 < 11; i26++) {
        int b_E_tmp;
        b_E_tmp = i26 + i18;
        b_E[b_E_tmp] += E[i26 + i20] * d6;
      }
      i20 += 11;
    }
    i18 += 11;
  }
  memset(&b_K[0], 0, 33U * (sizeof(double)));
  i23 = 0;
  i24 = 0;
  for (i25 = 0; i25 < 3; i25++) {
    int i27;
    i27 = 0;
    for (i28 = 0; i28 < 3; i28++) {
      double d7;
      d7 = R[i28 + i23];
      for (i30 = 0; i30 < 11; i30++) {
        int K_tmp;
        K_tmp = i30 + i24;
        b_K[K_tmp] += K[i30 + i27] * d7;
      }
      i27 += 11;
    }
    i23 += 3;
    i24 += 11;
  }
  memset(&c_E[0], 0, 121U * (sizeof(double)));
  memset(&c_K[0], 0, 121U * (sizeof(double)));
  for (i29 = 0; i29 < 11; i29++) {
    for (i31 = 0; i31 < 11; i31++) {
      double d8;
      d8 = E[i29 + (11 * i31)];
      for (i34 = 0; i34 < 11; i34++) {
        int c_E_tmp;
        c_E_tmp = i34 + (11 * i29);
        c_E[c_E_tmp] += b_E[i34 + (11 * i31)] * d8;
      }
    }
    for (i33 = 0; i33 < 3; i33++) {
      double d9;
      d9 = K[i29 + (11 * i33)];
      for (i36 = 0; i36 < 11; i36++) {
        int b_K_tmp;
        b_K_tmp = i36 + (11 * i29);
        c_K[b_K_tmp] += b_K[i36 + (11 * i33)] * d9;
      }
    }
  }
  for (i32 = 0; i32 < 121; i32++) {
    P_new[i32] = c_E[i32] + c_K[i32];
  }
  for (i35 = 0; i35 < 3; i35++) {
    double a_tmp;
    int b_a_tmp;
    int c_a_tmp;
    a_tmp = x[i35 + 1];
    c_a[3 * i35] = (a * ((double)iv[3 * i35])) + ((2.0 * x[1]) * a_tmp);
    b_a_tmp = (3 * i35) + 1;
    c_a[b_a_tmp] = (a * ((double)iv[b_a_tmp])) + ((2.0 * x[2]) * a_tmp);
    c_a_tmp = (3 * i35) + 2;
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
  for (i37 = 0; i37 < 9; i37++) {
    c_a[i37] -= b_dv[i37];
  }
  d10 = b[0];
  d11 = b[1];
  d12 = b[2];
  for (i38 = 0; i38 < 3; i38++) {
    b_y[i38] = y[i38] - (((c_a[i38] * d10) + (c_a[i38 + 3] * d11)) +
                         (c_a[i38 + 6] * d12));
  }
  d13 = b_y[0];
  d14 = b_y[1];
  d15 = b_y[2];
  for (i39 = 0; i39 < 11; i39++) {
    x_new[i39] =
        x[i39] + (((K[i39] * d13) + (K[i39 + 11] * d14)) + (K[i39 + 22] * d15));
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
  for (i = 0; i < 9; i++) {
    SD->pd->b_param.J[i] = dv[i];
    SD->pd->b_param.Jinv[i] = dv1[i];
  }
  SD->pd->b_param.altitude_initial = 420.0;
  SD->pd->b_param.c_aero = -0.035602020206989993;
  SD->pd->b_param.c_canard = 0.000649003220184504;
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
                              double dt_ctrl, const double where_it_is[2],
                              double pdyn, double delta_encoder,
                              struct0_T *ctrl_mem, double *u_motor,
                              double where_it_isnt[2], bool *w_status_ctrl) {
  static const double b_dv[5] = {-1.5707963267948966, 0.0, 0.0, 0.0, 0.0};
  static const signed char iv[5] = {7, 15, 25, 35, 45};
  static const signed char iv1[5] = {0, 0, 1, -1, 0};
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
    i = (int)iv[step_idx];
    if (b_time >= ((double)i)) {
      double d;
      double q;
      double r;
      double x;
      d = b_dv[step_idx];
      x = (((double)iv1[step_idx]) + ((b_time - ((double)i)) * d)) +
          3.1415926535897931;
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
  delta = delta_encoder / 2.0;
  pdyn_params = pdyn * SD->pd->param.c_canard;
  if (fabs(delta) < 0.005) {
    delta = 0.0;
  }
  delta_lp = (0.75 * ctrl_mem->delta_lp) + (0.25 * delta);
  w_dot_lp = (0.75 * ctrl_mem->w_dot_lp) +
             ((0.25 * (where_it_is[1] - ctrl_mem->w)) / dt_ctrl);
  r_idx_0 = pdyn_params * delta_lp;
  P[0] = ctrl_mem->P[0] + 1.0E-5;
  P[1] = ctrl_mem->P[1];
  P[2] = ctrl_mem->P[2];
  P[3] = ctrl_mem->P[3] + 1.0E-9;
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
  b_dv1[0] = 1.0 - (K[0] * r_idx_0);
  b_dv1[1] = 0.0 - (K[1] * r_idx_0);
  b_dv1[2] = 0.0 - (K[0] * pdyn_params);
  b_dv1[3] = 1.0 - (K[1] * pdyn_params);
  memset(&dv2[0], 0, (sizeof(double)) << 2);
  i1 = 0;
  d3 = b_dv1[0];
  d4 = b_dv1[1];
  d5 = b_dv1[2];
  d6 = b_dv1[3];
  for (i2 = 0; i2 < 2; i2++) {
    double d7;
    d7 = P[i1];
    dv2[i1] += d3 * d7;
    dv2[i1 + 1] += d4 * d7;
    d7 = P[i1 + 1];
    dv2[i1] += d5 * d7;
    dv2[i1 + 1] += d6 * d7;
    i1 += 2;
  }
  memset(&dv3[0], 0, (sizeof(double)) << 2);
  i3 = 0;
  d8 = dv2[0];
  d9 = dv2[1];
  d10 = dv2[2];
  d11 = dv2[3];
  for (i4 = 0; i4 < 2; i4++) {
    double d12;
    d12 = b_dv1[i4];
    dv3[i3] += d8 * d12;
    dv3[i3 + 1] += d9 * d12;
    b_K[i3] = K[0] * K[i4];
    d12 = b_dv1[i4 + 2];
    dv3[i3] += d10 * d12;
    dv3[i3 + 1] += d11 * d12;
    b_K[i3 + 1] = K[1] * K[i4];
    i3 += 2;
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
  *u_motor = fmin(fmax(((a * b_x) * deviation_idx_0) +
                           ((a * sqrt((2.0 * b_x) + (x_tmp + (blend * 20.0)))) *
                            deviation_idx_1),
                       -0.17453292519943295),
                  0.17453292519943295) *
             2.0;
  if (pdyn < 500.0) {
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
  double f_x[11];
  double g_x[11];
  double h_x[11];
  double i_x[11];
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
  double k_a;
  double k_expl_temp;
  double l_expl_temp;
  double m_expl_temp;
  double n_expl_temp;
  double o_expl_temp;
  double p_expl_temp;
  double t1_density;
  int b_i;
  int c_k;
  int f_k;
  int h_k;
  int i1;
  int i10;
  int i11;
  int i12;
  int i14;
  int i16;
  int i18;
  int i2;
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
  int i4;
  int i41;
  int i45;
  int i46;
  int i48;
  int i49;
  int i50;
  int i51;
  int i53;
  int i54;
  int i58;
  int i60;
  int i62;
  int i63;
  int i66;
  int i67;
  int i69;
  int i7;
  int i71;
  int i73;
  int i75;
  int i76;
  int i78;
  int i79;
  int i8;
  int i80;
  int i82;
  int i84;
  int i87;
  int i88;
  int i90;
  int i91;
  int i92;
  int i93;
  int i94;
  int i95;
  int i96;
  int i97;
  int i98;
  int j;
  signed char c_I[121];
  if (!((int)flight_phase)) {
    double ST[9];
    double d_a[9];
    double a[3];
    double a_norm;
    double b_a;
    double b_filtered;
    double c_a;
    double d7;
    double d8;
    double d9;
    double filtered;
    int i;
    int i6;
    int i9;
    if (sens_in->board_accel.status) {
      sens_filt->board_accel[0] = (0.0005 * sens_in->board_accel.meas[0]) +
                                  (0.9995 * sens_filt->board_accel[0]);
      sens_filt->board_accel[1] = (0.0005 * sens_in->board_accel.meas[1]) +
                                  (0.9995 * sens_filt->board_accel[1]);
      sens_filt->board_accel[2] = (0.0005 * sens_in->board_accel.meas[2]) +
                                  (0.9995 * sens_filt->board_accel[2]);
    }
    if (sens_in->board_gyro.status) {
      sens_filt->board_gyro[0] = (0.0005 * sens_in->board_gyro.meas[0]) +
                                 (0.9995 * sens_filt->board_gyro[0]);
      sens_filt->board_gyro[1] = (0.0005 * sens_in->board_gyro.meas[1]) +
                                 (0.9995 * sens_filt->board_gyro[1]);
      sens_filt->board_gyro[2] = (0.0005 * sens_in->board_gyro.meas[2]) +
                                 (0.9995 * sens_filt->board_gyro[2]);
    }
    if (sens_in->mti_accel.status) {
      sens_filt->mti_accel[0] = (0.0005 * sens_in->mti_accel.meas[0]) +
                                (0.9995 * sens_filt->mti_accel[0]);
      sens_filt->mti_accel[1] = (0.0005 * sens_in->mti_accel.meas[1]) +
                                (0.9995 * sens_filt->mti_accel[1]);
      sens_filt->mti_accel[2] = (0.0005 * sens_in->mti_accel.meas[2]) +
                                (0.9995 * sens_filt->mti_accel[2]);
    }
    if (sens_in->mti_gyro.status) {
      sens_filt->mti_gyro[0] = (0.0005 * sens_in->mti_gyro.meas[0]) +
                               (0.9995 * sens_filt->mti_gyro[0]);
      sens_filt->mti_gyro[1] = (0.0005 * sens_in->mti_gyro.meas[1]) +
                               (0.9995 * sens_filt->mti_gyro[1]);
      sens_filt->mti_gyro[2] = (0.0005 * sens_in->mti_gyro.meas[2]) +
                               (0.9995 * sens_filt->mti_gyro[2]);
    }
    if (sens_in->ad_accel.status) {
      sens_filt->ad_accel[0] = (0.0005 * sens_in->ad_accel.meas[0]) +
                               (0.9995 * sens_filt->ad_accel[0]);
      sens_filt->ad_accel[1] = (0.0005 * sens_in->ad_accel.meas[1]) +
                               (0.9995 * sens_filt->ad_accel[1]);
      sens_filt->ad_accel[2] = (0.0005 * sens_in->ad_accel.meas[2]) +
                               (0.9995 * sens_filt->ad_accel[2]);
    }
    if (sens_in->ad_gyro.status) {
      sens_filt->ad_gyro[0] = (0.0005 * sens_in->ad_gyro.meas[0]) +
                              (0.9995 * sens_filt->ad_gyro[0]);
      sens_filt->ad_gyro[1] = (0.0005 * sens_in->ad_gyro.meas[1]) +
                              (0.9995 * sens_filt->ad_gyro[1]);
      sens_filt->ad_gyro[2] = (0.0005 * sens_in->ad_gyro.meas[2]) +
                              (0.9995 * sens_filt->ad_gyro[2]);
    }
    filtered = sens_filt->board_baro;
    if (sens_in->board_baro.status) {
      filtered = (0.0005 * sens_in->board_baro.meas) +
                 (0.9995 * sens_filt->board_baro);
    }
    sens_filt->board_baro = filtered;
    if (sens_in->board_mag.status) {
      sens_filt->board_mag[0] = (0.0005 * sens_in->board_mag.meas[0]) +
                                (0.9995 * sens_filt->board_mag[0]);
      sens_filt->board_mag[1] = (0.0005 * sens_in->board_mag.meas[1]) +
                                (0.9995 * sens_filt->board_mag[1]);
      sens_filt->board_mag[2] = (0.0005 * sens_in->board_mag.meas[2]) +
                                (0.9995 * sens_filt->board_mag[2]);
    }
    b_filtered = sens_filt->mti_baro;
    if (sens_in->mti_baro.status) {
      b_filtered =
          (0.0005 * sens_in->mti_baro.meas) + (0.9995 * sens_filt->mti_baro);
    }
    sens_filt->mti_baro = b_filtered;
    if (sens_in->mti_mag.status) {
      sens_filt->mti_mag[0] = (0.0005 * sens_in->mti_mag.meas[0]) +
                              (0.9995 * sens_filt->mti_mag[0]);
      sens_filt->mti_mag[1] = (0.0005 * sens_in->mti_mag.meas[1]) +
                              (0.9995 * sens_filt->mti_mag[1]);
      sens_filt->mti_mag[2] = (0.0005 * sens_in->mti_mag.meas[2]) +
                              (0.9995 * sens_filt->mti_mag[2]);
    }
    a[0] = (sens_filt->board_accel[0] + sens_filt->mti_accel[0]) +
           sens_filt->ad_accel[0];
    a[1] = (sens_filt->board_accel[1] + sens_filt->mti_accel[1]) +
           sens_filt->ad_accel[1];
    a[2] = (sens_filt->board_accel[2] + sens_filt->mti_accel[2]) +
           sens_filt->ad_accel[2];
    a_norm = b_norm(a);
    if (a_norm < 1.0E-6) {
      q[0] = 1.0;
      q[1] = 0.0;
      q[2] = 0.0;
      q[3] = 0.0;
    } else {
      double b_absxk;
      double b_scale;
      double b_t;
      double qw;
      double qy;
      double qz;
      double y;
      qw = sqrt((0.5 * (a[0] / a_norm)) + 0.5);
      if (qw == 0.0) {
        qy = 1.0;
        qz = 0.0;
      } else {
        qy = (0.5 * (a[2] / a_norm)) / qw;
        qz = (-0.5 * (a[1] / a_norm)) / qw;
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
      q[0] = qw / y;
      q[1] = 0.0 / y;
      q[2] = qy / y;
      q[3] = qz / y;
    }
    x[0] = q[0];
    x[1] = q[1];
    x[2] = q[2];
    x[3] = q[3];
    x[10] = SD->pd->b_param.altitude_initial;
    x[4] = 0.0;
    x[7] = 0.0;
    bias->board_gyro[0] = sens_filt->board_gyro[0];
    bias->mti_gyro[0] = sens_filt->mti_gyro[0];
    bias->ad_gyro[0] = sens_filt->ad_gyro[0];
    x[5] = 0.0;
    x[8] = 0.0;
    bias->board_gyro[1] = sens_filt->board_gyro[1];
    bias->mti_gyro[1] = sens_filt->mti_gyro[1];
    bias->ad_gyro[1] = sens_filt->ad_gyro[1];
    x[6] = 0.0;
    x[9] = 0.0;
    bias->board_gyro[2] = sens_filt->board_gyro[2];
    bias->mti_gyro[2] = sens_filt->mti_gyro[2];
    bias->ad_gyro[2] = sens_filt->ad_gyro[2];
    b_a = (q[0] * q[0]) - (((q[1] * q[1]) + (q[2] * q[2])) + (q[3] * q[3]));
    c_a = 2.0 * q[0];
    i = 0;
    for (i1 = 0; i1 < 3; i1++) {
      double a_tmp;
      a_tmp = 2.0 * q[i1 + 1];
      d_a[i] = (b_a * b_b[i1]) + (a_tmp * q[1]);
      d_a[i + 1] = (b_a * b_b[i1 + 3]) + (a_tmp * q[2]);
      d_a[i + 2] = (b_a * b_b[i1 + 6]) + (a_tmp * q[3]);
      i += 3;
    }
    b_dv[0] = 0.0;
    b_dv[1] = c_a * (-q[3]);
    b_dv[2] = c_a * q[2];
    b_dv[3] = c_a * q[3];
    b_dv[4] = 0.0;
    b_dv[5] = c_a * (-q[1]);
    b_dv[6] = c_a * (-q[2]);
    b_dv[7] = c_a * q[1];
    b_dv[8] = 0.0;
    for (i3 = 0; i3 < 9; i3++) {
      ST[i3] = d_a[i3] - b_dv[i3];
    }
    bias->board_mag_earth[0] = 0.0;
    bias->board_mag_earth[1] = 0.0;
    bias->board_mag_earth[2] = 0.0;
    i6 = 0;
    for (i8 = 0; i8 < 3; i8++) {
      double d6;
      d6 = sens_filt->board_mag[i8];
      bias->board_mag_earth[0] += ST[i6] * d6;
      bias->board_mag_earth[1] += ST[i6 + 1] * d6;
      bias->board_mag_earth[2] += ST[i6 + 2] * d6;
      bias->mti_mag_earth[i8] = 0.0;
      i6 += 3;
    }
    i9 = 0;
    d7 = bias->mti_mag_earth[0];
    d8 = bias->mti_mag_earth[1];
    d9 = bias->mti_mag_earth[2];
    for (i10 = 0; i10 < 3; i10++) {
      double d10;
      d10 = sens_filt->mti_mag[i10];
      d7 += ST[i9] * d10;
      d8 += ST[i9 + 1] * d10;
      d9 += ST[i9 + 2] * d10;
      i9 += 3;
    }
    double t1_pressure;
    bias->mti_mag_earth[2] = d9;
    bias->mti_mag_earth[1] = d8;
    bias->mti_mag_earth[0] = d7;
    t1_pressure =
        airdata_atmos(SD->pd->b_param.altitude_initial, &e_expl_temp,
                      &t1_density, &f_expl_temp, &g_expl_temp, &h_expl_temp);
    bias->board_baro = filtered - t1_pressure;
    bias->mti_baro = b_filtered - t1_pressure;
    *w_status_nav = true;
  } else {
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
    double d_a[9];
    double dv4[9];
    double n_tilde[9];
    double skewed_exp_w_tmp[9];
    double w_exp_tilde[9];
    double b_dv1[4];
    double r_q_tmp[4];
    double C_total_a[3];
    double b_S[3];
    double c_r_q_tmp[3];
    double dn[3];
    double dv2[3];
    double e_x[3];
    double C_ad_w_idx_0;
    double C_total_a_tmp;
    double C_total_a_tmp_tmp;
    double absxk;
    double b;
    double b_C_total_a_tmp_tmp;
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
    double d13;
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
    double d27;
    double d29;
    double d3;
    double d30;
    double d32;
    double d4;
    double d57;
    double d58;
    double d59;
    double d_x;
    double dphi;
    double dphi_tmp;
    double e_a;
    double g_a;
    double h_a;
    double i_a;
    double j_a;
    double l_a;
    double n_a;
    double n_idx_0;
    double n_idx_1;
    double n_idx_2;
    double o_a;
    double p_a;
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
    int i74;
    int i77;
    int i81;
    int i85;
    int i86;
    int k;
    signed char b_I[16];
    signed char w_exp_tilde_tmp[9];
    d = 9.9999999999999981E+9 * ((double)sens_in->ad_gyro.status);
    C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_in->board_accel.status);
    b_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_in->mti_accel.status);
    c_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((double)sens_in->ad_accel.status);
    C_total_a_tmp =
        (C_total_a_tmp_tmp + b_C_total_a_tmp_tmp) + c_C_total_a_tmp_tmp;
    C_total_a[0] = C_total_a_tmp;
    d1 = 9.9999999999999981E+9 * ((double)sens_in->board_gyro.status);
    d2 = 9.9999999999999981E+9 * ((double)sens_in->mti_gyro.status);
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
      dn[0] = x[4] / dphi_tmp;
      dn[1] = x[5] / dphi_tmp;
      dn[2] = x[6] / dphi_tmp;
      n_idx_0 = x[4] / dphi_tmp;
      n_idx_1 = x[5] / dphi_tmp;
      n_idx_2 = x[6] / dphi_tmp;
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
    e_a = sin(b_dphi_tmp);
    b_x = cos(b_dphi_tmp);
    for (i2 = 0; i2 < 9; i2++) {
      w_exp_tilde_tmp[i2] = (signed char)0;
    }
    memset(&b_n_tilde[0], 0, 9U * (sizeof(double)));
    k = 0;
    b_k = 0;
    for (c_k = 0; c_k < 3; c_k++) {
      int i5;
      w_exp_tilde_tmp[k] = (signed char)1;
      i5 = 0;
      for (i7 = 0; i7 < 3; i7++) {
        double d5;
        d5 = n_tilde[i7 + b_k];
        b_n_tilde[b_k] += n_tilde[i5] * d5;
        b_n_tilde[b_k + 1] += n_tilde[i5 + 1] * d5;
        b_n_tilde[b_k + 2] += n_tilde[i5 + 2] * d5;
        i5 += 3;
      }
      k += 4;
      b_k += 3;
    }
    for (i4 = 0; i4 < 9; i4++) {
      w_exp_tilde[i4] = (((double)w_exp_tilde_tmp[i4]) - (e_a * n_tilde[i4])) +
                        ((1.0 - b_x) * b_n_tilde[i4]);
    }
    double f_a;
    f_a = b_norm(&x[7]);
    airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp,
                  &d_expl_temp);
    g_a = (0.5 * t1_density) * (f_a * f_a);
    h_a = SD->pd->c_param.c_aero * SD->pd->c_param.Cn_alpha;
    i_a = (x[0] * x[0]) - (((x[1] * x[1]) + (x[2] * x[2])) + (x[3] * x[3]));
    j_a = 2.0 * x[0];
    for (i11 = 0; i11 < 3; i11++) {
      double b_a_tmp;
      int c_a_tmp;
      int d_a_tmp;
      b_a_tmp = x[i11 + 1];
      d_a[3 * i11] = (i_a * b_b[3 * i11]) + ((2.0 * x[1]) * b_a_tmp);
      c_a_tmp = (3 * i11) + 1;
      d_a[c_a_tmp] = (i_a * b_b[c_a_tmp]) + ((2.0 * x[2]) * b_a_tmp);
      d_a_tmp = (3 * i11) + 2;
      d_a[d_a_tmp] = (i_a * b_b[d_a_tmp]) + ((2.0 * x[3]) * b_a_tmp);
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
    for (i12 = 0; i12 < 9; i12++) {
      S[i12] = d_a[i12] - b_dv[i12];
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
        double d12;
        d12 = SD->pd->c_param.J[i16 + i13];
        b_w_exp_tilde[i13] += w_exp_tilde[i15] * d12;
        b_w_exp_tilde[i13 + 1] += w_exp_tilde[i15 + 1] * d12;
        b_w_exp_tilde[i13 + 2] += w_exp_tilde[i15 + 2] * d12;
        i15 += 3;
      }
      double d11;
      d11 = x[i14 + 4];
      c_w_exp_tilde[0] += b_w_exp_tilde[i13] * d11;
      c_w_exp_tilde[1] += b_w_exp_tilde[i13 + 1] * d11;
      c_w_exp_tilde[2] += b_w_exp_tilde[i13 + 2] * d11;
      i13 += 3;
    }
    dv2[0] = 0.0;
    dv2[1] = g_a * (h_a * sin(atan2(x[9], x[7])));
    dv2[2] = g_a * (h_a * (-sin(atan2(x[8], x[7]))));
    memset(&dv3[0], 0, 3U * (sizeof(double)));
    memset(&b_dt[0], 0, 3U * (sizeof(double)));
    memset(&c_dt[0], 0, 3U * (sizeof(double)));
    i17 = 0;
    d13 = dv3[0];
    d14 = dv3[1];
    d15 = dv3[2];
    d16 = b_dt[0];
    d17 = b_dt[1];
    d18 = b_dt[2];
    d19 = x[7];
    d20 = x[8];
    d21 = x[9];
    d22 = c_dt[0];
    d23 = c_dt[1];
    d24 = c_dt[2];
    for (i18 = 0; i18 < 3; i18++) {
      double d25;
      double d26;
      double d28;
      double d31;
      double d33;
      double d35;
      d25 = SD->pd->c_param.Jinv[i17];
      d26 = c_w_exp_tilde[i18];
      d13 += d25 * d26;
      d28 = dv2[i18];
      d16 += (dt * d25) * d28;
      d31 = S[i17];
      d22 += (dt * d31) * SD->pd->c_param.g[i18];
      d33 = d31 * d19;
      d25 = SD->pd->c_param.Jinv[i17 + 1];
      d14 += d25 * d26;
      d17 += (dt * d25) * d28;
      d31 = S[i17 + 1];
      d23 += (dt * d31) * SD->pd->c_param.g[i18];
      d33 += d31 * d20;
      d25 = SD->pd->c_param.Jinv[i17 + 2];
      d15 += d25 * d26;
      d18 += (dt * d25) * d28;
      d31 = S[i17 + 2];
      d24 += (dt * d31) * SD->pd->c_param.g[i18];
      d33 += d31 * d21;
      d35 = C_total_a[i18];
      c_w_exp_tilde[i18] =
          (((w_exp_tilde[i18] * d19) + (w_exp_tilde[i18 + 3] * d20)) +
           (w_exp_tilde[i18 + 6] * d21)) +
          (dt *
           ((((C_total_a_tmp_tmp / d35) * sens_in->board_accel.meas[i18]) +
             ((b_C_total_a_tmp_tmp / d35) * sens_in->mti_accel.meas[i18])) +
            ((c_C_total_a_tmp_tmp / d35) * sens_in->ad_accel.meas[i18])));
      b_S[i18] = d33;
      i17 += 3;
    }
    memset(&c_q[0], 0, (sizeof(double)) << 2);
    i19 = 0;
    d27 = c_q[0];
    d29 = c_q[1];
    d30 = c_q[2];
    d32 = c_q[3];
    for (i20 = 0; i20 < 4; i20++) {
      double d34;
      d34 = b_dv1[i20];
      d27 += b_q[i19] * d34;
      d29 += b_q[i19 + 1] * d34;
      d30 += b_q[i19 + 2] * d34;
      d32 += b_q[i19 + 3] * d34;
      i19 += 4;
    }
    x_pred[0] = d27;
    x_pred[1] = d29;
    x_pred[2] = d30;
    x_pred[3] = d32;
    x_pred[4] = d13 + d16;
    x_pred[7] = c_w_exp_tilde[0] + d22;
    x_pred[5] = d14 + d17;
    x_pred[8] = c_w_exp_tilde[1] + d23;
    x_pred[6] = d15 + d18;
    x_pred[9] = c_w_exp_tilde[2] + d24;
    x_pred[10] = x[10] + (dt * b_S[0]);
    memset(&F[0], 0, 121U * (sizeof(double)));
    l_a = 0.5 * dt;
    W_dt[0] = 0.0;
    W_dt[4] = l_a * (-x[4]);
    W_dt[8] = l_a * (-x[5]);
    W_dt[12] = l_a * (-x[6]);
    W_dt[1] = l_a * x[4];
    W_dt[5] = 0.0;
    W_dt[9] = l_a * x[6];
    W_dt[13] = l_a * (-x[5]);
    W_dt[2] = l_a * x[5];
    W_dt[6] = l_a * (-x[6]);
    W_dt[10] = 0.0;
    W_dt[14] = l_a * x[4];
    W_dt[3] = l_a * x[6];
    W_dt[7] = l_a * x[5];
    W_dt[11] = l_a * (-x[4]);
    W_dt[15] = 0.0;
    memset(&c_b[0], 0, (sizeof(double)) << 4);
    i21 = 0;
    for (i22 = 0; i22 < 4; i22++) {
      int i23;
      i23 = 0;
      for (i25 = 0; i25 < 4; i25++) {
        double d36;
        d36 = W_dt[i25 + i21];
        c_b[i21] += W_dt[i23] * d36;
        c_b[i21 + 1] += W_dt[i23 + 1] * d36;
        c_b[i21 + 2] += W_dt[i23 + 2] * d36;
        c_b[i21 + 3] += W_dt[i23 + 3] * d36;
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
        double d37;
        d37 = c_b[i30 + e_k];
        b_W_dt[e_k] += W_dt[i28] * d37;
        d_b[e_k] += c_b[i28] * d37;
        b_W_dt[e_k + 1] += W_dt[i28 + 1] * d37;
        d_b[e_k + 1] += c_b[i28 + 1] * d37;
        b_W_dt[e_k + 2] += W_dt[i28 + 2] * d37;
        d_b[e_k + 2] += c_b[i28 + 2] * d37;
        b_W_dt[e_k + 3] += W_dt[i28 + 3] * d37;
        d_b[e_k + 3] += c_b[i28 + 3] * d37;
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
    e_a_tmp = l_a * q[0];
    m_a[0] = e_a_tmp;
    f_a_tmp = l_a * (-q[1]);
    m_a[4] = f_a_tmp;
    g_a_tmp = l_a * (-q[2]);
    m_a[8] = g_a_tmp;
    h_a_tmp = l_a * (-q[3]);
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
    i31 = 0;
    i32 = 0;
    for (i33 = 0; i33 < 3; i33++) {
      F[i31 + 44] = m_a[i32 + 4];
      F[i31 + 45] = m_a[i32 + 5];
      F[i31 + 46] = m_a[i32 + 6];
      F[i31 + 47] = m_a[i32 + 7];
      i31 += 11;
      i32 += 4;
    }
    n_a = (0.5 * SD->pd->d_param.c_aero) * SD->pd->d_param.Cn_alpha;
    airdata_atmos(x[10], &m_expl_temp, &t1_density, &n_expl_temp, &o_expl_temp,
                  &p_expl_temp);
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
        double d38;
        d38 = n_tilde[i38 + i34];
        b_n_tilde[i34] += n_tilde[i36] * d38;
        b_n_tilde[i34 + 1] += n_tilde[i36 + 1] * d38;
        b_n_tilde[i34 + 2] += n_tilde[i36 + 2] * d38;
        i36 += 3;
      }
      i34 += 3;
    }
    for (i37 = 0; i37 < 9; i37++) {
      w_exp_tilde[i37] =
          (((double)w_exp_tilde_tmp[i37]) - (e_a * n_tilde[i37])) +
          ((1.0 - b_x) * b_n_tilde[i37]);
    }
    memset(&b_dv[0], 0, 9U * (sizeof(double)));
    i39 = 0;
    i40 = 0;
    for (i41 = 0; i41 < 3; i41++) {
      int i43;
      i43 = 0;
      for (i45 = 0; i45 < 3; i45++) {
        double d39;
        d39 = w_exp_tilde[i45 + i39];
        b_dv[i39] += SD->pd->d_param.Jinv[i43] * d39;
        b_dv[i39 + 1] += SD->pd->d_param.Jinv[i43 + 1] * d39;
        b_dv[i39 + 2] += SD->pd->d_param.Jinv[i43 + 2] * d39;
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
        double d40;
        d40 = SD->pd->d_param.J[i48 + i42];
        F[i44 + 48] += b_dv[i47] * d40;
        F[i44 + 49] += b_dv[i47 + 1] * d40;
        F[i44 + 50] += b_dv[i47 + 2] * d40;
        i47 += 3;
      }
      i42 += 3;
      i44 += 11;
    }
    b_dv[1] = t1_density * (n_a * x[9]);
    b_dv[4] = 0.0;
    b_dv[7] = t1_density * (n_a * x[7]);
    b_dv[2] = t1_density * (n_a * (-x[8]));
    b_dv[5] = t1_density * (n_a * (-x[7]));
    b_dv[8] = 0.0;
    c_x = 0.0;
    for (i49 = 0; i49 < 3; i49++) {
      double d41;
      double d42;
      double d43;
      int F_tmp;
      b_dv[3 * i49] = 0.0;
      F_tmp = 11 * (i49 + 7);
      d41 = 0.0;
      d42 = 0.0;
      d43 = 0.0;
      for (i50 = 0; i50 < 3; i50++) {
        double d44;
        d44 = b_dv[i50 + (3 * i49)];
        d41 += (dt * SD->pd->d_param.Jinv[3 * i50]) * d44;
        d42 += (dt * SD->pd->d_param.Jinv[(3 * i50) + 1]) * d44;
        d43 += (dt * SD->pd->d_param.Jinv[(3 * i50) + 2]) * d44;
      }
      F[F_tmp + 6] = d43;
      F[F_tmp + 5] = d42;
      F[F_tmp + 4] = d41;
      c_x += x[i49 + 1] * SD->pd->d_param.g[i49];
    }
    d_x = x[0];
    e_x[0] = (x[2] * SD->pd->d_param.g[2]) - (SD->pd->d_param.g[1] * x[3]);
    e_x[1] = (SD->pd->d_param.g[0] * x[3]) - (x[1] * SD->pd->d_param.g[2]);
    e_x[2] = (x[1] * SD->pd->d_param.g[1]) - (SD->pd->d_param.g[0] * x[2]);
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
    o_a = 0.5 * (dt * dt);
    memset(&c_skewed_exp_w_tmp[0], 0, 9U * (sizeof(double)));
    memset(&b_dv[0], 0, 9U * (sizeof(double)));
    r_q_tmp[0] = x[0];
    b_r_q_tmp = 0.0;
    for (i51 = 0; i51 < 3; i51++) {
      double b_F_tmp;
      double d45;
      double d46;
      int c_F_tmp;
      int d_F_tmp;
      int e_F_tmp;
      F[i51 + 7] = dt * (2.0 * ((d_x * SD->pd->d_param.g[i51]) - e_x[i51]));
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
      d45 = c_skewed_exp_w_tmp[3 * i51];
      d46 = b_dv[3 * i51];
      for (i54 = 0; i54 < 3; i54++) {
        double d47;
        double d48;
        int b_skewed_exp_w_tmp_tmp;
        int i55;
        int skewed_exp_w_tmp_tmp;
        i55 = i54 + (3 * i51);
        d47 = b_skewed_exp_w_tmp[i55];
        d48 = skewed_exp_w_tmp[i55];
        d45 += skewed_exp_w_tmp[3 * i54] * d47;
        d46 += (2.0 * b_skewed_exp_w_tmp[3 * i54]) * d48;
        skewed_exp_w_tmp_tmp = (3 * i54) + 1;
        c_skewed_exp_w_tmp[d_F_tmp] +=
            skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d47;
        b_dv[d_F_tmp] += (2.0 * b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp]) * d48;
        b_skewed_exp_w_tmp_tmp = (3 * i54) + 2;
        c_skewed_exp_w_tmp[e_F_tmp] +=
            skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d47;
        b_dv[e_F_tmp] +=
            (2.0 * b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp]) * d48;
      }
      int f_F_tmp;
      int g_F_tmp;
      b_dv[3 * i51] = d46;
      c_skewed_exp_w_tmp[3 * i51] = d45;
      f_F_tmp = 11 * (i51 + 4);
      F[f_F_tmp + 7] = (dt * skewed_exp_w_tmp[3 * i51]) + (o_a * (d45 - d46));
      g_F_tmp = 11 * (i51 + 7);
      F[g_F_tmp + 7] = w_exp_tilde[3 * i51];
      F[f_F_tmp + 8] = (dt * skewed_exp_w_tmp[d_F_tmp]) +
                       (o_a * (c_skewed_exp_w_tmp[d_F_tmp] - b_dv[d_F_tmp]));
      F[g_F_tmp + 8] = w_exp_tilde[d_F_tmp];
      F[f_F_tmp + 9] = (dt * skewed_exp_w_tmp[e_F_tmp]) +
                       (o_a * (c_skewed_exp_w_tmp[e_F_tmp] - b_dv[e_F_tmp]));
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
    double q_a;
    F[10] = d_dt[0];
    F[21] = d_dt[3];
    F[32] = d_dt[6];
    F[43] = d_dt[9];
    p_a = (r_q_tmp[0] * r_q_tmp[0]) -
          (((r_q_tmp[1] * r_q_tmp[1]) + (r_q_tmp[2] * r_q_tmp[2])) +
           (r_q_tmp[3] * r_q_tmp[3]));
    q_a = 2.0 * r_q_tmp[0];
    b_dv[0] = 0.0;
    b_dv[3] = q_a * (-r_q_tmp[3]);
    b_dv[6] = q_a * r_q_tmp[2];
    b_dv[1] = q_a * r_q_tmp[3];
    b_dv[4] = 0.0;
    b_dv[7] = q_a * (-r_q_tmp[1]);
    b_dv[2] = q_a * (-r_q_tmp[2]);
    b_dv[5] = q_a * r_q_tmp[1];
    b_dv[8] = 0.0;
    i56 = 0;
    i57 = 0;
    for (i58 = 0; i58 < 3; i58++) {
      double l_a_tmp;
      l_a_tmp = r_q_tmp[i58 + 1];
      d_a[i56] = (p_a * b_b[i56]) + ((2.0 * r_q_tmp[1]) * l_a_tmp);
      d_a[i56 + 1] = (p_a * b_b[i56 + 1]) + ((2.0 * r_q_tmp[2]) * l_a_tmp);
      d_a[i56 + 2] = (p_a * b_b[i56 + 2]) + ((2.0 * r_q_tmp[3]) * l_a_tmp);
      F[i57 + 87] = dt * (d_a[i56] - b_dv[i56]);
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
        double d49;
        d49 = P[i63 + i59];
        for (i67 = 0; i67 < 11; i67++) {
          int h_F_tmp;
          h_F_tmp = i67 + i59;
          b_F[h_F_tmp] += F[i67 + i61] * d49;
        }
        i61 += 11;
      }
      i59 += 11;
    }
    for (i62 = 0; i62 < 11; i62++) {
      int i64;
      i64 = 0;
      for (i66 = 0; i66 < 11; i66++) {
        double d50;
        int i70;
        d50 = 0.0;
        i70 = 0;
        for (i71 = 0; i71 < 11; i71++) {
          d50 += b_F[i70 + i62] * F[i70 + i66];
          i70 += 11;
        }
        int P_pred_tmp;
        P_pred_tmp = i64 + i62;
        P_pred[P_pred_tmp] = d50 + Q[P_pred_tmp];
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
    mrdiv(&P_pred[44], b_P_pred, K);
    memset(&c_I[0], 0, 121U * (sizeof(signed char)));
    g_k = 0;
    for (h_k = 0; h_k < 11; h_k++) {
      c_I[g_k] = (signed char)1;
      g_k += 12;
    }
    i72 = 0;
    for (i73 = 0; i73 < 4; i73++) {
      for (i75 = 0; i75 < 11; i75++) {
        int E_tmp;
        E_tmp = i75 + i72;
        E[E_tmp] = (double)c_I[E_tmp];
      }
      i72 += 11;
    }
    i74 = 0;
    for (i76 = 0; i76 < 3; i76++) {
      for (i78 = 0; i78 < 11; i78++) {
        int b_E_tmp;
        b_E_tmp = i78 + i74;
        E[b_E_tmp + 44] = ((double)c_I[b_E_tmp + 44]) - K[b_E_tmp];
      }
      i74 += 11;
    }
    i77 = 0;
    for (i79 = 0; i79 < 4; i79++) {
      for (i80 = 0; i80 < 11; i80++) {
        int c_E_tmp;
        c_E_tmp = (i80 + i77) + 77;
        E[c_E_tmp] = (double)c_I[c_E_tmp];
      }
      i77 += 11;
    }
    memset(&b_E[0], 0, 121U * (sizeof(double)));
    i81 = 0;
    for (i82 = 0; i82 < 11; i82++) {
      int i83;
      i83 = 0;
      for (i84 = 0; i84 < 11; i84++) {
        double d51;
        d51 = P_pred[i84 + i81];
        for (i88 = 0; i88 < 11; i88++) {
          int d_E_tmp;
          d_E_tmp = i88 + i81;
          b_E[d_E_tmp] += E[i88 + i83] * d51;
        }
        i83 += 11;
      }
      i81 += 11;
    }
    memset(&b_K[0], 0, 33U * (sizeof(double)));
    i85 = 0;
    i86 = 0;
    for (i87 = 0; i87 < 3; i87++) {
      int i89;
      i89 = 0;
      for (i90 = 0; i90 < 3; i90++) {
        double d52;
        d52 = R[i90 + i85];
        for (i92 = 0; i92 < 11; i92++) {
          int K_tmp;
          K_tmp = i92 + i86;
          b_K[K_tmp] += K[i92 + i89] * d52;
        }
        i89 += 11;
      }
      i85 += 3;
      i86 += 11;
    }
    memset(&c_E[0], 0, 121U * (sizeof(double)));
    memset(&c_K[0], 0, 121U * (sizeof(double)));
    for (i91 = 0; i91 < 11; i91++) {
      for (i94 = 0; i94 < 11; i94++) {
        double d53;
        d53 = E[i91 + (11 * i94)];
        for (i96 = 0; i96 < 11; i96++) {
          int e_E_tmp;
          e_E_tmp = i96 + (11 * i91);
          c_E[e_E_tmp] += b_E[i96 + (11 * i94)] * d53;
        }
      }
      for (i95 = 0; i95 < 3; i95++) {
        double d56;
        d56 = K[i91 + (11 * i95)];
        for (i97 = 0; i97 < 11; i97++) {
          int b_K_tmp;
          b_K_tmp = i97 + (11 * i91);
          c_K[b_K_tmp] += b_K[i97 + (11 * i95)] * d56;
        }
      }
    }
    for (i93 = 0; i93 < 121; i93++) {
      P[i93] = c_E[i93] + c_K[i93];
    }
    double d54;
    double d55;
    d54 = d1 / d3;
    d55 = d2 / d3;
    d57 = ((((d1 / d4) * (sens_in->board_gyro.meas[0] - bias->board_gyro[0])) +
            ((d2 / d4) * (sens_in->mti_gyro.meas[0] - bias->mti_gyro[0]))) +
           (C_ad_w_idx_0 * (sens_in->ad_gyro.meas[0] - bias->ad_gyro[0]))) -
          x_pred[4];
    d58 = (((d54 * (sens_in->board_gyro.meas[1] - bias->board_gyro[1])) +
            (d55 * (sens_in->mti_gyro.meas[1] - bias->mti_gyro[1]))) +
           (d * (sens_in->ad_gyro.meas[1] - bias->ad_gyro[1]))) -
          x_pred[5];
    d59 = (((d54 * (sens_in->board_gyro.meas[2] - bias->board_gyro[2])) +
            (d55 * (sens_in->mti_gyro.meas[2] - bias->mti_gyro[2]))) +
           (d * (sens_in->ad_gyro.meas[2] - bias->ad_gyro[2]))) -
          x_pred[6];
    for (i98 = 0; i98 < 11; i98++) {
      x[i98] = x_pred[i98] +
               (((K[i98] * d57) + (K[i98 + 11] * d58)) + (K[i98 + 22] * d59));
    }
    c_scale = 3.3121686421112381E-170;
    c_absxk = fabs(x[0]);
    if (c_absxk > 3.3121686421112381E-170) {
      b_q_mag = 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / 3.3121686421112381E-170;
      b_q_mag = c_t * c_t;
    }
    c_absxk = fabs(x[1]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = ((b_q_mag * c_t) * c_t) + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    c_absxk = fabs(x[2]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = ((b_q_mag * c_t) * c_t) + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    c_absxk = fabs(x[3]);
    if (c_absxk > c_scale) {
      c_t = c_scale / c_absxk;
      b_q_mag = ((b_q_mag * c_t) * c_t) + 1.0;
      c_scale = c_absxk;
    } else {
      c_t = c_absxk / c_scale;
      b_q_mag += c_t * c_t;
    }
    b_q_mag = c_scale * sqrt(b_q_mag);
    x[0] /= b_q_mag;
    x[1] /= b_q_mag;
    x[2] /= b_q_mag;
    x[3] /= b_q_mag;
    if (sens_in->board_baro.status) {
      memcpy(&f_x[0], &x[0], 11U * (sizeof(double)));
      memcpy(&b_P[0], &P[0], 121U * (sizeof(double)));
      b_ekf_correct(f_x, b_P, sens_in->board_baro.meas, bias->board_baro, x, P);
    }
    if (sens_in->mti_baro.status) {
      memcpy(&g_x[0], &x[0], 11U * (sizeof(double)));
      memcpy(&c_P[0], &P[0], 121U * (sizeof(double)));
      b_ekf_correct(g_x, c_P, sens_in->mti_baro.meas, bias->mti_baro, x, P);
    }
    if (sens_in->board_mag.status) {
      memcpy(&h_x[0], &x[0], 11U * (sizeof(double)));
      memcpy(&d_P[0], &P[0], 121U * (sizeof(double)));
      ekf_correct(h_x, d_P, sens_in->board_mag.meas, bias->board_mag_earth, b_b,
                  x, P);
    }
    if (sens_in->mti_mag.status) {
      memcpy(&i_x[0], &x[0], 11U * (sizeof(double)));
      memcpy(&e_P[0], &P[0], 121U * (sizeof(double)));
      ekf_correct(i_x, e_P, sens_in->mti_mag.meas, bias->mti_mag_earth, b_b, x,
                  P);
    }
    *w_status_nav = false;
  }
  k_a = b_norm(&x[7]);
  airdata_atmos(x[10], &i_expl_temp, &t1_density, &j_expl_temp, &k_expl_temp,
                &l_expl_temp);
  *pdyn = (0.5 * t1_density) * (k_a * k_a);
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
