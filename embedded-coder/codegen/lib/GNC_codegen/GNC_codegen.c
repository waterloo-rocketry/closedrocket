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

static void b_eye(double b_I[9]);

static double b_norm(const double x[3]);

static double b_xnrm2(int n, const double x[11], int ix0);

static void b_xzlascl(double cfrom, double cto, double A[11]);

static void c_eye(double b_I[16]);

static double c_norm(const double x[4]);

static void controller_codegen_entry_init(GNC_codegenStackData *SD);

static void d_eye(double b_I[121]);

static double d_norm(const double x[121]);

static void diag(double d[4]);

static void dynamics_init(GNC_codegenStackData *SD);

static void dynamics_jacobian_init(GNC_codegenStackData *SD);

static void ekf_correct(const double x[11], const double P[121],
                        const double y[3], const double b[3], const double R[9],
                        double x_new[11], double P_new[121]);

static void eye(double b_I[4]);

static void inv(const double x[9], double y[9]);

static void pad_filter_init(GNC_codegenStackData *SD);

static void xaxpy(int n, double a, const double x[11], int ix0, double y[121],
                  int iy0);

static double xnrm2(int n, const double x[121], int ix0);

static double xrotg(double *a, double *b, double *s);

static void xzlascl(double cfrom, double cto, double A[121]);

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
  double b_dv[121];
  double c_E[121];
  double H[11];
  double K[11];
  double b_H[11];
  double b_P[11];
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
  int layer_idx_0;
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
  for (i = 0; i < 11; i++) {
    double d;
    d = b_H[i];
    for (i1 = 0; i1 < 11; i1++) {
      d += H[i1] * P[i1 + (11 * i)];
    }
    b_H[i] = d;
    c_H += d * H[i];
  }
  c_b = 1.0 / (c_H + 100.0);
  memset(&b_P[0], 0, 11U * (sizeof(double)));
  for (i2 = 0; i2 < 11; i2++) {
    for (i3 = 0; i3 < 11; i3++) {
      b_P[i3] += P[i3 + (11 * i2)] * H[i2];
    }
  }
  for (i4 = 0; i4 < 11; i4++) {
    K[i4] = b_P[i4] * c_b;
  }
  d_eye(b_dv);
  for (i5 = 0; i5 < 11; i5++) {
    for (i6 = 0; i6 < 11; i6++) {
      int E_tmp;
      E_tmp = i6 + (11 * i5);
      E[E_tmp] = b_dv[E_tmp] - (K[i6] * H[i5]);
    }
  }
  memset(&b_E[0], 0, 121U * (sizeof(double)));
  for (i7 = 0; i7 < 11; i7++) {
    for (i8 = 0; i8 < 11; i8++) {
      double d1;
      d1 = P[i8 + (11 * i7)];
      for (i10 = 0; i10 < 11; i10++) {
        int b_E_tmp;
        b_E_tmp = i10 + (11 * i7);
        b_E[b_E_tmp] += E[i10 + (11 * i8)] * d1;
      }
    }
  }
  memset(&c_E[0], 0, 121U * (sizeof(double)));
  for (i9 = 0; i9 < 11; i9++) {
    for (i11 = 0; i11 < 11; i11++) {
      double d2;
      d2 = E[i9 + (11 * i11)];
      for (i13 = 0; i13 < 11; i13++) {
        int c_E_tmp;
        c_E_tmp = i13 + (11 * i9);
        c_E[c_E_tmp] += b_E[i13 + (11 * i11)] * d2;
      }
      b_K[i11 + (11 * i9)] = (K[i11] * 100.0) * K[i9];
    }
  }
  for (i12 = 0; i12 < 121; i12++) {
    P_new[i12] = c_E[i12] + b_K[i12];
  }
  for (b_i = 0; b_i < 11; b_i++) {
    x_new[b_i] = x[b_i] + (K[b_i] * b_b);
  }
  double q_mag;
  q_mag = c_norm(&x_new[0]);
  x_new[0] /= q_mag;
  x_new[1] /= q_mag;
  x_new[2] /= q_mag;
  x_new[3] /= q_mag;
}

static void b_eye(double b_I[9]) {
  memset(&b_I[0], 0, 9U * (sizeof(double)));
  b_I[0] = 1.0;
  b_I[4] = 1.0;
  b_I[8] = 1.0;
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

static double b_xnrm2(int n, const double x[11], int ix0) {
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (k = ix0; k < kend; k++) {
        double absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = ((y * t) * t) + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

static void b_xzlascl(double cfrom, double cto, double A[11]) {
  double cfromc;
  double ctoc;
  int i;
  bool notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }
    for (i = 0; i < 11; i++) {
      A[i] *= mul;
    }
  }
}

static void c_eye(double b_I[16]) {
  memset(&b_I[0], 0, 16U * (sizeof(double)));
  b_I[0] = 1.0;
  b_I[5] = 1.0;
  b_I[10] = 1.0;
  b_I[15] = 1.0;
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

static void d_eye(double b_I[121]) {
  int k;
  memset(&b_I[0], 0, 121U * (sizeof(double)));
  for (k = 0; k < 11; k++) {
    b_I[k + (11 * k)] = 1.0;
  }
}

static double d_norm(const double x[121]) {
  double A[121];
  double S[11];
  double e[11];
  double s[11];
  double work[11];
  double a__1;
  double a__2;
  double a__3;
  double anrm;
  double b_f;
  double b_sn;
  double c_sn;
  double cscale;
  double d_sn;
  double f;
  double rt;
  double sn;
  double snorm;
  int b_jj;
  int b_k;
  int b_q;
  int c_ii;
  int c_jj;
  int c_k;
  int d_k;
  int e_k;
  int f_k;
  int g_k;
  int h_k;
  int i_k;
  int iter;
  int j_k;
  int jj;
  int k;
  int k_k;
  int m;
  int q;
  int qp1;
  bool doscale;
  memcpy(&A[0], &x[0], 121U * (sizeof(double)));
  memset(&s[0], 0, 11U * (sizeof(double)));
  memset(&e[0], 0, 11U * (sizeof(double)));
  memset(&work[0], 0, 11U * (sizeof(double)));
  doscale = false;
  anrm = 0.0;
  for (k = 0; k < 121; k++) {
    double absxk;
    absxk = fabs(x[k]);
    if (absxk > anrm) {
      anrm = absxk;
    }
  }
  cscale = anrm;
  if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
    doscale = true;
    cscale = 6.7178761075670888E-139;
    xzlascl(anrm, cscale, A);
  } else if (anrm > 1.4885657073574029E+138) {
    doscale = true;
    cscale = 1.4885657073574029E+138;
    xzlascl(anrm, cscale, A);
  }
  for (q = 0; q < 10; q++) {
    double nrm;
    int qq;
    int qq_tmp;
    bool apply_transform;
    qp1 = q + 2;
    qq_tmp = q + (11 * q);
    qq = qq_tmp + 1;
    apply_transform = false;
    nrm = xnrm2(11 - q, A, qq_tmp + 1);
    if (nrm > 0.0) {
      double d1;
      apply_transform = true;
      if (A[qq_tmp] < 0.0) {
        d1 = -nrm;
      } else {
        d1 = nrm;
      }
      s[q] = d1;
      if (fabs(d1) >= 1.0020841800044864E-292) {
        double a;
        int i1;
        a = 1.0 / d1;
        i1 = (qq_tmp - q) + 11;
        for (d_k = qq; d_k <= i1; d_k++) {
          A[d_k - 1] *= a;
        }
      } else {
        int i;
        i = (qq_tmp - q) + 11;
        for (b_k = qq; b_k <= i; b_k++) {
          A[b_k - 1] /= s[q];
        }
      }
      A[qq_tmp]++;
      s[q] = -s[q];
    } else {
      s[q] = 0.0;
    }
    for (jj = qp1; jj < 12; jj++) {
      int qjj;
      qjj = q + (11 * (jj - 1));
      if (apply_transform) {
        double b_d;
        double c_a;
        int n;
        n = 10 - q;
        b_d = 0.0;
        for (c_k = 0; c_k <= n; c_k++) {
          b_d += A[qq_tmp + c_k] * A[qjj + c_k];
        }
        c_a = -(b_d / A[qq_tmp]);
        if (c_a != 0.0) {
          int i2;
          i2 = 11 - q;
          for (g_k = 0; g_k < i2; g_k++) {
            int A_tmp;
            A_tmp = qjj + g_k;
            A[A_tmp] += c_a * A[qq_tmp + g_k];
          }
        }
      }
      e[jj - 1] = A[qjj];
    }
    if ((q + 1) <= 9) {
      nrm = b_xnrm2(10 - q, e, q + 2);
      if (nrm == 0.0) {
        e[q] = 0.0;
      } else {
        double b_a;
        if (e[q + 1] < 0.0) {
          e[q] = -nrm;
        } else {
          e[q] = nrm;
        }
        b_a = e[q];
        if (fabs(e[q]) >= 1.0020841800044864E-292) {
          double d_a;
          d_a = 1.0 / e[q];
          for (f_k = qp1; f_k < 12; f_k++) {
            e[f_k - 1] *= d_a;
          }
        } else {
          for (e_k = qp1; e_k < 12; e_k++) {
            e[e_k - 1] /= b_a;
          }
        }
        e[q + 1]++;
        e[q] = -e[q];
        for (c_ii = qp1; c_ii < 12; c_ii++) {
          work[c_ii - 1] = 0.0;
        }
        for (b_jj = qp1; b_jj < 12; b_jj++) {
          double d4;
          d4 = e[b_jj - 1];
          if (d4 != 0.0) {
            int i4;
            int ix;
            ix = q + (11 * (b_jj - 1));
            i4 = 10 - q;
            for (j_k = 0; j_k < i4; j_k++) {
              int work_tmp;
              work_tmp = (q + j_k) + 1;
              work[work_tmp] += d4 * A[(ix + j_k) + 1];
            }
          }
        }
        for (c_jj = qp1; c_jj < 12; c_jj++) {
          xaxpy(10 - q, (-e[c_jj - 1]) / e[q + 1], work, q + 2, A,
                (q + (11 * (c_jj - 1))) + 2);
        }
      }
    }
  }
  m = 9;
  s[10] = A[120];
  e[9] = A[119];
  e[10] = 0.0;
  iter = 0;
  snorm = 0.0;
  for (b_q = 0; b_q < 11; b_q++) {
    double d;
    double d2;
    d = s[b_q];
    d2 = d;
    if (d != 0.0) {
      rt = fabs(d);
      d2 = rt;
      s[b_q] = rt;
      if ((b_q + 1) < 11) {
        e[b_q] /= d / rt;
      }
    }
    if ((b_q + 1) < 11) {
      double d3;
      d3 = e[b_q];
      if (d3 != 0.0) {
        rt = fabs(d3);
        e[b_q] = rt;
        s[b_q + 1] *= rt / d3;
      }
    }
    snorm = fmax(snorm, fmax(fabs(d2), fabs(e[b_q])));
  }
  while (((m + 2) > 0) && (iter < 75)) {
    int c_q;
    int ii;
    int kase;
    ii = m + 1;
    int exitg1;
    do {
      exitg1 = 0;
      c_q = ii;
      if (ii == 0) {
        exitg1 = 1;
      } else {
        double ztest0;
        ztest0 = fabs(e[ii - 1]);
        if (((ztest0 <=
              (2.2204460492503131E-16 * (fabs(s[ii - 1]) + fabs(s[ii])))) ||
             (ztest0 <= 1.0020841800044864E-292)) ||
            ((iter > 20) && (ztest0 <= (2.2204460492503131E-16 * snorm)))) {
          e[ii - 1] = 0.0;
          exitg1 = 1;
        } else {
          ii--;
        }
      }
    } while (exitg1 == 0);
    if (ii == (m + 1)) {
      kase = 4;
    } else {
      int b_ii;
      int qs;
      bool exitg2;
      qs = m + 2;
      b_ii = m + 2;
      exitg2 = false;
      while ((!((int)exitg2)) && (b_ii >= ii)) {
        qs = b_ii;
        if (b_ii == ii) {
          exitg2 = true;
        } else {
          double test;
          double ztest;
          test = 0.0;
          if (b_ii < (m + 2)) {
            test = fabs(e[b_ii - 1]);
          }
          if (b_ii > (ii + 1)) {
            test += fabs(e[b_ii - 2]);
          }
          ztest = fabs(s[b_ii - 1]);
          if ((ztest <= (2.2204460492503131E-16 * test)) ||
              (ztest <= 1.0020841800044864E-292)) {
            s[b_ii - 1] = 0.0;
            exitg2 = true;
          } else {
            b_ii--;
          }
        }
      }
      if (qs == ii) {
        kase = 3;
      } else if (qs == (m + 2)) {
        kase = 1;
      } else {
        kase = 2;
        c_q = qs;
      }
    }
    switch (kase) {
    case 1: {
      int i3;
      f = e[m];
      e[m] = 0.0;
      i3 = m + 1;
      for (i_k = i3; i_k >= (c_q + 1); i_k--) {
        double cs;
        cs = xrotg(&s[i_k - 1], &f, &sn);
        if (i_k > (c_q + 1)) {
          double b_f_tmp;
          b_f_tmp = e[i_k - 2];
          f = (-sn) * b_f_tmp;
          e[i_k - 2] = b_f_tmp * cs;
        }
      }
    } break;
    case 2: {
      f = e[c_q - 1];
      e[c_q - 1] = 0.0;
      for (h_k = c_q + 1; h_k <= (m + 2); h_k++) {
        double b_cs;
        double f_tmp;
        a__1 = f;
        b_cs = xrotg(&s[h_k - 1], &a__1, &b_sn);
        f_tmp = e[h_k - 1];
        f = (-b_sn) * f_tmp;
        e[h_k - 1] = f_tmp * b_cs;
      }
    } break;
    case 3: {
      double b;
      double c;
      double emm1;
      double g;
      double scale;
      double scale_tmp;
      double shift;
      double sm;
      double smm1;
      double sqds;
      int mm1;
      mm1 = m + 1;
      scale_tmp = s[m + 1];
      scale = fmax(fmax(fmax(fmax(fabs(scale_tmp), fabs(s[m])), fabs(e[m])),
                        fabs(s[c_q])),
                   fabs(e[c_q]));
      sm = scale_tmp / scale;
      smm1 = s[m] / scale;
      emm1 = e[m] / scale;
      sqds = s[c_q] / scale;
      b = (((smm1 + sm) * (smm1 - sm)) + (emm1 * emm1)) / 2.0;
      c = sm * emm1;
      c *= c;
      if ((b != 0.0) || (c != 0.0)) {
        shift = sqrt((b * b) + c);
        if (b < 0.0) {
          shift = -shift;
        }
        shift = c / (b + shift);
      } else {
        shift = 0.0;
      }
      f = ((sqds + sm) * (sqds - sm)) + shift;
      g = sqds * (e[c_q] / scale);
      for (k_k = c_q + 1; k_k <= mm1; k_k++) {
        double c_cs;
        double c_f_tmp;
        double d_cs;
        double d_f_tmp;
        double e_f_tmp;
        b_f = f;
        a__2 = g;
        c_cs = xrotg(&b_f, &a__2, &c_sn);
        if (k_k > (c_q + 1)) {
          e[k_k - 2] = b_f;
        }
        c_f_tmp = e[k_k - 1];
        d_f_tmp = s[k_k - 1];
        e[k_k - 1] = (c_cs * c_f_tmp) - (c_sn * d_f_tmp);
        a__3 = c_sn * s[k_k];
        s[k_k] *= c_cs;
        s[k_k - 1] = (c_cs * d_f_tmp) + (c_sn * c_f_tmp);
        d_cs = xrotg(&s[k_k - 1], &a__3, &d_sn);
        e_f_tmp = e[k_k - 1];
        f = (d_cs * e_f_tmp) + (d_sn * s[k_k]);
        s[k_k] = ((-d_sn) * e_f_tmp) + (d_cs * s[k_k]);
        g = d_sn * e[k_k];
        e[k_k] *= d_cs;
      }
      e[m] = f;
      iter++;
    } break;
    default:
      if (s[c_q] < 0.0) {
        s[c_q] = -s[c_q];
      }
      qp1 = c_q + 1;
      while (((c_q + 1) < 11) && (s[c_q] < s[qp1])) {
        rt = s[c_q];
        s[c_q] = s[qp1];
        s[qp1] = rt;
        c_q = qp1;
        qp1++;
      }
      iter = 0;
      m--;
      break;
    }
  }
  memcpy(&S[0], &s[0], 11U * (sizeof(double)));
  if (doscale) {
    b_xzlascl(cscale, anrm, S);
  }
  return S[0];
}

static void diag(double d[4]) {
  d[1] = 0.0;
  d[2] = 0.0;
  d[0] = 1.0E-5;
  d[3] = 1.0E-9;
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
  double b_dv1[121];
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
  double b_a;
  double b_x;
  double d12;
  double d13;
  double d14;
  double d15;
  double d16;
  double d17;
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
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 11; i2++) {
      y_tmp[i2 + (11 * i1)] = H[i1 + (3 * i2)];
    }
  }
  memset(&b_H[0], 0, 33U * (sizeof(double)));
  for (i3 = 0; i3 < 11; i3++) {
    double d;
    int e_H_tmp;
    int f_H_tmp;
    d = b_H[3 * i3];
    e_H_tmp = (3 * i3) + 1;
    f_H_tmp = (3 * i3) + 2;
    for (i5 = 0; i5 < 11; i5++) {
      double d1;
      d1 = P[i5 + (11 * i3)];
      d += H[3 * i5] * d1;
      b_H[e_H_tmp] += H[(3 * i5) + 1] * d1;
      b_H[f_H_tmp] += H[(3 * i5) + 2] * d1;
    }
    b_H[3 * i3] = d;
  }
  memset(&b_P[0], 0, 33U * (sizeof(double)));
  for (i4 = 0; i4 < 3; i4++) {
    for (i6 = 0; i6 < 3; i6++) {
      double d2;
      d2 = 0.0;
      for (i8 = 0; i8 < 11; i8++) {
        d2 += b_H[i4 + (3 * i8)] * y_tmp[i8 + (11 * i6)];
      }
      int g_H_tmp;
      g_H_tmp = i4 + (3 * i6);
      c_H[g_H_tmp] = d2 + R[g_H_tmp];
    }
    for (i7 = 0; i7 < 11; i7++) {
      double d3;
      d3 = y_tmp[i7 + (11 * i4)];
      for (i10 = 0; i10 < 11; i10++) {
        int P_tmp;
        P_tmp = i10 + (11 * i4);
        b_P[P_tmp] += P[i10 + (11 * i7)] * d3;
      }
    }
  }
  inv(c_H, b_dv);
  memset(&K[0], 0, 33U * (sizeof(double)));
  for (i9 = 0; i9 < 3; i9++) {
    for (i11 = 0; i11 < 3; i11++) {
      double d4;
      d4 = b_dv[i11 + (3 * i9)];
      for (i13 = 0; i13 < 11; i13++) {
        int K_tmp;
        K_tmp = i13 + (11 * i9);
        K[K_tmp] += b_P[i13 + (11 * i11)] * d4;
      }
    }
  }
  d_eye(b_dv1);
  for (i12 = 0; i12 < 11; i12++) {
    double d5;
    double d6;
    double d7;
    d5 = K[i12];
    d6 = K[i12 + 11];
    d7 = K[i12 + 22];
    for (i15 = 0; i15 < 11; i15++) {
      int E_tmp;
      E_tmp = i12 + (11 * i15);
      E[E_tmp] = b_dv1[E_tmp] - (((d5 * H[3 * i15]) + (d6 * H[(3 * i15) + 1])) +
                                 (d7 * H[(3 * i15) + 2]));
    }
  }
  memset(&b_E[0], 0, 121U * (sizeof(double)));
  for (i14 = 0; i14 < 11; i14++) {
    for (i16 = 0; i16 < 11; i16++) {
      double d8;
      d8 = P[i16 + (11 * i14)];
      for (i18 = 0; i18 < 11; i18++) {
        int b_E_tmp;
        b_E_tmp = i18 + (11 * i14);
        b_E[b_E_tmp] += E[i18 + (11 * i16)] * d8;
      }
    }
  }
  memset(&b_K[0], 0, 33U * (sizeof(double)));
  for (i17 = 0; i17 < 3; i17++) {
    for (i19 = 0; i19 < 3; i19++) {
      double d9;
      d9 = R[i19 + (3 * i17)];
      for (i20 = 0; i20 < 11; i20++) {
        int b_K_tmp;
        b_K_tmp = i20 + (11 * i17);
        b_K[b_K_tmp] += K[i20 + (11 * i19)] * d9;
      }
    }
  }
  memset(&c_E[0], 0, 121U * (sizeof(double)));
  memset(&c_K[0], 0, 121U * (sizeof(double)));
  for (i21 = 0; i21 < 11; i21++) {
    for (i22 = 0; i22 < 11; i22++) {
      double d10;
      d10 = E[i21 + (11 * i22)];
      for (i25 = 0; i25 < 11; i25++) {
        int c_E_tmp;
        c_E_tmp = i25 + (11 * i21);
        c_E[c_E_tmp] += b_E[i25 + (11 * i22)] * d10;
      }
    }
    for (i24 = 0; i24 < 3; i24++) {
      double d11;
      d11 = K[i21 + (11 * i24)];
      for (i27 = 0; i27 < 11; i27++) {
        int c_K_tmp;
        c_K_tmp = i27 + (11 * i21);
        c_K[c_K_tmp] += b_K[i27 + (11 * i24)] * d11;
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
    c_a[3 * i26] = (a * ((double)iv[3 * i26])) + ((2.0 * x[1]) * a_tmp);
    b_a_tmp = (3 * i26) + 1;
    c_a[b_a_tmp] = (a * ((double)iv[b_a_tmp])) + ((2.0 * x[2]) * a_tmp);
    c_a_tmp = (3 * i26) + 2;
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
  for (i28 = 0; i28 < 9; i28++) {
    c_a[i28] -= b_dv[i28];
  }
  d12 = b[0];
  d13 = b[1];
  d14 = b[2];
  for (i29 = 0; i29 < 3; i29++) {
    b_y[i29] = y[i29] - (((c_a[i29] * d12) + (c_a[i29 + 3] * d13)) +
                         (c_a[i29 + 6] * d14));
  }
  d15 = b_y[0];
  d16 = b_y[1];
  d17 = b_y[2];
  for (i30 = 0; i30 < 11; i30++) {
    x_new[i30] =
        x[i30] + (((K[i30] * d15) + (K[i30 + 11] * d16)) + (K[i30 + 22] * d17));
  }
  double q_mag;
  q_mag = c_norm(&x_new[0]);
  x_new[0] /= q_mag;
  x_new[1] /= q_mag;
  x_new[2] /= q_mag;
  x_new[3] /= q_mag;
}

static void eye(double b_I[4]) {
  b_I[1] = 0.0;
  b_I[2] = 0.0;
  b_I[0] = 1.0;
  b_I[3] = 1.0;
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

static void xaxpy(int n, double a, const double x[11], int ix0, double y[121],
                  int iy0) {
  int k;
  if ((n >= 1) && (a != 0.0)) {
    for (k = 0; k < n; k++) {
      int i;
      i = (iy0 + k) - 1;
      y[i] += a * x[(ix0 + k) - 1];
    }
  }
}

static double xnrm2(int n, const double x[121], int ix0) {
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (k = ix0; k < kend; k++) {
        double absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = ((y * t) * t) + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

static double xrotg(double *a, double *b, double *s) {
  double absa;
  double absb;
  double b_c;
  double b_s;
  double c;
  double roe;
  double scale;
  roe = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    roe = *a;
  }
  scale = absa + absb;
  if (scale == 0.0) {
    b_s = 0.0;
    b_c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    double ads;
    double bds;
    double r;
    ads = absa / scale;
    bds = absb / scale;
    r = scale * sqrt((ads * ads) + (bds * bds));
    if (roe < 0.0) {
      r = -r;
    }
    b_c = (*a) / r;
    b_s = (*b) / r;
    if (absa > absb) {
      *b = b_s;
    } else if (b_c != 0.0) {
      *b = 1.0 / b_c;
    } else {
      *b = 1.0;
    }
    *a = r;
  }
  c = b_c;
  *s = b_s;
  return c;
}

static void xzlascl(double cfrom, double cto, double A[121]) {
  double cfromc;
  double ctoc;
  int i;
  bool notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }
    for (i = 0; i < 121; i++) {
      A[i] *= mul;
    }
  }
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
  double b_tmp;
  double b_tmp_tmp;
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
  pdyn_params = pdyn * SD->pd->param.c_canard;
  c_delta = delta / 2.0;
  if (fabs(c_delta) < 0.005) {
    c_delta = 0.0;
  }
  c_delta = (0.75 * ctrl_mem_in->d_old) + (0.25 * c_delta);
  w_dot = (0.75 * ctrl_mem_in->w_dot_old) +
          ((0.25 * (xR[1] - ctrl_mem_in->w_old)) / dt_ctrl);
  r_idx_0 = pdyn_params * c_delta;
  diag(b_dv);
  P[0] = ctrl_mem_in->P_minus[0] + b_dv[0];
  P[1] = ctrl_mem_in->P_minus[1] + b_dv[1];
  P[2] = ctrl_mem_in->P_minus[2] + b_dv[2];
  P[3] = ctrl_mem_in->P_minus[3] + b_dv[3];
  memset(&b_r[0], 0, (sizeof(double)) << 1);
  d = r_idx_0 * P[0];
  d1 = pdyn_params * P[3];
  c_r = (((b_r[0] + d) + (pdyn_params * P[1])) * r_idx_0) +
        (((b_r[1] + (r_idx_0 * P[2])) + d1) * pdyn_params);
  K[0] = (d + (P[2] * pdyn_params)) / (c_r + 1.0);
  K[1] = ((P[1] * r_idx_0) + d1) / (c_r + 1.0);
  b = w_dot - ((r_idx_0 * ctrl_mem_in->coeffs[0]) +
               (pdyn_params * ctrl_mem_in->coeffs[1]));
  eye(b_dv);
  ctrl_mem_out->coeffs[0] = ctrl_mem_in->coeffs[0] + (K[0] * b);
  b_dv1[0] = b_dv[0] - (K[0] * r_idx_0);
  b_dv1[1] = b_dv[1] - (K[1] * r_idx_0);
  ctrl_mem_out->coeffs[1] = ctrl_mem_in->coeffs[1] + (K[1] * b);
  b_dv1[2] = b_dv[2] - (K[0] * pdyn_params);
  b_dv1[3] = b_dv[3] - (K[1] * pdyn_params);
  memset(&b_dv[0], 0, (sizeof(double)) << 2);
  d2 = b_dv1[0];
  d3 = b_dv1[1];
  d4 = b_dv1[2];
  d5 = b_dv1[3];
  for (i = 0; i < 2; i++) {
    double d6;
    double d7;
    double d9;
    int i1;
    d6 = P[2 * i];
    d7 = b_dv[2 * i] + (d2 * d6);
    i1 = (2 * i) + 1;
    d9 = b_dv[i1] + (d3 * d6);
    d6 = P[i1];
    d7 += d4 * d6;
    b_dv[2 * i] = d7;
    d9 += d5 * d6;
    b_dv[i1] = d9;
  }
  memset(&dv2[0], 0, (sizeof(double)) << 2);
  d8 = b_dv[0];
  d10 = b_dv[1];
  d11 = b_dv[2];
  d12 = b_dv[3];
  for (i2 = 0; i2 < 2; i2++) {
    double d13;
    double d14;
    double d15;
    int i3;
    d13 = b_dv1[i2];
    d14 = dv2[2 * i2] + (d8 * d13);
    i3 = (2 * i2) + 1;
    d15 = dv2[i3] + (d10 * d13);
    b_K[2 * i2] = K[0] * K[i2];
    d13 = b_dv1[i2 + 2];
    d14 += d11 * d13;
    dv2[2 * i2] = d14;
    d15 += d12 * d13;
    dv2[i3] = d15;
    b_K[i3] = K[1] * K[i2];
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
  b_tmp_tmp = (1.0 - blend) * 5.0;
  b_tmp = sqrt(b_tmp_tmp);
  K[0] = a * b_tmp;
  *u = fmin(fmax(((K[0] * xR[0]) +
                  ((a * sqrt((2.0 * b_tmp) + (b_tmp_tmp + (blend * 20.0)))) *
                   xR[1])) +
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
                              struct2_T *sens_filt_ret, double *cov_norm,
                              struct6_T *airdata, double roll_state[2]) {
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
  double E[121];
  double F[121];
  double b_E[121];
  double b_F[121];
  double b_P_ret[121];
  double c_E[121];
  double c_K[121];
  double c_P_ret[121];
  double d_P_ret[121];
  double dv6[121];
  double e_P_ret[121];
  double f_P_ret[121];
  double g_P_ret[121];
  double h_P_ret[121];
  double i_P_ret[121];
  double K[33];
  double b_K[33];
  double aBuffer[16];
  double b_W_dt[16];
  double b_aBuffer[16];
  double c_W_dt[16];
  double d_W_dt[16];
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
  int i;
  int i1;
  int i10;
  int i11;
  int i12;
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
  int i7;
  int i8;
  int i9;
  memcpy(&P_ret[0], &P[0], 121U * (sizeof(double)));
  *bias_ret = *bias;
  *sens_filt_ret = *sens_filt;
  if (!((int)flight_phase)) {
    double ST[9];
    double h_a[9];
    double a[3];
    double board_baro_f;
    double d10;
    double d12;
    double d13;
    double d5;
    double d6;
    double d7;
    double i_a;
    double j_a;
    double mti_baro_f;
    double qw;
    double qy;
    double qz;
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
    q[0] = qw;
    q[1] = 0.0;
    q[2] = qy;
    q[3] = qz;
    d6 = c_norm(q);
    d7 = qw / d6;
    q[0] = d7;
    x_ret[0] = d7;
    d7 = 0.0 / d6;
    q[1] = d7;
    x_ret[1] = d7;
    d7 = qy / d6;
    q[2] = d7;
    x_ret[2] = d7;
    d7 = qz / d6;
    q[3] = d7;
    x_ret[3] = d7;
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
    i_a = (q[0] * q[0]) - (((q[1] * q[1]) + (q[2] * q[2])) + (d7 * d7));
    j_a = 2.0 * q[0];
    for (i4 = 0; i4 < 3; i4++) {
      double d_a_tmp;
      d_a_tmp = 2.0 * q[i4 + 1];
      h_a[3 * i4] = (i_a * b_b[i4]) + (d_a_tmp * q[1]);
      h_a[(3 * i4) + 1] = (i_a * b_b[i4 + 3]) + (d_a_tmp * q[2]);
      h_a[(3 * i4) + 2] = (i_a * b_b[i4 + 6]) + (d_a_tmp * d7);
    }
    b_dv[0] = 0.0;
    b_dv[1] = j_a * (-d7);
    b_dv[2] = j_a * q[2];
    b_dv[3] = j_a * d7;
    b_dv[4] = 0.0;
    b_dv[5] = j_a * (-q[1]);
    b_dv[6] = j_a * (-q[2]);
    b_dv[7] = j_a * q[1];
    b_dv[8] = 0.0;
    for (i6 = 0; i6 < 9; i6++) {
      ST[i6] = h_a[i6] - b_dv[i6];
    }
    bias_ret->board_mag_earth[0] = 0.0;
    bias_ret->board_mag_earth[1] = 0.0;
    bias_ret->board_mag_earth[2] = 0.0;
    for (i7 = 0; i7 < 3; i7++) {
      double d11;
      d11 = sens_filt_ret->board_mag_f[i7];
      bias_ret->board_mag_earth[0] += ST[3 * i7] * d11;
      bias_ret->board_mag_earth[1] += ST[(3 * i7) + 1] * d11;
      bias_ret->board_mag_earth[2] += ST[(3 * i7) + 2] * d11;
      bias_ret->mti_mag_earth[i7] = 0.0;
    }
    d10 = bias_ret->mti_mag_earth[0];
    d12 = bias_ret->mti_mag_earth[1];
    d13 = bias_ret->mti_mag_earth[2];
    for (i8 = 0; i8 < 3; i8++) {
      double d14;
      d14 = sens_filt_ret->mti_mag_f[i8];
      d10 += ST[3 * i8] * d14;
      d12 += ST[(3 * i8) + 1] * d14;
      d13 += ST[(3 * i8) + 2] * d14;
    }
    double t1_pressure;
    bias_ret->mti_mag_earth[2] = d13;
    bias_ret->mti_mag_earth[1] = d12;
    bias_ret->mti_mag_earth[0] = d10;
    t1_pressure =
        airdata_atmos(SD->pd->b_param.elevation, &e_expl_temp, &t1_density,
                      &f_expl_temp, &g_expl_temp, &h_expl_temp);
    bias_ret->board_baro = board_baro_f - t1_pressure;
    bias_ret->mti_baro = mti_baro_f - t1_pressure;
  } else {
    double P_pred[121];
    double W_dt[16];
    double b_q[16];
    double dv4[16];
    double l_a[16];
    double d_dt[12];
    double x_pred[11];
    double S[9];
    double b_P_pred[9];
    double b_skewed_exp_w_tmp[9];
    double dv5[9];
    double h_a[9];
    double n_tilde[9];
    double skewed_exp_w_tmp[9];
    double w_exp_tilde[9];
    double w_exp_tilde_tmp[9];
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
    double b;
    double b_C_total_a_tmp_tmp;
    double b_a;
    double b_dphi_tmp;
    double b_r_q_tmp;
    double b_x;
    double c_C_total_a_tmp_tmp;
    double c_x;
    double d;
    double d1;
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
    double d28;
    double d3;
    double d31;
    double d32;
    double d34;
    double d35;
    double d4;
    double d71;
    double d72;
    double d73;
    double d_a;
    double dphi;
    double dphi_tmp;
    double e_a;
    double f_a;
    double g_a;
    double k_a;
    double m_a;
    double n_a;
    double n_idx_0;
    double n_idx_1;
    double n_idx_2;
    double o_a;
    double q_mag;
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
    b_eye(w_exp_tilde_tmp);
    memset(&b_n_tilde[0], 0, 9U * (sizeof(double)));
    for (i = 0; i < 3; i++) {
      double d8;
      int b_n_tilde_tmp;
      int n_tilde_tmp;
      d8 = b_n_tilde[3 * i];
      n_tilde_tmp = (3 * i) + 1;
      b_n_tilde_tmp = (3 * i) + 2;
      for (i2 = 0; i2 < 3; i2++) {
        double d9;
        d9 = n_tilde[i2 + (3 * i)];
        d8 += n_tilde[3 * i2] * d9;
        b_n_tilde[n_tilde_tmp] += n_tilde[(3 * i2) + 1] * d9;
        b_n_tilde[b_n_tilde_tmp] += n_tilde[(3 * i2) + 2] * d9;
      }
      b_n_tilde[3 * i] = d8;
    }
    for (i1 = 0; i1 < 9; i1++) {
      w_exp_tilde[i1] = (w_exp_tilde_tmp[i1] - (b_a * n_tilde[i1])) +
                        ((1.0 - b_x) * b_n_tilde[i1]);
    }
    double c_a;
    c_a = b_norm(&x[7]);
    airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp,
                  &d_expl_temp);
    d_a = (0.5 * t1_density) * (c_a * c_a);
    e_a = SD->pd->c_param.c_aero * SD->pd->c_param.Cn_alpha;
    f_a = (x[0] * x[0]) - (((x[1] * x[1]) + (x[2] * x[2])) + (x[3] * x[3]));
    g_a = 2.0 * x[0];
    for (i3 = 0; i3 < 3; i3++) {
      double a_tmp;
      int b_a_tmp;
      int c_a_tmp;
      a_tmp = x[i3 + 1];
      h_a[3 * i3] = (f_a * b_b[3 * i3]) + ((2.0 * x[1]) * a_tmp);
      b_a_tmp = (3 * i3) + 1;
      h_a[b_a_tmp] = (f_a * b_b[b_a_tmp]) + ((2.0 * x[2]) * a_tmp);
      c_a_tmp = (3 * i3) + 2;
      h_a[c_a_tmp] = (f_a * b_b[c_a_tmp]) + ((2.0 * x[3]) * a_tmp);
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
    for (i5 = 0; i5 < 9; i5++) {
      S[i5] = h_a[i5] - b_dv[i5];
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
    for (i9 = 0; i9 < 3; i9++) {
      double d15;
      int b_w_exp_tilde_tmp;
      int c_w_exp_tilde_tmp;
      b_dv1[i9 + 1] = dn[i9] * b;
      d15 = b_w_exp_tilde[3 * i9];
      b_w_exp_tilde_tmp = (3 * i9) + 1;
      c_w_exp_tilde_tmp = (3 * i9) + 2;
      for (i10 = 0; i10 < 3; i10++) {
        double d16;
        d16 = SD->pd->c_param.J[i10 + (3 * i9)];
        d15 += w_exp_tilde[3 * i10] * d16;
        b_w_exp_tilde[b_w_exp_tilde_tmp] += w_exp_tilde[(3 * i10) + 1] * d16;
        b_w_exp_tilde[c_w_exp_tilde_tmp] += w_exp_tilde[(3 * i10) + 2] * d16;
      }
      double d_w_exp_tilde_tmp;
      b_w_exp_tilde[3 * i9] = d15;
      d_w_exp_tilde_tmp = x[i9 + 4];
      c_w_exp_tilde[0] += d15 * d_w_exp_tilde_tmp;
      c_w_exp_tilde[1] += b_w_exp_tilde[(3 * i9) + 1] * d_w_exp_tilde_tmp;
      c_w_exp_tilde[2] += b_w_exp_tilde[(3 * i9) + 2] * d_w_exp_tilde_tmp;
    }
    dv2[0] = 0.0;
    dv2[1] = d_a * (e_a * sin(atan2(x[9], x[7])));
    dv2[2] = d_a * (e_a * (-sin(atan2(x[8], x[7]))));
    memset(&dv3[0], 0, 3U * (sizeof(double)));
    memset(&b_dt[0], 0, 3U * (sizeof(double)));
    memset(&c_dt[0], 0, 3U * (sizeof(double)));
    d17 = dv3[0];
    d18 = dv3[1];
    d19 = dv3[2];
    d20 = b_dt[0];
    d21 = b_dt[1];
    d22 = b_dt[2];
    d23 = c_dt[0];
    d24 = c_dt[1];
    d25 = c_dt[2];
    d26 = x[7];
    d27 = x[8];
    d28 = x[9];
    for (i11 = 0; i11 < 3; i11++) {
      double d29;
      double d30;
      double d33;
      double d36;
      double d37;
      double d39;
      double d40;
      int i13;
      int i14;
      d29 = SD->pd->c_param.Jinv[3 * i11];
      d30 = c_w_exp_tilde[i11];
      d17 += d29 * d30;
      d33 = dv2[i11];
      d20 += (dt * d29) * d33;
      d36 = S[3 * i11];
      d37 = SD->pd->c_param.g[i11];
      d23 += (dt * d36) * d37;
      d39 = d36 * d26;
      i13 = (3 * i11) + 1;
      d29 = SD->pd->c_param.Jinv[i13];
      d18 += d29 * d30;
      d21 += (dt * d29) * d33;
      d36 = S[i13];
      d24 += (dt * d36) * d37;
      d39 += d36 * d27;
      i14 = (3 * i11) + 2;
      d29 = SD->pd->c_param.Jinv[i14];
      d19 += d29 * d30;
      d22 += (dt * d29) * d33;
      d36 = S[i14];
      d25 += (dt * d36) * d37;
      d39 += d36 * d28;
      d40 = C_total_a[i11];
      c_w_exp_tilde[i11] =
          (((w_exp_tilde[i11] * d26) + (w_exp_tilde[i11 + 3] * d27)) +
           (w_exp_tilde[i11 + 6] * d28)) +
          (dt *
           ((((C_total_a_tmp_tmp / d40) * sens_input->board_accel.meas[i11]) +
             ((b_C_total_a_tmp_tmp / d40) * sens_input->mti_accel.meas[i11])) +
            ((c_C_total_a_tmp_tmp / d40) * sens_input->ad_accel.meas[i11])));
      b_S[i11] = d39;
    }
    memset(&c_q[0], 0, (sizeof(double)) << 2);
    d31 = c_q[0];
    d32 = c_q[1];
    d34 = c_q[2];
    d35 = c_q[3];
    for (i12 = 0; i12 < 4; i12++) {
      double d38;
      d38 = b_dv1[i12];
      d31 += b_q[4 * i12] * d38;
      d32 += b_q[(4 * i12) + 1] * d38;
      d34 += b_q[(4 * i12) + 2] * d38;
      d35 += b_q[(4 * i12) + 3] * d38;
    }
    double W_dt_tmp;
    double b_W_dt_tmp;
    double c_W_dt_tmp;
    double d_W_dt_tmp;
    double e_W_dt_tmp;
    double f_W_dt_tmp;
    x_pred[0] = d31;
    x_pred[1] = d32;
    x_pred[2] = d34;
    x_pred[3] = d35;
    x_pred[4] = d17 + d20;
    x_pred[7] = c_w_exp_tilde[0] + d23;
    x_pred[5] = d18 + d21;
    x_pred[8] = c_w_exp_tilde[1] + d24;
    x_pred[6] = d19 + d22;
    x_pred[9] = c_w_exp_tilde[2] + d25;
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
    memset(&aBuffer[0], 0, (sizeof(double)) << 4);
    c_eye(dv4);
    memset(&b_W_dt[0], 0, (sizeof(double)) << 4);
    memset(&c_W_dt[0], 0, (sizeof(double)) << 4);
    memset(&d_W_dt[0], 0, (sizeof(double)) << 4);
    for (i15 = 0; i15 < 4; i15++) {
      double d41;
      double d42;
      double d43;
      double d47;
      int aBuffer_tmp;
      int b_aBuffer_tmp;
      int d_aBuffer_tmp;
      int g_W_dt_tmp;
      int h_W_dt_tmp;
      int i_W_dt_tmp;
      d41 = aBuffer[4 * i15];
      d42 = b_W_dt[4 * i15];
      d43 = c_W_dt[4 * i15];
      aBuffer_tmp = (4 * i15) + 1;
      b_aBuffer_tmp = (4 * i15) + 2;
      d_aBuffer_tmp = (4 * i15) + 3;
      for (i17 = 0; i17 < 4; i17++) {
        double d45;
        double g_aBuffer_tmp;
        double h_aBuffer_tmp;
        double i_aBuffer_tmp;
        double j_aBuffer_tmp;
        d45 = W_dt[i17 + (4 * i15)];
        g_aBuffer_tmp = W_dt[4 * i17] * d45;
        d41 += g_aBuffer_tmp;
        d42 += g_aBuffer_tmp;
        d43 += g_aBuffer_tmp;
        h_aBuffer_tmp = W_dt[(4 * i17) + 1] * d45;
        aBuffer[aBuffer_tmp] += h_aBuffer_tmp;
        b_W_dt[aBuffer_tmp] += h_aBuffer_tmp;
        c_W_dt[aBuffer_tmp] += h_aBuffer_tmp;
        i_aBuffer_tmp = W_dt[(4 * i17) + 2] * d45;
        aBuffer[b_aBuffer_tmp] += i_aBuffer_tmp;
        b_W_dt[b_aBuffer_tmp] += i_aBuffer_tmp;
        c_W_dt[b_aBuffer_tmp] += i_aBuffer_tmp;
        j_aBuffer_tmp = W_dt[(4 * i17) + 3] * d45;
        aBuffer[d_aBuffer_tmp] += j_aBuffer_tmp;
        b_W_dt[d_aBuffer_tmp] += j_aBuffer_tmp;
        c_W_dt[d_aBuffer_tmp] += j_aBuffer_tmp;
      }
      c_W_dt[4 * i15] = d43;
      b_W_dt[4 * i15] = d42;
      aBuffer[4 * i15] = d41;
      d47 = d_W_dt[4 * i15];
      g_W_dt_tmp = (4 * i15) + 1;
      h_W_dt_tmp = (4 * i15) + 2;
      i_W_dt_tmp = (4 * i15) + 3;
      for (i19 = 0; i19 < 4; i19++) {
        double d48;
        d48 = c_W_dt[i19 + (4 * i15)];
        d47 += W_dt[4 * i19] * d48;
        d_W_dt[g_W_dt_tmp] += W_dt[(4 * i19) + 1] * d48;
        d_W_dt[h_W_dt_tmp] += W_dt[(4 * i19) + 2] * d48;
        d_W_dt[i_W_dt_tmp] += W_dt[(4 * i19) + 3] * d48;
      }
      d_W_dt[4 * i15] = d47;
    }
    memset(&b_aBuffer[0], 0, (sizeof(double)) << 4);
    for (i16 = 0; i16 < 4; i16++) {
      double d44;
      int c_aBuffer_tmp;
      int e_aBuffer_tmp;
      int f_aBuffer_tmp;
      d44 = b_aBuffer[4 * i16];
      c_aBuffer_tmp = (4 * i16) + 1;
      e_aBuffer_tmp = (4 * i16) + 2;
      f_aBuffer_tmp = (4 * i16) + 3;
      for (i18 = 0; i18 < 4; i18++) {
        double d46;
        d46 = aBuffer[i18 + (4 * i16)];
        d44 += aBuffer[4 * i18] * d46;
        b_aBuffer[c_aBuffer_tmp] += aBuffer[(4 * i18) + 1] * d46;
        b_aBuffer[e_aBuffer_tmp] += aBuffer[(4 * i18) + 2] * d46;
        b_aBuffer[f_aBuffer_tmp] += aBuffer[(4 * i18) + 3] * d46;
      }
      int F_tmp;
      int b_F_tmp;
      int c_F_tmp;
      b_aBuffer[4 * i16] = d44;
      F[11 * i16] =
          (((dv4[4 * i16] + W_dt[4 * i16]) + (0.5 * b_W_dt[4 * i16])) +
           (0.16666666666666666 * d_W_dt[4 * i16])) +
          (0.041666666666666664 * d44);
      F_tmp = (4 * i16) + 1;
      F[(11 * i16) + 1] =
          (((dv4[F_tmp] + W_dt[F_tmp]) + (0.5 * b_W_dt[F_tmp])) +
           (0.16666666666666666 * d_W_dt[F_tmp])) +
          (0.041666666666666664 * b_aBuffer[F_tmp]);
      b_F_tmp = (4 * i16) + 2;
      F[(11 * i16) + 2] =
          (((dv4[b_F_tmp] + W_dt[b_F_tmp]) + (0.5 * b_W_dt[b_F_tmp])) +
           (0.16666666666666666 * d_W_dt[b_F_tmp])) +
          (0.041666666666666664 * b_aBuffer[b_F_tmp]);
      c_F_tmp = (4 * i16) + 3;
      F[(11 * i16) + 3] =
          (((dv4[c_F_tmp] + W_dt[c_F_tmp]) + (0.5 * b_W_dt[c_F_tmp])) +
           (0.16666666666666666 * d_W_dt[c_F_tmp])) +
          (0.041666666666666664 * b_aBuffer[c_F_tmp]);
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
    for (i20 = 0; i20 < 3; i20++) {
      int d_F_tmp;
      int e_F_tmp;
      d_F_tmp = 4 * (i20 + 1);
      e_F_tmp = 11 * (i20 + 4);
      F[e_F_tmp] = l_a[d_F_tmp];
      F[e_F_tmp + 1] = l_a[d_F_tmp + 1];
      F[e_F_tmp + 2] = l_a[d_F_tmp + 2];
      F[e_F_tmp + 3] = l_a[d_F_tmp + 3];
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
    for (i21 = 0; i21 < 3; i21++) {
      double d49;
      int c_n_tilde_tmp;
      int d_n_tilde_tmp;
      d49 = b_n_tilde[3 * i21];
      c_n_tilde_tmp = (3 * i21) + 1;
      d_n_tilde_tmp = (3 * i21) + 2;
      for (i23 = 0; i23 < 3; i23++) {
        double d50;
        d50 = n_tilde[i23 + (3 * i21)];
        d49 += n_tilde[3 * i23] * d50;
        b_n_tilde[c_n_tilde_tmp] += n_tilde[(3 * i23) + 1] * d50;
        b_n_tilde[d_n_tilde_tmp] += n_tilde[(3 * i23) + 2] * d50;
      }
      b_n_tilde[3 * i21] = d49;
    }
    for (i22 = 0; i22 < 9; i22++) {
      w_exp_tilde[i22] = (w_exp_tilde_tmp[i22] - (b_a * n_tilde[i22])) +
                         ((1.0 - b_x) * b_n_tilde[i22]);
    }
    memset(&b_dv[0], 0, 9U * (sizeof(double)));
    for (i24 = 0; i24 < 3; i24++) {
      double d51;
      int i26;
      int i27;
      d51 = b_dv[3 * i24];
      i26 = (3 * i24) + 1;
      i27 = (3 * i24) + 2;
      for (i29 = 0; i29 < 3; i29++) {
        double d53;
        d53 = w_exp_tilde[i29 + (3 * i24)];
        d51 += SD->pd->d_param.Jinv[3 * i29] * d53;
        b_dv[i26] += SD->pd->d_param.Jinv[(3 * i29) + 1] * d53;
        b_dv[i27] += SD->pd->d_param.Jinv[(3 * i29) + 2] * d53;
        F[(i29 + (11 * (i24 + 4))) + 4] = 0.0;
      }
      b_dv[3 * i24] = d51;
    }
    for (i25 = 0; i25 < 3; i25++) {
      int F_tmp_tmp;
      F_tmp_tmp = 11 * (i25 + 4);
      for (i28 = 0; i28 < 3; i28++) {
        double d52;
        d52 = SD->pd->d_param.J[i28 + (3 * i25)];
        F[F_tmp_tmp + 4] += b_dv[3 * i28] * d52;
        F[F_tmp_tmp + 5] += b_dv[(3 * i28) + 1] * d52;
        F[F_tmp_tmp + 6] += b_dv[(3 * i28) + 2] * d52;
      }
    }
    b_dv[1] = t1_density * (m_a * x[9]);
    b_dv[4] = 0.0;
    b_dv[7] = t1_density * (m_a * x[7]);
    b_dv[2] = t1_density * (m_a * (-x[8]));
    b_dv[5] = t1_density * (m_a * (-x[7]));
    b_dv[8] = 0.0;
    c_x = 0.0;
    for (i30 = 0; i30 < 3; i30++) {
      double d54;
      double d55;
      double d56;
      int f_F_tmp;
      b_dv[3 * i30] = 0.0;
      f_F_tmp = 11 * (i30 + 7);
      d54 = 0.0;
      d55 = 0.0;
      d56 = 0.0;
      for (i31 = 0; i31 < 3; i31++) {
        double d57;
        d57 = b_dv[i31 + (3 * i30)];
        d54 += (dt * SD->pd->d_param.Jinv[3 * i31]) * d57;
        d55 += (dt * SD->pd->d_param.Jinv[(3 * i31) + 1]) * d57;
        d56 += (dt * SD->pd->d_param.Jinv[(3 * i31) + 2]) * d57;
      }
      F[f_F_tmp + 6] = d56;
      F[f_F_tmp + 5] = d55;
      F[f_F_tmp + 4] = d54;
      c_x += x[i30 + 1] * SD->pd->d_param.g[i30];
    }
    d_x[0] = (x[2] * SD->pd->d_param.g[2]) - (SD->pd->d_param.g[1] * x[3]);
    d_x[1] = (SD->pd->d_param.g[0] * x[3]) - (x[1] * SD->pd->d_param.g[2]);
    d_x[2] = (x[1] * SD->pd->d_param.g[1]) - (SD->pd->d_param.g[0] * x[2]);
    dv5[0] = 0.0;
    dv5[3] = x[0] * (-SD->pd->d_param.g[2]);
    dv5[6] = x[0] * SD->pd->d_param.g[1];
    dv5[1] = x[0] * SD->pd->d_param.g[2];
    dv5[4] = 0.0;
    dv5[7] = x[0] * (-SD->pd->d_param.g[0]);
    dv5[2] = x[0] * (-SD->pd->d_param.g[1]);
    dv5[5] = x[0] * SD->pd->d_param.g[0];
    dv5[8] = 0.0;
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
    for (i32 = 0; i32 < 3; i32++) {
      double d58;
      double d59;
      double g_F_tmp;
      int h_F_tmp;
      int i_F_tmp;
      int j_F_tmp;
      F[i32 + 7] = dt * (2.0 * ((x[0] * SD->pd->d_param.g[i32]) - d_x[i32]));
      g_F_tmp = x[i32 + 1];
      h_F_tmp = 11 * (i32 + 1);
      F[h_F_tmp + 7] =
          dt *
          (2.0 * ((((c_x * b_b[3 * i32]) + (x[1] * SD->pd->d_param.g[i32])) -
                   (SD->pd->d_param.g[0] * g_F_tmp)) +
                  dv5[3 * i32]));
      i_F_tmp = (3 * i32) + 1;
      F[h_F_tmp + 8] =
          dt *
          (2.0 * ((((c_x * b_b[i_F_tmp]) + (x[2] * SD->pd->d_param.g[i32])) -
                   (SD->pd->d_param.g[1] * g_F_tmp)) +
                  dv5[i_F_tmp]));
      j_F_tmp = (3 * i32) + 2;
      F[h_F_tmp + 9] =
          dt *
          (2.0 * ((((c_x * b_b[j_F_tmp]) + (x[3] * SD->pd->d_param.g[i32])) -
                   (SD->pd->d_param.g[2] * g_F_tmp)) +
                  dv5[j_F_tmp]));
      d58 = c_skewed_exp_w_tmp[3 * i32];
      d59 = b_dv[3 * i32];
      for (i34 = 0; i34 < 3; i34++) {
        double d60;
        double d61;
        int b_skewed_exp_w_tmp_tmp;
        int i35;
        int skewed_exp_w_tmp_tmp;
        i35 = i34 + (3 * i32);
        d60 = b_skewed_exp_w_tmp[i35];
        d61 = skewed_exp_w_tmp[i35];
        d58 += skewed_exp_w_tmp[3 * i34] * d60;
        d59 += (2.0 * b_skewed_exp_w_tmp[3 * i34]) * d61;
        skewed_exp_w_tmp_tmp = (3 * i34) + 1;
        c_skewed_exp_w_tmp[i_F_tmp] +=
            skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d60;
        b_dv[i_F_tmp] += (2.0 * b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp]) * d61;
        b_skewed_exp_w_tmp_tmp = (3 * i34) + 2;
        c_skewed_exp_w_tmp[j_F_tmp] +=
            skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d60;
        b_dv[j_F_tmp] +=
            (2.0 * b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp]) * d61;
      }
      int k_F_tmp;
      int l_F_tmp;
      b_dv[3 * i32] = d59;
      c_skewed_exp_w_tmp[3 * i32] = d58;
      k_F_tmp = 11 * (i32 + 4);
      F[k_F_tmp + 7] = (dt * skewed_exp_w_tmp[3 * i32]) + (n_a * (d58 - d59));
      l_F_tmp = 11 * (i32 + 7);
      F[l_F_tmp + 7] = w_exp_tilde[3 * i32];
      F[k_F_tmp + 8] = (dt * skewed_exp_w_tmp[i_F_tmp]) +
                       (n_a * (c_skewed_exp_w_tmp[i_F_tmp] - b_dv[i_F_tmp]));
      F[l_F_tmp + 8] = w_exp_tilde[i_F_tmp];
      F[k_F_tmp + 9] = (dt * skewed_exp_w_tmp[j_F_tmp]) +
                       (n_a * (c_skewed_exp_w_tmp[j_F_tmp] - b_dv[j_F_tmp]));
      F[l_F_tmp + 9] = w_exp_tilde[j_F_tmp];
      r_q_tmp[i32 + 1] = -g_F_tmp;
      b_r_q_tmp += (-g_F_tmp) * x[i32 + 7];
    }
    c_r_q_tmp[0] = (r_q_tmp[2] * x[9]) - (r_q_tmp[3] * x[8]);
    c_r_q_tmp[1] = (r_q_tmp[3] * x[7]) - (r_q_tmp[1] * x[9]);
    c_r_q_tmp[2] = (r_q_tmp[1] * x[8]) - (r_q_tmp[2] * x[7]);
    for (i33 = 0; i33 < 3; i33++) {
      double b_dt_tmp;
      double dt_tmp;
      int c_dt_tmp;
      int d_dt_tmp;
      int e_dt_tmp;
      dt_tmp = x[i33 + 7];
      d_dt[i33] = dt * (2.0 * ((r_q_tmp[0] * dt_tmp) - c_r_q_tmp[i33]));
      b_dt_tmp = r_q_tmp[i33 + 1];
      c_dt_tmp = 3 * (i33 + 1);
      d_dt[c_dt_tmp] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[3 * i33]) + (r_q_tmp[1] * dt_tmp)) -
                        (x[7] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[3 * i33])));
      d_dt_tmp = (3 * i33) + 1;
      d_dt[c_dt_tmp + 1] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[d_dt_tmp]) + (r_q_tmp[2] * dt_tmp)) -
                        (x[8] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[d_dt_tmp])));
      e_dt_tmp = (3 * i33) + 2;
      d_dt[c_dt_tmp + 2] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[e_dt_tmp]) + (r_q_tmp[3] * dt_tmp)) -
                        (x[9] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[e_dt_tmp])));
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
    for (i36 = 0; i36 < 3; i36++) {
      F[(11 * (i36 + 7)) + 10] =
          dt *
          (((o_a * b_b[3 * i36]) + ((2.0 * r_q_tmp[1]) * r_q_tmp[i36 + 1])) -
           b_dv[3 * i36]);
    }
    F[120] = 1.0;
    memset(&b_F[0], 0, 121U * (sizeof(double)));
    for (i37 = 0; i37 < 11; i37++) {
      for (i38 = 0; i38 < 11; i38++) {
        double d62;
        d62 = P[i38 + (11 * i37)];
        for (i41 = 0; i41 < 11; i41++) {
          int m_F_tmp;
          m_F_tmp = i41 + (11 * i37);
          b_F[m_F_tmp] += F[i41 + (11 * i38)] * d62;
        }
      }
    }
    for (i39 = 0; i39 < 11; i39++) {
      for (i40 = 0; i40 < 11; i40++) {
        double d63;
        d63 = 0.0;
        for (i43 = 0; i43 < 11; i43++) {
          d63 += b_F[i39 + (11 * i43)] * F[i40 + (11 * i43)];
        }
        int c_P_pred_tmp;
        c_P_pred_tmp = i39 + (11 * i40);
        P_pred[c_P_pred_tmp] = d63 + Q[c_P_pred_tmp];
      }
    }
    for (i42 = 0; i42 < 3; i42++) {
      int P_pred_tmp;
      int b_P_pred_tmp;
      int d_P_pred_tmp;
      P_pred_tmp = 11 * (i42 + 4);
      b_P_pred[3 * i42] = P_pred[P_pred_tmp + 4] + R[3 * i42];
      b_P_pred_tmp = (3 * i42) + 1;
      b_P_pred[b_P_pred_tmp] = P_pred[P_pred_tmp + 5] + R[b_P_pred_tmp];
      d_P_pred_tmp = (3 * i42) + 2;
      b_P_pred[d_P_pred_tmp] = P_pred[P_pred_tmp + 6] + R[d_P_pred_tmp];
    }
    inv(b_P_pred, b_dv);
    memset(&K[0], 0, 33U * (sizeof(double)));
    for (i44 = 0; i44 < 3; i44++) {
      for (i45 = 0; i45 < 3; i45++) {
        double d64;
        d64 = b_dv[i45 + (3 * i44)];
        for (i46 = 0; i46 < 11; i46++) {
          int K_tmp;
          K_tmp = i46 + (11 * i44);
          K[K_tmp] += P_pred[i46 + (11 * (i45 + 4))] * d64;
        }
      }
    }
    d_eye(dv6);
    memcpy(&E[0], &dv6[0], 44U * (sizeof(double)));
    for (i47 = 0; i47 < 33; i47++) {
      E[i47 + 44] = dv6[i47 + 44] - K[i47];
    }
    memcpy(&E[77], &dv6[77], 44U * (sizeof(double)));
    memset(&b_E[0], 0, 121U * (sizeof(double)));
    for (i48 = 0; i48 < 11; i48++) {
      for (i49 = 0; i49 < 11; i49++) {
        double d65;
        d65 = P_pred[i49 + (11 * i48)];
        for (i51 = 0; i51 < 11; i51++) {
          int E_tmp;
          E_tmp = i51 + (11 * i48);
          b_E[E_tmp] += E[i51 + (11 * i49)] * d65;
        }
      }
    }
    memset(&b_K[0], 0, 33U * (sizeof(double)));
    for (i50 = 0; i50 < 3; i50++) {
      for (i52 = 0; i52 < 3; i52++) {
        double d66;
        d66 = R[i52 + (3 * i50)];
        for (i53 = 0; i53 < 11; i53++) {
          int b_K_tmp;
          b_K_tmp = i53 + (11 * i50);
          b_K[b_K_tmp] += K[i53 + (11 * i52)] * d66;
        }
      }
    }
    memset(&c_E[0], 0, 121U * (sizeof(double)));
    memset(&c_K[0], 0, 121U * (sizeof(double)));
    for (i54 = 0; i54 < 11; i54++) {
      for (i55 = 0; i55 < 11; i55++) {
        double d67;
        d67 = E[i54 + (11 * i55)];
        for (i58 = 0; i58 < 11; i58++) {
          int b_E_tmp;
          b_E_tmp = i58 + (11 * i54);
          c_E[b_E_tmp] += b_E[i58 + (11 * i55)] * d67;
        }
      }
      for (i57 = 0; i57 < 3; i57++) {
        double d69;
        d69 = K[i54 + (11 * i57)];
        for (i59 = 0; i59 < 11; i59++) {
          int c_K_tmp;
          c_K_tmp = i59 + (11 * i54);
          c_K[c_K_tmp] += b_K[i59 + (11 * i57)] * d69;
        }
      }
    }
    for (i56 = 0; i56 < 121; i56++) {
      P_ret[i56] = c_E[i56] + c_K[i56];
    }
    double d68;
    double d70;
    d68 = d1 / d3;
    d70 = d2 / d3;
    d71 =
        ((((d1 / d4) * (sens_input->board_gyro.meas[0] - bias->board_gyro[0])) +
          ((d2 / d4) * (sens_input->mti_gyro.meas[0] - bias->mti_gyro[0]))) +
         (C_ad_w_idx_0 * (sens_input->ad_gyro.meas[0] - bias->ad_gyro[0]))) -
        x_pred[4];
    d72 = (((d68 * (sens_input->board_gyro.meas[1] - bias->board_gyro[1])) +
            (d70 * (sens_input->mti_gyro.meas[1] - bias->mti_gyro[1]))) +
           (d * (sens_input->ad_gyro.meas[1] - bias->ad_gyro[1]))) -
          x_pred[5];
    d73 = (((d68 * (sens_input->board_gyro.meas[2] - bias->board_gyro[2])) +
            (d70 * (sens_input->mti_gyro.meas[2] - bias->mti_gyro[2]))) +
           (d * (sens_input->ad_gyro.meas[2] - bias->ad_gyro[2]))) -
          x_pred[6];
    for (i60 = 0; i60 < 11; i60++) {
      x_ret[i60] = x_pred[i60] + (((K[i60] * d71) + (K[i60 + 11] * d72)) +
                                  (K[i60 + 22] * d73));
    }
    double b_q_mag;
    b_q_mag = c_norm(&x_ret[0]);
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
  *cov_norm = d_norm(P);
  airdata->pressure = airdata_atmos(x[10], &airdata->temperature,
                                    &airdata->density, &airdata->sonic_speed,
                                    &airdata->mach, &airdata->dynamic_pressure);
  airspeed = b_norm(&x[7]);
  airdata->mach = airspeed / airdata->sonic_speed;
  airdata->dynamic_pressure = (0.5 * airdata->density) * (airspeed * airspeed);
  roll_state[0] =
      atan2(2.0 * ((x[2] * x[3]) + (x[0] * x[1])),
            (((x[0] * x[0]) - (x[1] * x[1])) - (x[2] * x[2])) + (x[3] * x[3]));
  roll_state[1] = x[4];
}
