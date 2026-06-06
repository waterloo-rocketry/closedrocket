#include "GNC_codegen.h"
#include "GNC_codegen_types.h"
#include <math.h>
#include <string.h>

static const real_T dv[9] = {0.46, 0.0, 0.0, 0.0, 49.5, 0.0, 0.0, 0.0, 49.5};

static const real_T dv1[9] = {2.1739130434782608,   0.0, 0.0, 0.0,
                              0.020202020202020204, 0.0, 0.0, 0.0,
                              0.020202020202020204};

static real_T airdata_atmos(real_T altitude, real_T *airdata_temperature,
                            real_T *airdata_density,
                            real_T *airdata_sonic_speed, real_T *airdata_mach,
                            real_T *airdata_dynamic_pressure);

static void b_ekf_correct(const real_T x[11], const real_T P[121], real_T y,
                          real_T b, real_T x_new[11], real_T P_new[121]);

static real_T b_norm(const real_T x[3]);

static real_T b_xnrm2(int32_T n, const real_T x[11], int32_T ix0);

static void b_xzlascl(real_T cfrom, real_T cto, real_T A[11]);

static real_T c_norm(const real_T x[121]);

static void controller_codegen_entry_init(GNC_codegenStackData *SD);

static void dynamics_init(GNC_codegenStackData *SD);

static void dynamics_jacobian_init(GNC_codegenStackData *SD);

static void ekf_correct(const real_T x[11], const real_T P[121],
                        const real_T y[3], const real_T b[3], const real_T R[9],
                        real_T x_new[11], real_T P_new[121]);

static void inv(const real_T x[9], real_T y[9]);

static void pad_filter_init(GNC_codegenStackData *SD);

static void xaxpy(int32_T n, real_T a, const real_T x[11], int32_T ix0,
                  real_T y[121], int32_T iy0);

static real_T xnrm2(int32_T n, const real_T x[121], int32_T ix0);

static real_T xrotg(real_T *a, real_T *b, real_T *s);

static void xzlascl(real_T cfrom, real_T cto, real_T A[121]);

static real_T airdata_atmos(real_T altitude, real_T *airdata_temperature,
                            real_T *airdata_density,
                            real_T *airdata_sonic_speed, real_T *airdata_mach,
                            real_T *airdata_dynamic_pressure) {
  real_T airdata_pressure;
  real_T layer_idx_1;
  real_T layer_idx_2;
  real_T layer_idx_3;
  real_T pressure;
  real_T temperature;
  int32_T layer_idx_0;
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
                                          (altitude - ((real_T)layer_idx_0))),
                                   9.81 / (287.0579 * layer_idx_3));
    }
  } else {
    pressure = layer_idx_1 * pow(1.0 - ((layer_idx_3 / layer_idx_2) *
                                        (altitude - ((real_T)layer_idx_0))),
                                 9.81 / (287.0579 * layer_idx_3));
  }
  temperature =
      layer_idx_2 - (layer_idx_3 * (altitude - ((real_T)layer_idx_0)));
  *airdata_density = pressure / (287.0579 * temperature);
  *airdata_sonic_speed = sqrt(401.88106 * temperature);
  airdata_pressure = pressure;
  *airdata_temperature = temperature;
  *airdata_mach = 0.0;
  *airdata_dynamic_pressure = 0.0;
  return airdata_pressure;
}

static void b_ekf_correct(const real_T x[11], const real_T P[121], real_T y,
                          real_T b, real_T x_new[11], real_T P_new[121]) {
  real_T E[121];
  real_T b_E[121];
  real_T b_K[121];
  real_T c_E[121];
  real_T H[11];
  real_T K[11];
  real_T b_H[11];
  real_T b_P[11];
  real_T absxk;
  real_T airdata_altitude_pressure;
  real_T altitude;
  real_T altitude_ratio;
  real_T b_b;
  real_T b_expl_temp;
  real_T c_H;
  real_T c_b;
  real_T c_expl_temp;
  real_T d_expl_temp;
  real_T e_expl_temp;
  real_T expl_temp;
  real_T layer_idx_1;
  real_T layer_idx_2;
  real_T layer_idx_3;
  real_T q_mag;
  real_T scale;
  real_T t;
  real_T t0_pressure;
  int32_T b_i;
  int32_T i;
  int32_T i1;
  int32_T i10;
  int32_T i11;
  int32_T i12;
  int32_T i13;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T i5;
  int32_T i6;
  int32_T i7;
  int32_T i8;
  int32_T i9;
  int32_T k;
  int32_T layer_idx_0;
  int8_T b_I[121];
  t0_pressure = airdata_atmos(x[10], &expl_temp, &b_expl_temp, &c_expl_temp,
                              &d_expl_temp, &e_expl_temp);
  b_b = y - (t0_pressure + b);
  memset(&H[0], 0, 11U * (sizeof(real_T)));
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
                     (altitude - ((real_T)layer_idx_0))),
              (9.81 / (287.0579 * layer_idx_3)) - 1.0);
    }
  } else {
    airdata_altitude_pressure =
        ((((-layer_idx_1) * 9.81) / (layer_idx_2 * 287.0579)) *
         (altitude_ratio * altitude_ratio)) *
        pow(1.0 - ((layer_idx_3 / layer_idx_2) *
                   (altitude - ((real_T)layer_idx_0))),
            (9.81 / (287.0579 * layer_idx_3)) - 1.0);
  }
  H[10] = airdata_altitude_pressure;
  memset(&b_H[0], 0, 11U * (sizeof(real_T)));
  c_H = 0.0;
  for (i = 0; i < 11; i++) {
    real_T d;
    d = b_H[i];
    for (i1 = 0; i1 < 11; i1++) {
      d += H[i1] * P[i1 + (11 * i)];
    }
    b_H[i] = d;
    c_H += d * H[i];
  }
  c_b = 1.0 / (c_H + 100.0);
  memset(&b_P[0], 0, 11U * (sizeof(real_T)));
  for (i2 = 0; i2 < 11; i2++) {
    for (i3 = 0; i3 < 11; i3++) {
      b_P[i3] += P[i3 + (11 * i2)] * H[i2];
    }
  }
  for (i4 = 0; i4 < 11; i4++) {
    K[i4] = b_P[i4] * c_b;
  }
  memset(&b_I[0], 0, 121U * (sizeof(int8_T)));
  for (k = 0; k < 11; k++) {
    b_I[k + (11 * k)] = (int8_T)1;
  }
  for (i5 = 0; i5 < 11; i5++) {
    for (i6 = 0; i6 < 11; i6++) {
      int32_T E_tmp;
      E_tmp = i6 + (11 * i5);
      E[E_tmp] = ((real_T)b_I[E_tmp]) - (K[i6] * H[i5]);
    }
  }
  memset(&b_E[0], 0, 121U * (sizeof(real_T)));
  for (i7 = 0; i7 < 11; i7++) {
    for (i8 = 0; i8 < 11; i8++) {
      real_T d1;
      d1 = P[i8 + (11 * i7)];
      for (i10 = 0; i10 < 11; i10++) {
        int32_T b_E_tmp;
        b_E_tmp = i10 + (11 * i7);
        b_E[b_E_tmp] += E[i10 + (11 * i8)] * d1;
      }
    }
  }
  memset(&c_E[0], 0, 121U * (sizeof(real_T)));
  for (i9 = 0; i9 < 11; i9++) {
    for (i11 = 0; i11 < 11; i11++) {
      real_T d2;
      d2 = E[i9 + (11 * i11)];
      for (i13 = 0; i13 < 11; i13++) {
        int32_T c_E_tmp;
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

static real_T b_norm(const real_T x[3]) {
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
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

static real_T b_xnrm2(int32_T n, const real_T x[11], int32_T ix0) {
  real_T y;
  int32_T k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (k = ix0; k < kend; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = ((y * t) * t) + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

static void b_xzlascl(real_T cfrom, real_T cto, real_T A[11]) {
  real_T cfromc;
  real_T ctoc;
  int32_T i;
  boolean_T notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    real_T cfrom1;
    real_T cto1;
    real_T mul;
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

static real_T c_norm(const real_T x[121]) {
  real_T A[121];
  real_T S[11];
  real_T e[11];
  real_T s[11];
  real_T work[11];
  real_T a__1;
  real_T a__2;
  real_T a__3;
  real_T anrm;
  real_T b_f;
  real_T b_sn;
  real_T c_sn;
  real_T cscale;
  real_T d_sn;
  real_T f;
  real_T rt;
  real_T sn;
  real_T snorm;
  int32_T b_jj;
  int32_T b_k;
  int32_T b_q;
  int32_T c_ii;
  int32_T c_jj;
  int32_T c_k;
  int32_T d_k;
  int32_T e_k;
  int32_T f_k;
  int32_T g_k;
  int32_T h_k;
  int32_T i_k;
  int32_T iter;
  int32_T j_k;
  int32_T jj;
  int32_T k;
  int32_T k_k;
  int32_T m;
  int32_T q;
  int32_T qp1;
  boolean_T doscale;
  memcpy(&A[0], &x[0], 121U * (sizeof(real_T)));
  memset(&s[0], 0, 11U * (sizeof(real_T)));
  memset(&e[0], 0, 11U * (sizeof(real_T)));
  memset(&work[0], 0, 11U * (sizeof(real_T)));
  doscale = false;
  anrm = 0.0;
  for (k = 0; k < 121; k++) {
    real_T absxk;
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
    real_T nrm;
    int32_T qq;
    int32_T qq_tmp;
    boolean_T apply_transform;
    qp1 = q + 2;
    qq_tmp = q + (11 * q);
    qq = qq_tmp + 1;
    apply_transform = false;
    nrm = xnrm2(11 - q, A, qq_tmp + 1);
    if (nrm > 0.0) {
      real_T d1;
      apply_transform = true;
      if (A[qq_tmp] < 0.0) {
        d1 = -nrm;
      } else {
        d1 = nrm;
      }
      s[q] = d1;
      if (fabs(d1) >= 1.0020841800044864E-292) {
        real_T a;
        int32_T i1;
        a = 1.0 / d1;
        i1 = (qq_tmp - q) + 11;
        for (d_k = qq; d_k <= i1; d_k++) {
          A[d_k - 1] *= a;
        }
      } else {
        int32_T i;
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
      int32_T qjj;
      qjj = q + (11 * (jj - 1));
      if (apply_transform) {
        real_T b_d;
        real_T c_a;
        int32_T n;
        n = 10 - q;
        b_d = 0.0;
        for (c_k = 0; c_k <= n; c_k++) {
          b_d += A[qq_tmp + c_k] * A[qjj + c_k];
        }
        c_a = -(b_d / A[qq_tmp]);
        if (c_a != 0.0) {
          int32_T i2;
          i2 = 11 - q;
          for (g_k = 0; g_k < i2; g_k++) {
            int32_T A_tmp;
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
        real_T b_a;
        if (e[q + 1] < 0.0) {
          e[q] = -nrm;
        } else {
          e[q] = nrm;
        }
        b_a = e[q];
        if (fabs(e[q]) >= 1.0020841800044864E-292) {
          real_T d_a;
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
          real_T d4;
          d4 = e[b_jj - 1];
          if (d4 != 0.0) {
            int32_T i4;
            int32_T ix;
            ix = q + (11 * (b_jj - 1));
            i4 = 10 - q;
            for (j_k = 0; j_k < i4; j_k++) {
              int32_T work_tmp;
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
    real_T d;
    real_T d2;
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
      real_T d3;
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
    int32_T c_q;
    int32_T ii;
    int32_T kase;
    ii = m + 1;
    int32_T exitg1;
    do {
      exitg1 = 0;
      c_q = ii;
      if (ii == 0) {
        exitg1 = 1;
      } else {
        real_T ztest0;
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
      int32_T b_ii;
      int32_T qs;
      boolean_T exitg2;
      qs = m + 2;
      b_ii = m + 2;
      exitg2 = false;
      while ((!((int32_T)exitg2)) && (b_ii >= ii)) {
        qs = b_ii;
        if (b_ii == ii) {
          exitg2 = true;
        } else {
          real_T test;
          real_T ztest;
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
      int32_T i3;
      f = e[m];
      e[m] = 0.0;
      i3 = m + 1;
      for (i_k = i3; i_k >= (c_q + 1); i_k--) {
        real_T cs;
        cs = xrotg(&s[i_k - 1], &f, &sn);
        if (i_k > (c_q + 1)) {
          real_T b_f_tmp;
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
        real_T b_cs;
        real_T f_tmp;
        a__1 = f;
        b_cs = xrotg(&s[h_k - 1], &a__1, &b_sn);
        f_tmp = e[h_k - 1];
        f = (-b_sn) * f_tmp;
        e[h_k - 1] = f_tmp * b_cs;
      }
    } break;
    case 3: {
      real_T b;
      real_T c;
      real_T emm1;
      real_T g;
      real_T scale;
      real_T scale_tmp;
      real_T shift;
      real_T sm;
      real_T smm1;
      real_T sqds;
      int32_T mm1;
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
        real_T c_cs;
        real_T c_f_tmp;
        real_T d_cs;
        real_T d_f_tmp;
        real_T e_f_tmp;
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
  memcpy(&S[0], &s[0], 11U * (sizeof(real_T)));
  if (doscale) {
    b_xzlascl(cscale, anrm, S);
  }
  return S[0];
}

static void controller_codegen_entry_init(GNC_codegenStackData *SD) {
  int32_T i;
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
  int32_T i;
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
  int32_T i;
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

static void ekf_correct(const real_T x[11], const real_T P[121],
                        const real_T y[3], const real_T b[3], const real_T R[9],
                        real_T x_new[11], real_T P_new[121]) {
  static const int8_T iv[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  real_T E[121];
  real_T b_E[121];
  real_T c_E[121];
  real_T c_K[121];
  real_T H[33];
  real_T K[33];
  real_T b_H[33];
  real_T b_K[33];
  real_T b_P[33];
  real_T y_tmp[33];
  real_T b_dv[9];
  real_T c_H[9];
  real_T c_a[9];
  real_T b_y[3];
  real_T c_x[3];
  real_T a;
  real_T absxk;
  real_T b_a;
  real_T b_x;
  real_T d12;
  real_T d13;
  real_T d14;
  real_T d15;
  real_T d16;
  real_T d17;
  real_T q_mag;
  real_T scale;
  real_T t;
  int32_T i;
  int32_T i1;
  int32_T i10;
  int32_T i11;
  int32_T i12;
  int32_T i13;
  int32_T i14;
  int32_T i15;
  int32_T i16;
  int32_T i17;
  int32_T i18;
  int32_T i19;
  int32_T i2;
  int32_T i20;
  int32_T i21;
  int32_T i22;
  int32_T i23;
  int32_T i24;
  int32_T i25;
  int32_T i26;
  int32_T i27;
  int32_T i28;
  int32_T i29;
  int32_T i3;
  int32_T i30;
  int32_T i4;
  int32_T i5;
  int32_T i6;
  int32_T i7;
  int32_T i8;
  int32_T i9;
  int32_T k;
  int8_T b_I[121];
  a = (x[0] * x[0]) - (((x[1] * x[1]) + (x[2] * x[2])) + (x[3] * x[3]));
  b_a = 2.0 * x[0];
  memset(&H[0], 0, 33U * (sizeof(real_T)));
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
    real_T H_tmp;
    int32_T b_H_tmp;
    int32_T c_H_tmp;
    int32_T d_H_tmp;
    H[i] = 2.0 * ((x[0] * b[i]) - c_x[i]);
    H_tmp = x[i + 1];
    b_H_tmp = 3 * (i + 1);
    H[b_H_tmp] =
        2.0 *
        ((((b_x * ((real_T)iv[3 * i])) + (x[1] * b[i])) - (b[0] * H_tmp)) +
         b_dv[3 * i]);
    c_H_tmp = (3 * i) + 1;
    H[b_H_tmp + 1] =
        2.0 *
        ((((b_x * ((real_T)iv[c_H_tmp])) + (x[2] * b[i])) - (b[1] * H_tmp)) +
         b_dv[c_H_tmp]);
    d_H_tmp = (3 * i) + 2;
    H[b_H_tmp + 2] =
        2.0 *
        ((((b_x * ((real_T)iv[d_H_tmp])) + (x[3] * b[i])) - (b[2] * H_tmp)) +
         b_dv[d_H_tmp]);
  }
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 11; i2++) {
      y_tmp[i2 + (11 * i1)] = H[i1 + (3 * i2)];
    }
  }
  memset(&b_H[0], 0, 33U * (sizeof(real_T)));
  for (i3 = 0; i3 < 11; i3++) {
    real_T d;
    int32_T e_H_tmp;
    int32_T f_H_tmp;
    d = b_H[3 * i3];
    e_H_tmp = (3 * i3) + 1;
    f_H_tmp = (3 * i3) + 2;
    for (i5 = 0; i5 < 11; i5++) {
      real_T d1;
      d1 = P[i5 + (11 * i3)];
      d += H[3 * i5] * d1;
      b_H[e_H_tmp] += H[(3 * i5) + 1] * d1;
      b_H[f_H_tmp] += H[(3 * i5) + 2] * d1;
    }
    b_H[3 * i3] = d;
  }
  memset(&b_P[0], 0, 33U * (sizeof(real_T)));
  for (i4 = 0; i4 < 3; i4++) {
    for (i6 = 0; i6 < 3; i6++) {
      real_T d2;
      d2 = 0.0;
      for (i8 = 0; i8 < 11; i8++) {
        d2 += b_H[i4 + (3 * i8)] * y_tmp[i8 + (11 * i6)];
      }
      int32_T g_H_tmp;
      g_H_tmp = i4 + (3 * i6);
      c_H[g_H_tmp] = d2 + R[g_H_tmp];
    }
    for (i7 = 0; i7 < 11; i7++) {
      real_T d3;
      d3 = y_tmp[i7 + (11 * i4)];
      for (i10 = 0; i10 < 11; i10++) {
        int32_T P_tmp;
        P_tmp = i10 + (11 * i4);
        b_P[P_tmp] += P[i10 + (11 * i7)] * d3;
      }
    }
  }
  inv(c_H, b_dv);
  memset(&K[0], 0, 33U * (sizeof(real_T)));
  for (i9 = 0; i9 < 3; i9++) {
    for (i11 = 0; i11 < 3; i11++) {
      real_T d4;
      d4 = b_dv[i11 + (3 * i9)];
      for (i12 = 0; i12 < 11; i12++) {
        int32_T K_tmp;
        K_tmp = i12 + (11 * i9);
        K[K_tmp] += b_P[i12 + (11 * i11)] * d4;
      }
    }
  }
  memset(&b_I[0], 0, 121U * (sizeof(int8_T)));
  for (k = 0; k < 11; k++) {
    b_I[k + (11 * k)] = (int8_T)1;
  }
  for (i13 = 0; i13 < 11; i13++) {
    real_T d5;
    real_T d6;
    real_T d7;
    d5 = K[i13];
    d6 = K[i13 + 11];
    d7 = K[i13 + 22];
    for (i15 = 0; i15 < 11; i15++) {
      int32_T E_tmp;
      E_tmp = i13 + (11 * i15);
      E[E_tmp] = ((real_T)b_I[E_tmp]) -
                 (((d5 * H[3 * i15]) + (d6 * H[(3 * i15) + 1])) +
                  (d7 * H[(3 * i15) + 2]));
    }
  }
  memset(&b_E[0], 0, 121U * (sizeof(real_T)));
  for (i14 = 0; i14 < 11; i14++) {
    for (i16 = 0; i16 < 11; i16++) {
      real_T d8;
      d8 = P[i16 + (11 * i14)];
      for (i18 = 0; i18 < 11; i18++) {
        int32_T b_E_tmp;
        b_E_tmp = i18 + (11 * i14);
        b_E[b_E_tmp] += E[i18 + (11 * i16)] * d8;
      }
    }
  }
  memset(&b_K[0], 0, 33U * (sizeof(real_T)));
  for (i17 = 0; i17 < 3; i17++) {
    for (i19 = 0; i19 < 3; i19++) {
      real_T d9;
      d9 = R[i19 + (3 * i17)];
      for (i20 = 0; i20 < 11; i20++) {
        int32_T b_K_tmp;
        b_K_tmp = i20 + (11 * i17);
        b_K[b_K_tmp] += K[i20 + (11 * i19)] * d9;
      }
    }
  }
  memset(&c_E[0], 0, 121U * (sizeof(real_T)));
  memset(&c_K[0], 0, 121U * (sizeof(real_T)));
  for (i21 = 0; i21 < 11; i21++) {
    for (i22 = 0; i22 < 11; i22++) {
      real_T d10;
      d10 = E[i21 + (11 * i22)];
      for (i25 = 0; i25 < 11; i25++) {
        int32_T c_E_tmp;
        c_E_tmp = i25 + (11 * i21);
        c_E[c_E_tmp] += b_E[i25 + (11 * i22)] * d10;
      }
    }
    for (i24 = 0; i24 < 3; i24++) {
      real_T d11;
      d11 = K[i21 + (11 * i24)];
      for (i27 = 0; i27 < 11; i27++) {
        int32_T c_K_tmp;
        c_K_tmp = i27 + (11 * i21);
        c_K[c_K_tmp] += b_K[i27 + (11 * i24)] * d11;
      }
    }
  }
  for (i23 = 0; i23 < 121; i23++) {
    P_new[i23] = c_E[i23] + c_K[i23];
  }
  for (i26 = 0; i26 < 3; i26++) {
    real_T a_tmp;
    int32_T b_a_tmp;
    int32_T c_a_tmp;
    a_tmp = x[i26 + 1];
    c_a[3 * i26] = (a * ((real_T)iv[3 * i26])) + ((2.0 * x[1]) * a_tmp);
    b_a_tmp = (3 * i26) + 1;
    c_a[b_a_tmp] = (a * ((real_T)iv[b_a_tmp])) + ((2.0 * x[2]) * a_tmp);
    c_a_tmp = (3 * i26) + 2;
    c_a[c_a_tmp] = (a * ((real_T)iv[c_a_tmp])) + ((2.0 * x[3]) * a_tmp);
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

static void inv(const real_T x[9], real_T y[9]) {
  real_T b_x[9];
  real_T absx11;
  real_T absx21;
  real_T absx31;
  real_T t2;
  real_T t3;
  int32_T p1;
  int32_T p2;
  int32_T p3;
  memcpy(&b_x[0], &x[0], 9U * (sizeof(real_T)));
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
    real_T t1;
    int32_T itmp;
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
  int32_T i;
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

static void xaxpy(int32_T n, real_T a, const real_T x[11], int32_T ix0,
                  real_T y[121], int32_T iy0) {
  int32_T k;
  if ((n >= 1) && (a != 0.0)) {
    for (k = 0; k < n; k++) {
      int32_T i;
      i = (iy0 + k) - 1;
      y[i] += a * x[(ix0 + k) - 1];
    }
  }
}

static real_T xnrm2(int32_T n, const real_T x[121], int32_T ix0) {
  real_T y;
  int32_T k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (k = ix0; k < kend; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = ((y * t) * t) + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * sqrt(y);
    }
  }
  return y;
}

static real_T xrotg(real_T *a, real_T *b, real_T *s) {
  real_T absa;
  real_T absb;
  real_T b_c;
  real_T b_s;
  real_T c;
  real_T roe;
  real_T scale;
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
    real_T ads;
    real_T bds;
    real_T r;
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

static void xzlascl(real_T cfrom, real_T cto, real_T A[121]) {
  real_T cfromc;
  real_T ctoc;
  int32_T i;
  boolean_T notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    real_T cfrom1;
    real_T cto1;
    real_T mul;
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

void controller_codegen_entry(GNC_codegenStackData *SD, real_T b_time,
                              real_T dt_ctrl, const real_T xR[2], real_T pdyn,
                              real_T delta, const struct0_T *ctrl_mem_in,
                              real_T *u, real_T *r, struct0_T *ctrl_mem_out) {
  real_T P[4];
  real_T b_K[4];
  real_T b_dv[4];
  real_T b_dv1[4];
  real_T dv2[4];
  real_T K[2];
  real_T b_r[2];
  real_T L_delta;
  real_T a;
  real_T b;
  real_T blend;
  real_T c_delta;
  real_T c_r;
  real_T d;
  real_T d1;
  real_T d10;
  real_T d11;
  real_T d12;
  real_T d2;
  real_T d3;
  real_T d4;
  real_T d5;
  real_T d8;
  real_T pdyn_params;
  real_T r_idx_0;
  real_T w_dot;
  real_T x;
  real_T x_tmp;
  int32_T i;
  int32_T i2;
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
  memset(&b_r[0], 0, (sizeof(real_T)) << 1);
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
  memset(&b_dv1[0], 0, (sizeof(real_T)) << 2);
  d2 = b_dv[0];
  d3 = b_dv[1];
  d4 = b_dv[2];
  d5 = b_dv[3];
  for (i = 0; i < 2; i++) {
    real_T d6;
    real_T d7;
    real_T d9;
    int32_T i1;
    d6 = P[2 * i];
    d7 = b_dv1[2 * i] + (d2 * d6);
    i1 = (2 * i) + 1;
    d9 = b_dv1[i1] + (d3 * d6);
    d6 = P[i1];
    d7 += d4 * d6;
    b_dv1[2 * i] = d7;
    d9 += d5 * d6;
    b_dv1[i1] = d9;
  }
  memset(&dv2[0], 0, (sizeof(real_T)) << 2);
  d8 = b_dv1[0];
  d10 = b_dv1[1];
  d11 = b_dv1[2];
  d12 = b_dv1[3];
  for (i2 = 0; i2 < 2; i2++) {
    real_T d13;
    real_T d14;
    real_T d15;
    int32_T i3;
    d13 = b_dv[i2];
    d14 = dv2[2 * i2] + (d8 * d13);
    i3 = (2 * i2) + 1;
    d15 = dv2[i3] + (d10 * d13);
    b_K[2 * i2] = K[0] * K[i2];
    d13 = b_dv[i2 + 2];
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

void navigation_codegen_entry(GNC_codegenStackData *SD, real_T dt,
                              boolean_T flight_phase, const real_T x[11],
                              const real_T P[121], const struct1_T *bias,
                              const struct2_T *sens_filt,
                              const struct3_T *sens_input, real_T x_ret[11],
                              real_T P_ret[121], struct1_T *bias_ret,
                              struct2_T *sens_filt_ret, real_T *cov_norm,
                              struct6_T *airdata, real_T roll_state[2]) {
  static const real_T Q[121] = {
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
  static const real_T R[9] = {1.0E-9, 0.0, 0.0, 0.0,   1.0E-9,
                              0.0,    0.0, 0.0, 1.0E-9};
  static const real_T b_b[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  real_T F[121];
  real_T b_E[121];
  real_T b_F[121];
  real_T b_P_ret[121];
  real_T c_E[121];
  real_T c_K[121];
  real_T c_P_ret[121];
  real_T d_P_ret[121];
  real_T e_P_ret[121];
  real_T f_P_ret[121];
  real_T g_P_ret[121];
  real_T h_P_ret[121];
  real_T i_P_ret[121];
  real_T K[33];
  real_T b_K[33];
  real_T b_W_dt[16];
  real_T c_b[16];
  real_T d_b[16];
  real_T b_x_ret[11];
  real_T c_x_ret[11];
  real_T d_x_ret[11];
  real_T e_x_ret[11];
  real_T f_x_ret[11];
  real_T g_x_ret[11];
  real_T h_x_ret[11];
  real_T i_x_ret[11];
  real_T b_dv[9];
  real_T b_n_tilde[9];
  real_T b_w_exp_tilde[9];
  real_T c_skewed_exp_w_tmp[9];
  real_T c_q[4];
  real_T q[4];
  real_T b_dt[3];
  real_T c_dt[3];
  real_T c_w_exp_tilde[3];
  real_T dv3[3];
  real_T airspeed;
  real_T b_expl_temp;
  real_T c_expl_temp;
  real_T d_expl_temp;
  real_T e_expl_temp;
  real_T expl_temp;
  real_T f_expl_temp;
  real_T g_expl_temp;
  real_T h_expl_temp;
  real_T i_expl_temp;
  real_T j_expl_temp;
  real_T k_expl_temp;
  real_T l_expl_temp;
  real_T t1_density;
  int32_T b_k;
  int32_T c_k;
  int32_T i;
  int32_T i1;
  int32_T i10;
  int32_T i11;
  int32_T i12;
  int32_T i15;
  int32_T i16;
  int32_T i17;
  int32_T i18;
  int32_T i19;
  int32_T i2;
  int32_T i20;
  int32_T i21;
  int32_T i22;
  int32_T i23;
  int32_T i24;
  int32_T i27;
  int32_T i28;
  int32_T i29;
  int32_T i3;
  int32_T i30;
  int32_T i31;
  int32_T i32;
  int32_T i33;
  int32_T i35;
  int32_T i36;
  int32_T i37;
  int32_T i38;
  int32_T i39;
  int32_T i4;
  int32_T i40;
  int32_T i41;
  int32_T i42;
  int32_T i43;
  int32_T i44;
  int32_T i45;
  int32_T i46;
  int32_T i47;
  int32_T i48;
  int32_T i49;
  int32_T i5;
  int32_T i50;
  int32_T i51;
  int32_T i52;
  int32_T i53;
  int32_T i54;
  int32_T i55;
  int32_T i56;
  int32_T i57;
  int32_T i58;
  int32_T i59;
  int32_T i6;
  int32_T i60;
  int32_T i61;
  int32_T i7;
  int32_T i8;
  int32_T i9;
  int32_T k;
  int8_T c_I[121];
  memcpy(&P_ret[0], &P[0], 121U * (sizeof(real_T)));
  *bias_ret = *bias;
  *sens_filt_ret = *sens_filt;
  if (!((int32_T)flight_phase)) {
    real_T ST[9];
    real_T h_a[9];
    real_T a[3];
    real_T b_absxk;
    real_T b_scale;
    real_T b_t;
    real_T board_baro_f;
    real_T d11;
    real_T d12;
    real_T d5;
    real_T d6;
    real_T d9;
    real_T d_a;
    real_T f_a;
    real_T mti_baro_f;
    real_T qw;
    real_T qy;
    real_T qz;
    real_T y;
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
    f_a = 2.0 * q[0];
    for (i3 = 0; i3 < 3; i3++) {
      real_T a_tmp;
      a_tmp = 2.0 * q[i3 + 1];
      h_a[3 * i3] = (d_a * b_b[i3]) + (a_tmp * q[1]);
      h_a[(3 * i3) + 1] = (d_a * b_b[i3 + 3]) + (a_tmp * q[2]);
      h_a[(3 * i3) + 2] = (d_a * b_b[i3 + 6]) + (a_tmp * d6);
    }
    b_dv[0] = 0.0;
    b_dv[1] = f_a * (-d6);
    b_dv[2] = f_a * q[2];
    b_dv[3] = f_a * d6;
    b_dv[4] = 0.0;
    b_dv[5] = f_a * (-q[1]);
    b_dv[6] = f_a * (-q[2]);
    b_dv[7] = f_a * q[1];
    b_dv[8] = 0.0;
    for (i5 = 0; i5 < 9; i5++) {
      ST[i5] = h_a[i5] - b_dv[i5];
    }
    bias_ret->board_mag_earth[0] = 0.0;
    bias_ret->board_mag_earth[1] = 0.0;
    bias_ret->board_mag_earth[2] = 0.0;
    for (i7 = 0; i7 < 3; i7++) {
      real_T d10;
      d10 = sens_filt_ret->board_mag_f[i7];
      bias_ret->board_mag_earth[0] += ST[3 * i7] * d10;
      bias_ret->board_mag_earth[1] += ST[(3 * i7) + 1] * d10;
      bias_ret->board_mag_earth[2] += ST[(3 * i7) + 2] * d10;
      bias_ret->mti_mag_earth[i7] = 0.0;
    }
    d9 = bias_ret->mti_mag_earth[0];
    d11 = bias_ret->mti_mag_earth[1];
    d12 = bias_ret->mti_mag_earth[2];
    for (i8 = 0; i8 < 3; i8++) {
      real_T d13;
      d13 = sens_filt_ret->mti_mag_f[i8];
      d9 += ST[3 * i8] * d13;
      d11 += ST[(3 * i8) + 1] * d13;
      d12 += ST[(3 * i8) + 2] * d13;
    }
    real_T t1_pressure;
    bias_ret->mti_mag_earth[2] = d12;
    bias_ret->mti_mag_earth[1] = d11;
    bias_ret->mti_mag_earth[0] = d9;
    t1_pressure =
        airdata_atmos(SD->pd->b_param.elevation, &e_expl_temp, &t1_density,
                      &f_expl_temp, &g_expl_temp, &h_expl_temp);
    bias_ret->board_baro = board_baro_f - t1_pressure;
    bias_ret->mti_baro = mti_baro_f - t1_pressure;
  } else {
    real_T E[121];
    real_T P_pred[121];
    real_T W_dt[16];
    real_T b_q[16];
    real_T l_a[16];
    real_T d_dt[12];
    real_T x_pred[11];
    real_T S[9];
    real_T b_P_pred[9];
    real_T b_skewed_exp_w_tmp[9];
    real_T dv4[9];
    real_T h_a[9];
    real_T n_tilde[9];
    real_T skewed_exp_w_tmp[9];
    real_T w_exp_tilde[9];
    real_T b_dv1[4];
    real_T r_q_tmp[4];
    real_T C_total_a[3];
    real_T b_S[3];
    real_T c_r_q_tmp[3];
    real_T d_x[3];
    real_T dn[3];
    real_T dv2[3];
    real_T C_ad_w_idx_0;
    real_T C_total_a_tmp;
    real_T C_total_a_tmp_tmp;
    real_T absxk;
    real_T b;
    real_T b_C_total_a_tmp_tmp;
    real_T b_a;
    real_T b_dphi_tmp;
    real_T b_q_mag;
    real_T b_r_q_tmp;
    real_T b_x;
    real_T c_C_total_a_tmp_tmp;
    real_T c_absxk;
    real_T c_scale;
    real_T c_t;
    real_T c_x;
    real_T d;
    real_T d1;
    real_T d16;
    real_T d17;
    real_T d18;
    real_T d19;
    real_T d2;
    real_T d20;
    real_T d21;
    real_T d22;
    real_T d23;
    real_T d24;
    real_T d25;
    real_T d26;
    real_T d27;
    real_T d3;
    real_T d30;
    real_T d31;
    real_T d33;
    real_T d34;
    real_T d4;
    real_T d67;
    real_T d68;
    real_T d69;
    real_T dphi;
    real_T dphi_tmp;
    real_T e_a;
    real_T g_a;
    real_T i_a;
    real_T j_a;
    real_T k_a;
    real_T m_a;
    real_T n_a;
    real_T n_idx_0;
    real_T n_idx_1;
    real_T n_idx_2;
    real_T o_a;
    real_T q_mag;
    real_T scale;
    real_T t;
    int8_T b_I[16];
    int8_T w_exp_tilde_tmp[9];
    d = 9.9999999999999981E+9 * ((real_T)sens_input->ad_gyro.status);
    C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((real_T)sens_input->board_accel.status);
    b_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((real_T)sens_input->mti_accel.status);
    c_C_total_a_tmp_tmp =
        1.0000000000000002E+14 * ((real_T)sens_input->ad_accel.status);
    C_total_a_tmp =
        (C_total_a_tmp_tmp + b_C_total_a_tmp_tmp) + c_C_total_a_tmp_tmp;
    C_total_a[0] = C_total_a_tmp;
    d1 = 9.9999999999999981E+9 * ((real_T)sens_input->board_gyro.status);
    d2 = 9.9999999999999981E+9 * ((real_T)sens_input->mti_gyro.status);
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
      real_T b_dn_tmp;
      real_T c_dn_tmp;
      real_T dn_tmp;
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
      w_exp_tilde_tmp[i] = (int8_T)0;
    }
    memset(&b_n_tilde[0], 0, 9U * (sizeof(real_T)));
    for (k = 0; k < 3; k++) {
      real_T d7;
      int32_T b_n_tilde_tmp;
      int32_T n_tilde_tmp;
      w_exp_tilde_tmp[k + (3 * k)] = (int8_T)1;
      d7 = b_n_tilde[3 * k];
      n_tilde_tmp = (3 * k) + 1;
      b_n_tilde_tmp = (3 * k) + 2;
      for (i2 = 0; i2 < 3; i2++) {
        real_T d8;
        d8 = n_tilde[i2 + (3 * k)];
        d7 += n_tilde[3 * i2] * d8;
        b_n_tilde[n_tilde_tmp] += n_tilde[(3 * i2) + 1] * d8;
        b_n_tilde[b_n_tilde_tmp] += n_tilde[(3 * i2) + 2] * d8;
      }
      b_n_tilde[3 * k] = d7;
    }
    for (i1 = 0; i1 < 9; i1++) {
      w_exp_tilde[i1] = (((real_T)w_exp_tilde_tmp[i1]) - (b_a * n_tilde[i1])) +
                        ((1.0 - b_x) * b_n_tilde[i1]);
    }
    real_T c_a;
    c_a = b_norm(&x[7]);
    airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp,
                  &d_expl_temp);
    e_a = (0.5 * t1_density) * (c_a * c_a);
    g_a = SD->pd->c_param.c_aero * SD->pd->c_param.Cn_alpha;
    i_a = (x[0] * x[0]) - (((x[1] * x[1]) + (x[2] * x[2])) + (x[3] * x[3]));
    j_a = 2.0 * x[0];
    for (i4 = 0; i4 < 3; i4++) {
      real_T b_a_tmp;
      int32_T c_a_tmp;
      int32_T d_a_tmp;
      b_a_tmp = x[i4 + 1];
      h_a[3 * i4] = (i_a * b_b[3 * i4]) + ((2.0 * x[1]) * b_a_tmp);
      c_a_tmp = (3 * i4) + 1;
      h_a[c_a_tmp] = (i_a * b_b[c_a_tmp]) + ((2.0 * x[2]) * b_a_tmp);
      d_a_tmp = (3 * i4) + 2;
      h_a[d_a_tmp] = (i_a * b_b[d_a_tmp]) + ((2.0 * x[3]) * b_a_tmp);
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
    memset(&b_w_exp_tilde[0], 0, 9U * (sizeof(real_T)));
    memset(&c_w_exp_tilde[0], 0, 3U * (sizeof(real_T)));
    for (i9 = 0; i9 < 3; i9++) {
      real_T d14;
      int32_T b_w_exp_tilde_tmp;
      int32_T c_w_exp_tilde_tmp;
      b_dv1[i9 + 1] = dn[i9] * b;
      d14 = b_w_exp_tilde[3 * i9];
      b_w_exp_tilde_tmp = (3 * i9) + 1;
      c_w_exp_tilde_tmp = (3 * i9) + 2;
      for (i10 = 0; i10 < 3; i10++) {
        real_T d15;
        d15 = SD->pd->c_param.J[i10 + (3 * i9)];
        d14 += w_exp_tilde[3 * i10] * d15;
        b_w_exp_tilde[b_w_exp_tilde_tmp] += w_exp_tilde[(3 * i10) + 1] * d15;
        b_w_exp_tilde[c_w_exp_tilde_tmp] += w_exp_tilde[(3 * i10) + 2] * d15;
      }
      real_T d_w_exp_tilde_tmp;
      b_w_exp_tilde[3 * i9] = d14;
      d_w_exp_tilde_tmp = x[i9 + 4];
      c_w_exp_tilde[0] += d14 * d_w_exp_tilde_tmp;
      c_w_exp_tilde[1] += b_w_exp_tilde[(3 * i9) + 1] * d_w_exp_tilde_tmp;
      c_w_exp_tilde[2] += b_w_exp_tilde[(3 * i9) + 2] * d_w_exp_tilde_tmp;
    }
    dv2[0] = 0.0;
    dv2[1] = e_a * (g_a * sin(atan2(x[9], x[7])));
    dv2[2] = e_a * (g_a * (-sin(atan2(x[8], x[7]))));
    memset(&dv3[0], 0, 3U * (sizeof(real_T)));
    memset(&b_dt[0], 0, 3U * (sizeof(real_T)));
    memset(&c_dt[0], 0, 3U * (sizeof(real_T)));
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
      real_T d28;
      real_T d29;
      real_T d32;
      real_T d35;
      real_T d36;
      real_T d38;
      real_T d39;
      int32_T i13;
      int32_T i14;
      d28 = SD->pd->c_param.Jinv[3 * i11];
      d29 = c_w_exp_tilde[i11];
      d16 += d28 * d29;
      d32 = dv2[i11];
      d19 += (dt * d28) * d32;
      d35 = S[3 * i11];
      d36 = SD->pd->c_param.g[i11];
      d22 += (dt * d35) * d36;
      d38 = d35 * d25;
      i13 = (3 * i11) + 1;
      d28 = SD->pd->c_param.Jinv[i13];
      d17 += d28 * d29;
      d20 += (dt * d28) * d32;
      d35 = S[i13];
      d23 += (dt * d35) * d36;
      d38 += d35 * d26;
      i14 = (3 * i11) + 2;
      d28 = SD->pd->c_param.Jinv[i14];
      d18 += d28 * d29;
      d21 += (dt * d28) * d32;
      d35 = S[i14];
      d24 += (dt * d35) * d36;
      d38 += d35 * d27;
      d39 = C_total_a[i11];
      c_w_exp_tilde[i11] =
          (((w_exp_tilde[i11] * d25) + (w_exp_tilde[i11 + 3] * d26)) +
           (w_exp_tilde[i11 + 6] * d27)) +
          (dt *
           ((((C_total_a_tmp_tmp / d39) * sens_input->board_accel.meas[i11]) +
             ((b_C_total_a_tmp_tmp / d39) * sens_input->mti_accel.meas[i11])) +
            ((c_C_total_a_tmp_tmp / d39) * sens_input->ad_accel.meas[i11])));
      b_S[i11] = d38;
    }
    memset(&c_q[0], 0, (sizeof(real_T)) << 2);
    d30 = c_q[0];
    d31 = c_q[1];
    d33 = c_q[2];
    d34 = c_q[3];
    for (i12 = 0; i12 < 4; i12++) {
      real_T d37;
      d37 = b_dv1[i12];
      d30 += b_q[4 * i12] * d37;
      d31 += b_q[(4 * i12) + 1] * d37;
      d33 += b_q[(4 * i12) + 2] * d37;
      d34 += b_q[(4 * i12) + 3] * d37;
    }
    real_T W_dt_tmp;
    real_T b_W_dt_tmp;
    real_T c_W_dt_tmp;
    real_T d_W_dt_tmp;
    real_T e_W_dt_tmp;
    real_T f_W_dt_tmp;
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
    x_pred[10] = x[10] + (dt * b_S[0]);
    memset(&F[0], 0, 121U * (sizeof(real_T)));
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
    memset(&c_b[0], 0, (sizeof(real_T)) << 4);
    for (i15 = 0; i15 < 4; i15++) {
      real_T d40;
      int32_T b_b_tmp;
      int32_T b_tmp;
      int32_T c_b_tmp;
      d40 = c_b[4 * i15];
      b_tmp = (4 * i15) + 1;
      b_b_tmp = (4 * i15) + 2;
      c_b_tmp = (4 * i15) + 3;
      for (i17 = 0; i17 < 4; i17++) {
        real_T d41;
        d41 = W_dt[i17 + (4 * i15)];
        d40 += W_dt[4 * i17] * d41;
        c_b[b_tmp] += W_dt[(4 * i17) + 1] * d41;
        c_b[b_b_tmp] += W_dt[(4 * i17) + 2] * d41;
        c_b[c_b_tmp] += W_dt[(4 * i17) + 3] * d41;
      }
      c_b[4 * i15] = d40;
    }
    for (i16 = 0; i16 < 16; i16++) {
      b_I[i16] = (int8_T)0;
    }
    memset(&b_W_dt[0], 0, (sizeof(real_T)) << 4);
    memset(&d_b[0], 0, (sizeof(real_T)) << 4);
    for (b_k = 0; b_k < 4; b_k++) {
      real_T d42;
      real_T d43;
      int32_T g_W_dt_tmp;
      int32_T h_W_dt_tmp;
      int32_T i_W_dt_tmp;
      b_I[b_k + (4 * b_k)] = (int8_T)1;
      d42 = b_W_dt[4 * b_k];
      d43 = d_b[4 * b_k];
      g_W_dt_tmp = (4 * b_k) + 1;
      h_W_dt_tmp = (4 * b_k) + 2;
      i_W_dt_tmp = (4 * b_k) + 3;
      for (i18 = 0; i18 < 4; i18++) {
        real_T d44;
        int32_T j_W_dt_tmp;
        int32_T k_W_dt_tmp;
        int32_T l_W_dt_tmp;
        d44 = c_b[i18 + (4 * b_k)];
        d42 += W_dt[4 * i18] * d44;
        d43 += c_b[4 * i18] * d44;
        j_W_dt_tmp = (4 * i18) + 1;
        b_W_dt[g_W_dt_tmp] += W_dt[j_W_dt_tmp] * d44;
        d_b[g_W_dt_tmp] += c_b[j_W_dt_tmp] * d44;
        k_W_dt_tmp = (4 * i18) + 2;
        b_W_dt[h_W_dt_tmp] += W_dt[k_W_dt_tmp] * d44;
        d_b[h_W_dt_tmp] += c_b[k_W_dt_tmp] * d44;
        l_W_dt_tmp = (4 * i18) + 3;
        b_W_dt[i_W_dt_tmp] += W_dt[l_W_dt_tmp] * d44;
        d_b[i_W_dt_tmp] += c_b[l_W_dt_tmp] * d44;
      }
      int32_T F_tmp;
      int32_T b_F_tmp;
      int32_T c_F_tmp;
      d_b[4 * b_k] = d43;
      b_W_dt[4 * b_k] = d42;
      F[11 * b_k] =
          (((((real_T)b_I[4 * b_k]) + W_dt[4 * b_k]) + (0.5 * c_b[4 * b_k])) +
           (0.16666666666666666 * d42)) +
          (0.041666666666666664 * d43);
      F_tmp = (4 * b_k) + 1;
      F[(11 * b_k) + 1] =
          (((((real_T)b_I[F_tmp]) + W_dt[F_tmp]) + (0.5 * c_b[F_tmp])) +
           (0.16666666666666666 * b_W_dt[F_tmp])) +
          (0.041666666666666664 * d_b[F_tmp]);
      b_F_tmp = (4 * b_k) + 2;
      F[(11 * b_k) + 2] =
          (((((real_T)b_I[b_F_tmp]) + W_dt[b_F_tmp]) + (0.5 * c_b[b_F_tmp])) +
           (0.16666666666666666 * b_W_dt[b_F_tmp])) +
          (0.041666666666666664 * d_b[b_F_tmp]);
      c_F_tmp = (4 * b_k) + 3;
      F[(11 * b_k) + 3] =
          (((((real_T)b_I[c_F_tmp]) + W_dt[c_F_tmp]) + (0.5 * c_b[c_F_tmp])) +
           (0.16666666666666666 * b_W_dt[c_F_tmp])) +
          (0.041666666666666664 * d_b[c_F_tmp]);
    }
    real_T e_a_tmp;
    real_T f_a_tmp;
    real_T g_a_tmp;
    real_T h_a_tmp;
    real_T i_a_tmp;
    real_T j_a_tmp;
    real_T k_a_tmp;
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
    for (i19 = 0; i19 < 3; i19++) {
      int32_T d_F_tmp;
      int32_T e_F_tmp;
      d_F_tmp = 4 * (i19 + 1);
      e_F_tmp = 11 * (i19 + 4);
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
    memset(&b_n_tilde[0], 0, 9U * (sizeof(real_T)));
    for (i20 = 0; i20 < 3; i20++) {
      real_T d45;
      int32_T c_n_tilde_tmp;
      int32_T d_n_tilde_tmp;
      d45 = b_n_tilde[3 * i20];
      c_n_tilde_tmp = (3 * i20) + 1;
      d_n_tilde_tmp = (3 * i20) + 2;
      for (i22 = 0; i22 < 3; i22++) {
        real_T d46;
        d46 = n_tilde[i22 + (3 * i20)];
        d45 += n_tilde[3 * i22] * d46;
        b_n_tilde[c_n_tilde_tmp] += n_tilde[(3 * i22) + 1] * d46;
        b_n_tilde[d_n_tilde_tmp] += n_tilde[(3 * i22) + 2] * d46;
      }
      b_n_tilde[3 * i20] = d45;
    }
    for (i21 = 0; i21 < 9; i21++) {
      w_exp_tilde[i21] =
          (((real_T)w_exp_tilde_tmp[i21]) - (b_a * n_tilde[i21])) +
          ((1.0 - b_x) * b_n_tilde[i21]);
    }
    memset(&b_dv[0], 0, 9U * (sizeof(real_T)));
    for (i23 = 0; i23 < 3; i23++) {
      real_T d47;
      int32_T i25;
      int32_T i26;
      d47 = b_dv[3 * i23];
      i25 = (3 * i23) + 1;
      i26 = (3 * i23) + 2;
      for (i28 = 0; i28 < 3; i28++) {
        real_T d49;
        d49 = w_exp_tilde[i28 + (3 * i23)];
        d47 += SD->pd->d_param.Jinv[3 * i28] * d49;
        b_dv[i25] += SD->pd->d_param.Jinv[(3 * i28) + 1] * d49;
        b_dv[i26] += SD->pd->d_param.Jinv[(3 * i28) + 2] * d49;
        F[(i28 + (11 * (i23 + 4))) + 4] = 0.0;
      }
      b_dv[3 * i23] = d47;
    }
    for (i24 = 0; i24 < 3; i24++) {
      int32_T F_tmp_tmp;
      F_tmp_tmp = 11 * (i24 + 4);
      for (i27 = 0; i27 < 3; i27++) {
        real_T d48;
        d48 = SD->pd->d_param.J[i27 + (3 * i24)];
        F[F_tmp_tmp + 4] += b_dv[3 * i27] * d48;
        F[F_tmp_tmp + 5] += b_dv[(3 * i27) + 1] * d48;
        F[F_tmp_tmp + 6] += b_dv[(3 * i27) + 2] * d48;
      }
    }
    b_dv[1] = t1_density * (m_a * x[9]);
    b_dv[4] = 0.0;
    b_dv[7] = t1_density * (m_a * x[7]);
    b_dv[2] = t1_density * (m_a * (-x[8]));
    b_dv[5] = t1_density * (m_a * (-x[7]));
    b_dv[8] = 0.0;
    c_x = 0.0;
    for (i29 = 0; i29 < 3; i29++) {
      real_T d50;
      real_T d51;
      real_T d52;
      int32_T f_F_tmp;
      b_dv[3 * i29] = 0.0;
      f_F_tmp = 11 * (i29 + 7);
      d50 = 0.0;
      d51 = 0.0;
      d52 = 0.0;
      for (i30 = 0; i30 < 3; i30++) {
        real_T d53;
        d53 = b_dv[i30 + (3 * i29)];
        d50 += (dt * SD->pd->d_param.Jinv[3 * i30]) * d53;
        d51 += (dt * SD->pd->d_param.Jinv[(3 * i30) + 1]) * d53;
        d52 += (dt * SD->pd->d_param.Jinv[(3 * i30) + 2]) * d53;
      }
      F[f_F_tmp + 6] = d52;
      F[f_F_tmp + 5] = d51;
      F[f_F_tmp + 4] = d50;
      c_x += x[i29 + 1] * SD->pd->d_param.g[i29];
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
    memset(&c_skewed_exp_w_tmp[0], 0, 9U * (sizeof(real_T)));
    memset(&b_dv[0], 0, 9U * (sizeof(real_T)));
    r_q_tmp[0] = x[0];
    b_r_q_tmp = 0.0;
    for (i31 = 0; i31 < 3; i31++) {
      real_T d54;
      real_T d55;
      real_T g_F_tmp;
      int32_T h_F_tmp;
      int32_T i_F_tmp;
      int32_T j_F_tmp;
      F[i31 + 7] = dt * (2.0 * ((x[0] * SD->pd->d_param.g[i31]) - d_x[i31]));
      g_F_tmp = x[i31 + 1];
      h_F_tmp = 11 * (i31 + 1);
      F[h_F_tmp + 7] =
          dt *
          (2.0 * ((((c_x * b_b[3 * i31]) + (x[1] * SD->pd->d_param.g[i31])) -
                   (SD->pd->d_param.g[0] * g_F_tmp)) +
                  dv4[3 * i31]));
      i_F_tmp = (3 * i31) + 1;
      F[h_F_tmp + 8] =
          dt *
          (2.0 * ((((c_x * b_b[i_F_tmp]) + (x[2] * SD->pd->d_param.g[i31])) -
                   (SD->pd->d_param.g[1] * g_F_tmp)) +
                  dv4[i_F_tmp]));
      j_F_tmp = (3 * i31) + 2;
      F[h_F_tmp + 9] =
          dt *
          (2.0 * ((((c_x * b_b[j_F_tmp]) + (x[3] * SD->pd->d_param.g[i31])) -
                   (SD->pd->d_param.g[2] * g_F_tmp)) +
                  dv4[j_F_tmp]));
      d54 = c_skewed_exp_w_tmp[3 * i31];
      d55 = b_dv[3 * i31];
      for (i33 = 0; i33 < 3; i33++) {
        real_T d56;
        real_T d57;
        int32_T b_skewed_exp_w_tmp_tmp;
        int32_T i34;
        int32_T skewed_exp_w_tmp_tmp;
        i34 = i33 + (3 * i31);
        d56 = b_skewed_exp_w_tmp[i34];
        d57 = skewed_exp_w_tmp[i34];
        d54 += skewed_exp_w_tmp[3 * i33] * d56;
        d55 += (2.0 * b_skewed_exp_w_tmp[3 * i33]) * d57;
        skewed_exp_w_tmp_tmp = (3 * i33) + 1;
        c_skewed_exp_w_tmp[i_F_tmp] +=
            skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d56;
        b_dv[i_F_tmp] += (2.0 * b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp]) * d57;
        b_skewed_exp_w_tmp_tmp = (3 * i33) + 2;
        c_skewed_exp_w_tmp[j_F_tmp] +=
            skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d56;
        b_dv[j_F_tmp] +=
            (2.0 * b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp]) * d57;
      }
      int32_T k_F_tmp;
      int32_T l_F_tmp;
      b_dv[3 * i31] = d55;
      c_skewed_exp_w_tmp[3 * i31] = d54;
      k_F_tmp = 11 * (i31 + 4);
      F[k_F_tmp + 7] = (dt * skewed_exp_w_tmp[3 * i31]) + (n_a * (d54 - d55));
      l_F_tmp = 11 * (i31 + 7);
      F[l_F_tmp + 7] = w_exp_tilde[3 * i31];
      F[k_F_tmp + 8] = (dt * skewed_exp_w_tmp[i_F_tmp]) +
                       (n_a * (c_skewed_exp_w_tmp[i_F_tmp] - b_dv[i_F_tmp]));
      F[l_F_tmp + 8] = w_exp_tilde[i_F_tmp];
      F[k_F_tmp + 9] = (dt * skewed_exp_w_tmp[j_F_tmp]) +
                       (n_a * (c_skewed_exp_w_tmp[j_F_tmp] - b_dv[j_F_tmp]));
      F[l_F_tmp + 9] = w_exp_tilde[j_F_tmp];
      r_q_tmp[i31 + 1] = -g_F_tmp;
      b_r_q_tmp += (-g_F_tmp) * x[i31 + 7];
    }
    c_r_q_tmp[0] = (r_q_tmp[2] * x[9]) - (r_q_tmp[3] * x[8]);
    c_r_q_tmp[1] = (r_q_tmp[3] * x[7]) - (r_q_tmp[1] * x[9]);
    c_r_q_tmp[2] = (r_q_tmp[1] * x[8]) - (r_q_tmp[2] * x[7]);
    for (i32 = 0; i32 < 3; i32++) {
      real_T b_dt_tmp;
      real_T dt_tmp;
      int32_T c_dt_tmp;
      int32_T d_dt_tmp;
      int32_T e_dt_tmp;
      dt_tmp = x[i32 + 7];
      d_dt[i32] = dt * (2.0 * ((r_q_tmp[0] * dt_tmp) - c_r_q_tmp[i32]));
      b_dt_tmp = r_q_tmp[i32 + 1];
      c_dt_tmp = 3 * (i32 + 1);
      d_dt[c_dt_tmp] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[3 * i32]) + (r_q_tmp[1] * dt_tmp)) -
                        (x[7] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[3 * i32])));
      d_dt_tmp = (3 * i32) + 1;
      d_dt[c_dt_tmp + 1] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[d_dt_tmp]) + (r_q_tmp[2] * dt_tmp)) -
                        (x[8] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[d_dt_tmp])));
      e_dt_tmp = (3 * i32) + 2;
      d_dt[c_dt_tmp + 2] =
          dt * (2.0 * ((((b_r_q_tmp * b_b[e_dt_tmp]) + (r_q_tmp[3] * dt_tmp)) -
                        (x[9] * b_dt_tmp)) +
                       (r_q_tmp[0] * skewed_exp_w_tmp[e_dt_tmp])));
    }
    real_T p_a;
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
    for (i35 = 0; i35 < 3; i35++) {
      F[(11 * (i35 + 7)) + 10] =
          dt *
          (((o_a * b_b[3 * i35]) + ((2.0 * r_q_tmp[1]) * r_q_tmp[i35 + 1])) -
           b_dv[3 * i35]);
    }
    F[120] = 1.0;
    memset(&b_F[0], 0, 121U * (sizeof(real_T)));
    for (i36 = 0; i36 < 11; i36++) {
      for (i37 = 0; i37 < 11; i37++) {
        real_T d58;
        d58 = P[i37 + (11 * i36)];
        for (i40 = 0; i40 < 11; i40++) {
          int32_T m_F_tmp;
          m_F_tmp = i40 + (11 * i36);
          b_F[m_F_tmp] += F[i40 + (11 * i37)] * d58;
        }
      }
    }
    for (i38 = 0; i38 < 11; i38++) {
      for (i39 = 0; i39 < 11; i39++) {
        real_T d59;
        d59 = 0.0;
        for (i42 = 0; i42 < 11; i42++) {
          d59 += b_F[i38 + (11 * i42)] * F[i39 + (11 * i42)];
        }
        int32_T c_P_pred_tmp;
        c_P_pred_tmp = i38 + (11 * i39);
        P_pred[c_P_pred_tmp] = d59 + Q[c_P_pred_tmp];
      }
    }
    for (i41 = 0; i41 < 3; i41++) {
      int32_T P_pred_tmp;
      int32_T b_P_pred_tmp;
      int32_T d_P_pred_tmp;
      P_pred_tmp = 11 * (i41 + 4);
      b_P_pred[3 * i41] = P_pred[P_pred_tmp + 4] + R[3 * i41];
      b_P_pred_tmp = (3 * i41) + 1;
      b_P_pred[b_P_pred_tmp] = P_pred[P_pred_tmp + 5] + R[b_P_pred_tmp];
      d_P_pred_tmp = (3 * i41) + 2;
      b_P_pred[d_P_pred_tmp] = P_pred[P_pred_tmp + 6] + R[d_P_pred_tmp];
    }
    inv(b_P_pred, b_dv);
    memset(&K[0], 0, 33U * (sizeof(real_T)));
    for (i43 = 0; i43 < 3; i43++) {
      for (i44 = 0; i44 < 3; i44++) {
        real_T d60;
        d60 = b_dv[i44 + (3 * i43)];
        for (i45 = 0; i45 < 11; i45++) {
          int32_T K_tmp;
          K_tmp = i45 + (11 * i43);
          K[K_tmp] += P_pred[i45 + (11 * (i44 + 4))] * d60;
        }
      }
    }
    memset(&c_I[0], 0, 121U * (sizeof(int8_T)));
    for (c_k = 0; c_k < 11; c_k++) {
      c_I[c_k + (11 * c_k)] = (int8_T)1;
    }
    for (i46 = 0; i46 < 44; i46++) {
      E[i46] = (real_T)c_I[i46];
    }
    for (i47 = 0; i47 < 33; i47++) {
      E[i47 + 44] = ((real_T)c_I[i47 + 44]) - K[i47];
    }
    for (i48 = 0; i48 < 44; i48++) {
      E[i48 + 77] = (real_T)c_I[i48 + 77];
    }
    memset(&b_E[0], 0, 121U * (sizeof(real_T)));
    for (i49 = 0; i49 < 11; i49++) {
      for (i50 = 0; i50 < 11; i50++) {
        real_T d61;
        d61 = P_pred[i50 + (11 * i49)];
        for (i52 = 0; i52 < 11; i52++) {
          int32_T E_tmp;
          E_tmp = i52 + (11 * i49);
          b_E[E_tmp] += E[i52 + (11 * i50)] * d61;
        }
      }
    }
    memset(&b_K[0], 0, 33U * (sizeof(real_T)));
    for (i51 = 0; i51 < 3; i51++) {
      for (i53 = 0; i53 < 3; i53++) {
        real_T d62;
        d62 = R[i53 + (3 * i51)];
        for (i54 = 0; i54 < 11; i54++) {
          int32_T b_K_tmp;
          b_K_tmp = i54 + (11 * i51);
          b_K[b_K_tmp] += K[i54 + (11 * i53)] * d62;
        }
      }
    }
    memset(&c_E[0], 0, 121U * (sizeof(real_T)));
    memset(&c_K[0], 0, 121U * (sizeof(real_T)));
    for (i55 = 0; i55 < 11; i55++) {
      for (i56 = 0; i56 < 11; i56++) {
        real_T d63;
        d63 = E[i55 + (11 * i56)];
        for (i59 = 0; i59 < 11; i59++) {
          int32_T b_E_tmp;
          b_E_tmp = i59 + (11 * i55);
          c_E[b_E_tmp] += b_E[i59 + (11 * i56)] * d63;
        }
      }
      for (i58 = 0; i58 < 3; i58++) {
        real_T d65;
        d65 = K[i55 + (11 * i58)];
        for (i60 = 0; i60 < 11; i60++) {
          int32_T c_K_tmp;
          c_K_tmp = i60 + (11 * i55);
          c_K[c_K_tmp] += b_K[i60 + (11 * i58)] * d65;
        }
      }
    }
    for (i57 = 0; i57 < 121; i57++) {
      P_ret[i57] = c_E[i57] + c_K[i57];
    }
    real_T d64;
    real_T d66;
    d64 = d1 / d3;
    d66 = d2 / d3;
    d67 =
        ((((d1 / d4) * (sens_input->board_gyro.meas[0] - bias->board_gyro[0])) +
          ((d2 / d4) * (sens_input->mti_gyro.meas[0] - bias->mti_gyro[0]))) +
         (C_ad_w_idx_0 * (sens_input->ad_gyro.meas[0] - bias->ad_gyro[0]))) -
        x_pred[4];
    d68 = (((d64 * (sens_input->board_gyro.meas[1] - bias->board_gyro[1])) +
            (d66 * (sens_input->mti_gyro.meas[1] - bias->mti_gyro[1]))) +
           (d * (sens_input->ad_gyro.meas[1] - bias->ad_gyro[1]))) -
          x_pred[5];
    d69 = (((d64 * (sens_input->board_gyro.meas[2] - bias->board_gyro[2])) +
            (d66 * (sens_input->mti_gyro.meas[2] - bias->mti_gyro[2]))) +
           (d * (sens_input->ad_gyro.meas[2] - bias->ad_gyro[2]))) -
          x_pred[6];
    for (i61 = 0; i61 < 11; i61++) {
      x_ret[i61] = x_pred[i61] + (((K[i61] * d67) + (K[i61 + 11] * d68)) +
                                  (K[i61 + 22] * d69));
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
      memcpy(&b_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&b_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      memcpy(&e_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&d_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      b_ekf_correct(e_x_ret, d_P_ret, sens_input->board_baro.meas,
                    bias->board_baro, x_ret, P_ret);
    }
    if (sens_input->mti_baro.status) {
      memcpy(&c_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&c_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      memcpy(&g_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&f_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      b_ekf_correct(g_x_ret, f_P_ret, sens_input->mti_baro.meas, bias->mti_baro,
                    x_ret, P_ret);
    }
    if (sens_input->board_mag.status) {
      memcpy(&d_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&e_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      memcpy(&h_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&h_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      ekf_correct(h_x_ret, h_P_ret, sens_input->board_mag.meas,
                  bias->board_mag_earth, b_b, x_ret, P_ret);
    }
    if (sens_input->mti_mag.status) {
      memcpy(&f_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&g_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      memcpy(&i_x_ret[0], &x_ret[0], 11U * (sizeof(real_T)));
      memcpy(&i_P_ret[0], &P_ret[0], 121U * (sizeof(real_T)));
      ekf_correct(i_x_ret, i_P_ret, sens_input->mti_mag.meas,
                  bias->mti_mag_earth, b_b, x_ret, P_ret);
    }
  }
  *cov_norm = c_norm(P);
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
