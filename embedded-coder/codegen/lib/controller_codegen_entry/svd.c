/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: svd.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
 */

/* Include Files */
#include "svd.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include "xrotg.h"
#include "xzlangeM.h"
#include "xzlascl.h"
#include <emmintrin.h>
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * Arguments    : const double A[121]
 *                double U[11]
 * Return Type  : void
 */
void svd(const double A[121], double U[11])
{
  double b_A[121];
  double e[11];
  double work[11];
  double anrm;
  double b;
  double cscale;
  double emm1;
  double nrm;
  double rt;
  double sm;
  double snorm;
  int jj;
  int k;
  int kase;
  int mm1;
  int q;
  int qjj;
  int qp1;
  int qq;
  int qq_tmp;
  boolean_T doscale;
  memcpy(&b_A[0], &A[0], 121U * sizeof(double));
  memset(&U[0], 0, 11U * sizeof(double));
  memset(&e[0], 0, 11U * sizeof(double));
  memset(&work[0], 0, 11U * sizeof(double));
  doscale = false;
  anrm = xzlangeM(A);
  cscale = anrm;
  if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
    doscale = true;
    cscale = 6.7178761075670888E-139;
    xzlascl(anrm, cscale, b_A);
  } else if (anrm > 1.4885657073574029E+138) {
    doscale = true;
    cscale = 1.4885657073574029E+138;
    xzlascl(anrm, cscale, b_A);
  }
  for (q = 0; q < 10; q++) {
    __m128d b_r;
    boolean_T apply_transform;
    qp1 = q + 2;
    qq_tmp = q + 11 * q;
    qq = qq_tmp + 1;
    apply_transform = false;
    nrm = xnrm2(11 - q, b_A, qq_tmp + 1);
    if (nrm > 0.0) {
      apply_transform = true;
      if (b_A[qq_tmp] < 0.0) {
        nrm = -nrm;
      }
      U[q] = nrm;
      if (fabs(nrm) >= 1.0020841800044864E-292) {
        nrm = 1.0 / nrm;
        kase = (qq_tmp - q) + 11;
        mm1 = ((((kase - qq_tmp) / 2) << 1) + qq_tmp) + 1;
        qjj = mm1 - 2;
        for (k = qq; k <= qjj; k += 2) {
          b_r = _mm_loadu_pd(&b_A[k - 1]);
          _mm_storeu_pd(&b_A[k - 1], _mm_mul_pd(_mm_set1_pd(nrm), b_r));
        }
        for (k = mm1; k <= kase; k++) {
          b_A[k - 1] *= nrm;
        }
      } else {
        kase = (qq_tmp - q) + 11;
        mm1 = ((((kase - qq_tmp) / 2) << 1) + qq_tmp) + 1;
        qjj = mm1 - 2;
        for (k = qq; k <= qjj; k += 2) {
          b_r = _mm_loadu_pd(&b_A[k - 1]);
          _mm_storeu_pd(&b_A[k - 1], _mm_div_pd(b_r, _mm_set1_pd(U[q])));
        }
        for (k = mm1; k <= kase; k++) {
          b_A[k - 1] /= U[q];
        }
      }
      b_A[qq_tmp]++;
      U[q] = -U[q];
    } else {
      U[q] = 0.0;
    }
    for (jj = qp1; jj < 12; jj++) {
      qjj = q + 11 * (jj - 1);
      if (apply_transform) {
        nrm = 0.0;
        kase = 10 - q;
        for (k = 0; k <= kase; k++) {
          nrm += b_A[qq_tmp + k] * b_A[qjj + k];
        }
        nrm = -(nrm / b_A[qq_tmp]);
        if (!(nrm == 0.0)) {
          kase = 11 - q;
          for (k = 0; k < kase; k++) {
            mm1 = qjj + k;
            b_A[mm1] += nrm * b_A[qq_tmp + k];
          }
        }
      }
      e[jj - 1] = b_A[qjj];
    }
    if (q + 1 <= 9) {
      nrm = b_xnrm2(10 - q, e, q + 2);
      if (nrm == 0.0) {
        e[q] = 0.0;
      } else {
        __m128d r1;
        if (e[q + 1] < 0.0) {
          e[q] = -nrm;
        } else {
          e[q] = nrm;
        }
        nrm = e[q];
        if (fabs(e[q]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[q];
          kase = ((((10 - q) / 2) << 1) + q) + 2;
          mm1 = kase - 2;
          for (k = qp1; k <= mm1; k += 2) {
            b_r = _mm_loadu_pd(&e[k - 1]);
            _mm_storeu_pd(&e[k - 1], _mm_mul_pd(_mm_set1_pd(nrm), b_r));
          }
          for (k = kase; k < 12; k++) {
            e[k - 1] *= nrm;
          }
        } else {
          kase = ((((10 - q) / 2) << 1) + q) + 2;
          mm1 = kase - 2;
          for (k = qp1; k <= mm1; k += 2) {
            b_r = _mm_loadu_pd(&e[k - 1]);
            _mm_storeu_pd(&e[k - 1], _mm_div_pd(b_r, _mm_set1_pd(nrm)));
          }
          for (k = kase; k < 12; k++) {
            e[k - 1] /= nrm;
          }
        }
        e[q + 1]++;
        e[q] = -e[q];
        for (k = qp1; k < 12; k++) {
          work[k - 1] = 0.0;
        }
        for (k = qp1; k < 12; k++) {
          nrm = e[k - 1];
          if (!(nrm == 0.0)) {
            qjj = q + 11 * (k - 1);
            qq = 10 - q;
            qq_tmp = ((10 - q) / 2) << 1;
            kase = qq_tmp - 2;
            for (jj = 0; jj <= kase; jj += 2) {
              b_r = _mm_loadu_pd(&b_A[(qjj + jj) + 1]);
              mm1 = (q + jj) + 1;
              r1 = _mm_loadu_pd(&work[mm1]);
              _mm_storeu_pd(&work[mm1],
                            _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(nrm), b_r)));
            }
            for (jj = qq_tmp; jj < qq; jj++) {
              kase = (q + jj) + 1;
              work[kase] += nrm * b_A[(qjj + jj) + 1];
            }
          }
        }
        for (k = qp1; k < 12; k++) {
          nrm = -e[k - 1] / e[q + 1];
          if (!(nrm == 0.0)) {
            qjj = (q + 11 * (k - 1)) + 1;
            qq = 10 - q;
            qq_tmp = ((10 - q) / 2) << 1;
            kase = qq_tmp - 2;
            for (jj = 0; jj <= kase; jj += 2) {
              b_r = _mm_loadu_pd(&work[(q + jj) + 1]);
              mm1 = qjj + jj;
              r1 = _mm_loadu_pd(&b_A[mm1]);
              _mm_storeu_pd(&b_A[mm1],
                            _mm_add_pd(r1, _mm_mul_pd(_mm_set1_pd(nrm), b_r)));
            }
            for (jj = qq_tmp; jj < qq; jj++) {
              kase = qjj + jj;
              b_A[kase] += nrm * work[(q + jj) + 1];
            }
          }
        }
      }
    }
  }
  qq = 9;
  U[10] = b_A[120];
  e[9] = b_A[119];
  e[10] = 0.0;
  qjj = 0;
  snorm = 0.0;
  for (k = 0; k < 11; k++) {
    nrm = U[k];
    emm1 = nrm;
    if (nrm != 0.0) {
      rt = fabs(nrm);
      emm1 = rt;
      U[k] = rt;
      if (k + 1 < 11) {
        e[k] /= nrm / rt;
      }
    }
    if (k + 1 < 11) {
      nrm = e[k];
      if (nrm != 0.0) {
        rt = fabs(nrm);
        e[k] = rt;
        U[k + 1] *= rt / nrm;
      }
    }
    snorm = fmax(snorm, fmax(fabs(emm1), fabs(e[k])));
  }
  while ((qq + 2 > 0) && (qjj < 75)) {
    boolean_T exitg1;
    qq_tmp = qq + 1;
    exitg1 = false;
    while (!(exitg1 || (qq_tmp == 0))) {
      nrm = fabs(e[qq_tmp - 1]);
      if ((nrm <=
           2.2204460492503131E-16 * (fabs(U[qq_tmp - 1]) + fabs(U[qq_tmp]))) ||
          (nrm <= 1.0020841800044864E-292) ||
          ((qjj > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
        e[qq_tmp - 1] = 0.0;
        exitg1 = true;
      } else {
        qq_tmp--;
      }
    }
    if (qq_tmp == qq + 1) {
      kase = 4;
    } else {
      mm1 = qq + 2;
      kase = qq + 2;
      exitg1 = false;
      while ((!exitg1) && (kase >= qq_tmp)) {
        mm1 = kase;
        if (kase == qq_tmp) {
          exitg1 = true;
        } else {
          nrm = 0.0;
          if (kase < qq + 2) {
            nrm = fabs(e[kase - 1]);
          }
          if (kase > qq_tmp + 1) {
            nrm += fabs(e[kase - 2]);
          }
          emm1 = fabs(U[kase - 1]);
          if ((emm1 <= 2.2204460492503131E-16 * nrm) ||
              (emm1 <= 1.0020841800044864E-292)) {
            U[kase - 1] = 0.0;
            exitg1 = true;
          } else {
            kase--;
          }
        }
      }
      if (mm1 == qq_tmp) {
        kase = 3;
      } else if (mm1 == qq + 2) {
        kase = 1;
      } else {
        kase = 2;
        qq_tmp = mm1;
      }
    }
    switch (kase) {
    case 1:
      b = e[qq];
      e[qq] = 0.0;
      kase = qq + 1;
      for (k = kase; k >= qq_tmp + 1; k--) {
        emm1 = xrotg(&U[k - 1], &b, &nrm);
        if (k > qq_tmp + 1) {
          rt = e[k - 2];
          b = -nrm * rt;
          e[k - 2] = rt * emm1;
        }
      }
      break;
    case 2:
      b = e[qq_tmp - 1];
      e[qq_tmp - 1] = 0.0;
      for (k = qq_tmp + 1; k <= qq + 2; k++) {
        emm1 = xrotg(&U[k - 1], &b, &nrm);
        rt = e[k - 1];
        b = -nrm * rt;
        e[k - 1] = rt * emm1;
      }
      break;
    case 3: {
      double scale;
      double sqds;
      mm1 = qq + 1;
      nrm = U[qq + 1];
      scale = fmax(fmax(fmax(fmax(fabs(nrm), fabs(U[qq])), fabs(e[qq])),
                        fabs(U[qq_tmp])),
                   fabs(e[qq_tmp]));
      sm = nrm / scale;
      nrm = U[qq] / scale;
      emm1 = e[qq] / scale;
      sqds = U[qq_tmp] / scale;
      b = ((nrm + sm) * (nrm - sm) + emm1 * emm1) / 2.0;
      nrm = sm * emm1;
      nrm *= nrm;
      if ((b != 0.0) || (nrm != 0.0)) {
        rt = sqrt(b * b + nrm);
        if (b < 0.0) {
          rt = -rt;
        }
        rt = nrm / (b + rt);
      } else {
        rt = 0.0;
      }
      b = (sqds + sm) * (sqds - sm) + rt;
      nrm = sqds * (e[qq_tmp] / scale);
      for (k = qq_tmp + 1; k <= mm1; k++) {
        rt = xrotg(&b, &nrm, &sm);
        if (k > qq_tmp + 1) {
          e[k - 2] = b;
        }
        nrm = e[k - 1];
        emm1 = U[k - 1];
        e[k - 1] = rt * nrm - sm * emm1;
        b = sm * U[k];
        U[k] *= rt;
        U[k - 1] = rt * emm1 + sm * nrm;
        rt = xrotg(&U[k - 1], &b, &emm1);
        nrm = e[k - 1];
        b = rt * nrm + emm1 * U[k];
        U[k] = -emm1 * nrm + rt * U[k];
        nrm = emm1 * e[k];
        e[k] *= rt;
      }
      e[qq] = b;
      qjj++;
    } break;
    default:
      if (U[qq_tmp] < 0.0) {
        U[qq_tmp] = -U[qq_tmp];
      }
      qp1 = qq_tmp + 1;
      while ((qq_tmp + 1 < 11) && (U[qq_tmp] < U[qp1])) {
        rt = U[qq_tmp];
        U[qq_tmp] = U[qp1];
        U[qp1] = rt;
        qq_tmp = qp1;
        qp1++;
      }
      qjj = 0;
      qq--;
      break;
    }
  }
  if (doscale) {
    b_xzlascl(cscale, anrm, U);
  }
}

/*
 * File trailer for svd.c
 *
 * [EOF]
 */
