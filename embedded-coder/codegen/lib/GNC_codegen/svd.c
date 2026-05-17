/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: svd.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:56:01
 */

/* Include Files */
#include "svd.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include "xrotg.h"
#include "xzlangeM.h"
#include "xzlascl.h"
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
  int ix;
  int iy;
  int jj;
  int k;
  int m;
  int q;
  int qp1;
  int qq;
  int qq_tmp;
  bool doscale;
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
    bool apply_transform;
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
        ix = (qq_tmp - q) + 11;
        for (k = qq; k <= ix; k++) {
          b_A[k - 1] *= nrm;
        }
      } else {
        ix = (qq_tmp - q) + 11;
        for (k = qq; k <= ix; k++) {
          b_A[k - 1] /= U[q];
        }
      }
      b_A[qq_tmp]++;
      U[q] = -U[q];
    } else {
      U[q] = 0.0;
    }
    for (jj = qp1; jj < 12; jj++) {
      qq = q + 11 * (jj - 1);
      if (apply_transform) {
        nrm = 0.0;
        ix = 10 - q;
        for (k = 0; k <= ix; k++) {
          nrm += b_A[qq_tmp + k] * b_A[qq + k];
        }
        nrm = -(nrm / b_A[qq_tmp]);
        if (!(nrm == 0.0)) {
          ix = 11 - q;
          for (k = 0; k < ix; k++) {
            iy = qq + k;
            b_A[iy] += nrm * b_A[qq_tmp + k];
          }
        }
      }
      e[jj - 1] = b_A[qq];
    }
    if (q + 1 <= 9) {
      nrm = b_xnrm2(10 - q, e, q + 2);
      if (nrm == 0.0) {
        e[q] = 0.0;
      } else {
        if (e[q + 1] < 0.0) {
          e[q] = -nrm;
        } else {
          e[q] = nrm;
        }
        nrm = e[q];
        if (fabs(e[q]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[q];
          for (k = qp1; k < 12; k++) {
            e[k - 1] *= nrm;
          }
        } else {
          for (k = qp1; k < 12; k++) {
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
            ix = q + 11 * (k - 1);
            iy = 10 - q;
            for (jj = 0; jj < iy; jj++) {
              qq = (q + jj) + 1;
              work[qq] += nrm * b_A[(ix + jj) + 1];
            }
          }
        }
        for (k = qp1; k < 12; k++) {
          nrm = -e[k - 1] / e[q + 1];
          if (!(nrm == 0.0)) {
            iy = (q + 11 * (k - 1)) + 1;
            ix = 10 - q;
            for (jj = 0; jj < ix; jj++) {
              qq = iy + jj;
              b_A[qq] += nrm * work[(q + jj) + 1];
            }
          }
        }
      }
    }
  }
  m = 9;
  U[10] = b_A[120];
  e[9] = b_A[119];
  e[10] = 0.0;
  qq = 0;
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
  while ((m + 2 > 0) && (qq < 75)) {
    bool exitg1;
    qq_tmp = m + 1;
    exitg1 = false;
    while (!(exitg1 || (qq_tmp == 0))) {
      nrm = fabs(e[qq_tmp - 1]);
      if ((nrm <=
           2.2204460492503131E-16 * (fabs(U[qq_tmp - 1]) + fabs(U[qq_tmp]))) ||
          (nrm <= 1.0020841800044864E-292) ||
          ((qq > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
        e[qq_tmp - 1] = 0.0;
        exitg1 = true;
      } else {
        qq_tmp--;
      }
    }
    if (qq_tmp == m + 1) {
      ix = 4;
    } else {
      iy = m + 2;
      ix = m + 2;
      exitg1 = false;
      while ((!exitg1) && (ix >= qq_tmp)) {
        iy = ix;
        if (ix == qq_tmp) {
          exitg1 = true;
        } else {
          nrm = 0.0;
          if (ix < m + 2) {
            nrm = fabs(e[ix - 1]);
          }
          if (ix > qq_tmp + 1) {
            nrm += fabs(e[ix - 2]);
          }
          emm1 = fabs(U[ix - 1]);
          if ((emm1 <= 2.2204460492503131E-16 * nrm) ||
              (emm1 <= 1.0020841800044864E-292)) {
            U[ix - 1] = 0.0;
            exitg1 = true;
          } else {
            ix--;
          }
        }
      }
      if (iy == qq_tmp) {
        ix = 3;
      } else if (iy == m + 2) {
        ix = 1;
      } else {
        ix = 2;
        qq_tmp = iy;
      }
    }
    switch (ix) {
    case 1:
      b = e[m];
      e[m] = 0.0;
      ix = m + 1;
      for (k = ix; k >= qq_tmp + 1; k--) {
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
      for (k = qq_tmp + 1; k <= m + 2; k++) {
        emm1 = xrotg(&U[k - 1], &b, &nrm);
        rt = e[k - 1];
        b = -nrm * rt;
        e[k - 1] = rt * emm1;
      }
      break;
    case 3: {
      double scale;
      double sqds;
      ix = m + 1;
      nrm = U[m + 1];
      scale = fmax(
          fmax(fmax(fmax(fabs(nrm), fabs(U[m])), fabs(e[m])), fabs(U[qq_tmp])),
          fabs(e[qq_tmp]));
      sm = nrm / scale;
      nrm = U[m] / scale;
      emm1 = e[m] / scale;
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
      for (k = qq_tmp + 1; k <= ix; k++) {
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
      e[m] = b;
      qq++;
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
      qq = 0;
      m--;
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
