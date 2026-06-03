/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: GNC_codegen_types.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 23:24:33
 */

#ifndef GNC_CODEGEN_TYPES_H
#define GNC_CODEGEN_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  double board_gyro[3];
  double mti_gyro[3];
  double ad_gyro[3];
  double board_mag_earth[3];
  double mti_mag_earth[3];
  double board_baro;
  double mti_baro;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
typedef struct {
  double board_accel_f[3];
  double board_gyro_f[3];
  double mti_accel_f[3];
  double mti_gyro_f[3];
  double ad_accel_f[3];
  double ad_gyro_f[3];
  double board_baro_f;
  double board_mag_f[3];
  double mti_baro_f;
  double mti_mag_f[3];
} struct2_T;
#endif /* typedef_struct2_T */

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  double meas[3];
  bool status;
} struct4_T;
#endif /* typedef_struct4_T */

#ifndef typedef_struct5_T
#define typedef_struct5_T
typedef struct {
  double meas;
  bool status;
} struct5_T;
#endif /* typedef_struct5_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
typedef struct {
  struct4_T board_accel;
  struct4_T board_gyro;
  struct4_T mti_accel;
  struct4_T mti_gyro;
  struct4_T ad_accel;
  struct4_T ad_gyro;
  struct5_T board_baro;
  struct4_T board_mag;
  struct5_T mti_baro;
  struct4_T mti_mag;
} struct3_T;
#endif /* typedef_struct3_T */

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  double coeffs[2];
  double w_old;
  double P_minus[4];
  double d_old;
  double w_dot_old;
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  double J[9];
  double Jinv[9];
  double g[3];
} struct_T;
#endif /* typedef_struct_T */

#endif
/*
 * File trailer for GNC_codegen_types.h
 *
 * [EOF]
 */
