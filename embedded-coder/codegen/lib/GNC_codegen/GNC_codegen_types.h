#ifndef GNC_CODEGEN_TYPES_H
#define GNC_CODEGEN_TYPES_H

#include "rtwtypes.h"

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
#endif

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
#endif

#ifndef typedef_struct4_T
#define typedef_struct4_T
typedef struct {
  double meas[3];
  bool status;
} struct4_T;
#endif

#ifndef typedef_struct5_T
#define typedef_struct5_T
typedef struct {
  double meas;
  bool status;
} struct5_T;
#endif

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
#endif

#ifndef typedef_struct6_T
#define typedef_struct6_T
typedef struct {
  double pressure;
  double temperature;
  double density;
  double sonic_speed;
  double mach;
  double dynamic_pressure;
} struct6_T;
#endif

#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  double coeffs[2];
  double w_old;
  double P_minus[4];
  double d_old;
  double w_dot_old;
} struct0_T;
#endif

#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  double Cn_alpha;
  double J[9];
  double Jinv[9];
  double c_aero;
  double c_canard;
  double elevation;
  double g[3];
} struct_T;
#endif

#ifndef c_typedef_GNC_codegenPersistent
#define c_typedef_GNC_codegenPersistent
typedef struct {
  struct_T param;
  struct_T b_param;
  struct_T c_param;
  struct_T d_param;
} GNC_codegenPersistentData;
#endif

#ifndef typedef_GNC_codegenStackData
#define typedef_GNC_codegenStackData
typedef struct {
  GNC_codegenPersistentData *pd;
} GNC_codegenStackData;
#endif

#endif
