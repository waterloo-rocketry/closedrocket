/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * navigation_codegen_entry_types.h
 *
 * Code generation for function 'navigation_codegen_entry'
 *
 */

#ifndef NAVIGATION_CODEGEN_ENTRY_TYPES_H
#define NAVIGATION_CODEGEN_ENTRY_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct0_T
#define typedef_struct0_T
typedef struct {
  double meas[3];
  boolean_T status;
} struct0_T;
#endif /* typedef_struct0_T */

#ifndef typedef_struct1_T
#define typedef_struct1_T
typedef struct {
  double meas;
  boolean_T status;
} struct1_T;
#endif /* typedef_struct1_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
typedef struct {
  double pressure;
  double temperature;
  double density;
  double sonic_speed;
  double mach;
  double dynamic_pressure;
} struct3_T;
#endif /* typedef_struct3_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
typedef struct {
  double q[4];
  double w[3];
  double v[3];
  double alt;
  double x[11];
} struct2_T;
#endif /* typedef_struct2_T */

#ifndef typedef_b_struct_T
#define typedef_b_struct_T
typedef struct {
  double J[9];
  double Jinv[9];
  double g[3];
} b_struct_T;
#endif /* typedef_b_struct_T */

#endif
/* End of code generation (navigation_codegen_entry_types.h) */
