/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: airdata_atmos.h
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 02-Jun-2026 23:24:33
 */

#ifndef AIRDATA_ATMOS_H
#define AIRDATA_ATMOS_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
double airdata_atmos(double altitude, double *airdata_temperature,
                     double *airdata_density, double *airdata_sonic_speed,
                     double *airdata_mach, double *airdata_dynamic_pressure);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for airdata_atmos.h
 *
 * [EOF]
 */
