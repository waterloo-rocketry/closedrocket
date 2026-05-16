/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * airdata_atmos.h
 *
 * Code generation for function 'airdata_atmos'
 *
 */

#ifndef AIRDATA_ATMOS_H
#define AIRDATA_ATMOS_H

/* Include files */
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
/* End of code generation (airdata_atmos.h) */
