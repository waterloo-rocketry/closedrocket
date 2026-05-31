/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: airdata_atmos.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 15:50:44
 */

/* Include Files */
#include "airdata_atmos.h"
#include "mpower.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Definitions */
/*
 * computes atmospheric stanadard data from altitude, according to US standard
 * atmosphere atmosphere data: static pressure, temperature, density, local
 * speed of sound calculations found in Stengel 2004, pp. 30
 *
 * Arguments    : double altitude
 *                double *airdata_temperature
 *                double *airdata_density
 *                double *airdata_sonic_speed
 *                double *airdata_mach
 *                double *airdata_dynamic_pressure
 * Return Type  : double
 */
double airdata_atmos(double altitude, double *airdata_temperature,
                     double *airdata_density, double *airdata_sonic_speed,
                     double *airdata_mach, double *airdata_dynamic_pressure)
{
  double airdata_pressure;
  double layer_idx_1;
  double layer_idx_2;
  double layer_idx_3;
  int layer_idx_0;
  /*  parameters */
  /*  adiabatic index */
  /*  specific gas constant for air */
  /*  troposphere */
  /*  tropopause */
  /*  stratosphere */
  /*  stratosphere 2 */
  /*  base height, P_base, T_base, lapse rate; */
  /*  mean earth radius */
  /*  zero height gravity */
  /*  geopotential altitude */
  altitude = 6.356766E+6 * altitude / (6.356766E+6 - altitude);
  /*  select atmosphere behaviour from table */
  layer_idx_0 = 0;
  layer_idx_1 = 101325.0;
  layer_idx_2 = 288.15;
  layer_idx_3 = 0.0065;
  if (altitude > 11000.0) {
    if (altitude < 20000.0) {
      layer_idx_0 = 11000;
      layer_idx_1 = 22632.1;
      layer_idx_2 = 216.65;
      layer_idx_3 = 0.0;
    } else if (altitude < 32000.0) {
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
  }
  /*  base height */
  /*  base pressure */
  /*  base temperature */
  /*  temperature lapse rate */
  if (layer_idx_3 == 0.0) {
    airdata_pressure =
        layer_idx_1 * exp(-9.81 * (altitude - (double)layer_idx_0) /
                          (287.0579 * layer_idx_2));
  } else {
    airdata_pressure =
        layer_idx_1 * mpower(1.0 - layer_idx_3 / layer_idx_2 *
                                       (altitude - (double)layer_idx_0),
                             9.81 / (287.0579 * layer_idx_3));
  }
  layer_idx_1 = layer_idx_2 - layer_idx_3 * (altitude - (double)layer_idx_0);
  *airdata_density = airdata_pressure / (287.0579 * layer_idx_1);
  *airdata_sonic_speed = sqrt(401.88106 * layer_idx_1);
  /*  return values */
  /*  [Pa] */
  *airdata_temperature = layer_idx_1;
  /*  [K] */
  /*  [kg/m^3] */
  /*  [m/s] */
  *airdata_mach = 0.0;
  /*  need to initialize struct field before calling airdata_dynamic */
  *airdata_dynamic_pressure = 0.0;
  /*  need to initialize struct field before calling airdata_dynamic */
  return airdata_pressure;
}

/*
 * File trailer for airdata_atmos.c
 *
 * [EOF]
 */
