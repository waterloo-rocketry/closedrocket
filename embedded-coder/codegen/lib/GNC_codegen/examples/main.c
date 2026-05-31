/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 15:50:44
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include Files */
#include "main.h"
#include "GNC_codegen_initialize.h"
#include "GNC_codegen_terminate.h"
#include "GNC_codegen_types.h"
#include "controller_codegen_entry.h"
#include "navigation_codegen_entry.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void argInit_11x11_real_T(double result[121]);

static void argInit_11x1_real_T(double result[11]);

static void argInit_2x1_real_T(double result[2]);

static void argInit_2x2_real_T(double result[4]);

static void argInit_3x1_real_T(double result[3]);

static bool argInit_boolean_T(void);

static double argInit_real_T(void);

static void argInit_struct0_T(struct0_T *result);

static void argInit_struct1_T(struct1_T *result);

static struct2_T argInit_struct2_T(void);

static struct3_T argInit_struct3_T(void);

static struct4_T argInit_struct4_T(void);

/* Function Definitions */
/*
 * Arguments    : double result[121]
 * Return Type  : void
 */
static void argInit_11x11_real_T(double result[121])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 121; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[11]
 * Return Type  : void
 */
static void argInit_11x1_real_T(double result[11])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 11; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[2]
 * Return Type  : void
 */
static void argInit_2x1_real_T(double result[2])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[4]
 * Return Type  : void
 */
static void argInit_2x2_real_T(double result[4])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 4; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : bool
 */
static bool argInit_boolean_T(void)
{
  return false;
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : struct0_T *result
 * Return Type  : void
 */
static void argInit_struct0_T(struct0_T *result)
{
  double result_tmp;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_3x1_real_T(result->board_gyro);
  result_tmp = argInit_real_T();
  result->board_baro = result_tmp;
  result->mti_baro = result_tmp;
  result->mti_gyro[0] = result->board_gyro[0];
  result->ad_gyro[0] = result->board_gyro[0];
  result->board_mag_earth[0] = result->board_gyro[0];
  result->mti_mag_earth[0] = result->board_gyro[0];
  result->mti_gyro[1] = result->board_gyro[1];
  result->ad_gyro[1] = result->board_gyro[1];
  result->board_mag_earth[1] = result->board_gyro[1];
  result->mti_mag_earth[1] = result->board_gyro[1];
  result->mti_gyro[2] = result->board_gyro[2];
  result->ad_gyro[2] = result->board_gyro[2];
  result->board_mag_earth[2] = result->board_gyro[2];
  result->mti_mag_earth[2] = result->board_gyro[2];
}

/*
 * Arguments    : struct1_T *result
 * Return Type  : void
 */
static void argInit_struct1_T(struct1_T *result)
{
  double result_tmp;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_3x1_real_T(result->board_accel_f);
  result_tmp = argInit_real_T();
  result->board_baro_f = result_tmp;
  result->mti_baro_f = result_tmp;
  result->board_gyro_f[0] = result->board_accel_f[0];
  result->mti_accel_f[0] = result->board_accel_f[0];
  result->mti_gyro_f[0] = result->board_accel_f[0];
  result->ad_accel_f[0] = result->board_accel_f[0];
  result->ad_gyro_f[0] = result->board_accel_f[0];
  result->board_mag_f[0] = result->board_accel_f[0];
  result->mti_mag_f[0] = result->board_accel_f[0];
  result->board_gyro_f[1] = result->board_accel_f[1];
  result->mti_accel_f[1] = result->board_accel_f[1];
  result->mti_gyro_f[1] = result->board_accel_f[1];
  result->ad_accel_f[1] = result->board_accel_f[1];
  result->ad_gyro_f[1] = result->board_accel_f[1];
  result->board_mag_f[1] = result->board_accel_f[1];
  result->mti_mag_f[1] = result->board_accel_f[1];
  result->board_gyro_f[2] = result->board_accel_f[2];
  result->mti_accel_f[2] = result->board_accel_f[2];
  result->mti_gyro_f[2] = result->board_accel_f[2];
  result->ad_accel_f[2] = result->board_accel_f[2];
  result->ad_gyro_f[2] = result->board_accel_f[2];
  result->board_mag_f[2] = result->board_accel_f[2];
  result->mti_mag_f[2] = result->board_accel_f[2];
}

/*
 * Arguments    : void
 * Return Type  : struct2_T
 */
static struct2_T argInit_struct2_T(void)
{
  struct2_T result;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_3x1_real_T(result.meas);
  result.status = argInit_boolean_T();
  return result;
}

/*
 * Arguments    : void
 * Return Type  : struct3_T
 */
static struct3_T argInit_struct3_T(void)
{
  struct3_T result;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_3x1_real_T(result.meas);
  result.status = argInit_real_T();
  return result;
}

/*
 * Arguments    : void
 * Return Type  : struct4_T
 */
static struct4_T argInit_struct4_T(void)
{
  struct4_T result;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  result.meas = argInit_real_T();
  result.status = argInit_boolean_T();
  return result;
}

/*
 * Arguments    : int argc
 *                char **argv
 * Return Type  : int
 */
int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* Initialize the application.
You do not need to do this more than one time. */
  GNC_codegen_initialize();
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_controller_codegen_entry();
  main_navigation_codegen_entry();
  /* Terminate the application.
You do not need to do this more than one time. */
  GNC_codegen_terminate();
  return 0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_controller_codegen_entry(void)
{
  double P_minus_ret[4];
  double dv[4];
  double coeffs_ret[2];
  double xR_tmp[2];
  double b_r;
  double d_old_ret;
  double time_tmp;
  double u;
  double w_dot_old_ret;
  double w_old_ret;
  /* Initialize function 'controller_codegen_entry' input arguments. */
  time_tmp = argInit_real_T();
  /* Initialize function input argument 'xR'. */
  argInit_2x1_real_T(xR_tmp);
  /* Initialize function input argument 'coeffs'. */
  /* Initialize function input argument 'P_minus'. */
  /* Call the entry-point 'controller_codegen_entry'. */
  argInit_2x2_real_T(dv);
  controller_codegen_entry(time_tmp, time_tmp, xR_tmp, time_tmp, time_tmp,
                           time_tmp, xR_tmp, dv, time_tmp, time_tmp, &u, &b_r,
                           coeffs_ret, &w_old_ret, P_minus_ret, &d_old_ret,
                           &w_dot_old_ret);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_navigation_codegen_entry(void)
{
  struct0_T b_r;
  struct0_T b_ret;
  struct1_T r1;
  struct1_T sf_ret;
  struct2_T board_accel_tmp;
  struct3_T r2;
  struct4_T r3;
  struct4_T r4;
  double P_ret[121];
  double dv1[121];
  double dv[11];
  double x_ret[11];
  /* Initialize function 'navigation_codegen_entry' input arguments. */
  /* Initialize function input argument 'x'. */
  /* Initialize function input argument 'P'. */
  /* Initialize function input argument 'b'. */
  /* Initialize function input argument 'sf'. */
  /* Initialize function input argument 'board_accel'. */
  board_accel_tmp = argInit_struct2_T();
  /* Initialize function input argument 'board_gyro'. */
  /* Initialize function input argument 'mti_accel'. */
  /* Initialize function input argument 'mti_gyro'. */
  /* Initialize function input argument 'ad_accel'. */
  /* Initialize function input argument 'ad_gyro'. */
  /* Initialize function input argument 'board_baro'. */
  /* Initialize function input argument 'board_mag'. */
  /* Initialize function input argument 'mti_baro'. */
  /* Initialize function input argument 'mti_mag'. */
  /* Call the entry-point 'navigation_codegen_entry'. */
  argInit_11x1_real_T(dv);
  argInit_11x11_real_T(dv1);
  argInit_struct0_T(&b_r);
  argInit_struct1_T(&r1);
  r2 = argInit_struct3_T();
  r3 = argInit_struct4_T();
  r4 = argInit_struct4_T();
  navigation_codegen_entry(argInit_real_T(), argInit_boolean_T(), dv, dv1, &b_r,
                           &r1, &board_accel_tmp, &board_accel_tmp,
                           &board_accel_tmp, &board_accel_tmp, &board_accel_tmp,
                           &r2, &r3, &board_accel_tmp, &r4, &board_accel_tmp,
                           x_ret, P_ret, &b_ret, &sf_ret);
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
