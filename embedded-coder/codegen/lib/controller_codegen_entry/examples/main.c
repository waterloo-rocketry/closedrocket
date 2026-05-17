/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 16-May-2026 22:44:24
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
#include "controller_codegen_entry.h"
#include "controller_codegen_entry_initialize.h"
#include "controller_codegen_entry_terminate.h"
#include "controller_codegen_entry_types.h"
#include "navigation_codegen_entry.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void argInit_2x1_real_T(double result[2]);

static void argInit_3x1_real_T(double result[3]);

static boolean_T argInit_boolean_T(void);

static double argInit_real_T(void);

static struct0_T argInit_struct0_T(void);

static struct1_T argInit_struct1_T(void);

/* Function Definitions */
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
 * Return Type  : boolean_T
 */
static boolean_T argInit_boolean_T(void)
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
 * Arguments    : void
 * Return Type  : struct0_T
 */
static struct0_T argInit_struct0_T(void)
{
  struct0_T result;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_3x1_real_T(result.meas);
  result.status = argInit_boolean_T();
  return result;
}

/*
 * Arguments    : void
 * Return Type  : struct1_T
 */
static struct1_T argInit_struct1_T(void)
{
  struct1_T result;
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
  controller_codegen_entry_initialize();
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_controller_codegen_entry();
  main_navigation_codegen_entry();
  /* Terminate the application.
You do not need to do this more than one time. */
  controller_codegen_entry_terminate();
  return 0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_controller_codegen_entry(void)
{
  double dv[2];
  double C_l_delta;
  double b_r;
  double time_tmp;
  double u;
  /* Initialize function 'controller_codegen_entry' input arguments. */
  time_tmp = argInit_real_T();
  /* Initialize function input argument 'xR'. */
  /* Call the entry-point 'controller_codegen_entry'. */
  argInit_2x1_real_T(dv);
  controller_codegen_entry(time_tmp, dv, time_tmp, time_tmp, &u, &b_r,
                           &C_l_delta);
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_navigation_codegen_entry(void)
{
  struct0_T board_accel_tmp;
  struct1_T b_r;
  struct1_T r1;
  struct2_T state;
  struct3_T airdata;
  double roll_state[2];
  double cov_norm;
  /* Initialize function 'navigation_codegen_entry' input arguments. */
  /* Initialize function input argument 'board_accel'. */
  board_accel_tmp = argInit_struct0_T();
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
  b_r = argInit_struct1_T();
  r1 = argInit_struct1_T();
  navigation_codegen_entry(argInit_real_T(), argInit_boolean_T(),
                           &board_accel_tmp, &board_accel_tmp, &board_accel_tmp,
                           &board_accel_tmp, &board_accel_tmp, &board_accel_tmp,
                           &b_r, &board_accel_tmp, &r1, &board_accel_tmp,
                           &state, &cov_norm, &airdata, roll_state);
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
