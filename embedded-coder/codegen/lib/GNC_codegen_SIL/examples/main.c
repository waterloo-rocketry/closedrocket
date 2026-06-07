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

#include "main.h"
#include "GNC_codegen_SIL.h"
#include "GNC_codegen_SIL_types.h"

static GNC_codegen_SILStackData GNC_codegen_SILStackDataGlobal;

static void argInit_11x11_real_T(double result[121]);

static void argInit_11x1_real_T(double result[11]);

static void argInit_2x1_real_T(double result[2]);

static void argInit_2x2_real_T(double result[4]);

static void argInit_3x1_real_T(double result[3]);

static boolean_T argInit_boolean_T(void);

static double argInit_real_T(void);

static void argInit_struct0_T(struct0_T *result);

static void argInit_struct1_T(struct1_T *result);

static void argInit_struct2_T(struct2_T *result);

static void argInit_struct3_T(struct3_T *result);

static struct4_T argInit_struct4_T(void);

static struct5_T argInit_struct5_T(void);

static void argInit_11x11_real_T(double result[121]) {
  int i;

  for (i = 0; i < 121; i++) {

    result[i] = argInit_real_T();
  }
}

static void argInit_11x1_real_T(double result[11]) {
  int idx0;

  for (idx0 = 0; idx0 < 11; idx0++) {

    result[idx0] = argInit_real_T();
  }
}

static void argInit_2x1_real_T(double result[2]) {
  int idx0;

  for (idx0 = 0; idx0 < 2; idx0++) {

    result[idx0] = argInit_real_T();
  }
}

static void argInit_2x2_real_T(double result[4]) {
  int i;

  for (i = 0; i < 4; i++) {

    result[i] = argInit_real_T();
  }
}

static void argInit_3x1_real_T(double result[3]) {
  int idx0;

  for (idx0 = 0; idx0 < 3; idx0++) {

    result[idx0] = argInit_real_T();
  }
}

static boolean_T argInit_boolean_T(void) { return false; }

static double argInit_real_T(void) { return 0.0; }

static void argInit_struct0_T(struct0_T *result) {
  double result_tmp;

  result_tmp = argInit_real_T();
  argInit_2x1_real_T(result->coeffs);
  result->w_old = result_tmp;
  argInit_2x2_real_T(result->P_minus);
  result->d_old = result_tmp;
  result->w_dot_old = result_tmp;
}

static void argInit_struct1_T(struct1_T *result) {
  double result_tmp;

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

static void argInit_struct2_T(struct2_T *result) {
  double result_tmp;

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

static void argInit_struct3_T(struct3_T *result) {

  result->board_accel = argInit_struct4_T();
  result->board_baro = argInit_struct5_T();
  result->mti_baro = argInit_struct5_T();
  result->board_gyro = result->board_accel;
  result->mti_accel = result->board_accel;
  result->mti_gyro = result->board_accel;
  result->ad_accel = result->board_accel;
  result->ad_gyro = result->board_accel;
  result->board_mag = result->board_accel;
  result->mti_mag = result->board_accel;
}

static struct4_T argInit_struct4_T(void) {
  struct4_T result;

  argInit_3x1_real_T(result.meas);
  result.status = argInit_boolean_T();
  return result;
}

static struct5_T argInit_struct5_T(void) {
  struct5_T result;

  result.meas = argInit_real_T();
  result.status = argInit_boolean_T();
  return result;
}

int main(int argc, char **argv) {
  static GNC_codegen_SILPersistentData c_GNC_codegen_SILPersistentData;
  (void)argc;
  (void)argv;
  GNC_codegen_SILStackDataGlobal.pd = &c_GNC_codegen_SILPersistentData;

  GNC_codegen_SIL_initialize(&GNC_codegen_SILStackDataGlobal);

  main_controller_codegen_entry();
  main_navigation_codegen_entry();

  GNC_codegen_SIL_terminate();
  return 0;
}

void main_controller_codegen_entry(void) {
  struct0_T ctrl_mem_out;
  struct0_T r;
  double b_dv[2];
  double b_r;
  double time_tmp;
  double u;

  time_tmp = argInit_real_T();

  argInit_2x1_real_T(b_dv);
  argInit_struct0_T(&r);
  controller_codegen_entry(&GNC_codegen_SILStackDataGlobal, time_tmp, time_tmp,
                           b_dv, time_tmp, time_tmp, &r, &u, &b_r,
                           &ctrl_mem_out);
}

void main_navigation_codegen_entry(void) {
  struct1_T bias_ret;
  struct1_T r;
  struct2_T r1;
  struct2_T sens_filt_ret;
  struct3_T r2;
  struct6_T airdata;
  double P_ret[121];
  double b_dv1[121];
  double b_dv[11];
  double x_ret[11];
  double roll_state[2];

  argInit_11x1_real_T(b_dv);
  argInit_11x11_real_T(b_dv1);
  argInit_struct1_T(&r);
  argInit_struct2_T(&r1);
  argInit_struct3_T(&r2);
  navigation_codegen_entry(&GNC_codegen_SILStackDataGlobal, argInit_real_T(),
                           argInit_boolean_T(), b_dv, b_dv1, &r, &r1, &r2,
                           x_ret, P_ret, &bias_ret, &sens_filt_ret, &airdata,
                           roll_state);
}
