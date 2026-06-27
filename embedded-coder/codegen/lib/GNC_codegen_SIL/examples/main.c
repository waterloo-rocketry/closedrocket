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
  result->w = result_tmp;
  argInit_2x2_real_T(result->P);
  result->delta_lp = result_tmp;
  result->w_dot_lp = result_tmp;
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

  argInit_3x1_real_T(result->board_accel);
  result_tmp = argInit_real_T();
  result->board_baro = result_tmp;
  result->mti_baro = result_tmp;
  result->board_gyro[0] = result->board_accel[0];
  result->mti_accel[0] = result->board_accel[0];
  result->mti_gyro[0] = result->board_accel[0];
  result->ad_accel[0] = result->board_accel[0];
  result->ad_gyro[0] = result->board_accel[0];
  result->board_mag[0] = result->board_accel[0];
  result->mti_mag[0] = result->board_accel[0];
  result->board_gyro[1] = result->board_accel[1];
  result->mti_accel[1] = result->board_accel[1];
  result->mti_gyro[1] = result->board_accel[1];
  result->ad_accel[1] = result->board_accel[1];
  result->ad_gyro[1] = result->board_accel[1];
  result->board_mag[1] = result->board_accel[1];
  result->mti_mag[1] = result->board_accel[1];
  result->board_gyro[2] = result->board_accel[2];
  result->mti_accel[2] = result->board_accel[2];
  result->mti_gyro[2] = result->board_accel[2];
  result->ad_accel[2] = result->board_accel[2];
  result->ad_gyro[2] = result->board_accel[2];
  result->board_mag[2] = result->board_accel[2];
  result->mti_mag[2] = result->board_accel[2];
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
  struct0_T ctrl_mem;
  double b_dv[2];
  double where_it_isnt[2];
  double time_tmp;
  double u_motor;
  boolean_T w_status_ctrl;

  time_tmp = argInit_real_T();

  argInit_struct0_T(&ctrl_mem);
  argInit_2x1_real_T(b_dv);
  controller_codegen_entry(&GNC_codegen_SILStackDataGlobal, time_tmp, time_tmp,
                           b_dv, time_tmp, time_tmp, &ctrl_mem, &u_motor,
                           where_it_isnt, &w_status_ctrl);
}

void main_navigation_codegen_entry(void) {
  struct1_T bias;
  struct2_T sens_filt;
  struct3_T r;
  double P[121];
  double x[11];
  double roll_state[2];
  double cov_norm;
  double pdyn;
  boolean_T w_status_nav;

  argInit_11x1_real_T(x);
  argInit_11x11_real_T(P);
  argInit_struct1_T(&bias);
  argInit_struct2_T(&sens_filt);
  argInit_struct3_T(&r);
  navigation_codegen_entry(&GNC_codegen_SILStackDataGlobal, argInit_real_T(),
                           argInit_boolean_T(), x, P, &bias, &sens_filt, &r,
                           &cov_norm, roll_state, &pdyn, &w_status_nav);
}
