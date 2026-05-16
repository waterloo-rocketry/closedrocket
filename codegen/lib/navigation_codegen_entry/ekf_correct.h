/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * ekf_correct.h
 *
 * Code generation for function 'ekf_correct'
 *
 */

#ifndef EKF_CORRECT_H
#define EKF_CORRECT_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_ekf_correct(const double x[11], const double b_P[121], double y,
                   double b, double x_new[11], double P_new[121]);

void ekf_correct(const double x[11], const double b_P[121], const double y[3],
                 const double b[3], const double R[9], double x_new[11],
                 double P_new[121]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (ekf_correct.h) */
