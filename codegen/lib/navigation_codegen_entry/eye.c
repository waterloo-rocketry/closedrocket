/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eye.c
 *
 * Code generation for function 'eye'
 *
 */

/* Include files */
#include "eye.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void b_eye(double b_I[16])
{
  memset(&b_I[0], 0, 16U * sizeof(double));
  b_I[0] = 1.0;
  b_I[5] = 1.0;
  b_I[10] = 1.0;
  b_I[15] = 1.0;
}

void c_eye(double b_I[121])
{
  int k;
  memset(&b_I[0], 0, 121U * sizeof(double));
  for (k = 0; k < 11; k++) {
    b_I[k + 11 * k] = 1.0;
  }
}

void eye(double b_I[9])
{
  memset(&b_I[0], 0, 9U * sizeof(double));
  b_I[0] = 1.0;
  b_I[4] = 1.0;
  b_I[8] = 1.0;
}

/* End of code generation (eye.c) */
