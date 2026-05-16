/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * controller_estimator.c
 *
 * Code generation for function 'controller_estimator'
 *
 */

/* Include files */
#include "controller_estimator.h"
#include "controller_codegen_entry_data.h"
#include <emmintrin.h>

/* Function Definitions */
void controller_estimator_init(void)
{
  __m128d r;
  double dv[2];
  w_old_not_empty = false;
  P_minus_not_empty = false;
  t = -0.01;
  /*  for /(time - t) */
  dv[0] = 0.0;
  dv[1] = 1.0;
  r = _mm_loadu_pd(&dv[0]);
  _mm_storeu_pd(&c[0],
                _mm_add_pd(_mm_set1_pd(2.0), _mm_mul_pd(_mm_set1_pd(-2.0), r)));
  /*  initial coefficient guess */
  d_old = 0.0;
  w_dot_old = 0.0;
}

/* End of code generation (controller_estimator.c) */
