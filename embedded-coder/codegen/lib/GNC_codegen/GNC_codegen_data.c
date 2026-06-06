/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: GNC_codegen_data.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 05-Jun-2026 20:31:45
 */

/* Include Files */
#include "GNC_codegen_data.h"
#include "GNC_codegen_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
struct_T b_param;

struct_T c_param;

struct_T d_param;

const struct_T r = {
    10.0,                                             /* Cn_alpha */
    {0.46, 0.0, 0.0, 0.0, 49.5, 0.0, 0.0, 0.0, 49.5}, /* J */
    {2.1739130434782608, 0.0, 0.0, 0.0, 0.020202020202020204, 0.0, 0.0, 0.0,
     0.020202020202020204}, /* Jinv */
    -0.016182736457722724,  /* c_aero */
    0.00061367415999999994, /* c_canard */
    420.0,                  /* elevation */
    {-9.81, 0.0, 0.0}       /* g */
};

bool isInitialized_GNC_codegen = false;

/*
 * File trailer for GNC_codegen_data.c
 *
 * [EOF]
 */
