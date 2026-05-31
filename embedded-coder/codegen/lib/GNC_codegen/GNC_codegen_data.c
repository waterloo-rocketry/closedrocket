/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: GNC_codegen_data.c
 *
 * MATLAB Coder version            : 25.2
 * C/C++ source code generated on  : 31-May-2026 15:50:44
 */

/* Include Files */
#include "GNC_codegen_data.h"
#include "GNC_codegen_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
struct_T param;

struct_T b_param;

const struct_T r = {
    {0.46, 0.0, 0.0, 0.0, 49.5, 0.0, 0.0, 0.0, 49.5}, /* J */
    {2.1739130434782608, 0.0, 0.0, 0.0, 0.020202020202020204, 0.0, 0.0, 0.0,
     0.020202020202020204}, /* Jinv */
    {-9.81, 0.0, 0.0}       /* g */
};

bool isInitialized_GNC_codegen = false;

/*
 * File trailer for GNC_codegen_data.c
 *
 * [EOF]
 */
