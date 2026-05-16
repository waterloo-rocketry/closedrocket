/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * controller_codegen_entry.h
 *
 * Code generation for function 'controller_codegen_entry'
 *
 */

#ifndef CONTROLLER_CODEGEN_ENTRY_H
#define CONTROLLER_CODEGEN_ENTRY_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void controller_codegen_entry(double b_time, const double xR[2],
                                     double pdyn, double delta, double *u,
                                     double *r, double *C_l_delta);

void controller_codegen_entry_init(void);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (controller_codegen_entry.h) */
