/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * controller_codegen_entry_initialize.c
 *
 * Code generation for function 'controller_codegen_entry_initialize'
 *
 */

/* Include files */
#include "controller_codegen_entry_initialize.h"
#include "controller_codegen_entry.h"
#include "controller_codegen_entry_data.h"
#include "controller_estimator.h"

/* Function Definitions */
void controller_codegen_entry_initialize(void)
{
  controller_codegen_entry_init();
  controller_estimator_init();
  isInitialized_controller_codegen_entry = true;
}

/* End of code generation (controller_codegen_entry_initialize.c) */
