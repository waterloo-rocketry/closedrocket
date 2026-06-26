/*
 * File: coder_assumptions.h
 *
 * Abstract: Coder assumptions header file
 */

#ifndef CODER_ASSUMPTIONS_H
#define CODER_ASSUMPTIONS_H

/* include model specific checks */
#include "GNC_codegen_SIL_ca.h"

/* global results variable mapping for static code */
#define CA_Expected                    CA_GNC_codegen_SIL_Exp
#define CA_Actual                      CA_GNC_codegen_SIL_Act
#define CA_TestResults                 CA_GNC_codegen_SIL_Res
#define CA_PortableWordSizes_Results   CA_GNC_codegen_SIL_PWSRes

/* entry point function mapping for static code */
#define CA_Run_Tests                   GNC_codegen_SIL_caRunTests
#define CA_MODEL_SPECIFIC_CHECKS_ENABLED 0
#endif                                 /* CODER_ASSUMPTIONS_H */
