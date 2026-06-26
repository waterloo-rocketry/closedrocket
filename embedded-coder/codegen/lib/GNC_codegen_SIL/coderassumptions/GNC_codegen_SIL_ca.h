/*
 * File: GNC_codegen_SIL_ca.h
 *
 * Abstract: Tests assumptions in the generated code.
 */

#ifndef GNC_CODEGEN_SIL_CA_H
#define GNC_CODEGEN_SIL_CA_H

/* preprocessor validation checks */
#include "GNC_codegen_SIL_ca_preproc.h"
#include "coder_assumptions_hwimpl.h"

/* variables holding test results */
extern CA_ChecksTestResults CA_GNC_codegen_SIL_Res;
extern CA_PWS_TestResults CA_GNC_codegen_SIL_PWSRes;

/* variables holding "expected" and "actual" hardware implementation */
extern const CA_Checks CA_GNC_codegen_SIL_Exp;
extern CA_Checks CA_GNC_codegen_SIL_Act;
extern const int numberOfImportedTypes;

/* entry point function to run tests */
void GNC_codegen_SIL_caRunTests(void);

#endif                                 /* GNC_CODEGEN_SIL_CA_H */
