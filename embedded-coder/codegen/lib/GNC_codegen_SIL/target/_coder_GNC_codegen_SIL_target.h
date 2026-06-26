#ifndef _CODER_GNC_CODEGEN_SIL_TARGET_H
#define _CODER_GNC_CODEGEN_SIL_TARGET_H

#include "GNC_codegen_SIL_types.h"
#include "rtwtypes.h"
#include "xil_interface_common.h"
#include "xil_target_interface.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void XILTarget_initialize(unsigned int fcnId);

extern void XILTarget_terminate(unsigned int fcnId);

extern XIL_PROCESSDATA_ERROR_CODE
xilTarget_controller_codegen_entry(unsigned int fcnId);

extern XIL_PROCESSDATA_ERROR_CODE
xilTarget_navigation_codegen_entry(unsigned int fcnId);

#ifdef __cplusplus
}
#endif

#endif
