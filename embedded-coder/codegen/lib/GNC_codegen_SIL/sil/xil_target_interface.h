/*
 * File: xil_target_interface.h
 *
 * Abstract: Target interface definitions.
 */

#ifndef XIL_TARGET_INTERFACE_H
#define XIL_TARGET_INTERFACE_H
#include "rtwtypes.h"
#include "xil_interface_common.h"
#include "xil_data_stream.h"

void xilPreEntryPoint(uint32_T entryPointID);
XIL_INTERFACE_ERROR_CODE xilPostEntryPoint(uint32_T entryPointID);
XIL_INTERFACE_ERROR_CODE xilPostTerminatePoint(uint32_T entryPointID);

#endif                                 /* XIL_TARGET_INTERFACE_H */
