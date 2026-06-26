/*
 * File: xil_host_interface.h
 *
 * Abstract: Host interface definitions.
 */

#ifndef XIL_HOST_INTERFACE_H
#define XIL_HOST_INTERFACE_H
#include "tmwtypes.h"
#include "coder/connectivity/memunit/MemUnitTransformerType.h"
#ifdef __cplusplus

extern "C"
{

#endif

  boolean_T xilReadData(uint8_T *dst, size_t numElements, mem_unit_type_T
                        dataTypeId);
  boolean_T xilWriteData(const uint8_T *src, size_t numElements, mem_unit_type_T
    dataTypeId);
  void xilPreEntryPointHost(uint32_T taskID);
  void xilEntryPointHost(uint32_T taskID);
  void xilPostEntryPointHost(uint32_T taskID);
  void xilInitializeHost(void (*cleanUp)(void));

#ifdef __cplusplus

}

#endif
#endif                                 /* XIL_HOST_INTERFACE_H */
