/*
 * File: xil_interface.c
 *
 * SIL generated interface for code: "GNC_codegen_SIL"
 *
 */

#include "GNC_codegen_SIL.h"
#include "xil_interface.h"
#include "xil_instrumentation.h"

/* Functions with a C call interface */
#ifdef __cplusplus

extern "C"
{

#endif

#include "xil_data_stream.h"
#include "codeinstr_data_stream.h"
#ifdef __cplusplus

}

#endif

#include "_coder_GNC_codegen_SIL_target.h"
#include "xil_target_interface.h"

/* In-the-Loop Interface functions - see xil_interface.h */
XIL_INTERFACE_ERROR_CODE xilProcessParams(uint32_T xilFcnId)
{
  /* Single In-the-Loop Component */
  if (xilFcnId != 0) {
    return XIL_INTERFACE_UNKNOWN_FCNID;
  }

  return XIL_INTERFACE_SUCCESS;
}

void xilUploadCodeInstrData(void * pData, uint32_T numMemUnits, uint32_T
  sectionId)
{
  /* Send code instrumentation data to host */
  if (codeInstrWriteData((MemUnit_T *) &numMemUnits, sizeof(numMemUnits)) !=
      XIL_DATA_STREAM_SUCCESS) {
    for (;;) ;
  }

  if (codeInstrWriteData((MemUnit_T *) &sectionId, sizeof(uint32_T)) !=
      XIL_DATA_STREAM_SUCCESS) {
    for (;;) ;
  }

  if (codeInstrWriteData((MemUnit_T *) pData, numMemUnits) !=
      XIL_DATA_STREAM_SUCCESS) {
    for (;;) ;
  }
}

XIL_INTERFACE_ERROR_CODE xilGetDataTypeInfo(void)
{
  {
    /* send response id code */
    MemUnit_T memUnitData = XIL_RESPONSE_TYPE_SIZE;
    if (xilWriteData(&memUnitData, sizeof(memUnitData)) !=
        XIL_DATA_STREAM_SUCCESS) {
      return XIL_INTERFACE_COMMS_FAILURE;
    }

    /* send type id */
    memUnitData = 6;
    if (xilWriteData(&memUnitData, sizeof(memUnitData)) !=
        XIL_DATA_STREAM_SUCCESS) {
      return XIL_INTERFACE_COMMS_FAILURE;
    }

    /* PIL_DOUBLE_SIZE should only be already defined for MathWorks testing */
#ifndef PIL_DOUBLE_SIZE
#define PIL_DOUBLE_SIZE                sizeof(double)
#endif

    /* send size in bytes */
    memUnitData = (MemUnit_T) PIL_DOUBLE_SIZE;

#ifndef HOST_WORD_ADDRESSABLE_TESTING

    /* convert MemUnits to bytes */
    memUnitData *= MEM_UNIT_BYTES;

#endif

    if (xilWriteData(&memUnitData, sizeof(memUnitData)) !=
        XIL_DATA_STREAM_SUCCESS) {
      return XIL_INTERFACE_COMMS_FAILURE;
    }
  }

  return XIL_INTERFACE_SUCCESS;
}

XIL_INTERFACE_ERROR_CODE xilInitialize(uint32_T xilFcnId)
{
  XIL_INTERFACE_ERROR_CODE errorCode = XIL_INTERFACE_SUCCESS;

  /* Invoke Initialization Function */
  if (xilFcnId == 0) {
    XILTarget_initialize(1U);
  } else {
    errorCode = XIL_INTERFACE_UNKNOWN_FCNID;
  }                                    /* if: } */

  return errorCode;
}

XIL_INTERFACE_ERROR_CODE xilPause(uint32_T xilFcnId)
{
  XIL_INTERFACE_ERROR_CODE errorCode = XIL_INTERFACE_SUCCESS;
  if (xilFcnId == 0) {
    /* Nothing to do */
  } else {
    errorCode = XIL_INTERFACE_UNKNOWN_FCNID;
  }                                    /* if */

  return errorCode;
}

XIL_INTERFACE_ERROR_CODE xilSystemInitialize(uint32_T xilFcnId)
{
  XIL_INTERFACE_ERROR_CODE errorCode = XIL_INTERFACE_SUCCESS;

  /* Single In-the-Loop Component */
  if (xilFcnId == 0) {
    /* No Function to Call */
    /* No Function to Call */
  } else {
    errorCode = XIL_INTERFACE_UNKNOWN_FCNID;
  }

  return errorCode;
}

XIL_INTERFACE_ERROR_CODE xilSystemReset(uint32_T xilFcnId)
{
  XIL_INTERFACE_ERROR_CODE errorCode = XIL_INTERFACE_SUCCESS;

  /* Single In-the-Loop Component */
  if (xilFcnId == 0) {
    /* No Function to Call */
  } else {
    errorCode = XIL_INTERFACE_UNKNOWN_FCNID;
  }

  return errorCode;
}

XIL_INTERFACE_ERROR_CODE xilGetHostToTargetData(uint32_T xilFcnId,
  XIL_COMMAND_TYPE_ENUM xilCommandType, uint32_T xilCommandIdx, XILIOData
  ** xilIOData)
{
  XIL_INTERFACE_ERROR_CODE errorCode = XIL_INTERFACE_SUCCESS;
  *xilIOData = 0;

  /* Single In-the-Loop Component */
  if (xilFcnId != 0) {
    errorCode = XIL_INTERFACE_UNKNOWN_FCNID;
    return errorCode;
  }

  errorCode = XIL_INTERFACE_UNKNOWN_TID;
  UNUSED_PARAMETER(xilCommandType);
  UNUSED_PARAMETER(xilCommandIdx);
  return errorCode;
}

XIL_INTERFACE_ERROR_CODE xilGetTargetToHostPreData(uint32_T xilFcnId,
  XIL_COMMAND_TYPE_ENUM xilCommandType, uint32_T xilCommandIdx, XILIOData
  ** xilIOData, MemUnit_T responseId, uint32_T serverFcnId)
{
  XIL_INTERFACE_ERROR_CODE errorCode = XIL_INTERFACE_SUCCESS;
  *xilIOData = 0;
  if (xilFcnId != 0) {
    errorCode = XIL_INTERFACE_UNKNOWN_FCNID;
    return errorCode;
  }                                    /* if */

  errorCode = XIL_INTERFACE_UNKNOWN_TID;
  UNUSED_PARAMETER(xilCommandType);
  UNUSED_PARAMETER(xilCommandIdx);
  UNUSED_PARAMETER(responseId);
  UNUSED_PARAMETER(serverFcnId);
  return errorCode;
}

XIL_INTERFACE_ERROR_CODE xilOutput(uint32_T xilFcnId, uint32_T xilTID)
{
  /* Single In-the-Loop Component */
  if (xilFcnId != 0) {
    return XIL_INTERFACE_UNKNOWN_FCNID;
  }

  switch (xilTID) {
   case 1:
    if (xilTarget_controller_codegen_entry(2U) != XIL_PROCESSDATA_SUCCESS ) {
      return XIL_INTERFACE_COMMS_FAILURE;
    }                                  /* if */
    break;

   case 2:
    if (xilTarget_navigation_codegen_entry(3U) != XIL_PROCESSDATA_SUCCESS ) {
      return XIL_INTERFACE_COMMS_FAILURE;
    }                                  /* if */
    break;

   default:
    return XIL_INTERFACE_UNKNOWN_TID;
  }                                    /* switch */

  return XIL_INTERFACE_SUCCESS;
}

XIL_INTERFACE_ERROR_CODE xilUpdate(uint32_T xilFcnId, uint32_T xilTID)
{
  /* Single In-the-Loop Component */
  if (xilFcnId != 0) {
    return XIL_INTERFACE_UNKNOWN_FCNID;
  }

  /* No Update Function */
  UNUSED_PARAMETER(xilTID);
  return XIL_INTERFACE_SUCCESS;
}

XIL_INTERFACE_ERROR_CODE xilGetTargetToHostData(uint32_T xilFcnId,
  XIL_COMMAND_TYPE_ENUM xilCommandType, uint32_T xilCommandIdx, XILIOData
  ** xilIOData, MemUnit_T responseId, uint32_T serverFcnId)
{
  XIL_INTERFACE_ERROR_CODE errorCode = XIL_INTERFACE_SUCCESS;

  /* Single In-the-Loop Component */
  *xilIOData = 0;
  if (xilFcnId != 0) {
    errorCode = XIL_INTERFACE_UNKNOWN_FCNID;
    return errorCode;
  }

  errorCode = XIL_INTERFACE_UNKNOWN_TID;
  UNUSED_PARAMETER(xilCommandType);
  UNUSED_PARAMETER(xilCommandIdx);
  UNUSED_PARAMETER(responseId);
  UNUSED_PARAMETER(serverFcnId);
  return errorCode;
}

XIL_INTERFACE_ERROR_CODE xilTerminate(uint32_T xilFcnId)
{
  if (xilFcnId != 0) {
    return XIL_INTERFACE_UNKNOWN_FCNID;
  }                                    /* if */

  /* Invoke any terminate Function */
  XILTarget_terminate(4U);
  return XIL_INTERFACE_SUCCESS;
}

XIL_INTERFACE_ERROR_CODE xilInitTargetData(void)
{
  return XIL_INTERFACE_SUCCESS;
}

XIL_INTERFACE_ERROR_CODE xilEnable(uint32_T xilFcnId, uint32_T xilTID)
{
  /* Single In-the-Loop Component */
  if (xilFcnId != 0) {
    return XIL_INTERFACE_UNKNOWN_FCNID;
  }

  UNUSED_PARAMETER(xilTID);

  /* No Enable Function - this function should never be called */
  return XIL_INTERFACE_UNKNOWN_TID;
}

XIL_INTERFACE_ERROR_CODE xilDisable(uint32_T xilFcnId, uint32_T xilTID)
{
  /* Single In-the-Loop Component */
  if (xilFcnId != 0) {
    return XIL_INTERFACE_UNKNOWN_FCNID;
  }

  UNUSED_PARAMETER(xilTID);

  /* No Disable Function - this function should never be called */
  return XIL_INTERFACE_UNKNOWN_TID;
}

void xilPreEntryPoint(uint32_T entryPointID)
{
  taskTimeStart_GNC_codegen_SIL(entryPointID);
}

XIL_INTERFACE_ERROR_CODE xilPostEntryPoint(uint32_T entryPointID)
{
  MemUnit_T responseId = XIL_RESPONSE_OUTPUT_DATA;
  taskTimeEnd_GNC_codegen_SIL(entryPointID);
  if (xilWriteData(&responseId, sizeof(responseId)) != XIL_DATA_STREAM_SUCCESS)
  {
    return XIL_INTERFACE_COMMS_FAILURE;
  }                                    /* if */

  return XIL_INTERFACE_SUCCESS;
}

XIL_INTERFACE_ERROR_CODE xilPostTerminatePoint(uint32_T entryPointID)
{
  taskTimeEnd_GNC_codegen_SIL(entryPointID);
  return XIL_INTERFACE_SUCCESS;
}
