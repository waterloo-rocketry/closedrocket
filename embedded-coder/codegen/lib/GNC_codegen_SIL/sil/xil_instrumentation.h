/*
 * File: xil_instrumentation.h
 *
 * Code generated for instrumentation.
 *
 */

/* Functions with a C call interface */
#ifdef __cplusplus

extern "C"
{

#endif

#include "host_timer_x86.h"
#ifdef __cplusplus

}

#endif

#include "rtwtypes.h"

/* Upload code instrumentation data point */
void xilUploadCodeInstrData(
  void* pData, uint32_T numMemUnits, uint32_T sectionId);

/* Uploads data */
void xilUploadProfilingData(uint32_T sectionId);

/* Pause/restart the timer while running code associated with storing and uploading the data. */
void xilProfilingTimerFreeze(void);
void xilProfilingTimerUnFreeze(void);

/* Code instrumentation method(s) for model GNC_codegen_SIL */
void profileStart_GNC_codegen_SIL(uint32_T sectionId);
void profileEnd_GNC_codegen_SIL(uint32_T sectionId);

/* Code instrumentation method(s) for model GNC_codegen_SIL */
void taskTimeStart_GNC_codegen_SIL(uint32_T sectionId);
void taskTimeEnd_GNC_codegen_SIL(uint32_T sectionId);
