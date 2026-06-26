#ifndef GNC_CODEGEN_SIL_H
#define GNC_CODEGEN_SIL_H

#include "GNC_codegen_SIL_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void GNC_codegen_SIL_initialize(GNC_codegen_SILStackData *b_SD);

extern void GNC_codegen_SIL_terminate(void);

extern void controller_codegen_entry(GNC_codegen_SILStackData *b_SD,
                                     double b_time, double dt_ctrl,
                                     const double xR[2], double pdyn,
                                     double delta_encoder, struct0_T *ctrl_mem,
                                     double *u_motor, double *r,
                                     boolean_T *w_status_ctrl);

extern void navigation_codegen_entry(GNC_codegen_SILStackData *b_SD, double dt,
                                     boolean_T flight_phase, double x[11],
                                     double P[121], struct1_T *bias,
                                     struct2_T *sens_filt,
                                     const struct3_T *sens_in, double *cov_norm,
                                     double roll_state[2], double *pdyn,
                                     boolean_T *w_status_nav);

#ifdef __cplusplus
}
#endif

#endif
