#ifndef GNC_CODEGEN_H
#define GNC_CODEGEN_H

#include "GNC_codegen_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

extern void GNC_codegen_initialize(GNC_codegenStackData *SD);

extern void GNC_codegen_terminate(void);

extern void controller_codegen_entry(GNC_codegenStackData *SD, double b_time,
                                     double dt_ctrl,
                                     const double where_it_is[2], double pdyn,
                                     double delta_encoder, struct0_T *ctrl_mem,
                                     double *u_motor, double where_it_isnt[2],
                                     bool *w_status_ctrl);

extern void navigation_codegen_entry(GNC_codegenStackData *SD, double dt,
                                     bool flight_phase, double x[11],
                                     double P[121], struct1_T *bias,
                                     struct2_T *sens_filt,
                                     const struct3_T *sens_in, double *cov_norm,
                                     double roll_state[2], double *pdyn,
                                     bool *w_status_nav);

#ifdef __cplusplus
}
#endif

#endif
