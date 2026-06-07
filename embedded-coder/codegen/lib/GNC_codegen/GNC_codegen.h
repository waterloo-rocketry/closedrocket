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
                                     double dt_ctrl, const double xR[2],
                                     double pdyn, double delta,
                                     const struct0_T *ctrl_mem_in, double *u,
                                     double *r, struct0_T *ctrl_mem_out);

extern void navigation_codegen_entry(
    GNC_codegenStackData *SD, double dt, bool flight_phase, const double x[11],
    const double P[121], const struct1_T *bias, const struct2_T *sens_filt,
    const struct3_T *sens_input, double x_ret[11], double P_ret[121],
    struct1_T *bias_ret, struct2_T *sens_filt_ret, struct6_T *airdata,
    double roll_state[2]);

#ifdef __cplusplus
}
#endif

#endif
