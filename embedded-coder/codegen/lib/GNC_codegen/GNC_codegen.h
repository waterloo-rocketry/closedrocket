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

extern void controller_codegen_entry(GNC_codegenStackData *SD, real_T b_time,
                                     real_T dt_ctrl, const real_T xR[2],
                                     real_T pdyn, real_T delta,
                                     const struct0_T *ctrl_mem_in, real_T *u,
                                     real_T *r, struct0_T *ctrl_mem_out);

extern void navigation_codegen_entry(
    GNC_codegenStackData *SD, real_T dt, boolean_T flight_phase,
    const real_T x[11], const real_T P[121], const struct1_T *bias,
    const struct2_T *sens_filt, const struct3_T *sens_input, real_T x_ret[11],
    real_T P_ret[121], struct1_T *bias_ret, struct2_T *sens_filt_ret,
    real_T *cov_norm, struct6_T *airdata, real_T roll_state[2]);

#ifdef __cplusplus
}
#endif

#endif
