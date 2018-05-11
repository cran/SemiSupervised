#include "general.h"
#include "joint_harmonic.h"
#include "LAE.h"
#include "Anchor_Reg.h"
#include "s4pm.h"

#include <R_ext/Rdynload.h>
static const R_CallMethodDef CEntries[] = {
  {"LAE", (DL_FUNC) &LAE, 6},
  {"lgraph", (DL_FUNC) &lgraph, 4},
  {"AREG_CV", (DL_FUNC) &AREG_CV, 12},
  {"ARIDGE", (DL_FUNC) &ARIDGE, 13},
  {"cv_jtharm_fit", (DL_FUNC) &cv_jtharm_fit, 8},
  {"jt_harm_fit", (DL_FUNC) &jt_harm_fit, 8},
  {"cv_s4pm_fit", (DL_FUNC) &cv_s4pm_fit, 10},
  {"s4pm_fit", (DL_FUNC) &s4pm_fit, 11},
  {NULL, NULL, 0}
};

#include <Rversion.h>
void R_init_SemiSupervised(DllInfo*);
