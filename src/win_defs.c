#include "win_defs.h"
void R_init_SemiSupervised(DllInfo *dll){
  R_registerRoutines(dll, NULL, CEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
