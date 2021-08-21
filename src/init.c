#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(kde_dir_vmf)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lcv_dir_vmf)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(lscv_dir_vmf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"kde_dir_vmf",  (DL_FUNC) &F77_NAME(kde_dir_vmf),   8},
    {"lcv_dir_vmf",  (DL_FUNC) &F77_NAME(lcv_dir_vmf),   7},
    {"lscv_dir_vmf", (DL_FUNC) &F77_NAME(lscv_dir_vmf), 10},
    {NULL, NULL, 0}
};

void R_init_DirStats(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
