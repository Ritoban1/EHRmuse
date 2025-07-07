#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the external C functions
extern void Update(double *beta, double *Sigma, double *alpha, double *mu, double *sigma, double *y, double *X, double *Z, 
                   double *T, int *ID, int *N, int *m, int *p, int *q, int *link, int *type);

extern void Update_homo(double *beta, double *alpha, double *mu, double *sigma, double *y, double *X, 
                        double *T, int *ID, int *N, int *m, int *p, int *link, int *type);

// Create the registration table
static const R_CMethodDef CEntries[] = {
    {"Update", (DL_FUNC) &Update, 16},
    {"Update_homo", (DL_FUNC) &Update_homo, 13},
    {NULL, NULL, 0}
};

// Register routines
void R_init_EHRmuse(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
