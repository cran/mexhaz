#include <stdlib.h> /* for NULL */
#include <Rinternals.h> /* for SEXP */
#include <R_ext/Rdynload.h>

extern SEXP DeltaBs0R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP paramt, SEXP matk, SEXP varcov, SEXP grad);
extern SEXP DeltaBs1R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP paramt, SEXP matk, SEXP totk, SEXP varcov, SEXP grad);
extern SEXP DeltaBs23R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP paramt, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP varcov, SEXP grad);
extern SEXP DeltaNsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP paramt, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2, SEXP varcov, SEXP grad);
extern SEXP DeltaWeibR(SEXP x, SEXP nph, SEXP fixobs, SEXP paramt, SEXP varcov, SEXP grad);
extern SEXP FrailtyAdapt(SEXP nodes, SEXP nodessquare, SEXP logweights, SEXP clust, SEXP clustd, SEXP expect, SEXP betal, SEXP betaL, SEXP A, SEXP var, SEXP muhatcond);
extern SEXP FrailtyAdaptL(SEXP nodes, SEXP nodessquare, SEXP logweights, SEXP clust, SEXP clustd, SEXP expect, SEXP betal, SEXP betaL0, SEXP betaL, SEXP A0, SEXP A, SEXP var, SEXP mh0, SEXP muhatcond);
extern SEXP HazardBs0C(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HazardBs0L(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HazardBs0R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HazardBs1C(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk, SEXP totk);
extern SEXP HazardBs1L(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk, SEXP totk);
extern SEXP HazardBs1R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk, SEXP totk);
extern SEXP HazardBs23C(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HazardBs23L(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HazardBs23R(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HazardNsC(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HazardNsL(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HazardNsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);

const static R_CallMethodDef R_CallDef[] = {
    {"DeltaBs0R",    (DL_FUNC) &DeltaBs0R,     8},
    {"DeltaBs1R",    (DL_FUNC) &DeltaBs1R,     9},
    {"DeltaBs23R",   (DL_FUNC) &DeltaBs23R,   12}, 
    {"DeltaNsR",     (DL_FUNC) &DeltaNsR,     15}, 
    {"DeltaWeibR",   (DL_FUNC) &DeltaWeibR,    6},
    {"FrailtyAdapt", (DL_FUNC) &FrailtyAdapt, 11},
    {"FrailtyAdaptL",(DL_FUNC) &FrailtyAdaptL,14},
    {"HazardBs0C",   (DL_FUNC) &HazardBs0C,    9},
    {"HazardBs0L",   (DL_FUNC) &HazardBs0L,    9},
    {"HazardBs0R",   (DL_FUNC) &HazardBs0R,    7},
    {"HazardBs1C",   (DL_FUNC) &HazardBs1C,   10},
    {"HazardBs1L",   (DL_FUNC) &HazardBs1L,   10},
    {"HazardBs1R",   (DL_FUNC) &HazardBs1R,    8},
    {"HazardBs23C",  (DL_FUNC) &HazardBs23C,  13},
    {"HazardBs23L",  (DL_FUNC) &HazardBs23L,  13},
    {"HazardBs23R",  (DL_FUNC) &HazardBs23R,  11},
    {"HazardNsC",    (DL_FUNC) &HazardNsC,    16},
    {"HazardNsL",    (DL_FUNC) &HazardNsL,    16},
    {"HazardNsR",    (DL_FUNC) &HazardNsR,    14},
    {NULL, NULL, 0}
};

void R_init_mexhaz(DllInfo *dll){
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
