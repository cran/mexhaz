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
extern SEXP HazardWeibC(SEXP x0, SEXP x, SEXP nph, SEXP fixobs, SEXP param, SEXP paramf);
extern SEXP HazardWeibL(SEXP x0, SEXP x, SEXP nph, SEXP fixobs, SEXP param, SEXP paramf);
extern SEXP HazardWeibR(SEXP x, SEXP nph, SEXP fixobs, SEXP param, SEXP paramf);
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
extern SEXP HGH_BsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HGH_BsRx(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP lambdaobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HGH_NsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HGH_NsRx(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP lambdaobs, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HGH_PwR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HGH_PwRx(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP lambdaobs, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HGH_WbR(SEXP x, SEXP nph, SEXP fixobs, SEXP param, SEXP paramf);
extern SEXP HGH_WbRx(SEXP x, SEXP nph, SEXP fixobs, SEXP lambdaobs, SEXP param, SEXP paramf);
extern SEXP HGHAggr_BsL(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HGHAggr_BsLx(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HGHAggr_BsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HGHAggr_BsRx(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk);
extern SEXP HGHAggr_NsL(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HGHAggr_NsLx(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HGHAggr_NsR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HGHAggr_NsRx(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP deg, SEXP n, SEXP lw, SEXP matk, SEXP totk, SEXP intk, SEXP nsadj1, SEXP nsadj2);
extern SEXP HGHAggr_Pois(SEXP fixobs, SEXP statobs, SEXP offobs, SEXP nbyclust, SEXP paramf);
extern SEXP HGHAggr_PwL(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HGHAggr_PwLx(SEXP x0, SEXP x, SEXP nph, SEXP timecat0, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HGHAggr_PwR(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HGHAggr_PwRx(SEXP x, SEXP nph, SEXP timecat, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf, SEXP matk);
extern SEXP HGHAggr_WbL(SEXP x0, SEXP x, SEXP nph, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf);
extern SEXP HGHAggr_WbLx(SEXP x0, SEXP x, SEXP nph, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf);
extern SEXP HGHAggr_WbR(SEXP x, SEXP nph, SEXP fixobs, SEXP statobs, SEXP nbyclust, SEXP param, SEXP paramf);
extern SEXP HGHAggr_WbRx(SEXP x, SEXP nph, SEXP fixobs, SEXP statobs, SEXP lambdaobs, SEXP nbyclust, SEXP param, SEXP paramf);
extern SEXP GaussProcNPH(SEXP vecnb, SEXP pDtau, SEXP lrowHXS, SEXP maxSXt, SEXP cst, SEXP idxpDx, SEXP keep);
extern SEXP GaussProcLIN(SEXP vecnb, SEXP event, SEXP pDtau, SEXP lrowHXS, SEXP maxSXt, SEXP cst, SEXP idxpDx, SEXP ordDx, SEXP nbval, SEXP keep);

const static R_CallMethodDef R_CallDef[] = {
    {"DeltaBs0R",    (DL_FUNC) &DeltaBs0R,     8},
    {"DeltaBs1R",    (DL_FUNC) &DeltaBs1R,     9},
    {"DeltaBs23R",   (DL_FUNC) &DeltaBs23R,   12}, 
    {"DeltaNsR",     (DL_FUNC) &DeltaNsR,     15}, 
    {"DeltaWeibR",   (DL_FUNC) &DeltaWeibR,    6},
    {"FrailtyAdapt", (DL_FUNC) &FrailtyAdapt, 11},
    {"FrailtyAdaptL",(DL_FUNC) &FrailtyAdaptL,14},
    {"HazardWeibC",  (DL_FUNC) &HazardWeibC,   6},
    {"HazardWeibL",  (DL_FUNC) &HazardWeibL,   6},
    {"HazardWeibR",  (DL_FUNC) &HazardWeibR,   5},
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
    {"HGH_BsR",      (DL_FUNC) &HGH_BsR,      11},
    {"HGH_BsRx",     (DL_FUNC) &HGH_BsRx,     12},
    {"HGH_NsR",      (DL_FUNC) &HGH_NsR,      14},
    {"HGH_NsRx",     (DL_FUNC) &HGH_NsRx,     15},
    {"HGH_PwR",      (DL_FUNC) &HGH_PwR,       7},
    {"HGH_PwRx",     (DL_FUNC) &HGH_PwRx,      8},
    {"HGH_WbR",      (DL_FUNC) &HGH_WbR,       5},
    {"HGH_WbRx",     (DL_FUNC) &HGH_WbRx,      6},
    {"HGHAggr_BsL",  (DL_FUNC) &HGHAggr_BsL,  15},
    {"HGHAggr_BsLx", (DL_FUNC) &HGHAggr_BsLx, 16},
    {"HGHAggr_BsR",  (DL_FUNC) &HGHAggr_BsR,  13},
    {"HGHAggr_BsRx", (DL_FUNC) &HGHAggr_BsRx, 14},
    {"HGHAggr_NsL",  (DL_FUNC) &HGHAggr_NsL,  18},
    {"HGHAggr_NsLx", (DL_FUNC) &HGHAggr_NsLx, 19},
    {"HGHAggr_NsR",  (DL_FUNC) &HGHAggr_NsR,  16},
    {"HGHAggr_NsRx", (DL_FUNC) &HGHAggr_NsRx, 17},
    {"HGHAggr_Pois", (DL_FUNC) &HGHAggr_Pois,  5},    
    {"HGHAggr_PwL",  (DL_FUNC) &HGHAggr_PwL,  11},
    {"HGHAggr_PwLx", (DL_FUNC) &HGHAggr_PwLx, 12},
    {"HGHAggr_PwR",  (DL_FUNC) &HGHAggr_PwR,   9},
    {"HGHAggr_PwRx", (DL_FUNC) &HGHAggr_PwRx, 10},
    {"HGHAggr_WbL",  (DL_FUNC) &HGHAggr_WbL,   8},
    {"HGHAggr_WbLx", (DL_FUNC) &HGHAggr_WbLx,  9},
    {"HGHAggr_WbR",  (DL_FUNC) &HGHAggr_WbR,   7},
    {"HGHAggr_WbRx", (DL_FUNC) &HGHAggr_WbRx,  8},
    {"GaussProcNPH", (DL_FUNC) &GaussProcNPH,  7},
    {"GaussProcLIN", (DL_FUNC) &GaussProcLIN, 10},
    {NULL, NULL, 0}
};

void R_init_mexhaz(DllInfo *dll){
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
