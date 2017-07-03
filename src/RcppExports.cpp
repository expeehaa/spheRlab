// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// slplotpolygon
void slplotpolygon(PIR pir, NumericVector lon, NumericVector lat, Function polygon, String colfill, String colborder, double borderlwd, int borderlty, bool ignorevisibility, bool removeidenticalneighbours, bool refineboundary, double refineboundaryprecision);
RcppExport SEXP spheRlab_slplotpolygon(SEXP pirSEXP, SEXP lonSEXP, SEXP latSEXP, SEXP polygonSEXP, SEXP colfillSEXP, SEXP colborderSEXP, SEXP borderlwdSEXP, SEXP borderltySEXP, SEXP ignorevisibilitySEXP, SEXP removeidenticalneighboursSEXP, SEXP refineboundarySEXP, SEXP refineboundaryprecisionSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< PIR >::type pir(pirSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lon(lonSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lat(latSEXP);
    Rcpp::traits::input_parameter< Function >::type polygon(polygonSEXP);
    Rcpp::traits::input_parameter< String >::type colfill(colfillSEXP);
    Rcpp::traits::input_parameter< String >::type colborder(colborderSEXP);
    Rcpp::traits::input_parameter< double >::type borderlwd(borderlwdSEXP);
    Rcpp::traits::input_parameter< int >::type borderlty(borderltySEXP);
    Rcpp::traits::input_parameter< bool >::type ignorevisibility(ignorevisibilitySEXP);
    Rcpp::traits::input_parameter< bool >::type removeidenticalneighbours(removeidenticalneighboursSEXP);
    Rcpp::traits::input_parameter< bool >::type refineboundary(refineboundarySEXP);
    Rcpp::traits::input_parameter< double >::type refineboundaryprecision(refineboundaryprecisionSEXP);
    slplotpolygon(pir, lon, lat, polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
    return R_NilValue;
END_RCPP
}
// slplotfieldloop
void slplotfieldloop(List plotinitres, NumericMatrix lonv, NumericMatrix latv, List colbar, IntegerVector colind, String colfill, String colborder, double borderlwd, int borderlty, int threads);
RcppExport SEXP spheRlab_slplotfieldloop(SEXP plotinitresSEXP, SEXP lonvSEXP, SEXP latvSEXP, SEXP colbarSEXP, SEXP colindSEXP, SEXP colfillSEXP, SEXP colborderSEXP, SEXP borderlwdSEXP, SEXP borderltySEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type plotinitres(plotinitresSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lonv(lonvSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type latv(latvSEXP);
    Rcpp::traits::input_parameter< List >::type colbar(colbarSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type colind(colindSEXP);
    Rcpp::traits::input_parameter< String >::type colfill(colfillSEXP);
    Rcpp::traits::input_parameter< String >::type colborder(colborderSEXP);
    Rcpp::traits::input_parameter< double >::type borderlwd(borderlwdSEXP);
    Rcpp::traits::input_parameter< int >::type borderlty(borderltySEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    slplotfieldloop(plotinitres, lonv, latv, colbar, colind, colfill, colborder, borderlwd, borderlty, threads);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"spheRlab_slplotpolygon", (DL_FUNC) &spheRlab_slplotpolygon, 12},
    {"spheRlab_slplotfieldloop", (DL_FUNC) &spheRlab_slplotfieldloop, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_spheRlab(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
