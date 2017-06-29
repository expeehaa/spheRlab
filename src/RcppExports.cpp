// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// slplotpolygon
void slplotpolygon(List plotinitres, NumericVector lon, NumericVector lat, String colfill, String colborder, double borderlwd, int borderlty, bool ignorevisibility, bool removeidenticalneighbours, bool refineboundary, double refineboundaryprecision);
RcppExport SEXP spheRlab_slplotpolygon(SEXP plotinitresSEXP, SEXP lonSEXP, SEXP latSEXP, SEXP colfillSEXP, SEXP colborderSEXP, SEXP borderlwdSEXP, SEXP borderltySEXP, SEXP ignorevisibilitySEXP, SEXP removeidenticalneighboursSEXP, SEXP refineboundarySEXP, SEXP refineboundaryprecisionSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type plotinitres(plotinitresSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lon(lonSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lat(latSEXP);
    Rcpp::traits::input_parameter< String >::type colfill(colfillSEXP);
    Rcpp::traits::input_parameter< String >::type colborder(colborderSEXP);
    Rcpp::traits::input_parameter< double >::type borderlwd(borderlwdSEXP);
    Rcpp::traits::input_parameter< int >::type borderlty(borderltySEXP);
    Rcpp::traits::input_parameter< bool >::type ignorevisibility(ignorevisibilitySEXP);
    Rcpp::traits::input_parameter< bool >::type removeidenticalneighbours(removeidenticalneighboursSEXP);
    Rcpp::traits::input_parameter< bool >::type refineboundary(refineboundarySEXP);
    Rcpp::traits::input_parameter< double >::type refineboundaryprecision(refineboundaryprecisionSEXP);
    slplotpolygon(plotinitres, lon, lat, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
    return R_NilValue;
END_RCPP
}
// slplotfieldloop
void slplotfieldloop(List plotinitres, NumericVector num, NumericMatrix lonv, NumericMatrix latv, List colbar, List colbarbreaks, String colfill, String colborder, bool colbarbreakslog, double borderlwd, int borderlty);
RcppExport SEXP spheRlab_slplotfieldloop(SEXP plotinitresSEXP, SEXP numSEXP, SEXP lonvSEXP, SEXP latvSEXP, SEXP colbarSEXP, SEXP colbarbreaksSEXP, SEXP colfillSEXP, SEXP colborderSEXP, SEXP colbarbreakslogSEXP, SEXP borderlwdSEXP, SEXP borderltySEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type plotinitres(plotinitresSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type num(numSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type lonv(lonvSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type latv(latvSEXP);
    Rcpp::traits::input_parameter< List >::type colbar(colbarSEXP);
    Rcpp::traits::input_parameter< List >::type colbarbreaks(colbarbreaksSEXP);
    Rcpp::traits::input_parameter< String >::type colfill(colfillSEXP);
    Rcpp::traits::input_parameter< String >::type colborder(colborderSEXP);
    Rcpp::traits::input_parameter< bool >::type colbarbreakslog(colbarbreakslogSEXP);
    Rcpp::traits::input_parameter< double >::type borderlwd(borderlwdSEXP);
    Rcpp::traits::input_parameter< int >::type borderlty(borderltySEXP);
    slplotfieldloop(plotinitres, num, lonv, latv, colbar, colbarbreaks, colfill, colborder, colbarbreakslog, borderlwd, borderlty);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"spheRlab_slplotpolygon", (DL_FUNC) &spheRlab_slplotpolygon, 11},
    {"spheRlab_slplotfieldloop", (DL_FUNC) &spheRlab_slplotfieldloop, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_spheRlab(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}