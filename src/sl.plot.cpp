// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
bool OPENMPEXISTS = true;
#endif
#ifndef _OPENMP
bool OPENMPEXISTS = false;
#endif
using namespace Rcpp;


IntegerVector which(LogicalVector lv){
  IntegerVector iv = IntegerVector::create();
  for(int i = 0; i < lv.length(); i++){
    if(lv[i]) iv[iv.length()] = i;
  }
  return iv;
}

template <int VECTYPE>
Vector<VECTYPE> shiftleft(Vector<VECTYPE> vector){
  Vector<VECTYPE> out(vector.length());
  for(int i = 1; i < vector.length(); i++){
    out[i-1] = vector[i];
  }
  out[vector.length()-1] = vector[0];
  return out;
}

IntegerVector renumber(int lower, int upper, int step = 1, ...){
  IntegerVector out = IntegerVector::create();
  for(int i = lower; i <= upper; i = i + step) out[out.length()] = i;
  return out;
}

template<int VECTYPE>
Vector<VECTYPE> removeEntries(Vector<VECTYPE> vec, LogicalVector logivec){
  while(logivec.length() < vec.length()) logivec[logivec.length()] = false;
  Vector<VECTYPE> out(sum(logivec));
  for(int i = 0; i < vec.length(); i++){
    if(logivec[i]) out[out.length()] = vec[i];
  }
  return out;
}

/// [[Rcpp::export("sl.plot.polygon")]]
// [[Rcpp::export]]
void slplotpolygon(List plotinitres, NumericVector lon, NumericVector lat, String colfill = "black", String colborder = "black", double borderlwd = 0.01, int borderlty = 1, bool ignorevisibility = false, bool removeidenticalneighbours = true, bool refineboundary = true, double refineboundaryprecision = 1){
  
  Environment env = Environment::global_env();
  
  Rcout << "Got global env! OpenMP? " << OPENMPEXISTS;
  int L = lon.length();
  if(L != lat.length()) {
    Rcpp::warning("lon and lat vector do not have the same length!");
    return;
  }
  Rcout << "length: " << L;
  
  if(removeidenticalneighbours){
    Function sllonlatidentical = env.get("sl.lonlat.identical");
    LogicalVector keep = !(sllonlatidentical(lon, lat, shiftleft(lon), shiftleft(lat)));
    lon = removeEntries(lon, keep);
    lat = removeEntries(lat, keep);
    L = lon.length();
  }
  
  Rcout << "length after 'rin': " << L;
  
  Rcpp::String projectionbuf = plotinitres["projection"];
  const char* projection = projectionbuf.get_cstring();
  
  Rcout << "projection: " << projection;
  
  if(!strcmp(projection, "platon") || !strcmp(projection, "3D")){
    int npir = 1;
    while(Rcpp::is<List>(plotinitres[npir])){
      slplotpolygon(plotinitres[npir], lon, lat, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
      npir++;
    }
  }
  
  Rcout << ignorevisibility;
  
  Function slvisshiftrot = env["sl.vis.shift.rot"];
  List vsrres = slvisshiftrot(plotinitres, lon, lat);
  LogicalVector visible = vsrres["visible"];
  if(ignorevisibility) visible.fill(true);
  int vissum = sum(visible);
  if(vissum == 0) return;
  Rcout << vissum;
  NumericVector x = vsrres["x"], y = vsrres["y"], rotlon = vsrres["rot.lon"], rotlat = vsrres["rot.lat"];
  
  
  bool vispartial = vissum < L ? true : false;
  Rcout << "going into if ? " << (vispartial && strcmp(projection, "lonlat"));
  if(vispartial && strcmp(projection, "lonlat")){
    LogicalVector visibleext = visible;
    visibleext[visibleext.length()] = visible[1];
    NumericVector xnew, ynew, rotlonnew, rotlatnew;
    int i = 0;
    while(!(!visibleext[i] && visibleext[i+1])) i++;
    do {
      xnew[xnew.length()] = x[i];
      ynew[ynew.length()] = y[i];
      rotlonnew[rotlonnew.length()] = rotlon[i];
      rotlatnew[rotlatnew.length()] = rotlat[i];
      i = (i % L) + 1;
    }
    while(!(!visibleext[i] && visibleext[i+1]));
    xnew[xnew.length()] = x[i];
    ynew[ynew.length()] = y[i];
    rotlonnew[rotlonnew.length()] = rotlon[i];
    rotlatnew[rotlatnew.length()] = rotlat[i];
    i = (i % L) + 1;
    xnew[xnew.length()] = x[i];
    ynew[ynew.length()] = y[i];
    rotlonnew[rotlonnew.length()] = rotlon[i];
    rotlatnew[rotlatnew.length()] = rotlat[i];
    x = xnew;
    y = ynew;
    rotlon = rotlonnew;
    rotlat = rotlatnew;
    L = x.length();
  }
  
  Rcout << "length: " << L;
  
  double xshift = plotinitres["xshift"], yshift = plotinitres["yshift"];
  
  if(!strcmp(projection, "lonlat")){
    NumericVector lonrange = plotinitres["lonlat.lonrange"], latrange = plotinitres["lonlat.latrange"];
    if(vispartial){
      //point(s) out of west boundary
      if(min(x) < lonrange[0]){
        ListOf<IntegerVector> inds = ((Function)env["sl.segment"])(x>lonrange[0], true);
        for(int i = 1; i < inds.size(); i++)
          slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        if(x[0] < lonrange[0]){
          if(x[1]-x[0]>180){
            y[0] = y[1] + (lonrange[1] - x[1])/(x[0] + 360 - x[1]) * (y[0] - y[1]);
            x[0] = lonrange[1];
          }
          else{
            y[0] = y[1] + (lonrange[0] - x[1])/(x[0] - x[1]) * (y[0] - y[1]);
            x[0] = lonrange[0];
          }
        }
        if(x[L-1] < lonrange[1]){
          if(x[L-2]-x[L-1]>180){
            y[L-1] = y[L-2] + (lonrange[1] - x[L-2])/(x[L-1] + 360 - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = lonrange[1];
          }
          else{
            y[L-1] = y[L-2] + (lonrange[0] - x[L-2])/(x[L-1] - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = lonrange[0];
          }
        }
      }
      //point(s) out of east boundary
      if(max(x) > lonrange[1]){
        ListOf<IntegerVector> inds = ((Function)env["sl.segment"])(x<lonrange[1], true);
        for(int i = 1; i < inds.size(); i++)
          slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        if (x[0] > lonrange[1]) {
          if (x[0]-x[1]>180) {
            y[0] = y[1] + (lonrange[0] - x[1])/(x[0]-360 - x[1]) * (y[0] - y[1]);
            x[0] = lonrange[0];
          } else {
            y[0] = y[1] + (lonrange[1] - x[1])/(x[0] - x[1]) * (y[0] - y[1]);
            x[0] = lonrange[0];
          }
        }
        if (x[L-1] > lonrange[1]) {
          if (x[L-1]-x[L-2]>180) {
            y[L-1] = y[L-2] + (lonrange[0] - x[L-2])/(x[L-1]-360 - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = lonrange[0];
          } else {
            y[L-1] = y[L-2] + (lonrange[1] - x[L-2])/(x[L-1] - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = lonrange[1];
          }
        }
      }
      //point(s) out of south boundary
      if(min(y) < latrange[0]){
        ListOf<IntegerVector> inds = ((Function)env["sl.segment"])(y>latrange[1], true);
        for(int i = 1; i < inds.size(); i++)
          slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        if (y[0] < latrange[0]) {
          if (abs(x[0]-x[1])>180) {
            x[0] += x[0]-x[1]>180 ? -360 : 360;
            x[0] = x[1] + (latrange[0] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            if (x[0]<lonrange[0]) x[0] += 360;
            if (x[0]>lonrange[1]) x[0] -= 360;
            y[0] = latrange[0];
          } 
          else {
            x[0] = x[1] + (latrange[0] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            y[0] = latrange[0];
          }
        }
        if (y[L-1] < latrange[0]) {
          if (abs(x[L-1]-x[L-2])>180) {
            x[L-1] += x[L-1]-x[L-2]>180 ? -360 : 360;
            x[L-1] = x[L-2] + (latrange[0] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            if (x[L-1]<lonrange[0]) x[L-1] += 360;
            if (x[L-1]>lonrange[1]) x[L-1] -= 360;
            y[L-1] = latrange[0];
          }
          else {
            x[L-1] = x[L-2] + (latrange[0] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            y[L-1] = latrange[0];
          }
        }
      }
      //point(s) out of north boundary
      if(max(y) > latrange[1]){
        ListOf<IntegerVector> inds = ((Function)env["sl.segment"])(y<latrange[1], true);
        for(int i = 1; i < inds.size(); i++)
          slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        if (y[0] > latrange[1]) {
          if (abs(x[0]-x[1])>180) {
            x[0] += x[0]-x[1]>180 ? -360 : 360;
            x[0] = x[1] + (latrange[1] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            if (x[0]<lonrange[0]) x[1] += 360;
            if (x[0]>lonrange[1]) x[1] -= 360;
            y[0] = latrange[1];
          } 
          else {
            x[0] = x[1] + (latrange[1] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            y[0] = latrange[1];
          }
        }
        if (y[L-1] > latrange[1]) {
          if (abs(x[L-1]-x[L-2])>180) {
            x[L-1] += x[L-1]-x[L-2]>180 ? -360 : 360;
            x[L-1] = x[L-1] + (latrange[1] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            if (x[L-1]<lonrange[0]) x[L-1] += 360;
            if (x[L-1]>lonrange[1]) x[L-1] -= 360;
            y[L-1] = latrange[1];
          } 
          else {
            x[L-1] = x[L-2] + (latrange[1] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            y[L-1] = latrange[1];
          }
        }
      }
    }
    
    if(max(x) - min(x) > 180 && max(abs(x - shiftleft(x))) > 180){
      IntegerVector l2r = which(shiftleft(x) - x > 180), r2l = which(shiftleft(x)-x<-180);
      int Nlr = l2r.length(), Nrl = r2l.length();
      if(Nlr != Nrl) {
        Rcpp::warning("This nasty polygon can not be plotted; it might be circular, crossing the lonlat boundary an uneven number of times, that is, it may contain one or the other pole. Consider splitting the polygon into better behaving pieces.");
        return;
      }
      
      if(Nlr > 1){
        if(l2r[0] > r2l[0]) r2l = shiftleft(r2l);
        for(int i = 1; i < Nlr; i++){
          IntegerVector right = r2l[i]%L+1 > l2r[i] ? renumber(l2r[i], r2l[i]%L+1) : (IntegerVector)(((Function)env["c"])(renumber(l2r[i],L), renumber(1, r2l[i]%L+1)));
          IntegerVector left = l2r[i%Nlr+1]%L+1 > r2l[i] ? renumber(r2l[i], l2r[i%Nlr+1]%L+1) : (IntegerVector)(((Function)env["c"])(renumber(r2l[i], L), renumber(1,l2r[i%Nlr+1]%L+1)));
          slplotpolygon(plotinitres, x[left], y[left], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
          slplotpolygon(plotinitres, x[right], y[right], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        }
      }
      
      IntegerVector right = r2l[1]%L+1 > l2r[1] ? renumber(l2r[1], r2l[1]%L+1) : (IntegerVector)(((Function)env["c"])(renumber(l2r[1],L), renumber(1, r2l[1]%L+1)));
      IntegerVector left = l2r[1%Nlr+1]%L+1 > r2l[1] ? renumber(r2l[1], l2r[1%Nlr+1]%L+1) : (IntegerVector)(((Function)env["c"])(renumber(r2l[1], L), renumber(1,l2r[1%Nlr+1]%L+1)));
      NumericVector xr = x[right], yr = y[right], xl = x[left], yl = y[left];
      int Lr = right.length(), Ll = left.length();
      
      if(xr[0] < lonrange[1]) xr[0] += 360;
      yr[0] = yr[1] + (lonrange[1] - xr[1])/(xr[0] - xr[1]) * (yr[0] - yr[1]);
      xr[0] = lonrange[1];
      
      if(xr[Lr-1] < lonrange[1]) xr[Lr-1] += 360;
      yr[Lr - 1] = yr[Lr-2] + (lonrange[1] - xr[Lr-2])/(xr[Lr-1] - xr[Lr-2]) * (yr[Lr-1] - yr[Lr-2]);
      xr[Lr-1] = lonrange[1];
      
      if (xl[0] > lonrange[0]) xl[0] -= 360;
      yl[0] = yl[1] + (lonrange[0] - xl[1])/(xl[0] - xl[1]) * (yl[0] - yl[1]);
      xl[0] = lonrange[0];
      
      if (xl[Ll - 1] > lonrange[0]) xl[Ll-1] -= 360;
      yl[Ll - 1] = yl[Ll-2] + (lonrange[0] - xl[Ll-2])/(xl[Ll-1] - xl[Ll-2]) * (yl[Ll-1] - yl[Ll-2]);
      xl[Ll-1] = lonrange[0];
      ((Function)env["polygon"])(xr+xshift, yr+yshift, colfill, borderlwd, borderlty, colborder);
      ((Function)env["polygon"])(xl+xshift, yl+yshift, colfill, borderlwd, borderlty, colborder);
    }
    else {
      ((Function)env["polygon"])(x+xshift, y+yshift, colfill, borderlwd, borderlty, colborder);
    }
  }
  else{
    Rcpp::stop("Other projections than 'lonlat' not yet implemented!");
  }
}

/// [[Rcpp::export(".sl.plot.field.loop")]]
// [[Rcpp::export]]
void slplotfieldloop(List plotinitres, NumericVector num, NumericMatrix lonv, NumericMatrix latv, List colbar, List colbarbreaks, String colfill = "black", String colborder = "black", bool colbarbreakslog = false, double borderlwd = 0.01, int borderlty = 1) {
  return;
}