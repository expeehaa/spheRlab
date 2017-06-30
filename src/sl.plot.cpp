// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

IntegerVector which(LogicalVector lv){
  IntegerVector iv(sum(lv));
  int counter = 0;
  for(int i = 0; i < lv.length(); i++){
    if(lv[i]) iv[counter++] = i;
  }
  return iv;
}

template <int VECTYPE>
Vector<VECTYPE> shiftleft(Vector<VECTYPE> vector){
  if(vector.length() <= 1) return vector;
  Vector<VECTYPE> out(vector.length());
  for(int i = 1; i < vector.length(); i++){
    out[i-1] = vector[i];
  }
  out[vector.length()-1] = vector[0];
  return out;
}

IntegerVector renumber(int lower, int upper, int step = 1){
  IntegerVector out(upper-lower+1);
  for(int i = lower; i <= upper; i = i + step) out[(i-lower)/step] = i;
  return out;
}

//template<int VECTYPE>
NumericVector removeEntries(NumericVector vec, LogicalVector logivec){
  while(logivec.length() < vec.length()) logivec[logivec.length()] = false;
  NumericVector out(sum(logivec));
  //Rcout << "logivec sum: " << sum(logivec) << '\n';
  int counter = 0;
  for(int i = 0; i < vec.length(); i++){
    //Rcout << logivec[i] << ", " << vec[i] << '\n';
    if(logivec[i]) out[counter++] = vec[i];
  }
  return out;
}

double abs(double v){
  return v < 0 ? -v : v;
}

template <int VECTYPE>
Vector<VECTYPE> cVecs(Vector<VECTYPE> vec1, Vector<VECTYPE> vec2){
  Vector<VECTYPE> out(vec1.length()+vec2.length());
  for(int i = 0; i < vec1.length(); i++){
    out[i] = vec1[i];
  }
  int l1 = vec1.length();
  for(int i = 0; i < vec2.length(); i++){
    out[i+l1] = vec2[i];
  }
  return out;
}

/// [[Rcpp::export("sl.plot.polygon")]]
// [[Rcpp::export]]
void slplotpolygon(List plotinitres, NumericVector lon, NumericVector lat, String colfill = "black", String colborder = "black", double borderlwd = 0.01, int borderlty = 1, bool ignorevisibility = false, bool removeidenticalneighbours = true, bool refineboundary = true, double refineboundaryprecision = 1){
  
  // Rcout << "polygon\n";
  
  Environment env("package:spheRlab");
  Environment envGraphics("package:graphics");
  Environment envBase("package:base");
  Environment envGlob = Environment::global_env();
  Function polygon = as<Function>(envGraphics["polygon"]);
  
  int L = lon.length();
  if(L != lat.length()) {
    Rcpp::warning("lon and lat vector do not have the same length!");
    return;
  }
  // Rcout << "length: " << L << '\n';
  // Rcout << removeidenticalneighbours << '\n';
  
  if(removeidenticalneighbours){
    // Rcout << wrap(lon) << ", " << wrap(lat) << '\n';
    Function sllonlatidentical = as<Function>(env["sl.lonlat.identical"]);
    // Rcout << "after wrap\n";
    LogicalVector keep = !as<LogicalVector>(sllonlatidentical(lon, lat, shiftleft(lon), shiftleft(lat)));
    lon = removeEntries(lon, keep);
    lat = removeEntries(lat, keep);
    L = lon.length();
  }
  
  // Rcout << "length after 'rin': " << L << '\n';
  
  const char* projection = as<String>(plotinitres["projection"]).get_cstring();
  
  // Rcout << "projection: " << projection << '\n';
  
  if(!strcmp(projection, "platon") || !strcmp(projection, "3D")){
    int npir = 1;
    while(Rcpp::is<List>(plotinitres[npir])){
      slplotpolygon(plotinitres[npir], lon, lat, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
      npir++;
    }
  }
  
  // Rcout << ignorevisibility << '\n';
  
  Function slvisshiftrot = Rcpp::as<Function>(env["sl.vis.shift.rot"]);
  List vsrres = slvisshiftrot(plotinitres, lon, lat);
  LogicalVector visible = vsrres["visible"];
  if(ignorevisibility) visible.fill(true);
  int vissum = sum(visible);
  if(vissum == 0) return;
  
  // Rcout << vissum << '\n';
  // Rcout << "vsrres: " << Rf_isNull(vsrres["rot.lon"]) << ", " << Rf_isNull(vsrres["rot.lat"]) << '\n';
  
  NumericVector x = !Rf_isNull(vsrres["x"]) ? vsrres["x"] : NumericVector::create(), 
    y = !Rf_isNull(vsrres["y"]) ? vsrres["y"] : NumericVector::create(),
    rotlon = !Rf_isNull(vsrres["rot.lon"]) ? vsrres["rot.lon"] : NumericVector::create(), 
    rotlat = !Rf_isNull(vsrres["rot.lat"]) ? vsrres["rot.lat"] : NumericVector::create();
  
  bool vispartial = vissum < L ? true : false;
  
  // Rcout << "going into if ? " << (vispartial && strcmp(projection, "lonlat")) << '\n';
  
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
  
  // Rcout << "length: " << L << '\n';
  
  double xshift = Rcpp::as<double>(plotinitres["xshift"]), yshift = Rcpp::as<double>(plotinitres["yshift"]);
  
  if(!strcmp(projection, "lonlat")){
    NumericVector lonrange = as<NumericVector>(plotinitres["lonlat.lonrange"]), latrange = as<NumericVector>(plotinitres["lonlat.latrange"]);
    if(vispartial){
      
      // Rcout << "vispartial" << '\n';
      //get function sl.segment
      Function slsegment = Rcpp::as<Function>(env["sl.segment"]);
      
      //point(s) out of west boundary
      // Rcout << "west" << '\n';
      
      if(min(x) < lonrange[0]){
        ListOf<IntegerVector> inds = as<ListOf<IntegerVector> >(slsegment(x>lonrange[0], _["extend"] = true));
        if(inds.size() > 1)
          for(int i = 1; i < inds.size(); i++)
            slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]-1];
        y = y[inds[0]-1];
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
      // Rcout << "east" << '\n';
      
      if(max(x) > lonrange[1]){
        ListOf<IntegerVector> inds = as<ListOf<IntegerVector> >(slsegment(x<lonrange[1], _["extend"] = true));
        if(inds.size() > 1)
          for(int i = 1; i < inds.size(); i++)
            slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]-1];
        y = y[inds[0]-1];
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
      // Rcout << "south" << '\n';
      
      if(min(y) < latrange[0]){
        ListOf<IntegerVector> inds = as<ListOf<IntegerVector> >(slsegment(y>latrange[1], _["extend"] = true));
        if(inds.size() > 1)
          for(int i = 1; i < inds.size(); i++)
            slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]-1];
        y = y[inds[0]-1];
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
      // Rcout << "north" << '\n';
      
      if(max(y) > latrange[1]){
        ListOf<IntegerVector> inds = as<ListOf<IntegerVector> >(slsegment(y<latrange[1], _["extend"] = true));
        if(inds.size() > 1){
          for(int i = 1; i < inds.size(); i++){
            slplotpolygon(plotinitres, x[inds[i]], y[inds[i]], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
          }
        }
        x = x[inds[0]-1];
        y = y[inds[0]-1];
        L = x.length();
        
        // Rcout << "x.length: " << L << '\n';
        
        // Rcout << (y[0] > latrange[1]) << '\n';
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
        // Rcout << (y[L-1] > latrange[1]) << '\n';
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
    
    // Rcout << "after vispartial" << '\n';
    if(max(x) - min(x) > 180 && max(abs(x - shiftleft(x))) > 180){
      IntegerVector l2r = which(shiftleft(x) - x > 180), r2l = which(shiftleft(x)-x<-180);
      int Nlr = l2r.length(), Nrl = r2l.length();
      if(Nlr != Nrl) {
        Rcpp::warning("This nasty polygon can not be plotted; it might be circular, crossing the lonlat boundary an uneven number of times, that is, it may contain one or the other pole. Consider splitting the polygon into better behaving pieces.");
        return;
      }
      // Rcout << "Nlr > 1 ? " << (Nlr > 1) << '\n';
      if(Nlr > 1){
        if(l2r[0] > r2l[0]) r2l = shiftleft(r2l);
        for (int i = 1; i< Nlr; i++) {
          IntegerVector right = r2l[i]%L+1 > l2r[i] ? renumber(l2r[i], r2l[i]%L+1) : cVecs(renumber(l2r[i], L), renumber(1, r2l[i]%L+1));
          IntegerVector left = l2r[i%Nlr]%L+1 > r2l[i] ? renumber(r2l[i], l2r[i%Nlr]%L+1) : cVecs(renumber(r2l[i], L), renumber(1, l2r[i%Nlr]%L+1));
          for(int o = 0; o < right.length(); o++)
            right[o] %= L;
          for(int o = 0; o < left.length(); o++)
            left[o] %= L;
          
          slplotpolygon(plotinitres, x[left], y[left], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
          slplotpolygon(plotinitres, x[right], y[right], colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        }
      }
      int i = 0;
      // Rcout << "right: " << (r2l[i]%L+1 > l2r[i]) << "\nleft: " << (l2r[i%Nlr]%L+1 > r2l[i]) << '\n';
      IntegerVector right = r2l[i]%L+1 > l2r[i] ? renumber(l2r[i], r2l[i]%L+1) : cVecs(renumber(l2r[i], L), renumber(1, r2l[i]%L+1));
      IntegerVector left = l2r[i%Nlr]%L+1 > r2l[i] ? renumber(r2l[i], l2r[i%Nlr]%L+1) : cVecs(renumber(r2l[i], L), renumber(1, l2r[i%Nlr]%L+1));
      for(int o = 0; o < right.length(); o++)
        right[o] %= L;
      for(int o = 0; o < left.length(); o++)
        left[o] %= L;
      
      // envGlob.assign("s2x", x);
      // envGlob.assign("s2y", y);
      // envGlob.assign("s2r2l", r2l);
      // envGlob.assign("s2l2r", l2r);
      // envGlob.assign("s2L", L);
      // envGlob.assign("s2right", right);
      // envGlob.assign("s2left", left);
      
      NumericVector xr = x[right], yr = y[right], xl = x[left], yl = y[left];
      int Lr = right.length(), Ll = left.length();
      
      // Rcout << "ranges1" << '\n';
      
      if(xr[0] < lonrange[1]) xr[0] += 360;
      yr[0] = yr[1] + (lonrange[1] - xr[1])/(xr[0] - xr[1]) * (yr[0] - yr[1]);
      xr[0] = lonrange[1];
      
      // Rcout << "ranges2" << '\n';
      
      if(xr[Lr-1] < lonrange[1]) xr[Lr-1] += 360;
      yr[Lr - 1] = yr[Lr-2] + (lonrange[1] - xr[Lr-2])/(xr[Lr-1] - xr[Lr-2]) * (yr[Lr-1] - yr[Lr-2]);
      xr[Lr-1] = lonrange[1];
      
      // Rcout << "ranges3" << '\n';
      
      if (xl[0] > lonrange[0]) xl[0] -= 360;
      yl[0] = yl[1] + (lonrange[0] - xl[1])/(xl[0] - xl[1]) * (yl[0] - yl[1]);
      xl[0] = lonrange[0];
      
      // Rcout << "ranges4" << '\n';
      
      if (xl[Ll - 1] > lonrange[0]) xl[Ll-1] -= 360;
      yl[Ll - 1] = yl[Ll-2] + (lonrange[0] - xl[Ll-2])/(xl[Ll-1] - xl[Ll-2]) * (yl[Ll-1] - yl[Ll-2]);
      xl[Ll-1] = lonrange[0];
      
      // Rcout << "plot polygons" << '\n';
      
      #pragma omp critical
      {
        polygon(_["x"] = xr+xshift, _["y"] = yr+yshift, _["col"] = colfill, _["lwd"] = borderlwd, _["lty"] = borderlty, _["border"] = colborder);
        polygon(_["x"] = xl+xshift, _["y"] = yl+yshift, _["col"] = colfill, _["lwd"] = borderlwd, _["lty"] = borderlty, _["border"] = colborder);
      }
    }
    else {
      #pragma omp critical
      {
        // Rcout << "just plot" << '\n';
        polygon(_["x"] = x+xshift, _["y"] = y+yshift, _["col"] = colfill, _["lwd"] = borderlwd, _["lty"] = borderlty, _["border"] = colborder);
      }
    }
  }
  else{
    Rcpp::stop("Other projections than 'lonlat' not yet implemented!");
  }
}

/// [[Rcpp::export(".sl.plot.field.loop")]]
// [[Rcpp::export]]
void slplotfieldloop(List plotinitres, NumericMatrix lonv, NumericMatrix latv, List colbar, IntegerVector colind, String colfill = "black", String colborder = "black", double borderlwd = 0.01, int borderlty = 1, int threads = 2) {
  // Rcout << "got into c++ for loop\n";
  // Rcout << colind.length() << ", " << max(colind) << ", " << colbar.size();
  #ifdef _OPENMP
  //setenv("OMP_STACKSIZE","100M",1);
  omp_set_num_threads(1);
  #pragma omp parallel for shared(plotinitres, lonv, latv, colbar, colind, colfill, colborder, borderlwd, borderlty) schedule(dynamic, 10)
  #endif
  for(int i = 0; i < lonv.rows(); i++){
    #ifdef _OPENMP
    Rcout << omp_get_thread_num() << ": " << i << '\n';
    #endif
    #ifndef _OPENMP
    Rcout << i << '\n';
    #endif
    String cbfill = colfill == "colbar" ? colbar[colind[i]-1] : colfill, cbborder = colborder == "colbar" ? colbar[colind[i]-1] : colborder;
    slplotpolygon(plotinitres, lonv(i, _), latv(i, _), cbfill, cbborder, borderlwd, borderlty);
  }
}