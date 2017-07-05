// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
//#define DEBUG
#define PLOTPOLYGONS
using namespace Rcpp;

class PIR{
public:
  String projection;
  NumericVector lonrange, latrange, xlim, ylim;
  int xshift, yshift;
  std::vector<PIR> pirlist;
  
  int is3D() {
    return (projection == "3D" || projection == "platon") ? 1 : (projection == "lonlat" || projection == "polar" || projection == "regpoly" ? 0 : -1);
  }
  
  PIR(const String& projection, const NumericVector& lonrange, const NumericVector& latrange, const NumericVector& xlim, const NumericVector& ylim, const int& xshift, const int& yshift) 
    : projection(projection), lonrange(lonrange), latrange(latrange), xlim(xlim), ylim(ylim), xshift(xshift), yshift(yshift){}
  
  PIR(const List& pir){
    projection = as<String>(pir["projection"]);
    
    if(is3D() == 1){
      for(int i = 1; Rcpp::is<List>(pir[i]); i++){
        pirlist.push_back(pir[i]);
      }
    }
    else if(is3D() == 0){
      lonrange = as<NumericVector>(pir["lonlat.lonrange"]);
      latrange = as<NumericVector>(pir["lonlat.latrange"]);
      xlim = as<NumericVector>(pir["xlim"]);
      ylim = as<NumericVector>(pir["ylim"]);
      xshift = as<int>(pir["xshift"]);
      yshift = as<int>(pir["yshift"]);
    }
    else{
      Rcpp::stop("projection type unknown!");
    }
  }
  
  PIR& operator=(const PIR &pir){
    projection = pir.projection;
    lonrange = pir.lonrange;
    latrange = pir.latrange;
    xlim = pir.xlim;
    ylim = pir.ylim;
    xshift = pir.xshift;
    yshift = pir.yshift;
    return *this;
  }
  
  PIR(const PIR& pir){
    projection = pir.projection;
    lonrange = pir.lonrange;
    latrange = pir.latrange;
    xlim = pir.xlim;
    ylim = pir.ylim;
    xshift = pir.xshift;
    yshift = pir.yshift;
  }
};

class VisShiftRot {
public:
  NumericVector x, y, rotlon, rotlat;
  LogicalVector visible;
  
  VisShiftRot(const NumericVector& x, const NumericVector& y, const NumericVector& rotlon, const NumericVector& rotlat, const LogicalVector& visible)
    : x(x), y(y), rotlon(rotlon), rotlat(rotlat), visible(visible){}
  VisShiftRot(const PIR& pir, const NumericVector& lon, const NumericVector& lat);
};

IntegerVector which(const LogicalVector& lv){
  IntegerVector iv(sum(lv));
  int counter = 0;
  for(int i = 0; i < lv.length(); i++){
    if(lv[i]) iv[counter++] = i;
  }
  return iv;
}

template <int VECTYPE>
Vector<VECTYPE> shiftleft(const Vector<VECTYPE>& vector){
  if(vector.length() <= 1) return vector;
  Vector<VECTYPE> out(vector.length());
  for(int i = 1; i < vector.length(); i++){
    out[i-1] = vector[i];
  }
  out[vector.length()-1] = vector[0];
  return out;
}

IntegerVector renumber(const int& lower, const int& upper, const int& step = 1){
  IntegerVector out(upper-lower+1);
  for(int i = lower; i <= upper; i = i + step) out[(i-lower)/step] = i;
  return out;
}

//template<int VECTYPE>
NumericVector removeEntries(const NumericVector& vec, LogicalVector& logivec){
  while(logivec.length() < vec.length()) logivec[logivec.length()] = false;
  NumericVector out(sum(logivec));
  int counter = 0;
  for(int i = 0; i < vec.length(); i++){
    if(logivec[i]) out[counter++] = vec[i];
  }
  return out;
}

double abs(double v){
  return v < 0 ? -v : v;
}

template <int VECTYPE>
Vector<VECTYPE> cVecs(const Vector<VECTYPE>& vec1, const Vector<VECTYPE>& vec2){
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

NumericVector cutoutVector(const NumericVector& vec, const double& minimum, const double& maximum){
  NumericVector ret(vec.begin(), vec.end());
  for(int i = 0; i < ret.length(); i++){
    while(ret[i] > maximum) ret[i] -= maximum - minimum;
    while(ret[i] <= minimum) ret[i] += maximum - minimum;
  }
  return ret;
}

IntegerVector modulo(const IntegerVector& vec, const int& n){
  IntegerVector ret(vec.length());
  for(int i = 0; i < vec.length(); i++)
    ret[i] = vec[i] % n;
  return ret;
}


List slsegment(const LogicalVector& logivec, const bool& extend = false, const bool& firstonly = false, const bool& segments = true){
  int L = logivec.length();
  IntegerVector f2t = which(!logivec & shiftleft(logivec));
  IntegerVector t2f = which(logivec & !shiftleft(logivec));
  
  if (extend) t2f = modulo(t2f, L) + 1;
  else f2t = modulo(f2t, L) + 1;
  int Nseg = t2f.length();

  if (Nseg > 1 && logivec[0]) t2f = shiftleft(t2f);
  
  
  if (segments) {
    IntegerVector seg = t2f[0] <= f2t[0] ? cVecs(renumber(f2t[0], L-1), renumber(0, t2f[0])) : renumber(f2t[0], t2f[0]-1);
    
    if (firstonly || Nseg == 1){
      List ret(1);
      ret[0] = seg;
      return ret;
    }
    else {
      List ret(Nseg);
      ret[0] = seg;
      for (int i = 1; i < Nseg; i++)
        ret[i] = t2f[i] <= f2t[i] ? cVecs(renumber(f2t[i], L-1), renumber(0, t2f[i])) : renumber(f2t[i], t2f[i]-1);
      return ret;
    }
  }
  else {
    List ret(2);
    if(firstonly){
      ret[0] = f2t[0] - 1;
      ret[1] = t2f[0] - 1;
    }
    else {
      ret[0] = f2t - 1;
      ret[1] = t2f - 1;
    }
    return ret;
  } 
}

LogicalVector sllonlatidentical(NumericVector lon1, NumericVector lat1, NumericVector lon2, NumericVector lat2, const bool& recycle = false, const double& tolerance = 0){
  lon1 = cutoutVector(lon1, -180, 180);
  lon2 = cutoutVector(lon2, -180, 180);
  
  IntegerVector Ls(4);
  Ls[0] = lon1.length();
  Ls[1] = lat1.length();
  Ls[2] = lon2.length();
  Ls[3] = lat2.length();
  int L = max(Ls);
  if(is_true(any(Ls < L))){
    if(recycle){
      lon1 = rep_len(lon1, L);
      lat1 = rep_len(lat1, L);
      lon2 = rep_len(lon2, L);
      lat2 = rep_len(lat2, L);
    }
    else{
      Rcpp::stop("input sizes do not match, consider setting 'recycle=TRUE'");
    }
  }
  
  return abs(lat1-lat2)<=tolerance & (abs(lon1-lon2)<=tolerance | abs(lat1+90)<=tolerance | abs(lat1-90)<=tolerance);
}

VisShiftRot::VisShiftRot(const PIR& pir, const NumericVector& lon, const NumericVector& lat){
  NumericVector rotlon(0), rotlat(0), x(0), y(0);
  LogicalVector visible(0);
  
  if (pir.projection == "lonlat") {
    NumericVector lon1 = cutoutVector(lon, -180, 180);
    visible = (lat>=pir.latrange[0] & lat<=pir.latrange[1] & lon1>=pir.lonrange[0] & lon1<=pir.lonrange[1]);
    x = lon1;
    y = lat;
  }
  else Rcpp::stop("Other projections than lonlat not yet implemented!");
  
  this->x = x;
  this->y = y;
  this->rotlon = rotlon;
  this->rotlat = rotlat;
  this->visible = visible;
}


void slplotpolygon(const PIR& plotinitres, const NumericVector& longitude, const NumericVector& latitude, const Function& polygon, const String& colfill = "black", const String& colborder = "black", const double borderlwd = 0.01, const int& borderlty = 1, const bool& ignorevisibility = false, const bool& removeidenticalneighbours = true, const bool& refineboundary = true, const double& refineboundaryprecision = 1){
  PIR pir(plotinitres);
  NumericVector lon = longitude, lat = latitude;
  
  #ifdef DEBUG
  #pragma omp critical
  Rcout << "polygon: " << lon << ", " << lat << ", " << colfill.get_cstring() << ", " << colborder.get_cstring() << ", " << borderlwd << ", " << borderlty << ", " << ignorevisibility << ", " << removeidenticalneighbours << ", " << refineboundary << ", " << refineboundaryprecision << '\n';
  #endif
  //
  int L = lon.length();
  if(L != lat.length()) {
    Rcpp::warning("lon and lat vector do not have the same length!");
    return;
  }
  #ifdef DEBUG
  #pragma omp critical
  Rcout << "length: " << L << '\n';
  Rcout << removeidenticalneighbours << '\n';
  #endif
  
  if(removeidenticalneighbours){
    LogicalVector keep = !sllonlatidentical(lon, lat, shiftleft(lon), shiftleft(lat));
    lon = removeEntries(lon, keep);
    lat = removeEntries(lat, keep);
    L = lon.length();
  }
  
  #ifdef DEBUG
  #pragma omp critical
  Rcout << "length after 'rin': " << L << '\n';
  Rcout << "projection: " << pir.projection.get_cstring() << '\n';
  #endif
  
  if(pir.is3D() == 1){
    for(int npir = 0; npir < pir.pirlist.size(); npir++){
      slplotpolygon(pir.pirlist[npir], lon, lat, polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
      npir++;
    }
  }
  
  #ifdef DEBUG
  #pragma omp critical
  Rcout << "ignore visibility" << ignorevisibility << '\n';
  #endif
  
  VisShiftRot vsrres(pir, lon, lat);
  LogicalVector visible = vsrres.visible;
  if(ignorevisibility) visible.fill(true);
  int vissum = sum(visible);
  if(vissum == 0) return;
  
  #ifdef DEBUG
  #pragma omp critical
  Rcout << vissum << '\n';
  #endif
  
  NumericVector x = vsrres.x, y = vsrres.y, rotlon = vsrres.rotlon, rotlat = vsrres.rotlat;
  
  bool vispartial = vissum < L ? true : false;
  
  #ifdef DEBUG
  #pragma omp critical
  Rcout << "going into if ? " << (vispartial && pir.projection != "lonlat") << '\n';
  #endif
  
  if(vispartial && pir.projection != "lonlat"){
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
  
  #ifdef DEBUG
  #pragma omp critical
  Rcout << "length: " << L << '\n';
  #endif
  
  if(pir.projection == "lonlat"){
    if(vispartial){
      
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "vispartial" << '\n';
      #endif
      
      //point(s) out of west boundary
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "west" << '\n';
      #endif
      
      if(min(x) < pir.lonrange[0]){
        ListOf<IntegerVector> inds = slsegment(x>pir.lonrange[0], true);
        if(inds.size() > 1)
          for(int i = 1; i < inds.size(); i++)
            slplotpolygon(pir, x[inds[i]], y[inds[i]], polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        if(x[0] < pir.lonrange[0]){
          if(x[1]-x[0]>180){
            y[0] = y[1] + (pir.lonrange[1] - x[1])/(x[0] + 360 - x[1]) * (y[0] - y[1]);
            x[0] = pir.lonrange[1];
          }
          else{
            y[0] = y[1] + (pir.lonrange[0] - x[1])/(x[0] - x[1]) * (y[0] - y[1]);
            x[0] = pir.lonrange[0];
          }
        }
        if(x[L-1] < pir.lonrange[1]){
          if(x[L-2]-x[L-1]>180){
            y[L-1] = y[L-2] + (pir.lonrange[1] - x[L-2])/(x[L-1] + 360 - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = pir.lonrange[1];
          }
          else{
            y[L-1] = y[L-2] + (pir.lonrange[0] - x[L-2])/(x[L-1] - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = pir.lonrange[0];
          }
        }
      }
      
      //point(s) out of east boundary
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "east" << '\n';
      #endif
      
      if(max(x) > pir.lonrange[1]){
        ListOf<IntegerVector> inds = slsegment(x<pir.lonrange[1], true);
        if(inds.size() > 1)
          for(int i = 1; i < inds.size(); i++)
            slplotpolygon(pir, x[inds[i]], y[inds[i]], polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        if (x[0] > pir.lonrange[1]) {
          if (x[0]-x[1]>180) {
            y[0] = y[1] + (pir.lonrange[0] - x[1])/(x[0]-360 - x[1]) * (y[0] - y[1]);
            x[0] = pir.lonrange[0];
          } else {
            y[0] = y[1] + (pir.lonrange[1] - x[1])/(x[0] - x[1]) * (y[0] - y[1]);
            x[0] = pir.lonrange[0];
          }
        }
        if (x[L-1] > pir.lonrange[1]) {
          if (x[L-1]-x[L-2]>180) {
            y[L-1] = y[L-2] + (pir.lonrange[0] - x[L-2])/(x[L-1]-360 - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = pir.lonrange[0];
          } else {
            y[L-1] = y[L-2] + (pir.lonrange[1] - x[L-2])/(x[L-1] - x[L-2]) * (y[L-1] - y[L-2]);
            x[L-1] = pir.lonrange[1];
          }
        }
      }
      
      //point(s) out of south boundary
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "south" << '\n';
      #endif
      
      if(min(y) < pir.latrange[0]){
        ListOf<IntegerVector> inds = slsegment(y>pir.latrange[1], true);
        if(inds.size() > 1)
          for(int i = 1; i < inds.size(); i++)
            slplotpolygon(pir, x[inds[i]], y[inds[i]], polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        if (y[0] < pir.latrange[0]) {
          if (abs(x[0]-x[1])>180) {
            x[0] += x[0]-x[1]>180 ? -360 : 360;
            x[0] = x[1] + (pir.latrange[0] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            if (x[0]<pir.lonrange[0]) x[0] += 360;
            if (x[0]>pir.lonrange[1]) x[0] -= 360;
            y[0] = pir.latrange[0];
          } 
          else {
            x[0] = x[1] + (pir.latrange[0] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            y[0] = pir.latrange[0];
          }
        }
        if (y[L-1] < pir.latrange[0]) {
          if (abs(x[L-1]-x[L-2])>180) {
            x[L-1] += x[L-1]-x[L-2]>180 ? -360 : 360;
            x[L-1] = x[L-2] + (pir.latrange[0] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            if (x[L-1]<pir.lonrange[0]) x[L-1] += 360;
            if (x[L-1]>pir.lonrange[1]) x[L-1] -= 360;
            y[L-1] = pir.latrange[0];
          }
          else {
            x[L-1] = x[L-2] + (pir.latrange[0] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            y[L-1] = pir.latrange[0];
          }
        }
      }
      
      //point(s) out of north boundary
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "north" << '\n';
      #endif
      
      if(max(y) > pir.latrange[1]){
        ListOf<IntegerVector> inds = slsegment(y<pir.latrange[1], true);
        
        if(inds.size() > 1){
          for(int i = 1; i < inds.size(); i++){
            slplotpolygon(pir, x[inds[i]], y[inds[i]], polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
          }
        }
        #ifdef DEBUG
        #pragma omp critical
        Rcout << "loop finished\n";
        Rcout << inds.size() << ", " << inds[0] << '\n';
        #endif
        
        x = x[inds[0]];
        y = y[inds[0]];
        L = x.length();
        
        #ifdef DEBUG
        #pragma omp critical
        Rcout << "x.length: " << L << '\n';
        Rcout << (y[0] > pir.latrange[1]) << '\n';
        #endif
        if (y[0] > pir.latrange[1]) {
          if (abs(x[0]-x[1])>180) {
            x[0] += x[0]-x[1]>180 ? -360 : 360;
            x[0] = x[1] + (pir.latrange[1] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            if (x[0]<pir.lonrange[0]) x[1] += 360;
            if (x[0]>pir.lonrange[1]) x[1] -= 360;
            y[0] = pir.latrange[1];
          } 
          else {
            x[0] = x[1] + (pir.latrange[1] - y[1])/(y[0] - y[1]) * (x[0] - x[1]);
            y[0] = pir.latrange[1];
          }
        }
        #ifdef DEBUG
        #pragma omp critical
        Rcout << (y[L-1] > pir.latrange[1]) << '\n';
        #endif
        if (y[L-1] > pir.latrange[1]) {
          if (abs(x[L-1]-x[L-2])>180) {
            x[L-1] += x[L-1]-x[L-2]>180 ? -360 : 360;
            x[L-1] = x[L-1] + (pir.latrange[1] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            if (x[L-1]<pir.lonrange[0]) x[L-1] += 360;
            if (x[L-1]>pir.lonrange[1]) x[L-1] -= 360;
            y[L-1] = pir.latrange[1];
          } 
          else {
            x[L-1] = x[L-2] + (pir.latrange[1] - y[L-2])/(y[L-1] - y[L-2]) * (x[L-1] - x[L-2]);
            y[L-1] = pir.latrange[1];
          }
        }
      }
    }
    
    #ifdef DEBUG
    #pragma omp critical
    Rcout << "after vispartial" << '\n';
    #endif
    if(max(x) - min(x) > 180 && max(abs(x - shiftleft(x))) > 180){
      IntegerVector l2r = which(shiftleft(x) - x > 180), r2l = which(shiftleft(x)-x<-180);
      int Nlr = l2r.length(), Nrl = r2l.length();
      if(Nlr != Nrl) {
        #pragma omp critical
        Rcpp::warning("This nasty polygon can not be plotted; it might be circular, crossing the lonlat boundary an uneven number of times, that is, it may contain one or the other pole. Consider splitting the polygon into better behaving pieces.");
        return;
      }
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "Nlr > 1 ? " << (Nlr > 1) << '\n';
      #endif
      if(Nlr > 1){
        if(l2r[0] > r2l[0]) r2l = shiftleft(r2l);
        for (int i = 1; i< Nlr; i++) {
          IntegerVector right = r2l[i]%L+1 > l2r[i] ? renumber(l2r[i], r2l[i]%L+1) : cVecs(renumber(l2r[i], L), renumber(1, r2l[i]%L+1));
          IntegerVector left = l2r[i%Nlr]%L+1 > r2l[i] ? renumber(r2l[i], l2r[i%Nlr]%L+1) : cVecs(renumber(r2l[i], L), renumber(1, l2r[i%Nlr]%L+1));
          for(int o = 0; o < right.length(); o++)
            right[o] %= L;
          for(int o = 0; o < left.length(); o++)
            left[o] %= L;
          
          slplotpolygon(pir, x[left], y[left], polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
          slplotpolygon(pir, x[right], y[right], polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
        }
      }
      int i = 0;
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "right: " << (r2l[i]%L+1 > l2r[i]) << "\nleft: " << (l2r[i%Nlr]%L+1 > r2l[i]) << '\n';
      #endif
      IntegerVector right = r2l[i]%L+1 > l2r[i] ? renumber(l2r[i], r2l[i]%L+1) : cVecs(renumber(l2r[i], L), renumber(1, r2l[i]%L+1));
      IntegerVector left = l2r[i%Nlr]%L+1 > r2l[i] ? renumber(r2l[i], l2r[i%Nlr]%L+1) : cVecs(renumber(r2l[i], L), renumber(1, l2r[i%Nlr]%L+1));
      for(int o = 0; o < right.length(); o++)
        right[o] %= L;
      for(int o = 0; o < left.length(); o++)
        left[o] %= L;
      
      NumericVector xr = x[right], yr = y[right], xl = x[left], yl = y[left];
      int Lr = right.length(), Ll = left.length();
      
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "ranges1" << '\n';
      #endif
      if(xr[0] < pir.lonrange[1]) xr[0] += 360;
      yr[0] = yr[1] + (pir.lonrange[1] - xr[1])/(xr[0] - xr[1]) * (yr[0] - yr[1]);
      xr[0] = pir.lonrange[1];
      
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "ranges2" << '\n';
      #endif
      if(xr[Lr-1] < pir.lonrange[1]) xr[Lr-1] += 360;
      yr[Lr - 1] = yr[Lr-2] + (pir.lonrange[1] - xr[Lr-2])/(xr[Lr-1] - xr[Lr-2]) * (yr[Lr-1] - yr[Lr-2]);
      xr[Lr-1] = pir.lonrange[1];
      
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "ranges3" << '\n';
      #endif
      if (xl[0] > pir.lonrange[0]) xl[0] -= 360;
      yl[0] = yl[1] + (pir.lonrange[0] - xl[1])/(xl[0] - xl[1]) * (yl[0] - yl[1]);
      xl[0] = pir.lonrange[0];
      
      #ifdef DEBUG
      #pragma omp critical
      Rcout << "ranges4" << '\n';
      #endif
      if (xl[Ll - 1] > pir.lonrange[0]) xl[Ll-1] -= 360;
      yl[Ll - 1] = yl[Ll-2] + (pir.lonrange[0] - xl[Ll-2])/(xl[Ll-1] - xl[Ll-2]) * (yl[Ll-1] - yl[Ll-2]);
      xl[Ll-1] = pir.lonrange[0];
      
      #ifdef PLOTPOLYGONS
      #pragma omp critical
      {
        #ifdef DEBUG
        Rcout << "plot polygons" << '\n';
        #endif
        
        
        polygon(_["x"] = xr+pir.xshift, _["y"] = yr+pir.yshift, _["col"] = colfill, _["lwd"] = borderlwd, _["lty"] = borderlty, _["border"] = colborder);
        polygon(_["x"] = xl+pir.xshift, _["y"] = yl+pir.yshift, _["col"] = colfill, _["lwd"] = borderlwd, _["lty"] = borderlty, _["border"] = colborder);
      }
      #endif
      
      return;
    }
    else {
      #ifdef PLOTPOLYGONS
      #pragma omp critical
      {
        #ifdef DEBUG
        Rcout << "just plot" << '\n';
        #endif
        polygon(_["x"] = x+pir.xshift, _["y"] = y+pir.yshift, _["col"] = colfill, _["lwd"] = borderlwd, _["lty"] = borderlty, _["border"] = colborder);
      }
      #endif
      return;
    }
  }
  else{
    Rcpp::stop("Other projections than 'lonlat' not yet implemented!");
  }
}

// [[Rcpp::export("sl.plot.polygon")]]
void slplotpolygon(List pir, NumericVector lon, NumericVector lat, Function polygon, String colfill = "black", String colborder = "black", double borderlwd = 0.01, int borderlty = 1, bool ignorevisibility = false, bool removeidenticalneighbours = true, bool refineboundary = true, double refineboundaryprecision = 1){
  PIR pir2(pir);
  slplotpolygon(pir2, lon, lat, polygon, colfill, colborder, borderlwd, borderlty, ignorevisibility, removeidenticalneighbours, refineboundary, refineboundaryprecision);
}

// [[Rcpp::export(".sl.plot.field.loop")]]
void slplotfieldloop(List plotinitres, NumericMatrix lonv, NumericMatrix latv, List colbar, IntegerVector colind, String colfill = "black", String colborder = "black", double borderlwd = 0.01, int borderlty = 1, int threads = 1) {
  #ifdef DEBUG
  Rcout << "got into c++ for loop\n" << colind.length() << ", " << max(colind) << ", " << colbar.size() << '\n';
  #endif
  
  Environment env("package:graphics");
  Function polygon = as<Function>(env["polygon"]);
  PIR pir(plotinitres);
  //NumericVector lon, lat;

  #ifdef _OPENMP
  //setenv("OMP_STACKSIZE","100M",1);
  omp_set_num_threads(threads);
  #pragma omp parallel for schedule(dynamic) ordered//shared(lonv, latv)
  #endif
  for(int i = 0; i < lonv.rows(); i++){
    #ifdef DEBUG
    #ifdef _OPENMP
    #pragma omp critical
    Rcout << omp_get_thread_num() << ": " << i << '\n';
    #endif
    #ifndef _OPENMP
    #pragma omp critical
    Rcout << i << '\n';
    #endif
    #endif
    
    NumericVector lon, lat;
    
    #pragma omp critical
    lon = lonv(i, _);
    #pragma omp critical
    lat = latv(i, _);
    
    slplotpolygon(pir, lon, lat, polygon, colfill == "colbar" ? colbar[colind[i]-1] : colfill, colborder == "colbar" ? colbar[colind[i]-1] : colborder, borderlwd, borderlty);
  }
#pragma omp barrier
}