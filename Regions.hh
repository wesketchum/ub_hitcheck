typedef enum ShortedRegion{
  kNominal,
  kShortedU,
  kShortedY,
  kUnknown  
} ShortedRegion_t;

bool isNominal(double y, double z){
  if( (z<54 || z>57) &&
      (z<93 || z>96) &&
      (z<101 || z>104) &&
      (z<246 || z>249) &&
      (z<290 || z>293) &&
      (z<347 || z>350) &&
      (z<399 || z>402) &&
      (z<414 || z>417) &&
      (z<700.9 || z>739.3) &&
      (z<808 || z>811) &&
      (z<821 || z>824) &&
      (z<875 || z>878) &&
      (y>z*0.577+14.77 || y<0.577*z-115.42 ) && 
      ( (y<-0.58*z+157 || y>-0.58*z+160) &&
	(y<-0.58*z+230 || y>-0.58*z+233) &&
	(y<-0.58*z+419.5 || y>-0.58*z+436) &&
	(y<-0.58*z+582 || y>-0.58*z+591) &&
	(y<-0.58*z+615 || y>-0.58*z+622) ) )
    return true;

  return false;
}

bool isShortedU(double y, double z){
  if( (y<z*0.577+14.77 && y>0.577*z+14.42)   ||
      (y<z*0.577+14.076 && y>z*0.577+7.84)   ||
      (y<z*0.577+7.49 && y>z*0.577+7.15)     ||
      (y<z*0.577+6.80 && y>z*0.577+3.68)     ||
      (y<z*0.577+0.91 && y>z*0.577+0.22)     ||
      (y<z*0.577-1.51 && y>z*0.577-2.55)     ||
      (y<z*0.577-3.25 && y>z*0.577-4.63)     ||
      (y<z*0.577-12.94 && y>z*0.577-21.60)   ||
      (y<z*0.577-24.72 && y>z*0.577-37.19)   ||
      (y<z*0.577-37.54 && y>z*0.577-50.70)   ||
      (y<z*0.577-56.25 && y>z*0.577-57.28)   ||
      (y<z*0.577-57.63 && y>z*0.577-63.17)   ||
      (y<z*0.577-63.52 && y>z*0.577-64.56)   ||
      (y<z*0.577-68.37 && y>z*0.577-76.68)   ||
      (y<z*0.577-77.03 && y>z*0.577-88.12)   ||
      (y<z*0.577-88.81 && y>z*0.577-90.19)   ||
      (y<z*0.577-90.54 && y>z*0.577-101.97)  ||
      (y<z*0.577-102.32 && y>z*0.577-108.9)  ||
      (y<z*0.577-109.25 && y>z*0.577-109.59) ||
      (y<z*0.577-109.93 && y>z*0.577-115.42) )
    return true;
  
  return false;
}

bool isShortedY(double, double z){
  if((z>700.9 && z<720.1) || (z>720.4 && z<724.6) || (z>724.9 && z<739.3))
    return true;
  
  return false;
}

ShortedRegion_t GetShortedRegionType(double y, double z){

  if(isNominal(y,z))
    return ShortedRegion_t::kNominal;
  if(isShortedU(y,z))
    return ShortedRegion_t::kShortedU;
  if(isShortedY(y,z))
    return ShortedRegion_t::kShortedY;
    
  return ShortedRegion_t::kUnknown;
}
