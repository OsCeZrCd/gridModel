#include "gridModel.hpp"
#include <math.h>
// dsFromRearranger
// Function that determines how softness changes when a rearranger
// is a vector distance site_dx away
// Parameters
// site_dx : std::vector<double>
// vector with the displacement between the current grid point and
// the rearranger in real units (can be 2D or 3D)

#define PI acos(-1)

double gridModel::dsFromRearranger(std::vector<double> site_dx)
{
  double dyz, dx, r, dS;
  double theta, cos2theta;
  double a;    // 0.0042 for T=0.05, 0.025 for T = 0.3
  double b;       // -2.8 for T=0.05, 0.0 for T = 0.3
  double c;
  double alpha;  // -0.80 for T=0.05, -5.1 for T = 0.3
  double flag;

  double a1, b1, c1, alpha1;

  double s_comp=0;
  double scale_factor = 1.0;

  if (T==0.05){

    // entao's fitting, same threshold, 0.1
    
    a = 0.002347; //0.07975;
    b = 1.0;   //0.02935;
    c = 0.0;   //-0.0027;
    s_comp = -2.905081 / 4168.0; // 1162.0;
    alpha = 0.6483; //0.6516;

    a1 = 0.0240; //0.00603874;
    alpha1 = -1.263; //-0.4973;
    b1 = 1.0; //1.3606;
    c1 = 0.0; //-0.1953;
   
  }

  else if (T==0.30){
    a = 0.025;
    b = 0.0;
    c = 0.0;
    alpha = -5.3;
    s_comp = 0;
  }


  //double s_comp = 0;

  // Computes:
  // 1. distance in yz (or y) plane (dyz)
  // 2. distance in x plane (dx)
  // 3. Euclidean distance (r)
  // 4. theta b/t x axis and vector

  dyz = 0.0;
  
  dx = site_dx[0] * lGrid;
  double dy = site_dx[1] * lGrid;
  double dz = site_dx[2] * lGrid;
  
  double cos_45 = std::cos(0.25 * PI);
  double sin_45 = std::sin(0.25 * PI);

  dx = site_dx[0] * cos_45 - site_dx[1] * sin_45;
  dy = site_dx[1] * sin_45 + site_dx[0] * cos_45;

  dyz = dy * dy + dz * dz;
  dyz = std::sqrt(dyz) * lGrid;
  
  dx *= lGrid;

  r = std::sqrt(dyz*dyz + dx*dx);
  
  theta = std::atan2(dyz, dx);
  
  if (theta > 0)
    theta = 0.5 * PI - theta;
  else
    theta = 0.5 * PI + theta;

  cos2theta = cos(2.0*(theta));

  double cos4theta = std::cos(4.0*(theta));

  // entao's fitting
  if (T == 0.05 && r>0 && r <=3.0){
    dS = scale_factor*a*pow(r,alpha)*(b*cos2theta+c*cos4theta);
    dS += s_comp;
  }
  else if (T == 0.05 && r>3.0 && r <= 10.0){
    dS = scale_factor*a1*pow(r,alpha1)*(b1*cos2theta+c1*cos4theta);
    dS += s_comp;
  }

  else 
    dS = 0;

  return dS;
}
