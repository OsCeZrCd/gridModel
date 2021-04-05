#include "gridModel.hpp"
#include <math.h>
// dsFromRearranger
// Function that determines how softness changes when a rearranger
// is a vector distance site_dx away
// Parameters
// site_dx : std::vector<double>
// vector with the displacement between the current grid point and
// the rearranger in real units (can be 2D or 3D)

double gridModel::dsFromRearranger(std::vector<double> site_dx)
{
  double dxy, dz, r, dS;
  double theta, cos2theta;
  double a;    // 0.0042 for T=0.05, 0.025 for T = 0.3
  double b;       // -2.8 for T=0.05, 0.0 for T = 0.3
  double alpha;  // -0.80 for T=0.05, -5.1 for T = 0.3

  if (T==0.05){
    a = 0.0042;
    b = -2.8;
    alpha = -0.80;
  }

  else if (T==0.30){
    a = 0.025;
    b = 0.0;
    alpha = -5.3;
  }

  // Computes:
  // 1. distance in xy (or x) plane (dxy)
  // 2. distance in z (or y) plane (dz)
  // 3. Euclidean distance (r)
  // 4. theta b/t z (or y) axis and vector

  dxy = 0.0;

  for (int dd=0; dd<d-1; dd++)
    dxy += site_dx[dd] * site_dx[dd];
  
  dxy = std::sqrt(dxy);
  dz = site_dx[d-1];
  r = std::sqrt(dxy*dxy + dz*dz);
  theta = std::atan2(dz, dxy);

  // Computes dS
  cos2theta = std::cos(2.0*(theta-0.785398));  // subtracts off 45
                                               // degrees b/c SB angle

  if (r > 0 && r < 6)
//    dS = a*pow(r,alpha)*(1.0+b*cos2theta);   // previous version of angular term
    dS = a*pow(r,alpha)*(b*cos2theta);         // new angular term, only accounts for unisotropic part
  
  else 
    dS = 0;

  return dS;
}
