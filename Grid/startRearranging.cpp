#include "gridModel.hpp"
#include <math.h>
#include <time.h>

// startRearranging function
//
// this is the function that determine when a grid point is going to start
// rearranging, from the grid points' strain and softness.
//
// Parmameters:
// e: in type of GeometryVectors, strain vector representing strain at the grid point.
// s: double, softness at the given grid point.

bool gridModel::startRearranging(int site, GeometryVector e, double s){
  double e2;
  
  double yieldStrain = deviatoricYieldStrain(site);
  //bool r_flag = false;  
  int r_flag = 0;
  double nu = 0.0;  // probability of rearranging when below yeild strain
  int rand2;
  double iRand;
  double A, B, C, D;
  A = 40.05;
  B = -25.57;
  C = 19.74;
  D = -11.27;
  
  double E = 1.0;
  double delta_F;
  double v0_G;
  if (T==0.05)
    E = 74.58;
  else if (T==0.30)
    E = 29.13;

  rand2 = rand();
  rand2 %= 10000; 
  // Calculate sum of strain components
  if (d==2)
    e2 = e.Modulus2();
  else
    e2 = e.Modulus2() + e.x[0]*e.x[3];

  if (e2 > (yieldStrain * yieldStrain)){
    r_flag = 1;
    s_induce_r += 1;
    //std::cout << yieldStrain <<std::endl;
  }

  else if (s>1.75155){
    if (T==0.05)
      nu = 0.00876925 * 1.0;
    else
      nu = 0.00876925 * 1.0;
    iRand = double(rand2) / 10000.0;
    if (nu > iRand){
      r_flag = 1;
    }
    else{
      r_flag = 0;
    }
  }

  else{

    delta_F = C + D * s - T * (A + B*s);
    
    v0_G = delta_F / (yieldStrain * yieldStrain);

    nu = exp(-(delta_F - v0_G*e2)/T);

    iRand = double (rand2) / 10000.0;
    
    if (nu > iRand){
      r_flag = 1;
    }
    else{
      r_flag = 0;
    }
  }

  return r_flag;

}
