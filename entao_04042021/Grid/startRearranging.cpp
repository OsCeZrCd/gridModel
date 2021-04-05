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
  bool r_flag = false;  
  double nu = 0.0;  // probability of rearranging when below yeild strain
  int rand2;
  double iRand;
  double A, B, C, D;
  A = 40.05;
  B = -25.57;
  C = 19.74;
  D = -11.27;
  
  double E = 1.0;
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
    r_flag = true;
    s_induce_r += 1;
  }
  
  else if (s>1.75155){
    nu = 0.00876925;
    iRand = double(rand2) / 10000.0;
    if (nu > iRand){
      r_flag = true;
    }
    else{
      r_flag = false;
    }
  }

  else{
    nu = exp((A+B*s) - (C+D*s-omega*e2*E)/T);
    
    iRand = double (rand2) / 10000.0;
    
    if (nu > iRand){
      r_flag = true;
    }
    else{
      r_flag = false;
    }
  }

  return r_flag;

}
