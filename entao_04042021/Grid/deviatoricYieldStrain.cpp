#include "gridModel.hpp"

double gridModel::deviatoricYieldStrain(int site){

  double s = alls[site];
  double meanYieldStrain;
  double k;

  if (T==0.05){
    meanYieldStrain = 0.0472 - 0.0045 * s;    
    k=1.84;
    if (s>2.0)
      meanYieldStrain = 0.0382;
  }
  else if (T==0.30){
    k = 2.02;
    if (s<=1.0)
      meanYieldStrain = -0.0102*s+0.0484;
    else
      meanYieldStrain = 0.0382;
  }                                         // both measured from MD simulation

  else{
    std::cout<<"Sorry, we do not support this temperature now."<<std::endl;
    exit(1);
  }

  //weibull distribution
  double lambda = meanYieldStrain / std::tgamma(1.0 + 1.0 / k);

  double yieldStrain = std::pow(-1.0 * std::log(1.0 - yieldStrainPx[site]), 1.0 / k) * lambda;
  return yieldStrain;
}
