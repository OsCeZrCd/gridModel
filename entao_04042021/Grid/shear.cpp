#include "gridModel.hpp"

// shear function, shears all the sites in the model

void gridModel::shear(){

#pragma omp parallel for schedule(static)     // design for parallel running
  for (int i=0; i<nSite; i++){
    this->alle[i].x[1] += 6.6e-4;
    if (T==0.05)
      this->alls[i] += 0.00413;        // 0.00413 for T=0.05
                                       // 0.00207 for T = 0.3
    else if (T==0.30)
      this->alls[i] += 0.00207;

    this->hasRearranged[i] = 0;
    this->rearrangingStep[i] = 0;
  }
}
