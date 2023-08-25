#include "gridModel.hpp"

// shear function, shears all the sites in the model

void gridModel::shear(){

#pragma omp parallel for schedule(static)     // design for parallel running
  for (int i=0; i<nSite; i++){
    
    this->alle[i].x[0] += (1.65e-4 * 1);
    this->alle[i].x[1] += (4.95e-4 * 1);
    this->alle[i].x[3] += (1.65e-4 * 1);

    if (T==0.05)
      this->alls[i] += 0.00413;
    else if (T==0.30)
      this->alls[i] += 0.00207;

    this->hasRearranged[i] = 0;
    this->rearrangingStep[i] = 0;
  }
}
