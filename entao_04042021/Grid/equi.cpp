#include "gridModel.hpp"

//  equi function, run without shearing for all the sites in the model
//  designed for test only, not used in Grid model simulation

void gridModel::equi(){
//   #pragma omp parallel for schedule(static)     // design for parallel running
  for (int i=0; i<nSite; i++){
    this->hasRearranged[i] = 0;
    this->rearrangingStep[i] = 0;
  }
}
