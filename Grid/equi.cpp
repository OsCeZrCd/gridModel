#include "gridModel.hpp"

void gridModel::equi(){
//   #pragma omp parallel for schedule(static)     // design for parallel running
  for (int i=0; i<nSite; i++){
    this->hasRearranged[i] = 0;
    this->rearrangingStep[i] = 0;
  }
}
