#include "gridModel.hpp"

void gridModel::allocate(){

  alle_buffer.resize(nSite);
  
  alle.resize(nSite);
  alls.resize(nSite);
  yieldStrainPx.resize(nSite);
  hasRearranged.resize(nSite);
  rearrangingStep.resize(nSite);
  
  alls_local.resize(nSite);
  
  rearrangingIntensity.resize(nSite);

  movingAverageTarget.resize(5);

  PR_bulk.resize(10);
  S_bulk.resize(10);

}
