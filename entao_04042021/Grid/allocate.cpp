#include "gridModel.hpp"

void gridModel::allocate(){

  alle.resize(nSite);
  alls.resize(nSite);
  alls_local.resize(nSite);
  yieldStrainPx.resize(nSite);
  hasRearranged.resize(nSite);
  rearrangingStep.resize(nSite);
  rearrangingIntensity.resize(nSite);
}
