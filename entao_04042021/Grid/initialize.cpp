#include "gridModel.hpp"

// Initializes initial strains and softnesses for model according to
// given distribution.
void gridModel::initialize()
{
  for (int site=0; site<nSite; site++)
  {
    for (int strain=0; strain<nStrain; strain++)  // nstrain is 2 for 2D,
                                                  // 5 for 3D
    {
      this->alle[site].x[strain] = this->eDistribution(this->rEngine);
    }
    this->alls[site] = this->sDistribution(this->rEngine);
    this->yieldStrainPx[site] = this->PxDistribution(this->rEngine);
    this->hasRearranged[site] = 0;
    this->rearrangingStep[site] = 0;
//    alls_local[site] = alls[site];      designed for restoring to local mean

  }
}

