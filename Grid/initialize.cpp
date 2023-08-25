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
      this->alle_buffer[site].x[strain] = alle[site].x[strain];
    }
    this->alls[site] = this->sDistribution(this->rEngine);
    alls_local[site] = alls[site];
    this->yieldStrainPx[site] = this->PxDistribution(this->rEngine);
    this->hasRearranged[site] = 0;
    this->rearrangingStep[site] = 0;
  }

  for (int j=0; j<9; j++){
    this->PR_bulk[j] = 0;
    this->S_bulk[j] = 0;
  }

  for (int i=0; i<5; i++)
    this->movingAverageTarget[i] = this->meanSoftness;
}

