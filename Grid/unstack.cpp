#include "gridModel.hpp"
// unstack function
// Takes a site index (site) and converts it to a site coordinate
// site_coor in (x, y), or (x, y, z) space
//
// Parameters:
// site: int, the index which runs from [0, nSite]
// site_coor: std::vector<int>
//   d dimensional vector that stores location of the site in the space.

void gridModel::unstack(int site, std::vector<int>& site_coor)
{
  if(d==2)
  {
    int x, y;
    y = site/nxGrid;
    x = site - y*nxGrid;
    site_coor[0] = x, site_coor[1] = y;
  }    
  else
  {
    int x, y, z;
    z = site / (nxGrid*nyGrid);
    y = (site - z*nxGrid*nyGrid) / nxGrid;
    x = site - z*nxGrid*nyGrid - y*nxGrid;
    site_coor[0] = x, site_coor[1] = y, site_coor[2] = z;
  
  }
}

