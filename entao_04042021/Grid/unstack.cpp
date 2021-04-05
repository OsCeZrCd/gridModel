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
    y = site/nGridPerSide;
    x = site - y*nGridPerSide;
    site_coor[0] = x, site_coor[1] = y;
  }    
  else
  {
    int x, y, z;
    z = site / (nGridPerSide*nGridPerSide);
    y = (site - z*nGridPerSide*nGridPerSide) / nGridPerSide;
    x = site - z*nGridPerSide*nGridPerSide - y*nGridPerSide;
    site_coor[0] = x, site_coor[1] = y, site_coor[2] = z;
  
  }
}

