#include "gridModel.hpp"
// stack function
// Take a site coordinate (site_coor) in (i, j) or (i, j, k) space
// and converts it to a site index (site)
//
// Parameters:
// site_coor: std::vector<int>
//   d dimensional vector that is the location of the site in (i, j) or
//   (i, j, k) space.
// site: int, site index which runs from [0, nSite]

void gridModel::stack(std::vector<int> site_coor, int& site)
{
  if (d==2)
    site = site_coor[0] +site_coor[1]*nxGrid;
  else
    site = site_coor[0] + site_coor[1]*nxGrid 
      + site_coor[2]*nxGrid*nyGrid;
}
