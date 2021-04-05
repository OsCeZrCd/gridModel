#include "gridModel.hpp"

// local_mean function, calculate the local mean softness of each site, average over -2 to +2

void gridModel::local_mean(){
  int count;
  
  std::vector<int> site_coor(d), site_coor_prime(d);
  int site_prime;
  for (int site=0; site<nSite; site++){
    unstack(site, site_coor);
    alls_local[site] = 0.0;

    for (int x=0; x<5; x++){    
      site_coor_prime[0] = site_coor[0] - x + 2;
      for (int y=0; y<5; y++){
        site_coor_prime[1] = site_coor[1] - y + 2;
        for (int z=0; z<5; z++){
          site_coor_prime[2] = site_coor[2] - z + 2;
          for (int k=0; k<3; k++){
            if(site_coor_prime[k] < 0)
              site_coor_prime[k] += nGridPerSide;
            else if (site_coor_prime[k] >= nGridPerSide)
              site_coor_prime[k] -= nGridPerSide;
          }
          stack(site_coor_prime, site_prime);
          alls_local[site] += alls[site_prime];
        }
      }
    }

    alls_local[site] /= 125.0;

  }
}
