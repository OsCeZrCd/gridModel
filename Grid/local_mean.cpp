#include "gridModel.hpp"

// shear function, shears all the sites in the model

void gridModel::local_mean(){
  int count;
//#pragma omp parallel for schedule(static)     // design for parallel running
  
  std::vector<int> site_coor(d), site_coor_prime(d);
  int site_prime;
  for (int site=0; site<nSite; site++){
    unstack(site, site_coor);
    alls_local[site] = 0.0;
    //count = 0;

    for (int x=0; x<5; x++){    
      site_coor_prime[0] = site_coor[0] - x + 2;
      for (int y=0; y<5; y++){
        site_coor_prime[1] = site_coor[1] - y + 2;
        for (int z=0; z<5; z++){
          site_coor_prime[2] = site_coor[2] - z + 2;
          
          if (site_coor_prime[0]<0)
            site_coor_prime[0] += nxGrid;
          else if (site_coor_prime[0] >= nxGrid)     
            site_coor_prime[0] -= nxGrid;
          if (site_coor_prime[1]<0)
            site_coor_prime[1] += nyGrid;
          else if (site_coor_prime[1] >= nyGrid)     
            site_coor_prime[1] -= nyGrid;
          if (site_coor_prime[2]<0)
            site_coor_prime[2] += nzGrid;
          else if (site_coor_prime[2] >= nzGrid)     
            site_coor_prime[2] -= nzGrid;
          
          stack(site_coor_prime, site_prime);
          alls_local[site] += alls[site_prime];
          //count += 1;
        }
      }
    }

    //alls_local[site] -= alls[site];
    alls_local[site] /= 125.0;
    //std::cout << count <<std::endl;
   // std::cout << alls_local[site] << std::endl;

  }
}
