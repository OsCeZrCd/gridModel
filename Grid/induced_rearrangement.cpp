#include "gridModel.hpp"
#include <math.h>
#include <time.h>

void gridModel::induced_rearrangement(int site){
      
  std::vector<int> site_coor(d), site_coor_shift(d);
  std::vector<int> site_coor_prime(d);
  
  unstack(site, site_coor);

  int site_shift;
  double dr, prob;
  int rand2;
  double iRand;

  for (int i=-3; i<3; i++){
    for (int j=-3; j<3; j++){
      for (int k=-3; k<3; k++){
        site_coor_shift[0] = site_coor[0] + i;
        site_coor_shift[1] = site_coor[1] + j;
        site_coor_shift[2] = site_coor[2] + k; 
         
        if (site_coor_shift[0]<0)
          site_coor_shift[0] += nxGrid;
        else if (site_coor_shift[0] >= nxGrid)    
          site_coor_shift[0] -= nxGrid;
        if (site_coor_shift[1]<0)
          site_coor_shift[1] += nyGrid;
        else if (site_coor_shift[1] >= nyGrid)     
          site_coor_shift[1] -= nyGrid;
        if (site_coor_shift[2]<0)
          site_coor_shift[2] += nzGrid;
        else if (site_coor_shift[2] >= nzGrid)     
          site_coor_shift[2] -= nzGrid;
        
        stack(site_coor_shift, site_shift);

        dr = std::sqrt(i*i + j*j + k*k);

        if (T==0.05){
          
          prob = 0.3 * exp(-dr/1.114);
        }
        else if (T==0.30){
          
          prob = 0.3*exp(-dr/0.748);
        }

        prob = 0.5 * prob;

        rand2 = rand();
        rand2 %= 10000;
        iRand = double(rand2) / 10000.0;
        
        if (iRand < prob && hasRearranged[site_shift]==0){
          hasRearranged[site_shift] = 0;
          rearrangingStep[site_shift] = 2;
        }
        
      }     
    }
  }
}
      
