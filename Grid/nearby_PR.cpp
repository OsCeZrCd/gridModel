#include "gridModel.hpp"
#include <math.h>
#include <time.h>

void gridModel::nearby_PR(int site, std::vector<double> PR_bulk){
      
  std::vector<int> site_coor(d), site_coor_shift(d);
  std::vector<int> site_coor_prime(d);
  
  unstack(site, site_coor);

  int site_shift;
  double dr, prob;
  int rand2;
  double iRand;
  double soft;


#pragma omp parallel for schedule(static)
  for (int i=-5; i<5; i++){
    for (int j=-5; j<5; j++){
      for (int k=-5; k<5; k++){
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
        soft = alls[site_shift];

        dr = std::sqrt(i*i + j*j + k*k);

        double alpha, beta;
        double lambda, c;
        double b1, b0;
        int bin_id;
        
        bin_id = int ((soft+3.0));
        if (bin_id<0)
          bin_id = 0;
        else if (bin_id>6)
          bin_id = 6; 

        if (T==0.05){ 
          alpha = exp((soft / -2.25) + 7.90);
          beta = 0.13 * soft - 1.91;
          prob = 1.0 * (alpha * exp(beta * dr) + 1) * PR_1[bin_id];//PR_bulk[bin_id]; //exp(0.83*soft-6.19);
        }

        else if (T==0.30){
          alpha = exp((soft / -2.92) + 3.94);
          beta = 0.10 * soft - 1.47;
          prob = 1.0 * (alpha * exp(beta * dr) + 1) * PR_2[bin_id];//PR_bulk[bin_id];//exp(0.65*soft-5.94); //PR_bulk[bin_id]; //exp(1.15*soft-6.57);

        }
  
        rand2 = rand();
        rand2 %= 10000;
        iRand = double(rand2) / 10000.0;

        //printf ("%lf %lf  ", iRand, prob);
        //fflush(stdout);
        
        if (iRand < prob && hasRearranged[site_shift]==0){
          hasRearranged[site_shift] = 0;
          rearrangingStep[site_shift] = 2;
        }
        
      }     
    }
  }
}
      
