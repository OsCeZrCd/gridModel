#include "gridModel.hpp"
#include <math.h>
#include <time.h>
bool gridModel::rearrange(std::string outputPrefix /*= ""*/)
{
  bool rearrangeHappened = false;
  int nStep = 0;
  s_induce_r = 0;

  // for each site:
  // 1) checks if site is about to rearrange
  // 2) if rearranged, calculate the strain distributed to 
  // the around grids: need to give out a value of 
  // how much strain is redistributed (intensity),
  // set to be 1 for current version

//  local_mean();   // calculate local mean softness, needed for restoring

#pragma omp parallel
  {
    srand(time(NULL));
    std::normal_distribution<double> noiseDistribution(0.0, 1.0);
    std::mt19937 threadEngine;
    std::vector<int> site_coor(d), site_coor_prime(d);
    std::vector<int> site_coor_shift(d);
    #pragma omp critical(random)
    {
      threadEngine.seed(rEngine());
    }
 
    int numRearrange;
    double local_rad[15];
    double num_count[15];

#pragma omp single
    {
      for (int site=0; site<nSite; site++){
        if (startRearranging(site, alle[site], alls[site])){   
          hasRearranged[site] = 1;  
          rearrangingStep[site] = 1;
          if (d==2){
            GeometryVector residual(this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine));
            rearrangingIntensity[site] = (alle[site]-residual) * (1.0 / rearrangeFrameLength);  
          }

          else{
            GeometryVector residual(this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine),  this->residualStrainDistribution(this->rEngine),  this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine));
            rearrangingIntensity[site] = (alle[site]-residual) * (1.0 / rearrangeFrameLength); 
          }
        }
      }
    
    }

    for (int r_id=0; r_id<15; r_id++){
      local_rad[r_id] = 0.0;
      num_count[r_id] = 0.0;
    }

    for (int site=0; site<nSite; site++){
      if (rearrangingStep[site] > 0){
        unstack(site, site_coor);

        for (int site_prime=0; site_prime<nSite; site_prime++){
          GeometryVector &e = alle[site_prime];

          unstack(site_prime, site_coor_prime);
          double dr = 0.0;
          for (int dd=0; dd<d; dd++){
            int dx = site_coor_prime[dd] - site_coor[dd];
            site_coor_shift[dd] = bufferCenter + dx;
            site_coor_shift[dd] = (site_coor_shift[dd] + nGridPerSide) % nGridPerSide; // PBC
            dr += double (dx * dx);
          }
          int site_shift;
          stack(site_coor_shift, site_shift);

          dr = std::sqrt(dr);
          
          //int bin_flag;
          //if (site == site_prime){
          //  local_rad[0] += alls[site_prime];
          //  num_count[0] += 1.0;
          //}
          //else{
          //  bin_flag = int (dr);
          //  if (bin_flag < 14){
          //    local_rad[bin_flag] += alls[site_prime];
          //    num_count[bin_flag] += 1.0;
          //  }
          //  else{
          //    local_rad[14] += alls[site_prime];
          //    num_count[14] += 1.0;
          //  }
          //}

          if (site == site_prime){
            local_rad[6] += alls[site_prime];
            num_count[6] += 1.0;
          }
          else if (dr < 2){
            local_rad[0] += alls[site_prime];
            num_count[0] += 1.0;
          }
          else if (dr < 3){
            local_rad[1] += alls[site_prime];
            num_count[1] += 1.0;
          }
          else if (dr < 4){
            local_rad[2] += alls[site_prime];
            num_count[2] += 1.0;
          }          
          else if (dr < 5){
            local_rad[3] += alls[site_prime];
            num_count[3] += 1.0;
          }
         // else if (dr < 6){
         //   local_rad[4] += alls[site_prime];
         //   num_count[4] += 1.0;
         // }
         // else if (dr < 7){
         //   local_rad[5] += alls[site_prime];
         //   num_count[5] += 1.0;
         // }
          else if (dr < 6){
            local_rad[4] += alls[site_prime];
            num_count[4] += 1.0;
          }
          else{
            local_rad[5] = c_soft; 
            num_count[5] = 1.0;
          }
        }
      }
    }
          
      for (int r_id=0; r_id<7; r_id++){
        local_rad[r_id] = local_rad[r_id] / num_count[r_id];     
      }

#pragma omp barrier

  //rearrangement affect other sites parameters

  // Updates:
  // 1) strains and softness around rearrangements
  // 2) number of rearrangements occuring currently
  // 3) number of steps a rearrangment has taken
 
   numRearrange = 0;
 
   for (int site=0; site<nSite; site++){
      if (rearrangingStep[site] > 0){
        unstack(site, site_coor);
#pragma omp for schedule(static)

        for (int site_prime=0; site_prime<nSite; site_prime++){
          GeometryVector &e = alle[site_prime];

          unstack(site_prime, site_coor_prime);
          double dr = 0.0;
          for (int dd=0; dd<d; dd++){
            int dx = site_coor_prime[dd] - site_coor[dd];
            site_coor_shift[dd] = bufferCenter + dx;
            site_coor_shift[dd] = (site_coor_shift[dd] + nGridPerSide) % nGridPerSide; // PBC
            dr += double (dx * dx);
          }
          int site_shift;
          stack(site_coor_shift, site_shift);

          for (int j=0; j<MaxDimension; j++){
              GeometryVector &de = dEBuffer[j][site_shift];
              e.AddFrom(rearrangingIntensity[site].x[j] * de);
          }

          double ds = dSBuffer[site_shift];
          double restore = 0.0;
          dr = std::sqrt(dr);

          double alpha, beta, cont;
          if (T==0.05){
            alpha = 0.1;   
            beta = -1.175;
            cont = 0;          // fitting results without constant

//            alpha = 0.112398;
//            beta = -1.59836;
//            cont = 0.006115;
                               // fitting results with constant            
          }

          else if (T==0.30){
            alpha = 0.165;
            beta = -0.237;
            cont = 0;           // fitting results without constant

//            alpha = 0.116996;
//            beta = -1.42444;
//            cont = 0.092151;
                                // fitting results with constant          
          }

          double softnessRestoringCoefficient = alpha * std::pow(dr, beta) + cont;
          // int bin_flag;
          // bin_flag = int(dr);
          //if (dr < 14){
          //  restore = softnessRestoringCoefficient * (local_rad[bin_flag] - alls[site_prime]);
          //}
          //else
          //  restore = softnessRestoringCoefficient * (local_rad[14] - alls[site_prime]);

          if (dr < 2){
            restore = softnessRestoringCoefficient * (local_rad[0] - alls[site_prime]); 
          }
          else if (dr < 3){
            restore = softnessRestoringCoefficient * (local_rad[1] - alls[site_prime]);
          }
          else if (dr < 4){
            restore = softnessRestoringCoefficient * (local_rad[2] - alls[site_prime]);
          }
          else if (dr < 5){
            restore = softnessRestoringCoefficient * (local_rad[3] - alls[site_prime]);
          }
          else if (dr < 6)
            restore = softnessRestoringCoefficient * (local_rad[4] - alls[site_prime]);
          else
            restore = softnessRestoringCoefficient * (local_rad[5] - alls[site_prime]);

          double harmonicDiffusion=0.0;
          if (T==0.05)
            harmonicDiffusion = noiseDistribution(threadEngine) * 
                                (sqrt(softnessRestoringCoefficient*(2-softnessRestoringCoefficient))*0.7197);
          
          else
            harmonicDiffusion = noiseDistribution(threadEngine) * 
                                (sqrt(softnessRestoringCoefficient*(2-softnessRestoringCoefficient))*0.7049);

          if (site != site_prime){
            alls[site_prime] += ds + restore + harmonicDiffusion;
          }
        }
        numRearrange++;
      }
    }
#pragma omp single
  
    {

  // Updates:
  // rearrangement site itself. This must be done separately 
  // from above to ensure rearranger strain --> 0
  // Why we want the rearanger strain --> 0?

      
      for (int site=0; site<nSite; site++){
        if (rearrangingStep[site] > 0){
          rearrangingStep[site] = 0;
          yieldStrainPx[site] = PxDistribution(rEngine);
         
          if (T==0.05){
            alls[site] = alls[site] - 0.2186 * (alls[site] - local_rad[6]);
            alls[site] += noiseDistribution(rEngine) * (sqrt(0.2186*(2-0.2186)) * 0.7197);
          }
          
          else{
            alls[site] = alls[site] - 0.2869 * (alls[site] - local_rad[6]);
            alls[site] += noiseDistribution(rEngine) * (sqrt(0.2869*(2-0.2869)) * 0.7049);
          }
        }
      }
      
    // If there has been 1 or more rearrangements, we say an "yielding" has occured
   
      if (numRearrange > 0){
        rearrangeHappened = true;
      // Outputs the number of rearrangements
        std::cout << s_induce_r << " " << numRearrange << " ";
      // Outputs the mean energy
        double sum = 0.0;
        for (auto &e : this->alle){
          if(d==2)
            sum += e.Modulus2();
          else
            sum += e.Modulus2() + e.x[0] * e.x[3];
        }
        std::cout << sum / alle.size() << " ";

      // Outputs the mean softness
        sum = 0.0;
        for (auto &s : this->alls){
          sum += s;
        }
        std::cout << sum / alls.size();
        std::cout << std::endl;
        c_soft = (sum/alls.size());

      }
    }
  }
 
  return rearrangeHappened;
}









