#include "gridModel.hpp"
#include <math.h>
#include <time.h>
bool gridModel::rearrange(std::string outputPrefix /*= ""*/)
{
  bool rearrangeHappened = false;
  int nStep = 0;
  s_induce_r = 0;
  double multi_fac = 1.0;
  
  // for each site:
  // 1) checks if site is about to rearrange
  // 2) if rearranged, calculate the strain distributed to 
  // the around grids: need to give out a value of 
  // how much strain is redistributed (intensity)


  /***/
  double delta_E;
  for (int i=0; i<nSite; i++){
    
    delta_E = std::sqrt( (alle[i].Modulus2() + alle[i].x[0] * alle[i].x[3]) )
      - std::sqrt( (alle_buffer[i].Modulus2() + alle_buffer[i].x[0] * alle_buffer[i].x[3]) ); 
    
    if (T == 0.05)
      alls[i] += 16.60 * 0.37 * delta_E + 11.05 * 0.76 * delta_E*delta_E;

    else
      alls[i] += 9.17 * 0.33 * delta_E + 13.50 * 0.78 * delta_E*delta_E;

    for (int j=0; j<nStrain; j++){
      alle_buffer[i].x[j] = alle[i].x[j];
    }
  }

  int bufferCenter[3];
  bufferCenter[0] = bufferCenterX;
  bufferCenter[1] = bufferCenterY;
  bufferCenter[2] = bufferCenterZ;
  int nGridPerSide[3];
  nGridPerSide[0] = nxGrid;
  nGridPerSide[1] = nyGrid;
  nGridPerSide[2] = nzGrid;


  local_mean();

 
#pragma omp parallel
  {
#pragma omp for schedule(static)
    for (int site=0; site<nSite; site++){
      int flag = startRearranging(site, alle[site], alls[site]);
      if (flag > 0){
        hasRearranged[site] = 1;  
        rearrangingStep[site] = flag;
      }
    }
  }

  if (step < 45){
    int bin_id;
    for (int site=0; site<nSite; site++){
      bin_id = int ((alls[site]+3.0));
      if (bin_id < 0)
        bin_id = 0;
      else if (bin_id > 6)
        bin_id = 6;
      S_bulk[bin_id] += 1.0;
      if (rearrangingStep[site] > 0)
        PR_bulk[bin_id] += 1.0;
    }

    for (int i=0; i<7; i++){
      PR_bulk[i] = PR_bulk[i] / S_bulk[i];
    }
  }

#pragma omp parallel
  {
#pragma omp for schedule(static)
    for (int site=0; site<nSite; site++){
      if (rearrangingStep[site] == 1){   
        //nearby_PR(site, PR_bulk);
        induced_rearrangement(site);
      }
    }
  }
  
 
  std::vector<int> site_coor_1(d), site_coor_prime_1(d);
  std::vector<int> site_coor_shift_1(d);

  double local_rad[6];
  double num_count[6];

  for (int r_id=0; r_id<6; r_id++){
    local_rad[r_id] = 0.0;
    num_count[r_id] = 0.0;
  }

  for (int site_1=0; site_1<nSite; site_1++){
    if (rearrangingStep[site_1] == 1){
      unstack(site_1, site_coor_1);
      for (int site_prime_1=0; site_prime_1<nSite; site_prime_1++){
        unstack(site_prime_1, site_coor_prime_1);
        double dr_1 = 0.0;
        for (int dd=0; dd<d; dd++){
          int dx_1 = site_coor_prime_1[dd] - site_coor_1[dd];
          site_coor_shift_1[dd] = bufferCenter[dd] + dx_1;
          site_coor_shift_1[dd] = (site_coor_shift_1[dd] + nGridPerSide[dd]) % nGridPerSide[dd]; // PBC
          if (dx_1 > 0.5 * nGridPerSide[dd])
            dx_1 -= nGridPerSide[dd];              
          else if (dx_1 < -0.5 * nGridPerSide[dd])
            dx_1 += nGridPerSide[dd];
          dr_1 += double (dx_1 * dx_1);
        }
        int site_shift_1;
        stack(site_coor_shift_1, site_shift_1);
        dr_1 = std::sqrt(dr_1);

        if (site_1 == site_prime_1){
          local_rad[0] += alls[site_prime_1];
          num_count[0] += 1.0;
        }            
        else if (dr_1 < 2){
          local_rad[1] += alls[site_prime_1];
          num_count[1] += 1.0;
        }
        else if (dr_1 < 3){
          local_rad[2] += alls[site_prime_1];
          num_count[2] += 1.0;
        }         
        else if (dr_1 < 4){
          local_rad[3] += alls[site_prime_1];
          num_count[3] += 1.0;
        }              
        else if (dr_1 < 5){
          local_rad[4] += alls[site_prime_1];
          num_count[4] += 1.0;
        }    
        else{
          local_rad[5] = c_soft; 
          num_count[5] = 1.0;
        }  
      }
    }
  }       
   
  for (int r_id=0; r_id<6; r_id++){
    local_rad[r_id] = local_rad[r_id] / num_count[r_id];    
//    printf ("%lf %lf\n", local_rad[r_id], num_count[r_id]);  
  }

  for (int r_id=0; r_id<5; r_id++){
    movingAverageTarget[r_id] = 0.1 * local_rad[r_id] + 0.9 * movingAverageTarget[r_id];
//    printf ("%lf ", movingAverageTarget[r_id]);  
  }

#pragma omp parallel
  {
    srand(time(NULL));
    //std::normal_distribution<double> noiseDistribution(0.0, 1.0);
    std::mt19937 threadEngine;
    std::vector<int> site_coor(d), site_coor_prime(d);
    std::vector<int> site_coor_shift(d);
    #pragma omp critical(random)
    {
      threadEngine.seed(rEngine());
    } 
    int numRearrange;
    double num_Re_count = 0;
#pragma omp single
    {
      for (int site=0; site<nSite; site++){
          
        if (hasRearranged[site] == 1){
          if (d==2){
            GeometryVector residual(this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine));
            rearrangingIntensity[site] = (alle[site]-residual) * (1.0 / rearrangeFrameLength);  
          }

          else{
            GeometryVector residual(this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine),  this->residualStrainDistribution(this->rEngine),  this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine));
            rearrangingIntensity[site] = (alle[site]-residual) * (1.0 / rearrangeFrameLength); 
          }
          num_Re_count += 1;
        }
      }
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
             site_coor_shift[dd] = bufferCenter[dd] + dx;
             while (site_coor_shift[dd] < 0)
               site_coor_shift[dd] += nGridPerSide[dd];
             while (site_coor_shift[dd] >= nGridPerSide[dd])
               site_coor_shift[dd] -= nGridPerSide[dd];

             //site_coor_shift[dd] = (site_coor_shift[dd] + nGridPerSide) % nGridPerSide; // PBC
              
             if (dx >= 0.5 * nGridPerSide[dd])
               dx -= nGridPerSide[dd];
             else if (dx < -0.5 * nGridPerSide[dd])
               dx += nGridPerSide[dd];
             dr += double (dx * dx);
           } 
           int site_shift;
           stack(site_coor_shift, site_shift);

           for (int j=0; j<MaxDimension; j++){
             GeometryVector &de = dEBuffer[j][site_shift];
             e.AddFrom(rearrangingIntensity[site].x[j] * de);
           }

           if (rearrangingStep[site] == 1){
             std::vector<double> site_dr(d);
             site_dr[0] = site_coor_shift[0] * lGrid;
             site_dr[1] = site_coor_shift[1] * lGrid;
             site_dr[2] = site_coor_shift[2] * lGrid;

             double ds = 0.0;
             double restore = 0.0;
             dr = std::sqrt(dr);

             double alpha, beta, constant;
             if (T==0.05){

             //EXPONENTIAL eta(r) = e^(alpha * rdist + beta) + c
             alpha = -0.69;
             beta = -1.82;
             constant = 0.0115;

             }

             else if (T==0.30){

             // eta(r) = e^(alpha * rdist + beta) + c
             alpha = -0.978;
             beta = -1.575;
             constant = 0.0591;

             }
             
             double softnessRestoringCoefficient = exp(alpha * dr + beta) + constant;

           if (dr < 2){
             restore = softnessRestoringCoefficient * (local_rad[1] - alls[site_prime]);  
           }
           else if (dr < 3){
             restore = softnessRestoringCoefficient * (local_rad[2] - alls[site_prime]);
           }
           else if (dr < 4){
             restore = softnessRestoringCoefficient * (local_rad[3] - alls[site_prime]); 
           }
           else if (dr < 5){
             restore = softnessRestoringCoefficient * (local_rad[4] - alls[site_prime]);
           }
           
           else
             restore = softnessRestoringCoefficient * (local_rad[5] - alls[site_prime]);

             double harmonicDiffusion=0.0;
             double stddev;
             if (T==0.05)
               stddev = std::sqrt(softnessRestoringCoefficient * (2.0 -
                     softnessRestoringCoefficient) * 0.7197 * 0.7197);
             else
               stddev = std::sqrt(softnessRestoringCoefficient * (2.0 -
                     softnessRestoringCoefficient) * 0.7049 * 0.7049);

             std::normal_distribution<double> noiseDistribution(0.0, stddev);
             harmonicDiffusion = noiseDistribution(threadEngine);

             restore *= multi_fac;
             if (site != site_prime){
               alls[site_prime] += (ds + restore + harmonicDiffusion) * 1.0;
               
             }
           }
         }
           
         if (rearrangingStep[site] == 1)
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
        if (rearrangingStep[site] == 1){
          rearrangingStep[site] = 0;
          yieldStrainPx[site] = PxDistribution(rEngine);
          for (int na=0; na<nStrain; na++){
            alle[site].x[na] = 0;
          }
          double stddev;

          if (T==0.05){

            // Fixxed origin fitting
            alls[site] = alls[site] - 0.2210 * (alls[site] - local_rad[0]) * multi_fac;
            stddev = std::sqrt(0.2210 * (2.0 - 0.2210) * 0.7197 * 0.7197);

            std::normal_distribution<double> noiseDistribution(0.0, stddev);
            alls[site] += noiseDistribution(rEngine);
          }
          
          else{

            alls[site] = alls[site] - 0.2412 * (alls[site] - local_rad[0]);
            stddev = std::sqrt(0.2412 * (2.0 - 0.2412) * 0.7049 * 0.7049);
                        
            std::normal_distribution<double> noiseDistribution(0.0, stddev);
            alls[site] += noiseDistribution(rEngine);
          }
        }
      }
      

    // If there has been 1 or more rearrangements, we say an "yielding" has occured
   
      if (numRearrange > 0){
        rearrangeHappened = true;
        
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
        //std::cout << ", mean s=" << sum / alls.size();
        
        std::cout << sum / alls.size();
        std::cout << std::endl;
        c_soft = (sum/alls.size());

        std::cout << PR_bulk[0] << " " << PR_bulk[1]<< " " << PR_bulk[2] << " "
          << PR_bulk[3] << " " << PR_bulk[4]<< " " << PR_bulk[5] << " "
          << PR_bulk[6] << std::endl; 

      }
    }
  
  }

  step++; 
//  printf ("%d\n", step);
//  fflush(stdout);
  
  return rearrangeHappened;
}









