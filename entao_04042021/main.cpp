#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>
#include <string>
#include <math.h>
#include "Grid/gridModel.hpp"
#include "GeometryVector/GeometryVector.h"
#include "Kiss_FFT/kiss_fft.h"
#include "Kiss_FFT/kiss_fftnd.h"

int main(){
  const int nGridPerSide = 20;        // num of grid points per side
  const double lenGridSide = 1.0;     // length of each grid side
  const int d = 3;                    // dimensionality of the space
  const double T = 0.30;
  int seed = 1;
  int rearrangeFrameLength = 1;
  double meanSoftness;                // for T=0.05, 
  double stdSoft;                     // sDistribution(0.0852, 0.7197)
                                      // for T=0.30
                                      // sDistribution(0.5751, 0.7049) 
  double c_soft;

  if (T==0.05){
    meanSoftness=0.0852;
    stdSoft=0.7197;
  }
  
  else if (T==0.30){
    meanSoftness=0.5751;
    stdSoft=0.7049;
  }
  
  c_soft = meanSoftness;

  std::normal_distribution<double> eDistribution(0.0, 0.00001);  // distribution of grid strains
  std::normal_distribution<double> sDistribution(meanSoftness, stdSoft);  // distribution of grid softness

  std::cout<<"distribution is:"<<sDistribution<<std::endl;
  gridModel model(nGridPerSide, lenGridSide, seed, d, T, meanSoftness, c_soft, rearrangeFrameLength, eDistribution, sDistribution);

  int numAvalanche = 0;
  int ii = 0;
  std::fstream strainFile("xyStrain.txt", std::fstream::out);
  std::fstream softFile("softness.txt", std::fstream::out);
  std::fstream sdistFile("softDist.txt", std::fstream::out);

  int equi_step = 0;

  int site;
  double totalExternalStrain = 0.0;
  while (totalExternalStrain < 8e-2){
    double sum=0.0;
    double count=0.0;

    for (auto &s : model.alls){
      sum += s;
      sdistFile << s;
      sdistFile << std::endl;
    }
    softFile << ii << ' '<<sum / model.alls.size();
    softFile << std::endl;

    double Strain = 6.6e-4;
    
    model.shear();
    totalExternalStrain += Strain;
    ii++;

    std::stringstream ss;
    ss << "yield_" << numAvalanche;

    bool yield = model.rearrange(ss.str());
    numAvalanche += yield;

    if (yield){
      std::cout << numAvalanche << "yielding so far.\n";
      if (ii%1==0){
        std::ostringstream ss1;
        ss1 << "xyRe_" << ii << ".cube";
        FILE *file = fopen(ss1.str().c_str(), "w");
//        FILE *file = fopen("xyRe_%d.cube", ii, "w");
        fprintf (file, "Grid Cube File\n");
        fprintf (file, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
        fprintf (file, "%d 0.0000 0.0000 0.0000\n", model.nSite);
        fprintf (file, "%d %lf 0.0000 0.0000\n", nGridPerSide, lenGridSide);
        fprintf (file, "%d 0.0000 %lf 0.0000\n", nGridPerSide, lenGridSide);
        fprintf (file, "%d 0.0000 0.0000 %lf\n", nGridPerSide, lenGridSide);

        int c_id = 1;
        for (int ix=0; ix<nGridPerSide; ix++){
          for (int iy=0; iy<nGridPerSide; iy++){
            for (int iz=0; iz<nGridPerSide; iz++){
              site = ix + iy*nGridPerSide + iz*nGridPerSide*nGridPerSide;
              c_id = (1+model.hasRearranged[site]);
              fprintf (file, "%d 0.0000 %d %d %d\n", c_id, ix, iy, iz);
            }
          }
        }

        for (int ix=0; ix<nGridPerSide; ix++){
          for (int iy=0; iy<nGridPerSide; iy++){
            for (int iz=0; iz<nGridPerSide; iz++){
              site = ix + iy*nGridPerSide + iz*nGridPerSide*nGridPerSide;
              fprintf (file, "%d ", int(model.hasRearranged[site]));
            }
            fprintf (file, "\n");
          }
        }

        fclose(file);
      }
    }

           
    auto outputStrainFunc = [&]() -> void {
      double sum = 0.0;
      for (auto s : model.alle){
        sum += s.x[1];  // xy direction
      }

      strainFile << totalExternalStrain << ' ' << sum / model.alle.size() << std::endl;
    };

    outputStrainFunc();
  
  }

}

