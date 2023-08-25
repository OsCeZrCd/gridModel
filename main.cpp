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
  const int nxGrid = 21;        // num of grid points per side
  const int nyGrid = 21;
  const int nzGrid = 21;
  const double lenGridSide = 1.0;    // length of each grid side
  const int d = 3;                   // dimensionality of the space
  const int nAvalanche = 100;        
  const double T = 0.05;
  int seed = 1;
  int rearrangeFrameLength = 1;
  double meanSoftness;        // for T=0.05, sDistribution(0.0852, 0.7197)
  double stdSoft;             // for T=0.3, sDistribution(0.5751, 0.7049)

  double c_soft;
  int step;

  double strain_para;

  if (T==0.05){
    meanSoftness=0.0852;
    stdSoft=0.7197;
    
    strain_para = 1.0;
  }
  
  else if (T==0.30){
    meanSoftness=0.5751;
    stdSoft=0.7049;
    
    strain_para = 1.0;
  }
  
  c_soft = meanSoftness;
  step = 0;

  std::normal_distribution<double> eDistribution(0.0, 0.00001);  // distribution of grid strains
  std::normal_distribution<double> sDistribution(meanSoftness, stdSoft);  // distribution of grid softness

  std::cout<<"distribution is:"<<sDistribution<<std::endl;
  gridModel model(nxGrid, nyGrid, nzGrid, lenGridSide, seed, d, T, meanSoftness, step, c_soft, rearrangeFrameLength, eDistribution, sDistribution);

  int numAvalanche = 0;
  int ii = 0;
  std::fstream strainFile("xyStrain.txt", std::fstream::out);
  std::fstream softFile("softness.txt", std::fstream::out);
  std::fstream sdistFile("softDist.txt", std::fstream::out);

  std::fstream dumpFile("dump.txt", std::fstream::out);

  int ds_count = 0;
  for (int i=0; i<model.nSite; i++){
    if (model.dSBuffer[i] != 0){
      ds_count += 1;
    }
  }
  std::cout<< "ds_count is: " <<  ds_count << std::endl;

  int equi_step = 0;

  int site;
  double totalExternalStrain = 0.0;
  while (totalExternalStrain < 7e-2){
    double sum=0.0;
    double count=0.0;

    for (auto &s : model.alls){
      sum += s;
      sdistFile << s;
      sdistFile << std::endl;
    }
    softFile << ii << ' '<<sum / model.alls.size();
    softFile << std::endl;

    double Strain = 6.6e-4 * 1.0;
    
    model.shear();
//    model.equi();
    totalExternalStrain += Strain;
    ii++;
  //  step++;

    std::stringstream ss;
    ss << "yield_" << numAvalanche;

    bool yield = model.rearrange(ss.str());
    int type;
    numAvalanche += yield; 

    dumpFile << "ITEM: TIMESTEP" << std::endl;
    dumpFile << ii << std::endl;
    dumpFile << "ITEM: NUMBER OF ATOMS" << std::endl;
    dumpFile << model.nSite << std::endl;
    dumpFile << "ITEM: BOX BOUNDS pp pp pp" << std::endl;
    dumpFile << "0 " << nxGrid << std::endl;
    dumpFile << "0 " << nyGrid << std::endl;
    dumpFile << "0 " << nzGrid << std::endl;
    dumpFile << "ITEM: ATOMS id type x y z soft e_xx e_xy e_xz e_yy e_yz et" << std::endl;

    for (int iz=0; iz<nzGrid; iz++){
      for (int iy=0; iy<nyGrid; iy++){
        for (int ix=0; ix<nxGrid; ix++){
          site = iz * nyGrid * nxGrid + iy * nxGrid + ix;
          if (int(model.hasRearranged[site]) > 0)
            type = 1;
          else
            type = int(model.rearrangingStep[site]); 
          dumpFile << site + 1 << ' ' 
            << type << ' '
            << ix << ' ' << iy << ' ' << iz << ' ' 
            << model.alls[site] <<' ' << model.alle[site].x[0] << ' '<< model.alle[site].x[1] << ' '
            << model.alle[site].x[2] << ' ' << model.alle[site].x[3] << ' ' << model.alle[site].x[4] << ' '
            << model.alle[site].Modulus2() + model.alle[site].x[0]*model.alle[site].x[3] << std::endl;
        }
      }
    }

    if (yield){
      std::cout << numAvalanche << "yielding so far.\n";
      if (ii%1==0){
        std::ostringstream ss1;
        ss1 << "xyRe_" << ii << ".cube";
        FILE *file = fopen(ss1.str().c_str(), "w");
        fprintf (file, "Grid Cube File\n");
        fprintf (file, "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
        fprintf (file, "%d 0.0000 0.0000 0.0000\n", model.nSite);
        fprintf (file, "%d %lf 0.0000 0.0000\n", nxGrid, lenGridSide);
        fprintf (file, "%d 0.0000 %lf 0.0000\n", nyGrid, lenGridSide);
        fprintf (file, "%d 0.0000 0.0000 %lf\n", nzGrid, lenGridSide);

        int c_id = 1;
        for (int ix=0; ix<nxGrid; ix++){
          for (int iy=0; iy<nyGrid; iy++){
            for (int iz=0; iz<nzGrid; iz++){
              site = ix + iy*nxGrid + iz*nxGrid*nyGrid;
              c_id = (1+model.hasRearranged[site]);
              fprintf (file, "%d 0.0000 %d %d %d\n", c_id, ix, iy, iz);
            }
          }
        }

        for (int ix=0; ix<nxGrid; ix++){
          for (int iy=0; iy<nyGrid; iy++){
            for (int iz=0; iz<nzGrid; iz++){
              site = ix + iy*nxGrid + iz*nxGrid*nyGrid;
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
        sum += s.x[1];  // xy
      }

      strainFile << totalExternalStrain * strain_para << ' ' << sum / model.alle.size() << std::endl;
    };

    outputStrainFunc();
  
  }


}

