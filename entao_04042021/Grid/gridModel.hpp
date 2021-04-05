#ifndef GRIDMODEL_HPP_
#define GRIDMODEL_HPP_

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <sstream>

#include "../GeometryVector/GeometryVector.h"
#include "../Kiss_FFT/kiss_fft.h"
#include "../Kiss_FFT/kiss_fftnd.h"
//#include <netcdf>

class gridModel
{
  public:
    int nGridPerSide;    // number of grid sites on each side
    double lGrid;        // Length of grid sites' side
    int d;               // Dimension
    int nSite;           // Number of grid sites total
    int nStrain;         // Number of independent strain tensor comps
    int bufferCenter;
    double omega=1.0;    // volume of the grid point
    double T;            // system temperature   
    double meanSoftness;
    int s_induce_r;   // stress induced rearrangements
    double c_soft;       // average softness at current step
    std::vector<double> dSBuffer, alls;
    std::vector<double> alls_local;

    // e is elastic strain
    std::vector<GeometryVector> dEBuffer[::MaxDimension], alle;  
    std::vector<GeometryVector> rearrangingIntensity;

    std::vector<char> rearrangingStep;
    std::vector<bool> hasRearranged;

    std::mt19937 rEngine;  // random number generator
    // distribution of strain initially    
    std::normal_distribution<double> eDistribution;

    std::normal_distribution<double> sDistribution;
    
    // 2nd generation new varaibles
    // distribution of elastic strain after rearrangement
    std::normal_distribution<double> residualStrainDistribution;

    std::uniform_real_distribution<double> PxDistribution;
    std::vector<double> yieldStrainPx;

    // length of a rearrangement in frames
    int rearrangeFrameLength;

    // Constructor: initialize variables

    gridModel(int nGirdPerSide, double lGrid, int seed, int d, double T,
        double meanSoftness, double c_soft, int rearrangeFrameLength,
        std::normal_distribution<double> eDist,
        std::normal_distribution<double> sDist): 
      rEngine(seed), residualStrainDistribution(0.0, 0.0),
      rearrangeFrameLength(rearrangeFrameLength), 
      nGridPerSide(nGirdPerSide), lGrid(lGrid),
      d(d), T(T), meanSoftness(meanSoftness), c_soft(c_soft), eDistribution(eDist), sDistribution(sDist)
    {
      nSite = 1;

      for (int dd=0; dd<d; dd++){
        nSite *= nGridPerSide;
        omega *= lGrid;
      }

      printf ("%d\n", nSite);
      nStrain = (d*d+d)/2 - 1;    // it is 2 for 2D, and 5 for 3D
     
      this -> allocate();
      this -> getBuffer();
      this -> initialize(); 
    }


    // initialize functions 
    void initialize();   // initialize strain and softness, 
                         // according to given distribution
    void unstack(int site, std::vector<int>& site_coor);
    void stack(std::vector<int> site_coor, int& site);  // this line is changed, no refer (&) for site_coor
    void getBuffer();    // initialize the strain change, 
                         // and softness change
    bool startRearranging(int site, GeometryVector e, double s);  // Given the E and S, determine whether the grid point is going to start rearranging
    double dsFromRearranger(std::vector<double> site_dx);  // Given the vector distance site_dx from the rearranger, calculate how the softness change
    void shear(); // shear all the sites in the grid model 
    void equi(); 
//    void updateStrainAndSoft(int site);  // update strain and softness 
    bool rearrange(std::string outputPrefix = "");

// new functions in 2nd generation    
    
    double deviatoricYieldStrain(int i); 
    void allocate();   
    void local_mean();   // calculate local softness, designed for restoring to local mean, 
                         // not used in current model

};
#endif
