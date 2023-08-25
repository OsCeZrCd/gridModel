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
    int nxGrid;    // number of grid sites on each side
    int nyGrid;
    int nzGrid;
    double lGrid;        // Length of grid sites' side
    int d;               // Dimension
    int nSite;           // Number of grid sites total
    int nStrain;         // Number of independent strain tensor comps
    int bufferCenterX;
    int bufferCenterY;
    int bufferCenterZ;
    double omega=1.0;    // volume of the grid point
    double T;            // system temperature   
    double meanSoftness;
    int s_induce_r;   // stress induced rearrangements
    int step;
    double c_soft;       // current average softness 
    
    double soft_const;  // soft constant for devatoric strain

    std::vector<double> dSBuffer, alls; 
    std::vector<double> alls_local;

    // e is elastic strain
    std::vector<GeometryVector> dEBuffer[::MaxDimension], alle, alle_buffer;  
    std::vector<GeometryVector> rearrangingIntensity;

    //std::vector<char> rearrangingStep;
    //std::vector<bool> hasRearranged;

    std::vector<int> rearrangingStep;
    std::vector<char> hasRearranged;


    std::mt19937 rEngine;  // random number generator
    // distribution of strain initially    
    std::normal_distribution<double> eDistribution;

    std::normal_distribution<double> sDistribution;
    
    // 2nd generation new varaibles
    // distribution of elastic strain after rearrangement
    std::normal_distribution<double> residualStrainDistribution;

    std::uniform_real_distribution<double> PxDistribution;
    std::vector<double> yieldStrainPx;

    std::vector<double> movingAverageTarget;
    std::vector<double> PR_bulk;
    std::vector<double> S_bulk;

    // length of a rearrangement in frames
    int rearrangeFrameLength;

    // Constructor: initialize variables

    gridModel(int nxGrid, int nyGrid, int nzGrid,  double lGrid, 
        int seed, int d, double T, double meanSoftness, double step, 
        double c_soft, int rearrangeFrameLength,
        std::normal_distribution<double> eDist,
        std::normal_distribution<double> sDist): 
      rEngine(seed), residualStrainDistribution(0.0, 0.00001),
      rearrangeFrameLength(rearrangeFrameLength), 
      nxGrid(nxGrid), nyGrid(nyGrid), nzGrid(nzGrid), lGrid(lGrid),
      d(d), T(T), meanSoftness(meanSoftness), step(step), c_soft(c_soft), eDistribution(eDist), sDistribution(sDist)
    {
      nSite = 1;
      printf ("%d %d %d\n", nxGrid, nyGrid, nzGrid);
      nSite = nxGrid * nyGrid * nzGrid;
      for (int dd=0; dd<d; dd++){
        //nSite *= nGridPerSide;
        omega *= lGrid;
      }

      printf ("%d\n", nSite);
      nStrain = (d*d+d)/2 - 1;    // it is 2 for 2D, and 5 for 3DKE
     
      this -> allocate();
      printf("h1\n");
      this -> getBuffer();
      printf("h2\n");
      this -> initialize();
      printf("h3\n");

      if (T==0.05)
        soft_const = 0.00413 / (4.95e-4 * 1);   //6.6e-4;  //3.884e-4;
      else
        soft_const = 0.00207 / (4.95e-4 * 1);   //4.010e-4;
    }


    // initialize functions 
    void initialize();   // initialize strain and softness, 
                         // according to given distribution
    void unstack(int site, std::vector<int>& site_coor);
    void stack(std::vector<int> site_coor, int& site);  // this line is changed, no refer (&) for site_coor
    void getBuffer();    // initialize the strain change, 
                         // and softness change
    bool startRearranging(int site, GeometryVector e, double s);  // Given the E and S, determine whether the grid point is going to start rearranging
    //double dsFromRearranger(std::vector<double> site_dx, const GeometryVector & rearrangingIntensity);  // Given the vector distance site_dx from the rearranger, calculate how the softness change
    double dsFromRearranger(std::vector<double> site_dx);

    void shear(); // shear all the sites in the grid model 
    void equi(); 
//    void updateStrainAndSoft(int site);  // update strain and softness 
    bool rearrange(std::string outputPrefix = "");

// new functions in 2nd generation    
    
    double deviatoricYieldStrain(int i); 
    void allocate();
    void local_mean();
    void induced_rearrangement(int site);
    //void nearby_PR(int site);//, std::vector<double> PR_bulk);
    void nearby_PR(int site, std::vector<double> PR_bulk);

    //    double ran2(void);

    std::vector<double> PR_1={1.60522534e-05, 4.36268517e-05, 1.18569159e-04, 3.22247537e-04, 
      8.75805108e-04, 2.38026516e-03, 6.46909021e-03,};
    
    std::vector<double> PR_2={5.89829970e-06, 1.81448963e-05, 5.58190122e-05, 1.71715620e-04,
     5.28247510e-04, 1.62504396e-03, 4.99911088e-03};
};
#endif
