#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <sstream>

#include "GeometryVector.h"
#include "kiss_fft.h"
#include "kiss_fftnd.h"
//#include "mgl2/mgl.h"

template <typename T>
void plot(const std::vector<T> &data, int nGridPerSide, std::string file)
{
    // mglData x(nGridPerSide, nGridPerSide);
    // for (int i = 0; i < nGridPerSide; i++)
    // {
    //     for (int j = 0; j < nGridPerSide; j++)
    //     {
    //         x.SetVal(mreal(data[i * nGridPerSide + j]), i, j);
    //         //std::cout<<data[i*nGridPerSide+j]<<" ";
    //     }
    //     //std::cout<<std::endl;
    // }
    // mglGraph gr;
    // gr.SetSize(1600, 1600);
    // //gr.Aspect(0.75, 1.0);
    // //gr.Colorbar(">kbcyr");
    // gr.Tile(x, "bkr");
    // gr.WritePNG((file + std::string(".png")).c_str());
}

class gridModel
{
public:
  int nGridPerSide;
  double lGrid;
  int d;
  int nSite;
  int bufferCenter;

  std::vector<double> dSBuffer, alls;
  std::vector<GeometryVector> dEBuffer, alle; //e is strain
  std::vector<char> rearrangingStep;
  std::vector<bool> hasRearranged;

  std::mt19937 rEngine;
  std::normal_distribution<double> eDistribution;
  std::normal_distribution<double> sDistribution;


  /********************************************************************
   * gridModel
   *
   * Initializes model variables.
   *
   * Parameters
   * ----------
   *  nGrid : int
   *    Number of grid points on each side of the box
   *  lGrid : double
   *    Length of the side of each grid point
   *  d : int
   *    Dimension of the box
   *  eDistribution : std::normal_distribution<double>
   *    distribution of initial strains (e)
   *  sDistribution : std::normal_distribution<double>
   *    distribution of initial softnesses (S)
   *******************************************************************/
  gridModel(int nGrid, double lGrid, int d, 
        std::normal_distribution<double> eDist, 
        std::normal_distribution<double> sDist) : rEngine(0), 
        nGridPerSide(nGrid), lGrid(lGrid), d(d), eDistribution(eDist),
        sDistribution(sDist)
  {
    nSite = 1;
    for(int dd=0; dd<d; dd++)
    {
      nSite *= nGrid;
    }
  }

  /********************************************************************
   * initialize
   *
   * Initializes initial strains and softnesses for model according to
   * given distribution.
   *
   *******************************************************************/
  void initialize()
  {   
    alle.resize(nSite);
    alls.resize(nSite);
    hasRearranged.resize(nSite);
    rearrangingStep.resize(nSite);
    for (int i = 0; i < nSite; i++)
    {   
      this->alle[i].x[0] = this->eDistribution(this->rEngine);
      this->alle[i].x[1] = this->eDistribution(this->rEngine);
      this->alls[i] = this->sDistribution(this->rEngine);
      this->hasRearranged[i] = 0;
      this->rearrangingStep[i] = 0;
    }   
    this->getBuffer();
  } 


  /********************************************************************
  * unstack
  *
  * Takes a site index (site) and converts it into a site coordinate
  * site_coor in i, j(, k) space.
  *
  * Parameters
  * ----------
  *  site : int
  *    site index which runs from [0, nSite]
  *  site_coor: std::vector<int>
  *    d dimensional vector that is the location of the site in i, j(, 
  *    k) space.
  *
  ********************************************************************/
  void unstack(int site, std::vector<int>& site_coor)
  {
    if(d==2)
    {
      int x, y;
      y = site/nGridPerSide;
      x = site - y*nGridPerSide;
      site_coor[0] = x, site_coor[1] = y;
    }
    else
    {
      int x, y, z;
      z = site / (nGridPerSide*nGridPerSide);
      y = (site - z*nGridPerSide*nGridPerSide) / nGridPerSide;
      x = site - z*nGridPerSide*nGridPerSide - y*nGridPerSide;
    }
  }

  /********************************************************************
  * stack
  *
  * Takes a site coordinate (site_coor) in i, j(, k) space and converts
  * it to a site index (site).
  *
  * Parameters
  * ----------
  *  site_coor: std::vector<int>
  *    d dimensional vector that is the location of the site in i, j(, 
  *    k) space.
  *  site : int
  *    site index which runs from [0, nSite]
  ********************************************************************/
  void stack(std::vector<int>& site_coor, int& site)
  {
    if(d==2)
    {
      site = site_coor[0]+site_coor[1]*nGridPerSide;
    }
    else
    {
      site = site_coor[0]+site_coor[1]*nGridPerSide;
      site += site_coor[2]*nGridPerSide*nGridPerSide;
    }
  }


  void getBuffer()
  {
    // Initializes change of strain (dE) and change of softness (dS)
    // buffers as well as bufferCenter
    bufferCenter = nGridPerSide / 2;
    dEBuffer.resize(nSite);
    dSBuffer.resize(nSite);

    // Calculates dSBuffer assuming rearranging grid point is at center
    //   - site_coor : i, j(, k) coordinate of site
    //   - site_dx : dx, dy(, dz) position difference from the grid 
    //         center
    std::vector<int> site_coor(d);
    std::vector<double> site_dx(d);
    for(int site = 0; site<nSite; site++)
    {
      unstack(site, site_coor);
      for(int dd=0; dd<d; dd++)
      {
        site_dx[dd] = (site_coor[dd]-bufferCenter) * lGrid;
      }
      dsFromRearranger(site_dx);
    }   

    //strain buffer calculated by Fourier transform
    kiss_fft_cpx *inbuf = new kiss_fft_cpx[nSite];
    kiss_fft_cpx *outbuf = new kiss_fft_cpx[nSite];

    //fill in inbuf
    for(int site=0; site<nSite; site++)
    {
      double fillBuff = -4.0;
      double q2 = 0.0;
      unstack(site, site_coor);
      for(int dd=0; dd<d; dd++)
      {
        int i = site_coor[dd];
        int ii = (i > bufferCenter) ? i - nGridPerSide : i;
        double L = nGridPerSide * lGrid;
        double q = 2.0 * M_PI * ii / L;
        q2 += q*q;
        fillBuff *= q*q;
      }
      for(int dd=0; dd<d; dd++)
      {
        fillBuff /= q2;
      }
      inbuf[site].r = fillBuff;
      inbuf[site].i = 0.0;
    }
    inbuf[0].r = 0;

    // Performs FFT in either 2D or 3D
    if(d==2)
    {
      int temp[2] = {nGridPerSide, nGridPerSide};
      kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
      kiss_fftnd(st, inbuf, outbuf);
      free(st);
    }
    else
    {
      int temp[3] = {nGridPerSide, nGridPerSide, nGridPerSide};
      kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
      kiss_fftnd(st, inbuf, outbuf);
      free(st);
    }

    // Fills in dEBuffer
    double factor = 0.02 / std::fabs(outbuf[0].r);
    int site_periodic;
    for(int site=0; site<nSite; site++)
    {
      unstack(site, site_coor);
      for(int dd=0; dd<d; dd++)
      {
        int i = site_coor[dd] - bufferCenter;
        int ii = (i < 0) ? i + nGridPerSide : i;
        site_coor[dd] = ii;
      }
      stack(site_coor, site_periodic);
      dEBuffer[site_periodic] = GeometryVector(outbuf[site].r*factor, 0.0);
    }
    
    // Deletes outbuf and inbuf to save space.
    delete[] outbuf;
    delete[] inbuf;

  }


    bool startRearranging(GeometryVector e, double s)
    {
        double yieldStrain = 0.07 - 0.01 * s;
        if (yieldStrain < 0.05)
            yieldStrain = 0.05;
        return e.Modulus2() > yieldStrain * yieldStrain;
    }

  /********************************************************************
   * dsFromRearranger
   *
   * Function that determines how softness changes when a rearranger
   * is a vector distance site_dx away
   *
   * Parameters
   * ----------
   *  site_dx : std::vector<double>
   *    vector with the displacement between the current grid point and
   *    the rearranger in real units (can be 2D or 3D)
   *
   *******************************************************************/ 
  double dsFromRearranger(std::vector<double> site_dx)
  {
    double dxy, dz, r, dS;

    // Computes:
    //   (1) Distance in xy (or x) plane (dxy)
    //   (2) Distance in z (or y) plane (dz)
    //   (3) Euclidean distance (r)
    for(int dd=0; dd<d-1; dd++)
    {
      dxy += site_dx[dd]*site_dx[dd];
    }
    dxy = std::sqrt(dxy);
    dz = site_dx[d-1];
    r = std::sqrt(dxy*dxy+dz*dz);
    
    // Computes change in softness
    if(r < 4.0)
    {
      dS = -0.03;
    }
    else if (r < 30)
    {
      dS = 1.0 / r / r / r;
      dS -= 0.16 / r / r * std::sin(2.0 * std::atan2(dz, dxy));
    }
    else
    {
      dS = 0.0;
    }

    return dS;
  }

  /********************************************************************
   * shear
   *
   * Function that shears all of the sites in the model.
   *
   *******************************************************************/ 
  void shear()
  {
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nSite; i++)
    {
      this->alle[i].x[0] += 1e-6;
      this->hasRearranged[i] = 0;
      this->rearrangingStep[i] = 0;
    }
  }

  /********************************************************************
   * deltaEnergy
   *
   * Function that determines the change in energy given that a grid
   * point rearranges.
   *
   * Parameters
   * ----------
   *  site : int
   *    site index of the rearranging grid point; site is in [0,nSite]
   *
   *******************************************************************/ 
  double deltaEnergy(int site)
  {
    // Initializes:
    //   - deltaEnergy (change in energy from making rearrangement)
    //   - site_coor (position in i, j(, k) space) of site
    //   - site_coor_prime (position in i, j(, k) space of other site)
    double deltaEnergy = 0.0;
    std::vector<int> site_coor, site_coor_prime;

    // For each site', calculates relative position from rearranger
    // then calculates the change in E given a change in strain.
    unstack(site, site_coor);
    #pragma omp parallel for reduction(+:deltaEnergy)
    for (int site_prime = 0; site_prime < nSite; site_prime++)
    {

      // Calculates the difference in site positions taking into account PBC
      unstack(site_prime, site_coor_prime);
      for(int dd=0; dd<d; dd++)
      {
        double dx = site_coor[dd] - site_coor_prime[dd];
        site_coor_prime[dd] = bufferCenter - dx;
        site_coor_prime[dd] = (site_coor_prime[dd]+nGridPerSide) % nGridPerSide;
      }
      stack(site_coor_prime, site_prime);

      // Calculates change in energy due to chane in strain at this site
      if(site!=site_prime)
      {
        GeometryVector &e = alle[site_prime];
        GeometryVector &de = dEBuffer[site_prime];
        deltaEnergy += (e + de).Modulus2() - e.Modulus2();
      }
    }
  }
 

  /********************************************************************
   * updateStrainAndSoft
   *
   * Function that updates the strain and softness of the system given
   * a rearrangement at grid point (site)
   *
   * Parameters
   * ----------
   *  site : int
   *    site index of the rearranging grid point; site is in [0,nSite]
   *
   *******************************************************************/ 
  void updateStrainAndSoft(int site)
  {

    int site_relative;
    double dx;
    std::vector<int> site_coor, site_coor_prime;

    // For each site', calculates relative position from rearranger
    unstack(site, site_coor);
    #pragma omp parallel for schedule(static)
    for(int site_prime=0; site_prime<nGridPerSide; site_prime++)
    {

      // Gets current strain of site_prime
      GeometryVector &e = alle[site_prime];

      // Gets relative site relative to rearranger
      unstack(site_prime, site_coor_prime);
      for(int dd=0; dd<d; dd++)
      {   
        dx = site_coor[dd] - site_coor_prime[dd];
        site_coor_prime[dd] = bufferCenter - dx; 
        site_coor_prime[dd] = (site_coor[dd]+nGridPerSide) % nGridPerSide;
      }
      stack(site_coor_prime, site_relative);

      // Gets change in strain for site positioned relative to
      // rearranger
      GeometryVector &de = dEBuffer[site_relative];

      // Updates strain and softness
      e.AddFrom(de);
      alls[site] += dSBuffer[site_relative];
    }

  }

  bool avalanche(std::string outputPrefix = "")
  {
    bool avalancheHappened = false;
    int nStep = 0;
    std::vector<int> site_coor(d);

    // Avalanche continues until no more particles are rearranging
    int numRearrange = 1;
    while (numRearrange > 0)
    {

      // For each site:
      //   (1) Checks if site has started to rearrange
      //   (2) If site has started to rearrange, checks if energy is
      //         increased. If it is, rearrangement is stopped.
      for (int site = 0; site < nSite; site++)
      {

        // If site starts to rearrange marks
        //   (1) The site has rearranged
        //   (2) This is the first step that the site has rearranged
        if (startRearranging(alle[site], alls[site]))
        {
          hasRearranged[site] = 1;
          rearrangingStep[site] = 1;
        }

        if(rearrangingStep[site] > 0)
        {   
          // Stop rearrangement if energy increases
          if (deltaEnergy(site) > 0)
          {   
            rearrangingStep[site] = 0;
          }   
        }
      }


      // Updates: 
      //   - strains and softnesses around rearrangements 
      //   - number of rearrangements occuring currently
      //   - number of steps a rearrangement has taken
      numRearrange = 0;
      for (int site = 0; site < nSite; site++)
      {
        if (rearrangingStep[site] > 0)
        {
          updateStrainAndSoft(site);
          numRearrange++;
          rearrangingStep[site]++;
        }
      }


      // Updates the rearrangement site itself
      for (int site = 0; site < nSite; site++)
      {
        if (rearrangingStep[site] > 0)
        {
          //carry out the rearrangement
          alle[site] = 0.0;
          alls[site] = sDistribution(rEngine);
        }
      }

      // If there has been 1 or more rearrangements, we say an 
      // "avalanche" has occured.
      if (numRearrange > 0)
      {
        avalancheHappened = true;

        // Outputs number of rearrangements
        std::cout << "num rearranger in this frame=" << numRearrange;

        // Outputs mean energy
        double sum = 0.0;
        for (auto &e : this->alle)
          sum += e.Modulus2();
        std::cout << ", mean energy=" << sum / alle.size();

        // Outputs mean softness
        sum = 0.0;
        for (auto &s : this->alls)
          sum += s;
        std::cout << ", mean s=" << sum / alls.size();
        std::cout << std::endl;

        // outputs something else
        if (outputPrefix != std::string(""))
        {
          std::stringstream ss;
          ss << outputPrefix << "_step_" << (nStep++);
          plot(this->rearrangingStep, nGridPerSide, ss.str());
        }
      }

    }

    return avalancheHappened;

  }


};

int main()
{

  // Sets: 
  //   - number of grid points per side (nGridPerSide)
  //   - length of each grid side (lenGridSide)
  //   - dimensionality of the space (d)
  //   - distribution of grid strains (eDistribution)
  //   - distribution of grid softnesses (sDistribution)
  const int nGridPerSide = 200;
  const double lenGridSide = 1.0;
  const int d = 2;
  const int nAvalanche = 100;
  std::normal_distribution<double> eDistribution(0.0, 0.01);
  std::normal_distribution<double> sDistribution(-2.0, 1.0);

  // Declares a grid model with:
  //   - nGridPerSide grid points on each side
  //   - lenGridSide grid side length
  //   - d dimensionality
  //   - eDistribution initial distribution of strains
  //   - sDistribution initial distribution of softnesses
  // Initializes model.
  gridModel model(nGridPerSide, lenGridSide, d, eDistribution, 
        sDistribution);




  /*
  model.initialize();

  // Shears until nAvalanche avalanches have occured.
  int numAvalanche = 0;
  while (numAvalanche < nAvalanche)
  {

    // Shears model
    model.shear();

    std::stringstream ss;
    ss << "avalanche_" << numAvalanche;

    bool avalanched = model.avalanche(ss.str());
    numAvalanche += avalanched;
    if (avalanched)
    {
      std::cout << numAvalanche << "avalanches so far.\n";
      plot(model.hasRearranged, nGridPerSide, ss.str());
    }
  }
  */
}
