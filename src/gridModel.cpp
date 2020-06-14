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
  int nGridPerSide;      // Number of grid sites on each side
  double lGrid;          // Length of grid sites' side
  int d;                 // Dimension
  int nSite;             // Number of grid sites total
  int nStrain;           // Number of independent strain tensor comps
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

    // Strain tensor is symetric => (d^2+d)/2 independent components
    // Incompressible assumption => -1
    nStrain = (d*d+d)/2-1;
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
    for (int site = 0; site < nSite; site++)
    { 
      for(int strain=0; strain<nStrain; strain++)
      {  
        this->alle[site].x[strain] = this->eDistribution(this->rEngine);
      }
      this->alls[site] = this->sDistribution(this->rEngine);
      this->hasRearranged[site] = 0;
      this->rearrangingStep[site] = 0;
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
      site_coor[0] = x, site_coor[1] = y, site_coor[2] = z;
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

  /********************************************************************
   * getBuffer
   *
   * Initializes the change of strain (dE) and change of softness (dS)
   * vectors.
   *
   ********************************************************************/ 
  void getBuffer()
  {

    // Declares variables for Fourier Transform
    int i, j, k;
    int site_periodic;
    double qx, qy, qz, q2;
    double factor;                     // Not quite sure what this is; should ask!
    double L = nGridPerSide * lGrid;

    //strain buffer calculated by Fourier transform
    kiss_fft_cpx *inbuf = new kiss_fft_cpx[nSite];
    kiss_fft_cpx *outbuf = new kiss_fft_cpx[nSite];

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


    // For each strain component, does Fourier Transform to obtain dEBuffer.
    for(int strain=0; strain<nStrain; strain++)
    {

      // Fill in inbuf of FT
      for(int site=0; site<nSite; site++)
      {
        unstack(site, site_coor);

        // Obtains Fourier components in x and y directions
        i = site_coor[0];
        i = (i > bufferCenter) ? i - nGridPerSide : i; // Enforces PBC
        qx = 2.0 * M_PI * i / L;
        j = site_coor[1];
        j = (j > bufferCenter) ? j - nGridPerSide : j;
        qy = 2.0 * M_PI * j / L;

        // Fourier components of shear stress due to local strain in 2D
        // and 3D from Picard, Elastic (2004)
        if(d==2)
        {
          q2 = qx*qx + qy*qy;

          if(strain==0)
          {
            // xy
            inbuf[site].r = -4.0*qx*qx*qy*qy/(q2*q2);
          }
          else
          {
            // xx
            inbuf[site].r = -2.0*qx*qy*(qx*qx - qy*qy)/(q2*q2);
          }
        }
        else
        {
          // Obtains Fourier component in z direction
          k = site_coor[2];
          k = (k > bufferCenter) ? k - nGridPerSide : k;
          qz = 2.0 * M_PI * k / L;
          q2 = qx*qx + qy*qy + qz*qz;

          if(strain==0)
          {
            // xy
            inbuf[site].r = (-4.0*qx*qx*qy*qy+qz*qz*q2)/(q2*q2);
          }
          else if(strain==1)
          {
            // yy
            inbuf[site].r = -2.0*qx*qy*(q2-2.0*qy*qy)/(q2*q2);
          }
          else if(strain==2)
          {
            // xz
            inbuf[site].r = qy*qz*(q2-4.0*qx*qx)/(q2*q2);
          }
          else if(strain==3)
          {
            // yz
            inbuf[site].r = qx*qz*(q2-4.0*qy*qy)/(q2*q2);
          }
          else
          {
            // zz
            inbuf[site].r = -4.0*qx*qy*qz*qz/(q2*q2);
          }
        }
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
      if(strain==0)
      {
        factor = 0.02 / std::fabs(outbuf[0].r);  // Don't understand this line...
      }
      for(int site=0; site<nSite; site++)
      {
        unstack(site, site_coor);
        for(int dd=0; dd<d; dd++)
        {
          i = site_coor[dd] - bufferCenter;   
          i = (i < 0) ? i + nGridPerSide : i;   // Enforces PBC
          site_coor[dd] = i;
        }
        stack(site_coor, site_periodic);

        dEBuffer[site_periodic].x[strain] = outbuf[site].r*factor;
      }
    }

    // Deletes outbuf and inbuf to save space.
    delete[] outbuf;
    delete[] inbuf;

  }

  /********************************************************************
   * startRearranging
   *
   * Function that determines when a grid point is going to start
   * rearranging from the grid points strain and softness.
   *
   * Parameters
   * ----------
   *  e : GeometryVector
   *    Strain vector representing strain at the grid point.
   *  s : double
   *    Softness at the given grid point.
   *******************************************************************/
  bool startRearranging(GeometryVector e, double s)
  {
    double e2;
    double yieldStrain = 0.07 - 0.01 * s;

    // Calculates yield strain
    if(yieldStrain < 0.05)
      yieldStrain = 0.05;

    // Calculates sum of strain components
    if(d==2)
    {
      e2 = e.Modulus2();
    }
    else
    {
      e2 = (e.Modulus2()+e.x[1]*e.x[4]);
    }

    return e2 > yieldStrain*yieldStrain;
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
    int  dx, site_prime_shift;
    std::vector<int> site_coor(d), site_coor_prime(d);
    std::vector<int> site_coor_shift(d);
    GeometryVector e_prime;

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
        dx = site_coor_prime[dd] - site_coor[dd];
        site_coor_shift[dd] = bufferCenter + dx;
        site_coor_shift[dd] = (site_coor_shift[dd]+nGridPerSide) 
              % nGridPerSide;
      }
      stack(site_coor_shift, site_prime_shift);

      // Calculates change in energy due to chane in strain at this site
      GeometryVector &e = alle[site_prime];
      GeometryVector &de = dEBuffer[site_prime_shift];
      e_prime = e + de;
      //std::cout << e.x[0] << ", " << e.x[1] << ", " << e.x[2] << ", ";
      //std::cout << e.x[3] << ", " << e.x[4] << "\n";
      //std::cout << de.x[0] << ", " << de.x[1] << ", " << de.x[2] << ", ";
      //std::cout << de.x[3] << ", " << de.x[4] << "\n";
      //std::cout << e_prime.x[0] << ", " << e_prime.x[1] << ", " << e_prime.x[2] << ", ";
      //std::cout << e_prime.x[3] << ", " << e_prime.x[4] << "\n";
      if(site!=site_prime)
      {
        // Non-rearranging site energy increases
        if(d==2)
        {
          deltaEnergy += (e_prime.Modulus2() - e.Modulus2());
        }
        else
        {
          deltaEnergy += ( 
                (e_prime.Modulus2()+e_prime.x[1]*e_prime.x[4]) - 
                (e.Modulus2()+e.x[1]*e.x[4]));
        }
      }
      else
      {
        // Rearranger strain goes to 0
        if(d==2)
        {
          deltaEnergy -= (e.Modulus2());
        }
        else
        {
          deltaEnergy -= (e.Modulus2()+e.x[1]*e.x[4]);
        }
      }
    }

    // Returns energy
    return deltaEnergy;
  }
 

  /********************************************************************
   * updateStrainAndSoft
   *
   * Function that updates the strain and softness of the system given
   * a rearrangement at grid point (site). Does not update the strain
   * or softness of the grid point itself.
   *
   * Parameters
   * ----------
   *  site : int
   *    site index of the rearranging grid point; site is in [0,nSite]
   *
   *******************************************************************/ 
  void updateStrainAndSoft(int site)
  {

    int site_shift;
    int dx;
    std::vector<int> site_coor(d), site_coor_prime(d);
    std::vector<int> site_coor_shift(d);

    // For each site', calculates relative position from rearranger
    unstack(site, site_coor);
    #pragma omp parallel for schedule(static)
    for(int site_prime=0; site_prime<nSite; site_prime++)
    {

      // Gets current strain of site_prime
      GeometryVector &e = alle[site_prime];

      // Gets relative site relative to rearranger
      unstack(site_prime, site_coor_prime);
      for(int dd=0; dd<d; dd++)
      {   
        dx = site_coor_prime[dd] - site_coor[dd];
        site_coor_shift[dd] = bufferCenter + dx; 
        site_coor_shift[dd] = (site_coor_shift[dd]+nGridPerSide) % 
              nGridPerSide; // Assumes PBC
      }
      stack(site_coor_shift, site_shift);

      if(site!=site_prime)
      {
        // Gets change in strain for site positioned relative to
        // rearranger
        GeometryVector &de = dEBuffer[site_shift];

        // Updates strain and softness
        e.AddFrom(de);
        alls[site_prime] += dSBuffer[site_shift];
      }
    }

  }

  bool avalanche(std::string outputPrefix = "")
  {
    bool avalancheHappened = false;
    int nStep = 0, numRearrange = 1;
    std::vector<int> site_coor(d);

    // Avalanche continues until no more particles are rearranging
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
        // std::cout << "1.1\n";
        if (startRearranging(alle[site], alls[site]))
        {
          hasRearranged[site] = 1;
          rearrangingStep[site] = 1;
        }

        if(rearrangingStep[site] > 0)
        {   
          // Stop rearrangement if energy increases
          if(deltaEnergy(site) > 0)
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

      // Updates the rearrangement site itself. This must be done
      // separately from above to ensure rearranger strain --> 0
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
          if(d==2)
          {
            sum += e.Modulus2();
          }
          else
          {
            sum += e.Modulus2() + e.x[1]*e.x[4];
          }
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
  const int nGridPerSide = 6;
  const double lenGridSide = 1.0;
  const int d = 3;
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
  model.initialize();

  // Check startRearranging
  //GeometryVector e;
  //e.x[0] = 0.069;
  //std::cout << "startRearranging(0.069,0) [FALSE] = " << model.startRearranging(e,0.0) << "\n";
  //e.x[0] = 0.049;
  //std::cout << "startRearranging(0.049,3.0) [FALSE] = " << model.startRearranging(e,3.0) << "\n";
  //e.x[0] = 0.071;
  //std::cout << "startRearranging(0.071,0) [TRUE] = " << model.startRearranging(e,0.0) << "\n";
  //e.x[0] = 0.051;
  //std::cout << "startRearranging(0.051,0) [TRUE] = " << model.startRearranging(e,3.0) << "\n";
  //e.x[0] = 0;
  //e.x[1] = 0.04099;
  //e.x[4] = 0.04099;
  //std::cout << "startRearranging(<0,0.04099,0,0,0.04099>,0) [TRUE] = " << model.startRearranging(e,0.0) << "\n";
  //e.x[4] = -0.0499;
  //std::cout << "startRearranging(<0,0.04099,0,0,-0.04099>,0) [FALSE] = " << model.startRearranging(e,0.0) << "\n";

  // Checks change in energy due to rearrangement
  //std::cout << "dE_0 = " << model.deltaEnergy(0) << "\n";
  //std::cout << "dE_1 = " << model.deltaEnergy(1) << "\n";
  //std::cout << "dE_790 = " << model.deltaEnergy(790) << "\n";
  //std::cout << "dE_10 = " << model.deltaEnergy(10) << "\n";
  //std::cout << "dE_100 = " << model.deltaEnergy(100) << "\n";
  // model.alle[0].x[0] = 1.0;
  // std::cout << "dE_0 alle[0].x[0] changed = " << model.deltaEnergy(0) << "\n";
  // std::cout << "dE_1 alle[0].x[0] changed = " << model.deltaEnergy(1) << "\n";
  //model.alle[0].x[0] = 0.0;

  // Checks updates of strain and softness
  //int site_orig = 770, site;
  //std::vector<int> site_coor(3);
  //std::cout << "Updates sofntess and strain around site " << site << "\n";
  //model.unstack(site_orig, site_coor);
  //std::cout << "site_coor = " << site_coor[0] << ", " << site_coor[1] << ", " << site_coor[2] << "\n";
  //std::cout << "Checks corect update to softness and strain applied to site_coor+<0, 1, 0>\n";
  //site_coor[0] = (site_coor[0] + 0 + nGridPerSide) % nGridPerSide;
  //site_coor[1] = (site_coor[1] + 1 + nGridPerSide) % nGridPerSide;
  //site_coor[2] = (site_coor[2] + 0 + nGridPerSide) % nGridPerSide;
  //model.stack(site_coor, site);
  //std::cout << "site = " << site << "\n";
  //std::cout << "site_coor = " << site_coor[0] << ", " << site_coor[1] << ", " << site_coor[2] << "\n";
  //e = model.alle[site];
  //std::cout << "e_{orig} = " << e.x[0] << ", " << e.x[1] << ", " << e.x[2] << ", ";
  //std::cout << e.x[3] << ", " << e.x[4] << "\n";
  //model.updateStrainAndSoft(site_orig);   // UPDATE OCCURS HERE!!
  //e = model.alle[site];
  //std::cout << "e_{new} = " << e.x[0] << ", " << e.x[1] << ", " << e.x[2] << ", ";
  //std::cout << e.x[3] << ", " << e.x[4] << "\n";
  //site_coor[0] = 5 + 0;
  //site_coor[1] = 5 + 1;
  //site_coor[2] = 5 + 0;
  //model.stack(site_coor, site);
  //std::cout << "de = " << model.dEBuffer[site].x[0] << ", " << model.dEBuffer[site].x[1];
  //std::cout << ", " << model.dEBuffer[site].x[2] << ", " << model.dEBuffer[site].x[3];
  //std::cout << model.dEBuffer[site].x[4] << "\n";


  //Shears until nAvalanche avalanches have occured.
  int numAvalanche = 0;
  int ii = 0;
  std::fstream strainFile("xyStrain.txt", std::fstream::out);
  while (numAvalanche < nAvalanche)
  {

    // Shears model
    model.shear();
    if(ii%1000==0)
      std::cout << "ii = " << ii << "\n";
    ii++;
    //numAvalanche += 1;

    std::stringstream ss;
    ss << "avalanche_" << numAvalanche;

    bool avalanched = model.avalanche(ss.str());

    numAvalanche += avalanched;
    if (avalanched)
    {
      std::cout << numAvalanche << "avalanches so far.\n";
      plot(model.hasRearranged, nGridPerSide, ss.str());
    }

    double sum=0.0;
    for(auto s : model.alle)
      sum+=s.x[0];
    strainFile<<sum/model.alle.size()<<std::endl;

  }
}
