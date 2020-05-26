#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

#include "GeometryVector.h"
#include "kiss_fft.h"
#include "kiss_fftnd.h"
#include "mgl2/mgl.h"

template <typename T>
void writebinary(const std::vector<T> &data, int nGridPerSide, std::string file)
{
  std::ofstream re_data;
  re_data.open ("data_hasRe.bin", std::ios::out | std::ios::binary | std::fstream::app);

  int nsum = nGridPerSide*nGridPerSide;
  int total_re = 0;
  for (int i = 0; i < nsum; i++)
  {
      if (data[i]==1){
          total_re +=1;
      }
  }
  re_data.write((char*)&total_re,sizeof(int));

  for (int i = 0; i < nsum; i++)
  {

        if (data[i]==1){
          int gindex=i;
          re_data.write((char*)&gindex,sizeof(int));
         }

  }
  std::cout<< "rearrangement written: " << total_re <<std::endl;
  re_data.close();

}

template <typename T>
void writebinary_scalar(const std::vector<T> &data, int nGridPerSide, std::string file)
{
  std::ofstream re_data;
  re_data.open ("data_Softness.bin", std::ios::out | std::ios::binary | std::fstream::app);

  double nsum = nGridPerSide*nGridPerSide;
  re_data.write((char*)&nsum,sizeof(double));

  for (int i = 0; i < nsum; i++)
  {

          double gindex=data[i];
          re_data.write((char*)&gindex,sizeof(double));

  }
  std::cout<< "softness written: " << nsum <<std::endl;
  re_data.close();
}


class gridModel
{
public:
    int nGridPerSide;
    double lGrid;
    int bufferCenter;
    double mean_softness;
    double dmean_softness;

    std::vector<double> dSBuffer, alls;
    std::vector<GeometryVector> dEBuffer, alle; //e is strain
    std::vector<char> rearrangingStep;
    std::vector<bool> hasRearranged;

    std::mt19937 rEngine;
    std::normal_distribution<double> sDistribution;
    std::normal_distribution<double> eDistribution;

    std::gamma_distribution<double> coeffDistribution;
    std::vector<double> yieldStrainCoeff;

    gridModel(int nGrid, double lGrid) : rEngine(0), eDistribution(0.0, 0.000001), sDistribution(-0.2, 2.0), coeffDistribution(1.5, 0.667), nGridPerSide(nGrid), lGrid(lGrid)
    {
    }

    bool startRearranging(GeometryVector e, double s, int i)
    {
        // double yieldStrain = 0.006 - 0.002 * s;
        // if (yieldStrain < 0.002)
        //     yieldStrain = 0.002;
        // return e.Modulus2() > yieldStrain * yieldStrain;

        double yieldStrain = 0.009 - 0.003 * s;
        if (yieldStrain < 0.003)
            yieldStrain = 0.003;
        yieldStrain*=yieldStrainCoeff[i];
        return e.Modulus2() > yieldStrain * yieldStrain;

        //return e.x[0] > yieldStrain;
    }
    // GeometryVector eFromRearranger(double dx, double dy, double r)
    // {
    //     double magnitude;
    //     if (r < 1.0)
    //         magnitude = 3e-2;
    //     else
    //         magnitude = 3e-2 / r / r;

    //     double theta = std::atan2(dy, dx);
    //     return magnitude * GeometryVector(std::cos(4 * theta), std::sin(4 * theta));
    //     //return magnitude * GeometryVector(std::cos(4 * theta), 0.0);
    // }
    double dsFromRearranger(double dx, double dy, double r)
    {
        if (r < 1.0)
             return 0.00;
        else if (r < 30)
            return 0.75*0.32 / r / r / r - 0.75*0.16 / r / r * (std::sin(2.0 * std::atan2(dy, dx)));
        else
            return 0.0;
    }
    void getBuffer()  //assume effect from all rearrangements are the same, so only need once
    {
        bufferCenter = nGridPerSide / 2;
        dEBuffer.resize(nGridPerSide * nGridPerSide); // vector of doubles
        dSBuffer.resize(nGridPerSide * nGridPerSide);
        for (int i = 0; i < nGridPerSide; i++)
            for (int j = 0; j < nGridPerSide; j++)
            {
                double dx = (i - bufferCenter) * lGrid;
                double dy = (j - bufferCenter) * lGrid;
                double r = std::sqrt(dx * dx + dy * dy);
                int index = i * nGridPerSide + j;
                // dEBuffer[index] = eFromRearranger(dx, dy, r);
                dSBuffer[index] = dsFromRearranger(dx, dy, r);
            }

        double factor;
        {
            //strain buffer calculated by Fourier transform
            int numPixels = nGridPerSide * nGridPerSide;
            kiss_fft_cpx *inbuf = new kiss_fft_cpx[numPixels];
            kiss_fft_cpx *outbuf = new kiss_fft_cpx[numPixels];
            //fill in inbuf
            for (int i = 0; i < nGridPerSide; i++)
                for (int j = 0; j < nGridPerSide; j++)
                {
                    int ii = (i > bufferCenter) ? i - nGridPerSide : i;
                    int jj = (j > bufferCenter) ? j - nGridPerSide : j;
                    double L = nGridPerSide * lGrid;
                    double pm = 2 * M_PI * ii / L;
                    double qn = 2 * M_PI * jj / L;
                    double q2 = pm * pm + qn * qn;

                    inbuf[i * nGridPerSide + j].r = -4 * pm * pm * qn * qn / q2 / q2;
                    inbuf[i * nGridPerSide + j].i = 0;
                }
            inbuf[0].r = 0;

            int temp[2] = {nGridPerSide, nGridPerSide};
            kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, 2, true, nullptr, nullptr);
            kiss_fftnd(st, inbuf, outbuf);
            //fill in dEBuffer
            factor = 0.005 / std::fabs(outbuf[0].r);

            for (int i = 0; i < nGridPerSide; i++)
                for (int j = 0; j < nGridPerSide; j++)
                {
                    int index = i * nGridPerSide + j;
                    int ii = i - bufferCenter;
                    while (ii < 0)
                        ii += nGridPerSide;
                    int jj = j - bufferCenter;
                    while (jj < 0)
                        jj += nGridPerSide;
                    int index2 = ii * nGridPerSide + jj;
                    dEBuffer[index2] = GeometryVector(outbuf[index].r * factor, 0.0);
                }

            free(st);
            delete[] outbuf;
            delete[] inbuf;
        }
        {
            //strain buffer calculated by Fourier transform, another component
            int numPixels = nGridPerSide * nGridPerSide;
            kiss_fft_cpx *inbuf = new kiss_fft_cpx[numPixels];
            kiss_fft_cpx *outbuf = new kiss_fft_cpx[numPixels];
            //fill in inbuf
            for (int i = 0; i < nGridPerSide; i++)
                for (int j = 0; j < nGridPerSide; j++)
                {
                    int ii = (i > bufferCenter) ? i - nGridPerSide : i;
                    int jj = (j > bufferCenter) ? j - nGridPerSide : j;
                    double L = nGridPerSide * lGrid;
                    double pm = 2 * M_PI * ii / L;
                    double qn = 2 * M_PI * jj / L;
                    double q2 = pm * pm + qn * qn;

                    inbuf[i * nGridPerSide + j].r = -2 * pm * qn * (pm * pm - qn * qn) / q2 / q2;
                    inbuf[i * nGridPerSide + j].i = 0;
                }
            inbuf[0].r = 0;

            int temp[2] = {nGridPerSide, nGridPerSide};
            kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, 2, true, nullptr, nullptr);
            kiss_fftnd(st, inbuf, outbuf);
            //fill in dEBuffer

            for (int i = 0; i < nGridPerSide; i++)
                for (int j = 0; j < nGridPerSide; j++)
                {
                    int index = i * nGridPerSide + j;
                    int ii = i - bufferCenter;
                    while (ii < 0)
                        ii += nGridPerSide;
                    int jj = j - bufferCenter;
                    while (jj < 0)
                        jj += nGridPerSide;
                    int index2 = ii * nGridPerSide + jj;
                    dEBuffer[index2].x[1] = outbuf[index].r * factor;
                }

            free(st);
            delete[] outbuf;
            delete[] inbuf;
        }

        //debug temp
        // for (int i = 0; i < nGridPerSide; i++)
        // {
        //     for (int j = 0; j < nGridPerSide; j++)
        //     {
        //         int index = i * nGridPerSide + j;
        //         std::cout << dEBuffer[index].x[0] << ' ';
        //     }
        //     std::cout << std::endl;
        // }
        // exit(0);
    }
    void initialize()
    {
        int nSite = nGridPerSide * nGridPerSide;
        alle.resize(nSite);
        alls.resize(nSite);
        hasRearranged.resize(nSite);
        rearrangingStep.resize(nSite);
        yieldStrainCoeff.resize(nSite);

        mean_softness = 0;
        dmean_softness = 2E-5;

        for (int i = 0; i < nSite; i++)
        {
            this->alle[i].x[0] = this->eDistribution(this->rEngine); //initialize strain with a distrbution e with a random number generator
            this->alle[i].x[1] = this->eDistribution(this->rEngine);
            this->alls[i] = this->sDistribution(this->rEngine); //initialize softness with a distrbution s
            this->yieldStrainCoeff[i] = this->coeffDistribution(this->rEngine);
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
        }
        this->getBuffer();
    }


    void shear()
    {

        double sum = 0.0;
        for (auto &s : this->alls)
            sum += s;
        double mean_s = sum / alls.size();

        int nSite = nGridPerSide * nGridPerSide;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < nSite; i++)
        {
            double E_increament = 1.0 + (this->alls[i] - mean_s) * 0.15;
            E_increament = std::max(0.0,E_increament);
            this->alle[i].x[0] += 1.0e-6 * E_increament;
            this->alle[i].x[1] += 0.0e-6;
            this->alls[i] += dmean_softness;
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
        }

        mean_softness = mean_softness + dmean_softness;

    }


    int avalanche(std::string outputPrefix = "")
    {
        bool avalancheHappened = false;
        int nSite = nGridPerSide * nGridPerSide;
        int nStep = 0;
        double deltaEnergy;
#pragma omp parallel
        {
            int numRearrange = 1;
            while (numRearrange > 0)
            {
#pragma omp for schedule(static)
                for (int i = 0; i < nSite; i++)
                    if (startRearranging(alle[i], alls[i],i))
                    {
                        // hasRearranged[i] = 1;
                        rearrangingStep[i] = 1;
                    }

                //stop rearrangements that increases energy
                for (int i = 0; i < nSite; i++)
                {
                    if (rearrangingStep[i] > 0)
                    {
                        //calculate energy difference
                        deltaEnergy = 0;
                        int rx = i / nGridPerSide;
                        int ry = i % nGridPerSide;
#pragma omp barrier
#pragma omp for schedule(static) reduction(+ \
                                           : deltaEnergy)
                   for (int x = 0; x < nGridPerSide; x++)
                   {
                       int xInBuffer = bufferCenter - rx + x;
                       // while (xInBuffer < 0)
                       //     xInBuffer += nGridPerSide;
                       // while (xInBuffer >= nGridPerSide)
                       //     xInBuffer -= nGridPerSide;
                       if (xInBuffer>=0 && xInBuffer<nGridPerSide)
                       {
                       for (int y = 0; y < nGridPerSide; y++)
                       {
                           int yInBuffer = bufferCenter - ry + y;
                           // while (yInBuffer < 0)
                           //     yInBuffer += nGridPerSide;
                           // while (yInBuffer >= nGridPerSide)
                           //     yInBuffer -= nGridPerSide;
                           if (yInBuffer>=0 && yInBuffer<nGridPerSide)
                           {
                           GeometryVector &e = alle[x * nGridPerSide + y];
                           GeometryVector &de = dEBuffer[xInBuffer * nGridPerSide + yInBuffer];

                           if (ry != y || rx != x)
                               deltaEnergy += (e + de).Modulus2() - e.Modulus2();
                           else
                               deltaEnergy -= e.Modulus2();
                           }
                       }
                       }
                   }
#pragma omp single
                        {
                            //stop if energy increases
                            // std::cout<<"de= "<<deltaEnergy<<' ';
                            if (deltaEnergy > 0)
                            {
                                //std::cout << "rearrangement at stage " << int(rearrangingStep[i]) << " refuted due to energy criteria\n";
                                rearrangingStep[i] = 0;
                            }
                        }
                    }
                }

#pragma omp barrier
                //rearrangement affect other sites parameters
                numRearrange = 0;
                for (int i = 0; i < nSite; i++)
                {
                    if (rearrangingStep[i] > 0)
                    {
                        //update softness and strain
                        int rx = i / nGridPerSide;
                        int ry = i % nGridPerSide;
#pragma omp for schedule(static)
                    for (int x = 0; x < nGridPerSide; x++)
                    {
                        int xInBuffer = bufferCenter - rx + x;
                        // while (xInBuffer < 0)
                        //     xInBuffer += nGridPerSide;
                        // while (xInBuffer >= nGridPerSide)
                        //     xInBuffer -= nGridPerSide;
                        if (xInBuffer>=0 && xInBuffer<nGridPerSide)
                        {
                        for (int y = 0; y < nGridPerSide; y++)
                        {
                            int yInBuffer = bufferCenter - ry + y;
                            // while (yInBuffer < 0)
                            //     yInBuffer += nGridPerSide;
                            // while (yInBuffer >= nGridPerSide)
                            //     yInBuffer -= nGridPerSide;
                            //alle[x * nGridPerSide + y] += dEBuffer[xInBuffer * nGridPerSide + yInBuffer];
                            if (yInBuffer>=0 && yInBuffer<nGridPerSide)
                            {
                            GeometryVector &e = alle[x * nGridPerSide + y];
                            GeometryVector &de = dEBuffer[xInBuffer * nGridPerSide + yInBuffer];
                            e.AddFrom(de);
                            alls[x * nGridPerSide + y] += dSBuffer[xInBuffer * nGridPerSide + yInBuffer];
                            }
                        }
                        }
                    }   numRearrange++;
                    }
                }
#pragma omp single
                {
                    for (int i = 0; i < nSite; i++)
                    {
                        if (rearrangingStep[i] > 0)
                        {
                            //carry out the rearrangement
                            rearrangingStep[i]++;

                            hasRearranged[i] = 1;
                            // if (rearrangingStep[i] > 4)
                            // {
                            //     rearrangingStep[i] = 0;
                            // }
                            alle[i] = 0.0;
                            alls[i] = sDistribution(rEngine) + mean_softness + 0.2;
                            yieldStrainCoeff[i] = coeffDistribution(rEngine);
                        }
                    }

                    std::cout << "Internal loop, number of rearrangement = " << numRearrange << std::endl;

                }
            }
        }

        int numRe_frame = 0;
        for (int i = 0; i < nSite; i++)
        {
            if (hasRearranged[i] == 1)
            {
                numRe_frame+=1;
            }

        }

        return numRe_frame;
    }




//
//
//     int avalanche(std::string outputPrefix = "")
//     {
//         bool avalancheHappened = false;
//         int nSite = nGridPerSide * nGridPerSide;
//         int nStep = 0;
//         double deltaEnergy;
//         // int loopnum = 0;
//         int numRearrange = 0;
// #pragma omp parallel
//         {
//             // int loopnum = 0;
//             //
//             // while (loopnum == 0)
//             // {
// #pragma omp for schedule(static)
//                 for (int i = 0; i < nSite; i++)
//                     if (startRearranging(alle[i], alls[i]))   // if yield strain cross rearrangement barrier
//                     {
//                         hasRearranged[i] = 1;
//                         rearrangingStep[i] = 1;
//                     }
//
//                 //stop rearrangements that increases energy
//                 for (int i = 0; i < nSite; i++)
//                 {
//                     if (rearrangingStep[i] > 0)
//                     {
//                         //calculate energy difference
//                         deltaEnergy = 0;
//                         int rx = i / nGridPerSide;
//                         int ry = i % nGridPerSide;
// #pragma omp barrier
// #pragma omp for schedule(static) reduction(+ \
//                                            : deltaEnergy)
//                         for (int x = 0; x < nGridPerSide; x++)
//                         {
//                             int xInBuffer = bufferCenter - rx + x;
//                             // while (xInBuffer < 0)
//                             //     xInBuffer += nGridPerSide;
//                             // while (xInBuffer >= nGridPerSide)
//                             //     xInBuffer -= nGridPerSide;
//                             if (xInBuffer>=0 && xInBuffer<nGridPerSide)
//                             {
//                             for (int y = 0; y < nGridPerSide; y++)
//                             {
//                                 int yInBuffer = bufferCenter - ry + y;
//                                 // while (yInBuffer < 0)
//                                 //     yInBuffer += nGridPerSide;
//                                 // while (yInBuffer >= nGridPerSide)
//                                 //     yInBuffer -= nGridPerSide;
//                                 if (yInBuffer>=0 && yInBuffer<nGridPerSide)
//                                 {
//                                 GeometryVector &e = alle[x * nGridPerSide + y];
//                                 GeometryVector &de = dEBuffer[xInBuffer * nGridPerSide + yInBuffer];
//
//                                 if (ry != y || rx != x)
//                                     deltaEnergy += (e + de).Modulus2() - e.Modulus2();
//                                 else
//                                     deltaEnergy -= e.Modulus2();
//                                 }
//                             }
//                             }
//                         }
// #pragma omp single
//                         {
//                             //stop if energy increases
//                             // std::cout<<"de= "<<deltaEnergy<<' ';
//                             if (deltaEnergy > 0)
//                             {
//                                 //std::cout << "rearrangement at stage " << int(rearrangingStep[i]) << " refuted due to energy criteria\n";
//                                 rearrangingStep[i] = 0;
//                                 hasRearranged[i] = 0;
//                             }
//                         }
//                     }
//                 }
//
// #pragma omp barrier
//                 //rearrangement affect other sites parameters
//                 numRearrange = 0;
//                 for (int i = 0; i < nSite; i++)
//                 {
//                     if (rearrangingStep[i] > 0)
//                     {
//                         //update softness and strain
//                         int rx = i / nGridPerSide;
//                         int ry = i % nGridPerSide;
// #pragma omp for schedule(static)
//                         for (int x = 0; x < nGridPerSide; x++)
//                         {
//                             int xInBuffer = bufferCenter - rx + x;
//                             // while (xInBuffer < 0)
//                             //     xInBuffer += nGridPerSide;
//                             // while (xInBuffer >= nGridPerSide)
//                             //     xInBuffer -= nGridPerSide;
//                             if (xInBuffer>=0 && xInBuffer<nGridPerSide)
//                             {
//                             for (int y = 0; y < nGridPerSide; y++)
//                             {
//                                 int yInBuffer = bufferCenter - ry + y;
//                                 // while (yInBuffer < 0)
//                                 //     yInBuffer += nGridPerSide;
//                                 // while (yInBuffer >= nGridPerSide)
//                                 //     yInBuffer -= nGridPerSide;
//                                 //alle[x * nGridPerSide + y] += dEBuffer[xInBuffer * nGridPerSide + yInBuffer];
//                                 if (yInBuffer>=0 && yInBuffer<nGridPerSide)
//                                 {
//                                 GeometryVector &e = alle[x * nGridPerSide + y];
//                                 GeometryVector &de = dEBuffer[xInBuffer * nGridPerSide + yInBuffer];
//                                 e.AddFrom(de);
//                                 alls[x * nGridPerSide + y] += dSBuffer[xInBuffer * nGridPerSide + yInBuffer];
//                                 }
//                             }
//                             }
//                         }
//
//                     }
//                 }
// #pragma omp single
//                 {
//                     for (int i = 0; i < nSite; i++)
//                     {
//                         if (rearrangingStep[i] > 0)
//                         {
//                             //carry out the rearrangement
//                             // rearrangingStep[i]++;
//                             // if (rearrangingStep[i] > 4)
//                             // {
//                             //     rearrangingStep[i] = 0;
//                             // }
//                             // finite lifetime
//                             double dissp_factor = 0;
//                             alle[i].x[0] = alle[i].x[0] *  dissp_factor;
//                             alle[i].x[1] = alle[i].x[1] *  dissp_factor;
//                             // alls[i] = alls[i] + (-0.2 * alls[i] + 0.1) * (1-dissp_factor);
//                             alls[i] = sDistribution(rEngine) + (0.4);
//                             numRearrange++;
//                             //alls[i] = sDistribution(rEngine);
//                             //rearrangingStep[i] = 0;
//                         }
//                     }
//
//                     // if (numRearrange > 0)
//                     // {
//                         // avalancheHappened = true;
//                         // std::cout << "num rearranger in this frame=" << numRearrange;
//                         // double sum = 0.0;
//                         // for (auto &e : this->alle)
//                         //     sum += e.Modulus2();
//                         // std::cout << ", mean energy=" << sum / alle.size();
//
//                         // double sum = 0.0;
//                         // for (auto &s : this->alls)
//                         //     sum += s;
//                         // std::cout << ", mean s=" << sum / alls.size();
//                         // std::cout << std::endl;
//
//                         // if (outputPrefix != std::string(""))
//                         // {
//                         //     std::stringstream ss;
//                         //     ss << outputPrefix << "_step_" << (nStep++);
//                         //     plot(this->rearrangingStep, nGridPerSide, ss.str());
//                         // }
//                     // }
//                 }
//                 // loopnum+=1;
//             // }
//         }
//         return numRearrange;
//     }
};

int main()
{
    const int nGridPerSide = 400;
    gridModel model(nGridPerSide, 1.0);
    model.initialize();
   int numAvalanche = 0;
//    while (numAvalanche < 100)
   int strainstep = 0;

   // std::vector<bool> hasRearranged_collect;
   // hasRearranged_collect.assign(0, nGridPerSide*nGridPerSide);
   double meanS = 0.0;

   while (strainstep<8000 && meanS<1.0)
    {
        //std::cout << "shearing\n";
        model.shear();
        //std::cout << "checking avalanche\n";

        std::stringstream ss;
        ss << "avalanche_" << numAvalanche;

        int numRe = model.avalanche(ss.str());



        // numAvalanche += avalanched;
        // if (avalanched)
        // {
            // std::cout << numAvalanche << "avalanches so far.\n";
          //  plot(model.hasRearranged, nGridPerSide, ss.str());
        // if (strainstep%25==0)
        // {

        double sum = 0.0;
        for (auto &s : model.alls)
           sum += s;

        meanS = sum / model.alls.size();

        //writebinary(model.hasRearranged, nGridPerSide, ss.str());
        writebinary<bool>(model.hasRearranged, nGridPerSide, ss.str());
        std::cout << "Currently at step: " << strainstep << ", Number of rearrangement:"<<  numRe;
        std::cout << ", mean s=" << sum / model.alls.size() << std::endl;
        // }

        strainstep+=1;

        if ((strainstep-1)%50==0)
        {
            writebinary_scalar<double>(model.alls, nGridPerSide, ss.str());
        }


        std::ofstream re_mean_data;
        re_mean_data.open ("data_mean.bin", std::ios::out | std::ios::binary | std::fstream::app);
        double nstep = strainstep;
        re_mean_data.write((char*)&nstep,sizeof(double));
        re_mean_data.write((char*)&meanS,sizeof(double));

        double etotal=0;
        for (auto &e : model.alle)
           etotal += e.Modulus2();

        re_mean_data.write((char*)&etotal,sizeof(double));

        re_mean_data.close();

        // }
    }
    // for (auto &s : model.alls)
    //     std::cout << s << std::endl;
}
