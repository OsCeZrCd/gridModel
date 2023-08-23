#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <ctime>

#include "GeometryVector.h"
#include "kiss_fft.h"
#include "kiss_fftnd.h"
#include "mgl2/mgl.h"
#include <netcdf>

template <typename T>
void plot(const std::vector<T> &data, int nGridPerSide, std::string file)
{
    // mglData x(nGridPerSide, nGridPerSide);
    // for (int i = 0; i < nGridPerSide; i++)
    // {
    //     for (int j = 0; j < nGridPerSide; j++)
    //     {
    //         mreal a = 0.0;
    //         if (data[i * nGridPerSide + j] > 0)
    //             a = 1.0;
    //         else if (i + 1 != nGridPerSide && data[(i + 1) * nGridPerSide + j] > 0)
    //             a = std::max(0.5, a);
    //         else if (i - 1 != -1 && data[(i - 1) * nGridPerSide + j] > 0)
    //             a = std::max(0.5, a);
    //         else if (j + 1 != nGridPerSide && data[i * nGridPerSide + j + 1] > 0)
    //             a = std::max(0.5, a);
    //         else if (j - 1 != -1 && data[i * nGridPerSide + j - 1] > 0)
    //             a = std::max(0.5, a);
    //         x.SetVal(a, i, j);
    //     }
    //     //std::cout<<std::endl;
    // }
    // mglGraph gr;
    // gr.SetSize(1600, 1600);
    // //gr.Aspect(0.75, 1.0);
    // //gr.Colorbar(">kbcyr");
    // gr.Tile(x, "bckyr");
    // gr.WritePNG((file + std::string(".png")).c_str());
}

const double meanSoftness = -0.04;
const double stdSoftness = 0.697;

// below not used
std::vector<double> numericalDsR =
    {
        0.90483741803596,
        1.10517091807565,
        1.349858807576,
        1.64872127070013,
        2.01375270747048,
        2.45960311115695,
        3.00416602394643,
        3.66929666761924,
        4.48168907033806,
        5.4739473917272,
        6.68589444227927,
        8.16616991256765,
        9.97418245481472,
        12.1824939607035,
        14.8797317248728,
        18.1741453694431,
        22.1979512814416,
        27.1126389206579,
        33.1154519586923,
};

std::vector<double> numericalDs =
    {
        0.115407163314733,
        -0.157866369001493,
        -0.122794150458982,
        -0.0111312471488491,
        -0.0123567435148076,
        -0.0195800437378993,
        -0.00140091760436869,
        -0.00346197619971859,
        0.000793079578640183,
        0.00106809091953078,
        0.00112736258290034,
        0.00110702362832669,
        0.000915153303260757,
        0.000483591475874019,
        0.000312333963820199,
        0.000132837059034092,
        4.22819533290796e-05,
        8.57444728501887e-06,
        1.17583919795709e-05,
};

// above not used

class gridModel
{
public:
    int nGridPerSide;
    double lGrid;
    int bufferCenter;
    double meanS_ins;
    double S_std_ins;
    double total_strain;

    std::vector<double> alls, alls_mean;

    //e is elastic strain
    //dEBuffer[0] and [1] correspond to
    // [0] response from converting an xy elastic strain to plastic strain
    // [1] response from converting an xx deviatoric elastic strain to plastic strain
    std::vector<GeometryVector> dEBuffer[::MaxDimension], alle, alle_0, allp;
    std::vector<double> Saverage, Saverage_step, SaverageNum;

    std::vector<int> rearrangingStep;
    std::vector<char> hasRearranged;

    std::mt19937 rEngine;
    std::mt19937 nebEngine;
    //distribution of softness after a rearrangement
    std::normal_distribution<double> sDistribution;
    //distribution of strain initially
    std::normal_distribution<double> eDistribution;
    //distribution of elastic strain after rearrangement
    std::normal_distribution<double> residualStrainDistribution;

    std::uniform_real_distribution<double> PxDistribution;
    std::uniform_real_distribution<double> Nebreprob;
    std::vector<double> yieldStrainPx;

    netCDF::NcFile dumpFile;
    netCDF::NcVar eVar, pVar, sVar, hasRearrangedVar;
    netCDF::NcVar coeffVar;

    //In this version of the program, length of a rearrangement in frames is determined from the intensity
    //int rearrangeFrameLength;

    gridModel(int nGrid, double lGrid, int seed) : rEngine(time(nullptr)),
                                                   nebEngine(time(nullptr)),
                                                   eDistribution(0.0, 0.0001),
                                                   residualStrainDistribution(0.0, 0.0),
                                                   sDistribution(meanSoftness, stdSoftness),
                                                   nGridPerSide(nGrid), lGrid(lGrid), meanS_ins(meanSoftness), S_std_ins(stdSoftness),total_strain(0)
    {
        this->allocate();
        this->getBuffer();
        this->initialize();
    }

    void openNewDumpFile(const std::string &filename)
    {
        dumpFile.open(filename, netCDF::NcFile::replace);
        int dim = 2;
        std::vector<netCDF::NcDim> dims, strainDims;

        netCDF::NcDim framedim = dumpFile.addDim("frames");
        dims.push_back(framedim);
        strainDims.push_back(framedim);
        for (int i = 0; i < dim; i++)
        {
            std::stringstream ss;
            ss << "dim_" << i;
            netCDF::NcDim temp = dumpFile.addDim(ss.str(), nGridPerSide);
            dims.push_back(temp);
            strainDims.push_back(temp);
        }
        strainDims.push_back(dumpFile.addDim("strainComponents", ::MaxDimension));
        eVar = dumpFile.addVar("strain", netCDF::ncDouble, strainDims);
        sVar = dumpFile.addVar("softness", netCDF::ncDouble, dims);
        hasRearrangedVar = dumpFile.addVar("hasRearranged", netCDF::ncByte, dims);
        coeffVar = dumpFile.addVar("yieldStrainPx", netCDF::ncDouble, dims);
        pVar = dumpFile.addVar("plastic", netCDF::ncDouble, strainDims);

        eVar.setCompression(true, true, 9);
        pVar.setCompression(true, true, 9);
        sVar.setCompression(true, true, 9);
        hasRearrangedVar.setCompression(true, true, 9);
        coeffVar.setCompression(true, true, 9);
    }
    void openExistingDumpFile(const std::string &filename)
    {
        dumpFile.open(filename, netCDF::NcFile::write);
        eVar = dumpFile.getVar("strain");
        pVar = dumpFile.getVar("plastic");
        sVar = dumpFile.getVar("softness");
        hasRearrangedVar = dumpFile.getVar("hasRearranged");
        coeffVar = dumpFile.getVar("yieldStrainPx");
    }
    void dump(bool writeStrain = false, bool writePlastic = false, bool writeSoftness = false, bool writeYieldStrainCoeff = false, bool writeHasRearranged = true)
    {
        int dim = 2;

        int framesAlreadyWritten = eVar.getDim(0).getSize();
        if (writeStrain)
        {
            std::vector<size_t> startp, countp;
            startp.push_back(framesAlreadyWritten); //start from the end of the previous frame
            for (int i = 0; i < dim; i++)
                startp.push_back(0);
            startp.push_back(0);
            countp.push_back(1); //write one frame
            for (int i = 0; i < dim; i++)
                countp.push_back(nGridPerSide);
            countp.push_back(::MaxDimension);
            eVar.putVar(startp, countp, alle.data());
        }

        if (writePlastic)
        {
            std::vector<size_t> startp, countp;
            startp.push_back(framesAlreadyWritten); //start from the end of the previous frame
            for (int i = 0; i < dim; i++)
                startp.push_back(0);
            startp.push_back(0);
            countp.push_back(1); //write one frame
            for (int i = 0; i < dim; i++)
                countp.push_back(nGridPerSide);
            countp.push_back(::MaxDimension);
            pVar.putVar(startp, countp, allp.data());
        }

        if (writeSoftness)
        {
            std::vector<size_t> startp, countp;
            startp.push_back(framesAlreadyWritten);
            for (int i = 0; i < dim; i++)
                startp.push_back(0);
            countp.push_back(1);
            for (int i = 0; i < dim; i++)
                countp.push_back(nGridPerSide);
            sVar.putVar(startp, countp, alls.data());
        }
        if (writeYieldStrainCoeff)
        {
            std::vector<size_t> startp, countp;
            startp.push_back(framesAlreadyWritten);
            for (int i = 0; i < dim; i++)
                startp.push_back(0);
            countp.push_back(1);
            for (int i = 0; i < dim; i++)
                countp.push_back(nGridPerSide);
            coeffVar.putVar(startp, countp, yieldStrainPx.data());
        }
        if (writeHasRearranged)
        {
            std::vector<size_t> startp, countp;
            startp.push_back(framesAlreadyWritten);
            for (int i = 0; i < dim; i++)
                startp.push_back(0);
            countp.push_back(1);
            for (int i = 0; i < dim; i++)
                countp.push_back(nGridPerSide);
            hasRearrangedVar.putVar(startp, countp, (void *)(hasRearranged.data()));
        }
    }

    double deviatoricYieldStrain(int i)
    {
        double s = alls[i];

        double meanYieldStrain = 0.01176 - 0.002451 * s;
        // double meanYieldStrain = 0.01249 - 0.002088 * s;
        //weibull distribution
        // double k = 1.6 + 0.4 * s ;
        double k = 2.1;
        // double k = 5.5;

        double lambda = meanYieldStrain / std::tgamma(1.0 + 1.0 / k);
        double yieldStrain = std::pow(-1.0 * std::log(1.0 - yieldStrainPx[i]), 1.0 / k) * lambda;
        yieldStrain = std::max(yieldStrain, 0.0);
        return yieldStrain;
    }

    bool startRearranging(int i)
    {
        double yieldStrain = deviatoricYieldStrain(i);
        auto e = alle[i];
        return e.Modulus2() > yieldStrain * yieldStrain;
    }
    double xyStrainDistanceToRarranging(int i)
    {
        double yieldStrain = deviatoricYieldStrain(i);
        auto e = alle[i];
        if (e.Modulus2() > yieldStrain * yieldStrain)
            return 0.0;
        else
            return std::sqrt(yieldStrain * yieldStrain - e.x[1] * e.x[1]) - e.x[0];
    }
    double minimumXyStrainDistanceToRarranging()
    {
        double minimum = std::numeric_limits<double>::max();
        for (int i = 0; i < this->alle.size(); i++)
            minimum = std::min(minimum, this->xyStrainDistanceToRarranging(i));
        return minimum;
    }

    double dsFromRearranger(double dx, double dy, double r, double s, const GeometryVector & rearrangingIntensity, std::mt19937 &engine, double meansneb, double emodold, double emodnew, int restep)
    {
      // const double angularContributionCoefficient=0.6044;
      const double angularContributionCoefficient=1.741*2.0;
      // const double angularContributionCoefficient=1.6;
      // if (r == 0.0)
          // return 0.0; // delta S of the rearranger is processed separately

      double meanContribution = 0.0;
      if (r==0)
      {
          // meanContribution = 0.0;
          // meanContribution = -0.05221 -  0.0000771860657566;
          meanContribution = 0.1125-0.0000235396432117;
      }
      // else if (r < 30)
      else
      {
          // meanContribution += 0.0;
          // meanContribution += 0.1717*std::exp(-r/0.6402) - 0.000186 +  0.00016645;
          // meanContribution += 0.05752 *std::pow(r, -3.1) -  0.000034202094;
          meanContribution += 0.1125 * std::exp(-r/0.6465)-0.0000235396432117;
          //contribution from volumetric strain
          meanContribution -= angularContributionCoefficient * rearrangingIntensity.x[0] / r / r * ( - 0.0 * std::cos(4.0 * std::atan2(dy, dx)) +  std::sin(2.0 * std::atan2(dy, dx))) / 3.0;
          meanContribution -= angularContributionCoefficient * rearrangingIntensity.x[1] / r / r * (  0.0 * std::cos(4.0 * std::atan2(dy, dx)) +   std::cos(2.0 * std::atan2(dy, dx))) / 3.0;
      }
      // else {
      //     meanContribution = 0.0;
      // }

      // const double alpha = 0.01, beta = -3.1;
      const double alpha = 0.0931;
      double restore = 0.0;
      double softnessRestoringCoefficient = 0;
      if (r > 0 && r < 15)
      {
          // softnessRestoringCoefficient = alpha * std::pow(r, beta);
          // softnessRestoringCoefficient = alpha * std::pow(r, beta);
          softnessRestoringCoefficient = alpha*std::exp(-r/0.6465);
          restore = softnessRestoringCoefficient * (meansneb - s);
      }

      if (r==0)
      {
         // restore = 0.0923 * (this->meanS_ins - s);
         restore = alpha * (meansneb - s);
         softnessRestoringCoefficient = alpha;
      }
      // else if (r < 20)
      // {
      //     double softnessRestoringCoefficient = -1e-5;
      //     restore = softnessRestoringCoefficient * (meanSoftness - s);
      // }

      double harmonicDiffusion = 0.0;
      if (r > 0 && r < 15)
      {
          double rcoefficient = std::sqrt(softnessRestoringCoefficient*(2.0-softnessRestoringCoefficient)*this->S_std_ins*this->S_std_ins);
          std::normal_distribution<double> noiseDistribution(0.0, rcoefficient);
          harmonicDiffusion = noiseDistribution(engine);
      }

      if (r == 0)
      {
          double rcoefficient = std::sqrt(softnessRestoringCoefficient*(2.0-softnessRestoringCoefficient)*this->S_std_ins*this->S_std_ins);
          std::normal_distribution<double> noiseDistribution(0.0, rcoefficient);
          harmonicDiffusion = noiseDistribution(engine);
      }

      // deviatoric strain effect
      double ds_devia = 0;
      // if (r>=0)
      // {
      double strain_energy_diff = (emodnew) - (emodold);
      ds_devia = strain_energy_diff * 339;
      // }

      double delta_s = 0;
      if (restep == 1)
      {
        delta_s = meanContribution + restore + harmonicDiffusion + ds_devia;
      } else if (restep == 2)
      {
        delta_s = ds_devia;
      }
      return delta_s;
    }

    void getBuffer()
    {
        bufferCenter = nGridPerSide / 2;

        for (int i = 0; i < MaxDimension; i++)
            dEBuffer[i].resize(nGridPerSide * nGridPerSide);

        double meanStrainDecrement = 1.0 / nGridPerSide / nGridPerSide;
        double factor[MaxDimension];

        {
            //strain buffer calculated by Fourier transform, xy to xy response
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
            //the interaction kernel at the rearranging site is -1.0 because
            //the rearranging site loses all of the redistributed strain
            factor[0] = (-1.0 + meanStrainDecrement) / outbuf[0].r;

            for (int i = 0; i < nGridPerSide; i++)
                for (int j = 0; j < nGridPerSide; j++)
                {
                    int index = i * nGridPerSide + j;
                    // int ii = i - bufferCenter;
                    // while (ii < 0)
                    //     ii += nGridPerSide;
                    // int jj = j - bufferCenter;
                    // while (jj < 0)
                    //     jj += nGridPerSide;

                    int ii = i + bufferCenter;
                    while (ii >= nGridPerSide)
                        ii -= nGridPerSide;
                    int jj = j + bufferCenter;
                    while (jj >= nGridPerSide)
                        jj -= nGridPerSide;

                    int index2 = ii * nGridPerSide + jj;
                    dEBuffer[0][index2] = GeometryVector(outbuf[index].r * factor[0] - meanStrainDecrement, 0.0);
                }

            free(st);
            delete[] outbuf;
            delete[] inbuf;
        }
        {
            //strain buffer calculated by Fourier transform, xx to xx response
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

                    double temp = (pm * pm - qn * qn);
                    inbuf[i * nGridPerSide + j].r = -1 * temp * temp / q2 / q2;
                    inbuf[i * nGridPerSide + j].i = 0;
                }
            inbuf[0].r = 0;

            int temp[2] = {nGridPerSide, nGridPerSide};
            kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, 2, true, nullptr, nullptr);
            kiss_fftnd(st, inbuf, outbuf);
            //fill in dEBuffer
            //the interaction kernel at the rearranging site is -1.0 because
            //the rearranging site loses all of the redistributed strain
            factor[1] = (-1.0 + meanStrainDecrement) / outbuf[0].r;

            for (int i = 0; i < nGridPerSide; i++)
                for (int j = 0; j < nGridPerSide; j++)
                {
                    int index = i * nGridPerSide + j;
                    // int ii = i - bufferCenter;
                    // while (ii < 0)
                    //     ii += nGridPerSide;
                    // int jj = j - bufferCenter;
                    // while (jj < 0)
                    //     jj += nGridPerSide;

                    int ii = i + bufferCenter;
                    while (ii >= nGridPerSide)
                        ii -= nGridPerSide;
                    int jj = j + bufferCenter;
                    while (jj >= nGridPerSide)
                        jj -= nGridPerSide;

                    int index2 = ii * nGridPerSide + jj;
                    dEBuffer[1][index2] = GeometryVector(0.0, outbuf[index].r * factor[1] - meanStrainDecrement);
                }

            free(st);
            delete[] outbuf;
            delete[] inbuf;
        }

        {
            //strain buffer calculated by Fourier transform, xy-to-xx and xx-to-xy component, which has the same formula
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
                    // int ii = i - bufferCenter;
                    // while (ii < 0)
                    //     ii += nGridPerSide;
                    // int jj = j - bufferCenter;
                    // while (jj < 0)
                    //     jj += nGridPerSide;
                    int ii = i + bufferCenter;
                    while (ii >= nGridPerSide)
                        ii -= nGridPerSide;
                    int jj = j + bufferCenter;
                    while (jj >= nGridPerSide)
                        jj -= nGridPerSide;

                    int index2 = ii * nGridPerSide + jj;
                    dEBuffer[0][index2].x[1] = outbuf[index].r * factor[0];
                    dEBuffer[1][index2].x[0] = outbuf[index].r * factor[1];
                }

            free(st);
            delete[] outbuf;
            delete[] inbuf;
        }

        //debug temp
        // std::cout << std::scientific;
        // std::cout << factor[0] << " " << factor[1] << std::endl;
        // std::vector<double> temp(nGridPerSide*nGridPerSide, 0.0);
        // for (int i = 0; i < nGridPerSide; i++)
        // {
        //     for (int j = 0; j < nGridPerSide; j++)
        //     {
        //         int index = i * nGridPerSide + j;
        //         temp[index]=dEBuffer[1][index].x[1];
        //     }
        // }
        // plot(temp, nGridPerSide, "e11");
        // exit(0);
    }
    void allocate()
    {
        int nSite = nGridPerSide * nGridPerSide;
        alle.resize(nSite);
        allp.resize(nSite);
        alle_0.resize(nSite);
        alls.resize(nSite);
        yieldStrainPx.resize(nSite);
        hasRearranged.resize(nSite);
        alls_mean.resize(nSite);
        rearrangingStep.resize(nSite);
        Saverage.resize(15);
        Saverage_step.resize(15);
        SaverageNum.resize(15);
    }
    void initialize()
    {
        int nSite = nGridPerSide * nGridPerSide;
        for (int i = 0; i < nSite; i++)
        {
            this->alle[i].x[0] = this->eDistribution(this->rEngine);
            this->alle[i].x[1] = this->eDistribution(this->rEngine);
            this->alle_0[i].x[0] = this->alle[i].x[0];
            this->alle_0[i].x[1] = this->alle[i].x[1];
            this->allp[i].x[0] = 0;
            this->allp[i].x[1] = 0;
            this->alls[i] = this->sDistribution(this->rEngine);
            this->alls_mean[i] = this->alls[i];
            this->yieldStrainPx[i] = this->PxDistribution(this->rEngine);
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
        }
        for (int i = 0; i < 15; i++)
        {
          this->Saverage[i] = meanSoftness;
          this->Saverage_step[i] = 0;
          this->SaverageNum[i] = 0;
        }
    }
    //read the dump file, initialize from the last frame where both softness and strain are recorded
    void initializeFromDumpFile(const std::string &filename)
    {
        int dim = 2;
        netCDF::NcFile inputFile(filename, netCDF::NcFile::read);
        netCDF::NcVar readEVar = inputFile.getVar("strain");
        netCDF::NcVar readSVar = inputFile.getVar("softness");
        netCDF::NcVar readCVar = inputFile.getVar("yieldStrainPx");

        int framesAlreadyWritten = readEVar.getDim(0).getSize();
        for (int framesToRead = framesAlreadyWritten - 1; framesToRead >= 0; framesToRead--)
        {
            {
                std::vector<size_t> startp, countp;
                startp.push_back(framesToRead); //start from the end of the previous frame
                for (int i = 0; i < dim; i++)
                    startp.push_back(0);
                startp.push_back(0);
                countp.push_back(1); //write one frame
                for (int i = 0; i < dim; i++)
                    countp.push_back(nGridPerSide);
                countp.push_back(::MaxDimension);
                readEVar.getVar(startp, countp, alle.data());
            }
            {
                std::vector<size_t> startp, countp;
                startp.push_back(framesToRead);
                for (int i = 0; i < dim; i++)
                    startp.push_back(0);
                countp.push_back(1);
                for (int i = 0; i < dim; i++)
                    countp.push_back(nGridPerSide);
                readSVar.getVar(startp, countp, alls.data());
            }
            {
                std::vector<size_t> startp, countp;
                startp.push_back(framesToRead);
                for (int i = 0; i < dim; i++)
                    startp.push_back(0);
                countp.push_back(1);
                for (int i = 0; i < dim; i++)
                    countp.push_back(nGridPerSide);
                readCVar.getVar(startp, countp, yieldStrainPx.data());
            }

            //netCDF use a fillValue to indicate non-written value
            //assume the zeroth element of alle and alls are fillValue
            //if any of the read data is not fillValue, frame data is valid, mission complete
            //otherwise, frame data is invalid, see if an earlier frame is valid by letting the loop continue;
            double &fillEValue = alle[0].x[0];
            double &fillSValue = alls[0];
            double &fillCValue = yieldStrainPx[0];
            for (int i = 1; i < nGridPerSide * nGridPerSide; i++)
            {
                if (alle[i].x[0] != fillEValue && alle[i].x[1] != fillEValue && alls[i] != fillSValue && yieldStrainPx[i] != fillCValue)
                {
                    std::cout << "initialized from dump file: " << filename << ", at frame " << framesToRead << std::endl;
                    return;
                }
            }
        }
        //should not reach here
        std::cerr << "gridModel::initializeFromDumpFile failed : found no frame with both softness and strain recorded.\n";
    }


    void initializeFromDumpFile_softness(const std::string &filename)
    {
        int dim = 2;
        netCDF::NcFile inputFile(filename, netCDF::NcFile::read);
        netCDF::NcVar readSVar = inputFile.getVar("softness");

      //  int framesAlreadyWritten = readSVar.getDim(0).getSize();
        int framesAlreadyWritten = 1;
        std::cout << "frames written" << framesAlreadyWritten << std::endl;

        for (int framesToRead = framesAlreadyWritten - 1; framesToRead >= 0; framesToRead--)
        {

            std::vector<size_t> startp, countp;
          //  startp.push_back(framesToRead);
            for (int i = 0; i < dim; i++)
                startp.push_back(0);
         //      countp.push_back(1);
            for (int i = 0; i < dim; i++)
                countp.push_back(nGridPerSide);
            readSVar.getVar(startp, countp, alls.data());

            //netCDF use a fillValue to indicate non-written value
            //assume the zeroth element of alle and alls are fillValue
            //if any of the read data is not fillValue, frame data is valid, mission complete
            //otherwise, frame data is invalid, see if an earlier frame is valid by letting the loop continue;
            double &fillSValue = alls[0];
            for (int i = 1; i < nGridPerSide * nGridPerSide; i++)
            {
                if (alls[i] != fillSValue)
                {
                    std::cout << "initialized correlated softness field from the dump file: " << filename << ", at frame " << framesToRead << std::endl;
                    return;
                }
            }
        }
        //should not reach here
        std::cerr << "gridModel::initializeFromDumpFile failed : found no frame with both softness recorded.\n";
    }

    //increase the first component of the strain tensor by the desired amount
    void shear(double strain = 1e-6)
    {
        int nSite = nGridPerSide * nGridPerSide;
        this->total_strain += strain;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < nSite; i++)
        {
            this->alle_0[i].x[0] = this->alle[i].x[0];
            this->alle_0[i].x[1] = this->alle[i].x[1];
            this->alle[i].x[0] += strain;
            // double strain_energy_diff = this->alle[i].x[0] * this->alle[i].x[0] + this->alle[i].x[1] * this->alle[i].x[1]  - (this->alle_0[i].x[0] * this->alle_0[i].x[0] + this->alle_0[i].x[1] * this->alle_0[i].x[1] );
            double strain_energy_diff = (this->alle[i].x[0] * this->alle[i].x[0] + this->alle[i].x[1] * this->alle[i].x[1]) - ((this->alle_0[i].x[0] * this->alle_0[i].x[0] + this->alle_0[i].x[1] * this->alle_0[i].x[1]));
            // double strain_energy_diff = std::sqrt(this->alle[i].x[0] * this->alle[i].x[0] + this->alle[i].x[1] * this->alle[i].x[1]) - std::sqrt((this->alle_0[i].x[0] * this->alle_0[i].x[0] + this->alle_0[i].x[1] * this->alle_0[i].x[1]));
            // this->alls[i] += strain_energy_diff * 1.094 + strain * 0.3223;
            this->alls[i] += strain_energy_diff * 339.0 + strain * 1.741*2.0/3.0; //0.427

            // this->alls[i] += strain * 1.5;
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
        }
    }
    bool avalanche(std::string outputPrefix = "")
    {
        bool avalancheHappened = false;
        int nSite = nGridPerSide * nGridPerSide;
        int nStep = 0;

        //if rearranging, the value is how much strain is redistributed per frame
        std::vector<GeometryVector> rearrangingIntensity;
        rearrangingIntensity.resize(nSite);
        std::vector<int> rearrangeFrameLength(nSite, 0);

       // back up previous strain
       // for (int i = 0; i < nSite; i++)
       //   {
       //     this->alle_0[i].x[0] = this->alle[i].x[0];
       //     this->alle_0[i].x[1] = this->alle[i].x[1];
       //   }

#pragma omp parallel
        {
            std::mt19937 threadEngine;
#pragma omp critical(random)
            {
                threadEngine.seed(rEngine());
            }
            int numRearrange = 1;
            while (numRearrange > 0 && numRearrange < 4000)
            {
#pragma omp single
                {
                    for (int i = 0; i < nSite; i++)
                      {
                        this->alle_0[i].x[0] = this->alle[i].x[0];
                        this->alle_0[i].x[1] = this->alle[i].x[1];
                        this->alls_mean[i] = this->alls[i];
                        // if (rearrangingStep[i] == 0 && startRearranging(i))
                        if (startRearranging(i))
                        {
                            rearrangingStep[i] = 1;
                            GeometryVector residual(this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine));
                            GeometryVector totalIntensity = (alle[i] - residual);
                            // rearrangeFrameLength[i] = std::max(int(std::ceil(std::sqrt(totalIntensity.Modulus2()) / 0.1)), 1);
                            rearrangeFrameLength[i] = 1;
                            rearrangingIntensity[i] = totalIntensity * (1.0 / rearrangeFrameLength[i]);
                            // GeometryVector &ep = allp[ind_neb];
                            for (int j = 0; j < MaxDimension; j++)
                            {
                                // ep.AddFrom(rearrangingIntensity[i].x[j]);
                                this->allp[i].x[j] = this->allp[i].x[j] + rearrangingIntensity[i].x[j];
                            }


                            int indx = i / nGridPerSide;
                            int indy = i % nGridPerSide;
                            for (int ix = -9; ix < 10 ; ix++)
                             {
                               for (int iy = -9; iy < 10; iy++)
                               {
                                 int nebx = indx + ix;
                                 int neby = indy + iy;

                                 if (nebx<0)
                                  {
                                   nebx = nebx + nGridPerSide;
                                  }
                                 if (neby<0)
                                  {
                                    neby = neby + nGridPerSide;
                                  }
                                 if (nebx>=nGridPerSide)
                                  {
                                    nebx = nebx - nGridPerSide;
                                  }
                                 if (neby>=nGridPerSide)
                                  {
                                    neby = neby - nGridPerSide;
                                  }

                                  int ind_neb = neby + nebx * nGridPerSide;
                                  double prob = this->Nebreprob(this->nebEngine);
                                  double r_neb = std::sqrt( (double)ix * (double)ix + (double)iy * (double)iy);
                                  // double coeffA = 0.6884 * std::exp( 0.311 * this->alls[ind_neb]);
                                  // double coeffA = 2.6884 * std::exp( 0.311 * this->alls[ind_neb]);
                                  // double coeffB = -0.09 * this->alls[ind_neb] + 1.002;
                                  // coeffB = std::max(coeffB, 0.0000001);
                                  double neb_prob = 1.0 * std::exp(- r_neb / 1.12);
                                  // neb_prob = neb_prob / 2.0;

                                 if (prob < neb_prob)
                                 {
                                   if (rearrangingStep[ind_neb] == 0)
                                   {
                                     rearrangingStep[ind_neb] = 2;
                                     GeometryVector residual(this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine));
                                     GeometryVector totalIntensity = (alle[ind_neb] - residual);
                                     // rearrangeFrameLength[i] = std::max(int(std::ceil(std::sqrt(totalIntensity.Modulus2()) / 0.1)), 1);
                                     rearrangeFrameLength[ind_neb] = 1;
                                     rearrangingIntensity[ind_neb] = totalIntensity * (1.0 / rearrangeFrameLength[ind_neb]);
                                     // GeometryVector &ep = allp[ind_neb];
                                     for (int j = 0; j < MaxDimension; j++)
                                     {
                                         // ep.AddFrom(rearrangingIntensity[ind_neb].x[j]);
                                         this->allp[ind_neb].x[j] = this->allp[ind_neb].x[j] + rearrangingIntensity[ind_neb].x[j];
                                     }
                                   }
                                 }
                               }
                             }
                        }
                      }

                    for (int i = 0; i < 15; i++)
                    {
                      this->Saverage_step[i] = 0;
                      this->SaverageNum[i] = 0;
                    }

                    double sum = 0.0;
                    for (auto &s : this->alls)
                        sum += s;
                    this->meanS_ins =  sum / alls.size();

                    double sum_std = 0.0;
                    for (auto &s : this->alls)
                        sum_std += (s-this->meanS_ins) * (s-this->meanS_ins);
                    this->S_std_ins = std::sqrt(sum_std / alls.size());

                }

#pragma omp barrier
                int numRearrange0 = 0;
                for (int i = 0; i < nSite; i++)
                {
                    if (rearrangingStep[i] == 1)
                    {
                        //update softness and strain
                        int rx = i / nGridPerSide;
                        int ry = i % nGridPerSide;
                #pragma omp for schedule(static)
                        for (int x = 0; x < nGridPerSide; x++)
                        {
                            int xInBuffer = bufferCenter - rx + x;
                            while (xInBuffer < 0)
                                xInBuffer += nGridPerSide;
                            while (xInBuffer >= nGridPerSide)
                                xInBuffer -= nGridPerSide;
                            for (int y = 0; y < nGridPerSide; y++)
                            {
                                int yInBuffer = bufferCenter - ry + y;
                                while (yInBuffer < 0)
                                    yInBuffer += nGridPerSide;
                                while (yInBuffer >= nGridPerSide)
                                    yInBuffer -= nGridPerSide;
                                //alle[x * nGridPerSide + y] += dEBuffer[xInBuffer * nGridPerSide + yInBuffer];

                                //softness has a restoring force
                                // double restore = 0.0;
                                double dx = (xInBuffer - bufferCenter) * lGrid;
                                double dy = (yInBuffer - bufferCenter) * lGrid;
                                double r = std::sqrt(dx * dx + dy * dy);
                                int idx = (int)std::floor(r);

                                if (idx < 15)
                                {
                                  this->Saverage_step[idx] += this->alls[x * nGridPerSide + y];
                                  this->SaverageNum[idx] += 1.0;
                                }

                            }
                        }
                        numRearrange0++;
                    }
                }

#pragma omp barrier

#pragma omp single
                  {
                      if (numRearrange0>0)
                      {
                          for (int i=0; i<15 ; i++)
                           {
                               // this->Saverage[i] = this->Saverage_step[i] / this->SaverageNum[i]; // this is an instantaneous measurement
                               if (numRearrange0 < 100.0 )
                               {
                                  this->Saverage[i] = (1.0-(double)numRearrange0/100.0) * this->Saverage[i] + ((double)numRearrange0/100.0) * this->Saverage_step[i] / this->SaverageNum[i];
                               } else
                               {
                                  this->Saverage[i] = this->Saverage_step[i] / this->SaverageNum[i];
                               }
                            // std::cout << "distance: " << i <<"Number of data: " << this->SaverageNum[i] << "Average Value: " << this->Saverage[i] << std::endl;
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
                            while (xInBuffer < 0)
                                xInBuffer += nGridPerSide;
                            while (xInBuffer >= nGridPerSide)
                                xInBuffer -= nGridPerSide;
                            for (int y = 0; y < nGridPerSide; y++)
                            {
                                int yInBuffer = bufferCenter - ry + y;
                                while (yInBuffer < 0)
                                    yInBuffer += nGridPerSide;
                                while (yInBuffer >= nGridPerSide)
                                    yInBuffer -= nGridPerSide;
                                //alle[x * nGridPerSide + y] += dEBuffer[xInBuffer * nGridPerSide + yInBuffer];
                                GeometryVector &e = alle[x * nGridPerSide + y];

                                double emod_old = alle[x * nGridPerSide + y].Modulus2();

                                for (int j = 0; j < MaxDimension; j++)
                                {
                                    GeometryVector &de = dEBuffer[j][xInBuffer * nGridPerSide + yInBuffer];
                                    e.AddFrom(rearrangingIntensity[i].x[j] * de);
                                }

                                double emod_new = alle[x * nGridPerSide + y].Modulus2();

                                //softness has a restoring force
                                double dx = (xInBuffer - bufferCenter) * lGrid;
                                double dy = (yInBuffer - bufferCenter) * lGrid;
                                double r = std::sqrt(dx * dx + dy * dy);
                                double meanS_neb = 0;
                                int idx = (int)std::floor(r);
                                if (idx < 15)
                                {
                                  meanS_neb = this->Saverage[idx];
                                }
                                else
                                {
                                  meanS_neb = this->meanS_ins;
                                }

                               //  double weight = 0;
                               //  double decayr = 1.0;
                               //  if (r<15)
                               //  {
                               //  int indx = x;
                               //  int indy = y;
                               //  for (int ix = -4; ix < 5 ; ix++)
                               //  {
                               //    for (int iy = -4; iy < 5; iy++)
                               //    {
                               //
                               //      int nebx = indx + ix;
                               //      int neby = indy + iy;
                               //
                               //      if (nebx<0)
                               //       {
                               //        nebx = nebx + nGridPerSide;
                               //       }
                               //      if (neby<0)
                               //       {
                               //         neby = neby + nGridPerSide;
                               //       }
                               //      if (nebx>=nGridPerSide)
                               //       {
                               //         nebx = nebx - nGridPerSide;
                               //       }
                               //      if (neby>=nGridPerSide)
                               //       {
                               //         neby = neby - nGridPerSide;
                               //       }
                               //
                               //      int ind_neb = neby + nebx * nGridPerSide;
                               //      meanS_neb += this->alls_mean[ind_neb] * std::exp(-r/decayr) ;
                               //      weight += std::exp(-r/decayr) ;
                               //
                               //    }
                               //  }
                               // }
                               //
                               // meanS_neb = (meanS_neb-this->alls_mean[x * nGridPerSide + y]) / (weight-1);

                              double ds = dsFromRearranger(dx, dy, r, alls[x * nGridPerSide + y], rearrangingIntensity[i], threadEngine, meanS_neb, emod_old, emod_new, rearrangingStep[i]);

                              // double strain_energy_diff = this->alle[x * nGridPerSide + y].x[0] * this->alle[x * nGridPerSide + y].x[0] - (this->alle_0[x * nGridPerSide + y].x[0] * this->alle_0[x * nGridPerSide + y].x[0]);
                              alls[x * nGridPerSide + y] += 1.0 * ds ;
                            }
                        }
                        numRearrange++;
                    }
                }
#pragma omp single
                {
                    for (int i = 0; i < nSite; i++)
                    {
                        // double strain_energy_diff = this->alle[i].x[0] * this->alle[i].x[0] - (this->alle_0[i].x[0] * this->alle_0[i].x[0]);
                        // double strain_energy_diff = this->alle[i].x[0] * this->alle[i].x[0] + this->alle[i].x[1] * this->alle[i].x[1] - (this->alle_0[i].x[0] * this->alle_0[i].x[0] + this->alle_0[i].x[1] * this->alle_0[i].x[1]);
                        // double strain_energy_diff = std::sqrt(this->alle[i].x[0] * this->alle[i].x[0] + this->alle[i].x[1] * this->alle[i].x[1]) - std::sqrt((this->alle_0[i].x[0] * this->alle_0[i].x[0] + this->alle_0[i].x[1] * this->alle_0[i].x[1]));
                        // // this->alls[i] += strain_energy_diff * 1.094;
                        // this->alls[i] += strain_energy_diff * 0.427;
                        // this->alle_0[i].x[0] = this->alle[i].x[0];
                        // this->alle_0[i].x[1] = this->alle[i].x[1];
                        if (rearrangingStep[i] > 0)
                        {
                            //rearrangement has a fixed number of steps
                            rearrangingStep[i]++;
                            if (rearrangingStep[i] > rearrangeFrameLength[i])
                            {
                                rearrangeFrameLength[i] = 0;

                                if (rearrangingStep[i]==2)
                                {
                                  hasRearranged[i] = 1;
                                } else if (rearrangingStep[i]==3)
                                {
                                  hasRearranged[i] = 2;
                                }

                                rearrangingStep[i] = 0;

                                //alls[i] = sDistribution(rEngine);
                                yieldStrainPx[i] = PxDistribution(rEngine);

                                //my simulation suggests this
                                // double dsCenter = std::min(-0.2 - 0.13 * alls[i], 0.25);
                                // alls[i] += dsCenter;
                            }
                        }
                    }

                    if (numRearrange > 0)
                    {
                        avalancheHappened = true;

                        std::cout << "Total strain" << this->total_strain <<", num rearranger in this frame=" << numRearrange;

                        double sum = 0.0;
                        for (auto &s : this->alls)
                            sum += s;
                        std::cout << ", mean s=" << sum / alls.size() << ", std s=" << this->S_std_ins;
                        std::cout << std::endl;

                        // if (sum / alls.size()>4.0)
                        // {
                        //   numRearrange = 0;
                        // }

                        if (outputPrefix != std::string(""))
                        {
                            std::stringstream ss;
                            ss << outputPrefix << "_step_" << (nStep++);
                            // plot(this->rearrangingStep, nGridPerSide, ss.str());
                        }
                    }
                }
            }
        }
        return avalancheHappened;
    }
};

inline bool fileExists(const std::string &name)
{
    std::ifstream f(name.c_str());
    return f.good();
}
int main()
{
    //make sure that the compiler does not add blank space in class GeometryVector
    //if the compiler does that, dumping may not work.
    assert(sizeof(GeometryVector) == MaxDimension * sizeof(double));

    const std::string ncFileName = "dump.nc";

    const int nGridPerSide = 119;
    // int seed, dumpLevel;
    // std::cin >> seed >> dumpLevel;
    int seed = 0;
    int dumpLevel = 10;
    gridModel model(nGridPerSide, 1.0, seed);
    if (fileExists(ncFileName))
    {
        model.initializeFromDumpFile(ncFileName);
        model.openExistingDumpFile(ncFileName);
    }
    else
    {
        model.openNewDumpFile(ncFileName);
        std::cout << "could not find netCDF dump file " << ncFileName << ", create a new one.\n";
        // const std::string ncFileName = "dump_sini.nc";
        // model.initializeFromDumpFile_softness(ncFileName);
    }

    int numAvalanche = 0;
    std::fstream strainFile("xyStrain.txt", std::fstream::out);
    double totalExternalStrain = 0.0;
    double meanS = 0.0;
    while (totalExternalStrain < 0.05 && !fileExists("stop.txt") && meanS < 4.0)
    {
        double strain = model.minimumXyStrainDistanceToRarranging() + 1e-10;
        // double strain = 0.000001;
        model.shear(strain);
        totalExternalStrain += strain;

        std::stringstream ss;
        ss << "avalanche_" << numAvalanche;

        bool avalanched;
        avalanched = model.avalanche("");

        double sumS_loop = 0.0;
        for (auto s: model.alls)
            sumS_loop += s;
        meanS = sumS_loop/model.alls.size();

        auto outputStrainFunc = [&]() -> void {
            double sum = 0.0;
            for (auto s : model.alle)
                sum += s.x[0];
            double sumS = 0.0;
            for (auto s: model.alls)
                sumS += s;

            if (avalanched)
            {
              strainFile << totalExternalStrain << ' ' << sum / model.alle.size() << ' ' << sumS/model.alls.size() << ' ' << 1.0 << std::endl;
            } else
            {
              strainFile << totalExternalStrain << ' ' << sum / model.alle.size() << ' ' << sumS/model.alls.size() << ' ' << 0.0 << std::endl;
            }
        };
        // outputStrainFunc();



        if (avalanched)
        {
            std::cout << numAvalanche << "avalanches so far.\n";
            if(dumpLevel>=10)
                model.dump(true, true, true, true, true);
            else if(dumpLevel>=5)
            {
                if(numAvalanche%100==0)
                    model.dump(true, true, true, true, true);
                else
                    model.dump(false, false, false, false, true);
            }
        }
        else
        {
            //the shear strain should be enough to trigger at least one rearrangement
            std::cerr << "Error in main : expected rearrangement did not occur\n";
            // exit(1);
        }

        numAvalanche += avalanched;

        outputStrainFunc();
    }
}
