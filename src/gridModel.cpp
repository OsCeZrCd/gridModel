#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>

#include "GeometryVector.h"
#include "kiss_fft.h"
#include "kiss_fftnd.h"
#include "mgl2/mgl.h"
#include <netcdf>

template <typename T>
void plot(const std::vector<T> &data, int nGridPerSide, std::string file)
{
    mglData x(nGridPerSide, nGridPerSide);
    for (int i = 0; i < nGridPerSide; i++)
    {
        for (int j = 0; j < nGridPerSide; j++)
        {
            mreal a = 0.0;
            if (data[i * nGridPerSide + j] > 0)
                a = 1.0;
            else if (i + 1 != nGridPerSide && data[(i + 1) * nGridPerSide + j] > 0)
                a = std::max(0.5, a);
            else if (i - 1 != -1 && data[(i - 1) * nGridPerSide + j] > 0)
                a = std::max(0.5, a);
            else if (j + 1 != nGridPerSide && data[i * nGridPerSide + j + 1] > 0)
                a = std::max(0.5, a);
            else if (j - 1 != -1 && data[i * nGridPerSide + j - 1] > 0)
                a = std::max(0.5, a);
            x.SetVal(a, i, j);
        }
        //std::cout<<std::endl;
    }
    mglGraph gr;
    gr.SetSize(1600, 1600);
    //gr.Aspect(0.75, 1.0);
    gr.SetRange('z', 0.0, 1.0);
    //workaround: I wanted kw and range: [0, 1], but I can only get range=[-1, 1], SetRange does not work
    //so I added something before kw
    gr.Tile(x, "bkw");
    //gr.Colorbar(">kw");
    gr.WritePNG((file + std::string(".png")).c_str());
}

const double meanSoftness = -1.660;
const double stdSoftness = 2.0;
const double dSoftnessDStrain2 = 101.3;
//const double dSoftnessDStrain2 = 0.0;

class gridModel
{
public:
    int nGridPerSide;
    double lGrid;
    int bufferCenter;

    std::vector<double> alls;

    //e is elastic strain
    //dEBuffer[0] and [1] correspond to
    // [0] response from converting an xy elastic strain to plastic strain
    // [1] response from converting an xx deviatoric elastic strain to plastic strain
    std::vector<GeometryVector> dEBuffer[::MaxDimension], alle;

    std::vector<int> rearrangingStep;
    std::vector<char> hasRearranged;

    std::mt19937 rEngine;
    //distribution of softness after a rearrangement
    std::normal_distribution<double> sDistribution;
    //distribution of strain initially
    std::normal_distribution<double> eDistribution;
    //distribution of elastic strain after rearrangement
    std::normal_distribution<double> residualStrainDistribution;

    std::uniform_real_distribution<double> PxDistribution;
    std::vector<double> yieldStrainPx;

    netCDF::NcFile dumpFile;
    netCDF::NcVar eVar, sVar, hasRearrangedVar;
    netCDF::NcVar coeffVar;

    std::vector<double> movingAverageTarget;

    //In this version of the program, length of a rearrangement in frames is determined from the intensity
    //int rearrangeFrameLength;

    double emaMeanShift;
    const double alpha = 0.6, beta = -3.1;

    double softnessChangeShift = 0.0;

    struct neighborRelease
    {
        int x, y;
        double strainReleaseProbability;
        neighborRelease(int x, int y, double i) : x(x), y(y), strainReleaseProbability(i)
        {
        }
    };
    std::vector<neighborRelease> neighborList;

    gridModel(int nGrid, double lGrid, int seed, double emaMeanShift) : rEngine(seed),
                                                                        eDistribution(0.0, 0.01),
                                                                        residualStrainDistribution(0.0, 0.0),
                                                                        sDistribution(meanSoftness, stdSoftness),
                                                                        nGridPerSide(nGrid), lGrid(lGrid),
                                                                        movingAverageTarget(11, meanSoftness),
                                                                        emaMeanShift(emaMeanShift)
    {
        this->allocate();
        this->getBuffer();
        this->initialize();

        neighborList.push_back(neighborRelease(0, 0, 1.0));
        double sumP = 1.0;
        for (int i = -10; i < 11; i++)
            for (int j = -10; j < 11; j++)
            {
                if (i != 0 || j != 0)
                {
                    double r = std::sqrt(double(i) * i + j * j) * lGrid;
                    double p = std::exp(-1.0 * r / 1.40179);
                    neighborList.push_back(neighborRelease(i, j, p));
                    sumP += p;
                }
            }
        std::cout << "sum strain release probability=" << sumP << std::endl;

        //calculate softnessChangeShift, which is a background softness change for every block assuming the rearranging intensity is 1.0
        double sumExpectedDs = 0.0;
        for (int i = 0; i < nGridPerSide; i++)
            for (int j = 0; j < nGridPerSide; j++)
            {
                double dx = (i - bufferCenter) * lGrid;
                double dy = (j - bufferCenter) * lGrid;
                double r = std::sqrt(dx * dx + dy * dy);
                if (r < 10)
                {
                    double softnessRestoringCoefficient = alpha * ((r > 0) ? std::pow(r, beta) : 1.0);
                    sumExpectedDs += softnessRestoringCoefficient * emaMeanShift;
                }
            }
        //debug temp
        std::cout << "sumExpectedDs=" << sumExpectedDs << std::endl;
        softnessChangeShift = ((-1.0) * sumExpectedDs) / nGridPerSide / nGridPerSide;
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
        eVar.setCompression(true, true, 9);
        sVar.setCompression(true, true, 9);
        hasRearrangedVar.setCompression(true, true, 9);
        coeffVar.setCompression(true, true, 9);
    }
    void openExistingDumpFile(const std::string &filename)
    {
        dumpFile.open(filename, netCDF::NcFile::write);
        eVar = dumpFile.getVar("strain");
        sVar = dumpFile.getVar("softness");
        hasRearrangedVar = dumpFile.getVar("hasRearranged");
        coeffVar = dumpFile.getVar("yieldStrainPx");
    }
    void dump(bool writeStrain = false, bool writeSoftness = false, bool writeYieldStrainCoeff = false, bool writeHasRearranged = true)
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

        double meanYieldStrain = 0.087261 - 0.0049821 * s;
        //weibull distribution
        double k = 2.0758 - 0.024151 * s + 0.0029429 * s * s;
        double lambda = meanYieldStrain / std::tgamma(1.0 + 1.0 / k);

        double yieldStrain = std::pow(-1.0 * std::log(1.0 - yieldStrainPx[i]), 1.0 / k) * lambda;
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

    double dsFromRearranger(double dx, double dy, double r, double s, const GeometryVector &rearrangingIntensity, std::mt19937 &engine)
    {
        //if (r == 0.0)
        //    return 0.0; // delta S of the rearranger is processed separately

        //double intensityModulus = std::sqrt(rearrangingIntensity.Modulus2());
        const double intensityModulus = 0.1; //strain release per frame, measured from particle simulation

        //deprecated parameter
        double softnessChangeShift2 = 0.0;
        //double softnessChangeShift2=(-1.0)*dSoftnessDStrain2*rearrangingIntensity.Modulus2()/nGridPerSide/nGridPerSide;

        double meanContribution = 0.0;
        if (r < 30 && r > 0)
        {
            //interpolate numericalDs
            /*auto iter = std::lower_bound(numericalDsR.begin(), numericalDsR.end(), r);
            if (iter == numericalDsR.begin() || iter == numericalDsR.end())
            {
                std::cerr << "Error in interpolation! check the range of arrays\n";
                return 0.0;
            }
            iter--;
            int index = iter - numericalDsR.begin();
            double fraction = (r - *iter) / (*(iter + 1) - *iter);
            meanContribution = numericalDs[index] + fraction * (numericalDs[index + 1] - numericalDs[index]);
            */

            //contribution from volumetric strain
            meanContribution -= 18.3 * rearrangingIntensity.x[0] / r / r * std::sin(2.0 * std::atan2(dy, dx));
            meanContribution -= 18.3 * rearrangingIntensity.x[1] / r / r * std::cos(2.0 * std::atan2(dy, dx));
        }
        meanContribution += intensityModulus * softnessChangeShift + softnessChangeShift2;

        double softnessRestoringCoefficient = alpha * ((r > 0) ? std::pow(r, beta) : 1.0);
        double eta = softnessRestoringCoefficient * intensityModulus;

        double restore = 0.0;
        if (r < 10)
        {
            int index = std::floor(r + 0.5);
            restore = eta * (movingAverageTarget[index] + emaMeanShift - s);
            movingAverageTarget[index] = 0.999 * movingAverageTarget[index] + 0.001 * s;
        }

        double harmonicDiffusion = 0.0;
        if (r < 10)
        {
            double stddev = std::sqrt(eta * (2.0 - eta)) * stdSoftness;
            std::normal_distribution<double> noiseDistribution(0.0, stddev);
            harmonicDiffusion = noiseDistribution(engine);
        }

        return meanContribution + restore + harmonicDiffusion;
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
        // std::vector<double> temp(nGridPerSide * nGridPerSide, 0.0);
        // for (int i = 0; i < nGridPerSide; i++)
        // {
        //     for (int j = 0; j < nGridPerSide; j++)
        //     {
        //         int index = i * nGridPerSide + j;
        //         temp[index] = dEBuffer[0][index].x[0];
        //     }
        // }
        // plot(temp, nGridPerSide, "e00");
        // for (int i = 0; i < nGridPerSide; i++)
        // {
        //     for (int j = 0; j < nGridPerSide; j++)
        //     {
        //         int index = i * nGridPerSide + j;
        //         temp[index] = dEBuffer[0][index].x[1];
        //     }
        // }
        // plot(temp, nGridPerSide, "e01");
        // for (int i = 0; i < nGridPerSide; i++)
        // {
        //     for (int j = 0; j < nGridPerSide; j++)
        //     {
        //         int index = i * nGridPerSide + j;
        //         temp[index] = dEBuffer[1][index].x[0];
        //     }
        // }
        // plot(temp, nGridPerSide, "e10");
        // for (int i = 0; i < nGridPerSide; i++)
        // {
        //     for (int j = 0; j < nGridPerSide; j++)
        //     {
        //         int index = i * nGridPerSide + j;
        //         temp[index] = dEBuffer[1][index].x[1];
        //     }
        // }
        // plot(temp, nGridPerSide, "e11");
        // exit(0);
    }
    void allocate()
    {
        int nSite = nGridPerSide * nGridPerSide;
        alle.resize(nSite);
        alls.resize(nSite);
        yieldStrainPx.resize(nSite);
        hasRearranged.resize(nSite);
        rearrangingStep.resize(nSite);
    }
    void initialize()
    {
        int nSite = nGridPerSide * nGridPerSide;
        for (int i = 0; i < nSite; i++)
        {
            this->alle[i].x[0] = this->eDistribution(this->rEngine);
            this->alle[i].x[1] = this->eDistribution(this->rEngine);
            this->alls[i] = this->sDistribution(this->rEngine);
            this->yieldStrainPx[i] = this->PxDistribution(this->rEngine);
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
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
    //increase the first component of the strain tensor by the desired amount
    void shear(double strain = 1e-6)
    {
        int nSite = nGridPerSide * nGridPerSide;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < nSite; i++)
        {
            double olde2 = this->alle[i].Modulus2();
            this->alle[i].x[0] += strain;
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
            this->alls[i] += (this->alle[i].Modulus2() - olde2) * dSoftnessDStrain2;
        }
    }

    double rearrangingIncentive(int i)
    {
        double yieldStrain = deviatoricYieldStrain(i);
        auto e = alle[i];
        return e.Modulus2() - yieldStrain * yieldStrain;
    }

    bool avalanche(std::string outputPrefix = "")
    {
        bool avalancheHappened = false;
        int nSite = nGridPerSide * nGridPerSide;
        int nStep = 0;
        netCDF::NcFile avalancheProcessDumpFile;
        netCDF::NcVar intensityVar;
        netCDF::NcVar softnessVar;
        netCDF::NcVar rearrangeCenterVar;

        if (outputPrefix != std::string(""))
        {
            avalancheProcessDumpFile.open(outputPrefix + ".nc", netCDF::NcFile::replace);
            int dim = 2;
            std::vector<netCDF::NcDim> strainDims;

            netCDF::NcDim framedim = avalancheProcessDumpFile.addDim("frames");
            strainDims.push_back(framedim);
            for (int i = 0; i < dim; i++)
            {
                std::stringstream ss;
                ss << "dim_" << i;
                netCDF::NcDim temp = avalancheProcessDumpFile.addDim(ss.str(), nGridPerSide);
                strainDims.push_back(temp);
            }

            rearrangeCenterVar = avalancheProcessDumpFile.addVar("rearrangeCenter", netCDF::ncByte, strainDims);
            rearrangeCenterVar.setCompression(true, true, 9);

            softnessVar = avalancheProcessDumpFile.addVar("softness", netCDF::ncDouble, strainDims);
            softnessVar.setCompression(true, true, 9);

            strainDims.push_back(avalancheProcessDumpFile.addDim("strainComponents", ::MaxDimension));

            intensityVar = avalancheProcessDumpFile.addVar("rearrangingIntensity", netCDF::ncDouble, strainDims);
            intensityVar.setCompression(true, true, 9);
        }

        //if rearranging, the value is how much strain is redistributed per frame
        std::vector<GeometryVector> rearrangingIntensity;
        rearrangingIntensity.resize(nSite);
        std::vector<int> rearrangeFrameLength(nSite, 0);
        //for strain release caused by a block itself yielding, set this to 1, and softness of all other sites will be updated
        //if strain release is just because a neighbor is rearranging, leave this as 0, softness update will be disabled
        std::vector<char> updateSoftness(nSite, 0);

        std::uniform_real_distribution<double> uDistribution(0.0, 1.0);

#pragma omp parallel
        {
            std::mt19937 threadEngine;
#pragma omp critical(random)
            {
                threadEngine.seed(rEngine());
            }
            int numRearrange = 1;
            while (numRearrange > 0)
            {
#pragma omp single
                {
                    //if a site is rearranging, do nothing
                    //otherwise, find out the site with the largest incentive, let it rearrange
                    int toRearrange = -1;
                    double maxIncentive = 0.0;
                    for (int i = 0; i < nSite; i++)
                        if (rearrangingStep[i] > 0)
                        {
                            toRearrange = -1;
                            break;
                        }
                        else
                        {
                            double incentive = rearrangingIncentive(i);
                            if (incentive > maxIncentive)
                            {
                                maxIncentive = incentive;
                                toRearrange = i;
                            }
                        }
                    if (toRearrange >= 0)
                    {
                        avalancheHappened = true;

                        updateSoftness[toRearrange] = 1;
                        std::vector<int> toReleaseStrain;
                        double sumTotalIntensity2 = 0.0;

                        //we call this 'neighborList', but it includes not only neighbors but also rearranger itself
                        for (auto n : neighborList)
                        {
                            if (uDistribution(rEngine) < n.strainReleaseProbability)
                            {
                                int x = toRearrange % nGridPerSide + n.x;
                                while (x < 0)
                                    x += nGridPerSide;
                                while (x >= nGridPerSide)
                                    x -= nGridPerSide;
                                int y = toRearrange / nGridPerSide + n.y;
                                while (y < 0)
                                    y += nGridPerSide;
                                while (y >= nGridPerSide)
                                    y -= nGridPerSide;
                                int toRearrange2 = y * nGridPerSide + x;

                                rearrangingStep[toRearrange2] = 1;
                                GeometryVector residual(this->residualStrainDistribution(this->rEngine), this->residualStrainDistribution(this->rEngine));
                                GeometryVector totalIntensity = alle[toRearrange2] - residual;
                                rearrangingIntensity[toRearrange2] = totalIntensity;
                                sumTotalIntensity2 += totalIntensity.Modulus2();
                                toReleaseStrain.push_back(toRearrange2);
                            }
                        }
                        int frameLength = std::max(int(std::ceil(std::sqrt(sumTotalIntensity2) / 0.1)), 1);

                        for (auto i : toReleaseStrain)
                        {
                            rearrangingIntensity[i].MultiplyFrom(1.0 / double(frameLength));
                            rearrangeFrameLength[i] = frameLength;
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
                                double olde2 = e.Modulus2();

                                for (int j = 0; j < MaxDimension; j++)
                                {
                                    GeometryVector &de = dEBuffer[j][xInBuffer * nGridPerSide + yInBuffer];
                                    e.AddFrom(rearrangingIntensity[i].x[j] * de);
                                }

                                if (updateSoftness[i])
                                {
                                    double dx = (xInBuffer - bufferCenter) * lGrid;
                                    double dy = (yInBuffer - bufferCenter) * lGrid;
                                    double r = std::sqrt(dx * dx + dy * dy);
                                    double ds = dsFromRearranger(dx, dy, r, alls[x * nGridPerSide + y], rearrangingIntensity[i], threadEngine);
                                    alls[x * nGridPerSide + y] += ds;
                                }
                                alls[x * nGridPerSide + y] += dSoftnessDStrain2 * (e.Modulus2() - olde2);
                            }
                        }
                        numRearrange++;
                    }
                }
#pragma omp single
                {
                    if (numRearrange > 0)
                    {
                        avalancheHappened = true;
                        //std::cout << "num rearranger in this frame=" << numRearrange << std::endl;

                        if (outputPrefix != std::string(""))
                        {
                            std::stringstream ss;
                            ss << outputPrefix << "_step_" << nStep;
                            //plot(this->rearrangingStep, nGridPerSide, ss.str());
                            int framesAlreadyWritten = intensityVar.getDim(0).getSize();
                            int dim = 2;
                            std::vector<size_t> startp, countp;
                            startp.push_back(framesAlreadyWritten); //start from the end of the previous frame
                            for (int i = 0; i < dim; i++)
                                startp.push_back(0);
                            countp.push_back(1); //write one frame
                            for (int i = 0; i < dim; i++)
                                countp.push_back(nGridPerSide);
                            softnessVar.putVar(startp, countp, alls.data());
                            rearrangeCenterVar.putVar(startp, countp, (void *)(updateSoftness.data()));

                            startp.push_back(0);
                            countp.push_back(::MaxDimension);
                            intensityVar.putVar(startp, countp, rearrangingIntensity.data());
                        }
                        nStep++;
                    }

                    for (int i = 0; i < nSite; i++)
                    {
                        if (rearrangingStep[i] > 0)
                        {
                            //rearrangement has a fixed number of steps
                            rearrangingStep[i]++;
                            if (rearrangingStep[i] > rearrangeFrameLength[i])
                            {
                                rearrangeFrameLength[i] = 0;
                                rearrangingIntensity[i] = GeometryVector(0.0, 0.0);
                                rearrangingStep[i] = 0;
                                hasRearranged[i] = 1;
                                updateSoftness[i] = 0;
                                yieldStrainPx[i] = PxDistribution(rEngine);
                            }
                        }
                    }
                }
            }
        }
        if (outputPrefix != std::string(""))
        {
            std::stringstream ss;
            ss << outputPrefix;
            //plot(this->hasRearranged, nGridPerSide, ss.str());
        }
        std::cout << "steps in this avalanche=" << nStep << std::endl;
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

    int nGridPerSide;
    int seed, dumpLevel;
    double emaMeanShift;
    double maxExternalStrain;
    std::cin >> seed >> dumpLevel >> emaMeanShift >> nGridPerSide >> maxExternalStrain;
    gridModel model(nGridPerSide, 1.0, seed, emaMeanShift);
    if (fileExists(ncFileName))
    {
        model.initializeFromDumpFile(ncFileName);
        model.openExistingDumpFile(ncFileName);
    }
    else
    {
        model.openNewDumpFile(ncFileName);
        std::cout << "could not find netCDF dump file " << ncFileName << ", create a new one.\n";
    }

    int numAvalanche = 0;
    std::fstream strainFile("xyStrain.txt", std::fstream::out);
    strainFile << std::setprecision(9);
    double totalExternalStrain = 0.0;
    double strainOverStep = 1e-10;
    while (totalExternalStrain < maxExternalStrain)
    {
        double strain = model.minimumXyStrainDistanceToRarranging() + strainOverStep;
        model.shear(strain);
        totalExternalStrain += strain;

        auto outputStrainFunc = [&]() -> void {
            double sum = 0.0, sum2 = 0.0;
            for (auto s : model.alle)
            {
                sum += s.x[0];
                sum2 += s.x[0] * s.x[0];
            }
            strainFile << totalExternalStrain << ' ' << sum / model.alle.size() << ' ' << sum2 / model.alle.size() << std::endl;
        };
        outputStrainFunc();

        bool avalanched;
        if (dumpLevel >= 15)
        {
            std::stringstream ss;
            ss << "avalanche_" << numAvalanche;
            avalanched = model.avalanche(ss.str());
        }
        else
            avalanched = model.avalanche("");

        if (avalanched)
        {
            std::cout << numAvalanche << "avalanches so far.\n";
            if (dumpLevel >= 10)
                model.dump(true, true, true, true);
            else if (dumpLevel >= 5)
            {
                if (numAvalanche % 100 == 0)
                    model.dump(true, true, true, true);
                else
                    model.dump(false, false, false, true);
            }
            numAvalanche++;
            outputStrainFunc();
        }
        else
        {
            //the shear strain should be enough to trigger at least one rearrangement
            strainOverStep *= 10;
            std::cout << "Error in main : expected rearrangement did not occur, increase strainOverStep to " << strainOverStep << "\n";
        }
    }
}
