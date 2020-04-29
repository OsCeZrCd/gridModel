#include <iostream>
#include <vector>
#include <random>
#include <cmath>

#include "GeometryVector.h"
#include "kiss_fft.h"
#include "kiss_fftnd.h"
#include "mgl2/mgl.h"

class gridModel
{
public:
    int nGridPerSide;
    double lGrid;
    int bufferCenter;

    std::vector<double> dSBuffer, alls;
    std::vector<GeometryVector> dEBuffer, alle; //e is strain
    std::vector<char> rearrangingStep;
    std::vector<bool> hasRearranged;

    std::mt19937 rEngine;
    std::normal_distribution<double> distriBution;

    gridModel(int nGrid, double lGrid) : rEngine(0), distriBution(-2.0, 2.0), nGridPerSide(nGrid), lGrid(lGrid)
    {
    }

    bool startRearranging(GeometryVector e, double s)
    {
        double yieldStrain = 0.07 - 0.01 * s;
        //return e.Modulus2() > yieldStrain*yieldStrain;
        return e.x[0] > yieldStrain;
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
        if (r < 4.0)
            return -0.1;
        else if (r < 30)
            return 1.0 / r / r / r;
        else
            return 0.0;
    }
    void getBuffer()
    {
        bufferCenter = nGridPerSide / 2;
        dEBuffer.resize(nGridPerSide * nGridPerSide);
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
            factor = 0.02 / std::fabs(outbuf[0].r);

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
        for (int i = 0; i < nSite; i++)
        {
            this->alle[i] = 0.0;
            this->alls[i] = this->distriBution(this->rEngine);
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
        }
        this->getBuffer();
    }
    void shear()
    {
        int nSite = nGridPerSide * nGridPerSide;
        for (int i = 0; i < nSite; i++)
        {
            this->alle[i].x[0] += 1e-5;
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
        }
    }
    bool avalanche()
    {
        bool avalancheHappened = false;
        int nSite = nGridPerSide * nGridPerSide;
#pragma omp parallel
        {
            int numRearrange = 1;
            while (numRearrange > 0)
            {
#pragma omp for schedule(static)
                for (int i = 0; i < nSite; i++)
                    if (startRearranging(alle[i], alls[i]))
                    {
                        hasRearranged[i] = 1;
                        rearrangingStep[i] = 1;
                    }
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
                                GeometryVector &de = dEBuffer[xInBuffer * nGridPerSide + yInBuffer];
                                e.AddFrom(de);
                                alls[x * nGridPerSide + y] += dSBuffer[xInBuffer * nGridPerSide + yInBuffer];
                            }
                        }
                        numRearrange++;
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
                            if (rearrangingStep[i] > 5)
                            {
                                rearrangingStep[i] = 0;
                                alle[i] = 0.0;
                                alls[i] = distriBution(rEngine);
                            }
                        }
                    }

                    if (numRearrange > 0)
                    {
                        avalancheHappened = true;
                        std::cout << "num rearranger in this frame=" << numRearrange << std::endl;
                    }
                }
            }
        }
        return avalancheHappened;
    }
};

template <typename T>
void plot(const std::vector<T> &data, int nGridPerSide, std::string file)
{
    mglData x(nGridPerSide, nGridPerSide);
    for (int i = 0; i < nGridPerSide; i++)
    {
        for (int j = 0; j < nGridPerSide; j++)
        {
            x.SetVal(mreal(data[i * nGridPerSide + j]), i, j);
            //std::cout<<data[i*nGridPerSide+j]<<" ";
        }
        //std::cout<<std::endl;
    }
    mglGraph gr;
    gr.SetSize(1600, 1600);
    //gr.Aspect(0.75, 1.0);
    //gr.Colorbar(">kbcyr");
    gr.Tile(x, "kbcyr");
    gr.WritePNG((file+std::string(".png")).c_str());
}
int main()
{
    const int nGridPerSide=300;
    gridModel model(nGridPerSide, 1.0);
    model.initialize();
    int numAvalanche = 0;
    while (numAvalanche < 5000)
    {
        //std::cout << "shearing\n";
        model.shear();
        //std::cout << "checking avalanche\n";
        bool avalanched = model.avalanche();
        numAvalanche += avalanched;
        if (avalanched)
        {
            std::cout << numAvalanche << "avalanches so far.\n";

            std::stringstream ss;
            ss<<"avalanche_"<<numAvalanche;
            plot(model.hasRearranged, nGridPerSide, ss.str());
        }
    }
    //for (auto &s : model.alls)
    //    std::cout << s << std::endl;
}