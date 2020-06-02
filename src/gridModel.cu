#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cuda_runtime.h>
#include <device_functions.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <curand.h>
#include <curand_kernel.h>
#include <string>
#include <time.h>
#include "GeometryVector.h"
#include "kiss_fft.h"
#include "kiss_fftnd.h"

void computeGridSize(int n, int blockSize, int &numBlocks, int &numThreads);
int iDivUp(int a, int b);


__global__ void energy_kernal(int nGridperside, double * E_xy, double * E_xx, double * E_buffer_xy, double * E_buffer_xx, int * pindex, double * energytotal, int num2check)
{

	int index = blockIdx.x*blockDim.x + threadIdx.x; // this is the parallel index, each thread correspond to a single grid location

	if (index >= nGridperside*nGridperside) return; // just to be safe

    int ix = index / nGridperside;
    int iy = index % nGridperside;
    int bufferCenter = nGridperside/2;

    double exx = E_xx[index];
    double exy = E_xy[index];

    double energy_old = exx * exx + exy * exy;

    for (int ip=0; ip < num2check ; ip++)   // loop through all the sites to check
    {
        int gid = pindex[ip];

    	int rx = gid / nGridperside;
        int ry = gid % nGridperside;

        int xInBuffer = bufferCenter - rx + ix;
    	while (xInBuffer < 0)
          xInBuffer += nGridperside;
    	while (xInBuffer >= nGridperside)
          xInBuffer -= nGridperside;

	    int yInBuffer = bufferCenter - ry + iy;
		while (yInBuffer < 0)
		  yInBuffer += nGridperside;
		while (yInBuffer >= nGridperside)
		  yInBuffer -= nGridperside;


    double dexx = E_buffer_xx[xInBuffer * nGridperside + yInBuffer];
    double dexy = E_buffer_xy[xInBuffer * nGridperside + yInBuffer];

    if (ry != iy || rx != ix)
    {
         double e_inc = (exx + dexx) * (exx + dexx) + (exy + dexy) * (exy + dexy) - energy_old;
//         energytotal[index] = e_inc;
         atomicAdd(&energytotal[ip], e_inc);   // a slow reduction method to ensure no core race
    }else
    {
    	 double e_dec = -energy_old;
//    	 energytotal[index] = e_dec;
         atomicAdd(&energytotal[ip], e_dec);
    }

    }
}

template <typename T>
void writebinary(const std::vector<T> &data, int nGridPerSide)
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
void writebinary_scalar(const std::vector<T> &data, int nGridPerSide)
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

void writebinary_test(double * data, int nGridPerSide)
{
  std::ofstream re_data;
  re_data.open ("data_test.bin", std::ios::out | std::ios::binary | std::fstream::app);

  double nsum = nGridPerSide*nGridPerSide;
  re_data.write((char*)&nsum,sizeof(double));

  for (int i = 0; i < nsum; i++)
  {

          double gindex=data[i];
//	      double gindex = i;
          re_data.write((char*)&gindex,sizeof(double));

  }
  std::cout<< "test written: " << nsum <<std::endl;
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
    double total_strain;

    std::vector<double> dSBuffer, alls;
    std::vector<GeometryVector> dEBuffer;
    std::vector<char> rearrangingStep;
    std::vector<bool> hasRearranged;

    std::mt19937 rEngine;
    std::normal_distribution<double> sDistribution;
    std::normal_distribution<double> eDistribution;
//    std::normal_distribution<double> deDistribution;
    std::uniform_real_distribution<double> deDistribution;
    std::normal_distribution<double> dsDistribution;


//    std::gamma_distribution<double> coeffDistribution;
    std::normal_distribution<double> coeffDistribution;
    std::vector<double> yieldStrainCoeff;
    std::vector<double> emittedStrainCoeff;

//    float *hPos;
//    float *dPos;

    double *dEbuffer_xy;
    double *dEbuffer_xx;
    double *dE_xy;
    double *dE_xx;
    double *denergy;

    double *hEbuffer_xy;
    double *hEbuffer_xx;
    double *hE_xy;
    double *hE_xx;
    double *henergy;




    gridModel(int nGrid, double lGrid) : rEngine(time(NULL)), eDistribution(0.0, 0.00001), deDistribution(0.0, 1.0), sDistribution(-0.2, 2.0), dsDistribution(0.0, 0.5), coeffDistribution(1.0, 0.14), nGridPerSide(nGrid), lGrid(lGrid)
    {
    }

    bool startRearranging(double e1, double e2, double s, int i)
    {
        double yieldStrain = 0.018 - 0.002 * s;
        if (yieldStrain < 0.014)
            yieldStrain = 0.014;
        yieldStrain*=yieldStrainCoeff[i];

//        double mod2 = e1*e1 + e2*e2;
          double mod2 = e1*e1;

//          double mod2_2 = sqrt( e1*e1);

//          if (mod2 >= yieldStrain)
//          {
//          	std::cout << "mod2:  " << mod2_2 << "   yield strain:  " << yieldStrain << "   pid:  " << i << std::endl;
//
//          }

          return sqrt(mod2) >= yieldStrain;
    }

    double checkRearranging(double e1, double e2, double s, int i)
    {
        double yieldStrain = 0.018 - 0.002 * s;
        if (yieldStrain < 0.014)
            yieldStrain = 0.014;
        yieldStrain*=yieldStrainCoeff[i];

//        double mod2 = sqrt( e1*e1 + e2*e2);

        double mod2 = sqrt( e1*e1);

        if (mod2 >= yieldStrain)
        {
        	std::cout << "mod2:  " << mod2 << "   yield strain:  " << yieldStrain << "   pid:  " << i << std::endl;
        	return 0.0;

        }

//        return yieldStrain - mod2;
        return yieldStrain - e1;
    }

    double dsFromRearranger(double dx, double dy, double r)
    {
        if (r < 1.0)
             return 0.00;
        else if (r < 30)
            return 1.0*0.52 / r / r / r - 1.0*0.16 / r / r * (std::sin(2.0 * std::atan2(dy, dx)));
//        	return  - 1.0*0.16 / r / r * (std::sin(2.0 * std::atan2(dy, dx)));
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
//            factor = 0.006 / std::fabs(outbuf[0].r);
            factor = 1.0 / std::fabs(outbuf[0].r);

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
                    dEBuffer[index2] = GeometryVector(outbuf[index].r * factor - 1.0/(double)numPixels, 0.0);
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
    }

    void initialize()
    {
        int nSite = nGridPerSide * nGridPerSide;
        total_strain = 0;
        alls.resize(nSite);
        hasRearranged.resize(nSite);
        rearrangingStep.resize(nSite);
        yieldStrainCoeff.resize(nSite);
        emittedStrainCoeff.resize(nSite);

        mean_softness = 0;
        dmean_softness = 1E-6;

        // initialize host arrays
        hEbuffer_xy = new double[nSite];
        hEbuffer_xx = new double[nSite];
        hE_xy = new double[nSite];
        hE_xx = new double[nSite];
        henergy = new double[1];

        memset(hEbuffer_xy, 0, nSite*sizeof(double));
        memset(hEbuffer_xx, 0, nSite*sizeof(double));
        memset(hE_xy, 0, nSite*sizeof(double));
        memset(hE_xx, 0, nSite*sizeof(double));
        memset(henergy, 0, sizeof(double));

        this->getBuffer();

        for (int i = 0; i < nSite; i++)
        {

            this->hE_xy[i] = this->eDistribution(this->rEngine);
            this->hE_xx[i] = this->eDistribution(this->rEngine);

            this->alls[i] = this->sDistribution(this->rEngine); //initialize softness with a distrbution s
            this->yieldStrainCoeff[i] = this->coeffDistribution(this->rEngine);
//            this->emittedStrainCoeff[i] = this->deDistribution(this->rEngine);
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;

            // this is a duplicate, probably don't need both
            this->hEbuffer_xy[i] = this->dEBuffer[i].x[0];
            this->hEbuffer_xx[i] = this->dEBuffer[i].x[1];

        }

        writebinary_test( hEbuffer_xy , nGridPerSide);
//        writebinary_test( hEbuffer_xx , nGridPerSide);


        // initialize and fill the device buffer
        size_t bytes_E = nSite*sizeof(double);

        cudaMalloc(&dEbuffer_xy, bytes_E);
        cudaMalloc(&dEbuffer_xx, bytes_E);
        cudaMalloc(&dE_xy, bytes_E);
        cudaMalloc(&dE_xx, bytes_E);
        cudaMalloc(&denergy, sizeof(double));

        cudaMemcpy(dEbuffer_xy, hEbuffer_xy, bytes_E, cudaMemcpyHostToDevice);
        cudaMemcpy(dEbuffer_xx, hEbuffer_xx, bytes_E, cudaMemcpyHostToDevice);
        cudaMemcpy(dE_xy, hE_xy, bytes_E, cudaMemcpyHostToDevice);
        cudaMemcpy(dE_xx, hE_xx, bytes_E, cudaMemcpyHostToDevice);
        cudaMemcpy(denergy, henergy, sizeof(double), cudaMemcpyHostToDevice);
    }


    void shear()
    {
        double sum = 0.0;
        for (auto &s : this->alls)
            sum += s;
        double mean_s = sum / alls.size();
        int nSite = nGridPerSide * nGridPerSide;
        double stepsize = 1.0;
        for (int i = 0; i < nSite; i++)
        {
        	double step_i = checkRearranging(hE_xy[i],hE_xx[i], alls[i],i);
        	stepsize = min(step_i,stepsize);
        }

        stepsize = max(0.0,stepsize);

        std::cout << "step size is: " << stepsize << std::endl;

        total_strain += stepsize;

#pragma omp parallel for schedule(static)
        for (int i = 0; i < nSite; i++)
        {
//            double E_increament = 1.0 + (this->alls[i] - mean_s) * 0.0;    // assign additional strain to soft particles
//            E_increament = std::max(0.0,E_increament);
        	double E_increament = 1.0 + 1E-10;
            this->hE_xy[i] += stepsize * E_increament;
            this->hE_xx[i] += 0.0e-6;

//            this->alls[i] += dmean_softness;
            this->hasRearranged[i] = 0;
            this->rearrangingStep[i] = 0;
        }

//        mean_softness = mean_softness + dmean_softness;
    }


    int avalanche(std::string outputPrefix = "")
    {
        int nSite = nGridPerSide * nGridPerSide;
//        double deltaEnergy;
        int numRearrange = 1;

#pragma omp parallel
        {
            while (numRearrange > 0)
            {
#pragma omp for schedule(static)
                for (int i = 0; i < nSite; i++)
                    if (startRearranging(hE_xy[i],hE_xx[i], alls[i],i))
                    {
                        rearrangingStep[i] = 1;
                        emittedStrainCoeff[i] = eDistribution(rEngine);
                    }
                    else
                    {
                    	rearrangingStep[i] = 0;
                    }
#pragma omp barrier
#pragma omp single
          {
////        	  // initiate a thrust vector to collect energy change for all potential site
//              thrust::host_vector<double> h_energy(nSite);
//              thrust::fill(h_energy.begin(),h_energy.end(), 0.0);
//              thrust::device_vector<double> d_energy = h_energy;
//              double *d_energy_pt = thrust::raw_pointer_cast( &d_energy[0] );
//
//              // initiate a thrust vector to collect the index of all potential site
//              thrust::host_vector<int> h_pindex(nSite);
//              thrust::fill(h_pindex.begin(),h_pindex.end(), 0);
//              thrust::device_vector<int> d_pindex = h_pindex;
//              int *d_pindex_pt = thrust::raw_pointer_cast( &d_pindex[0] );
//
//              // collect the index of all potential sites
//              int site_2_check=0;
//              for (int i = 0; i < nSite; i++)
//              {
//            	  if (rearrangingStep[i] > 0){
//            		  d_pindex[site_2_check] = i;
//            		  site_2_check+=1;
//            	  }
//
//              }
//              std::cout<< "total potential site checked: " << site_2_check << std::endl;
//
//              size_t bytes_E = nSite*sizeof(double);
//              cudaMemcpy(dE_xy, hE_xy, bytes_E, cudaMemcpyHostToDevice);
//              cudaMemcpy(dE_xx, hE_xx, bytes_E, cudaMemcpyHostToDevice);
//
//              //stop rearrangements that increases energy
//
//              // the idea is to assign each grid position to a gpu-thread, within each grid position, we loop and calculate energy change for all potential sites
//              int numBlocks, numThreads;
//              int nSite = nGridPerSide * nGridPerSide;
//              computeGridSize(nSite, 256, numBlocks, numThreads);
//
//              energy_kernal <<< numBlocks, numThreads >>> (nGridPerSide, dE_xy, dE_xx, dEbuffer_xy, dEbuffer_xx, d_pindex_pt, d_energy_pt, site_2_check);
//
//              cudaDeviceSynchronize();
//
//              for (int i = 0; i < site_2_check; i++)
//              {
////            	  std::cout << "energy change: " << d_energy[i] << std::endl;
//            	  if (d_energy[i]>0){
//            		  int index_grid = d_pindex[i];
//            		  rearrangingStep[index_grid] = 0;
//            	  }
//
//              }
            }


#pragma omp barrier
                //rearrangement affect other sites parameters
//                numRearrange = 0;
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
                        // if (xInBuffer>=0 && xInBuffer<nGridPerSide)
                        // {
                        for (int y = 0; y < nGridPerSide; y++)
                        {
                            int yInBuffer = bufferCenter - ry + y;
                             while (yInBuffer < 0)
                                 yInBuffer += nGridPerSide;
                             while (yInBuffer >= nGridPerSide)
                                 yInBuffer -= nGridPerSide;
//                            if (yInBuffer>=0 && yInBuffer<nGridPerSide)
//                            {

//                            hE_xx[x * nGridPerSide + y] +=  hEbuffer_xx[xInBuffer * nGridPerSide + yInBuffer] * emittedStrainCoeff[i];
                            hE_xy[x * nGridPerSide + y] +=  hEbuffer_xy[xInBuffer * nGridPerSide + yInBuffer] * ( hE_xy[i] - emittedStrainCoeff[i]) ;

                            alls[x * nGridPerSide + y] += dSBuffer[xInBuffer * nGridPerSide + yInBuffer];
//                            }
                        }
                        // }
                    }
//                    numRearrange++;
                    }
                }

                // copy new device strain array back to the geometry vector vector

#pragma omp single
                {
                	numRearrange=0; // this is temporary----------------------------------------------------

                    for (int i = 0; i < nSite; i++)
                    {
                        if (rearrangingStep[i] > 0)
                        {
                        	total_strain += (hE_xy[i] - emittedStrainCoeff[i]) / ( (double)nGridPerSide * (double)nGridPerSide  );
                            //carry out the rearrangement
                            rearrangingStep[i]++;
                            hasRearranged[i] = 1;
                            hE_xy[i] = emittedStrainCoeff[i];
//                            hE_xx[i] = eDistribution(rEngine);

                            alls[i] = dsDistribution(rEngine) + 0.25 - 0.4 * alls[i];
                            yieldStrainCoeff[i] = coeffDistribution(rEngine);
//                            emittedStrainCoeff[i] = deDistribution(rEngine);

                            numRearrange++;
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

//      std::cout << "got here" << std::endl;
        return numRe_frame;
    }
};

int main()
{
    const int nGridPerSide = 512;
    gridModel model(nGridPerSide, 1.0);
//    std::cout<<"before initialization"<<std::endl;
    model.initialize();
    int numAvalanche = 0;
    int strainstep = 0;
    double meanS = 0.0;

    std::cout<<"initialized"<<std::endl;
//    model.shear_cuda();
    int numRe = 0;
//
    while (strainstep<200000 && meanS<2.0)
    {


        double sum = 0.0;
        for (auto &s : model.alls)
           sum += s;

        meanS = sum / model.alls.size();
        std::ofstream re_mean_data;
        re_mean_data.open ("data_mean.bin", std::ios::out | std::ios::binary | std::fstream::app);
        double nstep = strainstep;
        re_mean_data.write((char*)&nstep,sizeof(double));
        re_mean_data.write((char*)&meanS,sizeof(double));
        double etotal=model.total_strain;
        double stotal = 0;
        for (int i = 0; i < nGridPerSide*nGridPerSide; i++)
         {
          stotal +=  model.hE_xy[i] ;
         }
        re_mean_data.write((char*)&stotal,sizeof(double));
        re_mean_data.write((char*)&etotal,sizeof(double));
        re_mean_data.close();
        std::cout << "Currently at step: " << strainstep << ",  total strain is: "<<   model.total_strain  <<", Number of rearrangement:"<<  numRe;
        std::cout << ", mean s=" << sum / model.alls.size() << ", total energy=" << etotal << std::endl;


        model.shear();

        std::stringstream ss;
        ss << "avalanche_" << numAvalanche;

        numRe = model.avalanche(ss.str());
//        int numRe = 0;


        writebinary<bool>(model.hasRearranged, nGridPerSide);

        strainstep+=1;
//        if ((strainstep-1)%500==0)
//        {
//            writebinary_scalar<double>(model.alls, nGridPerSide);
//        }



    }
}



//leave them alone, these are for determine how to assign threads
void computeGridSize(int n, int blockSize, int &numBlocks, int &numThreads)
{
    numThreads = std::min(blockSize, n);
    numBlocks = iDivUp(n, numThreads);
}

int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}
