#include "gridModel.hpp"
// getBuffer function
// Initializes the change of strain (dE) and change of softness (dS),
// both are vectors.

void gridModel::getBuffer()
{
  // Declares variables for Fourier Transform
  int i, j, k;
  int site_periodic;
  double qx, qy, qz, qx2, qy2, qz2, q2, q4;
  double L = nGridPerSide * lGrid;


  // Initializes change of strain (dE) and change of softness (dS)
  // buffers as well as bufferCenter
  bufferCenter = nGridPerSide / 2;
  for (int i=0; i<MaxDimension; i++)
    dEBuffer[i].resize(nSite);
  dSBuffer.resize(nSite);

  // Calculates dSBuffer assuming rearranging grid point is at center
  std::vector<int> center(d);
  std::vector<int> site_coor(d);    // 2d or 3d of site coordinates
  std::vector<double> site_dx(d);   // dr (in x, y, or z), 
                                    // position difference from the grid center
  
  for (int site=0; site<nSite; site++){
    unstack(site, site_coor);
    for (int dd=0; dd<d; dd++){
      site_dx[dd] = (site_coor[dd] - bufferCenter) * lGrid;
    }
    dSBuffer[site] = dsFromRearranger(site_dx);      // generate softness change
  }

  // set the ds for the center grid equals to 0
  if (d==2)
    dSBuffer[bufferCenter * nGridPerSide + bufferCenter] = 0.0;
  else
    dSBuffer[bufferCenter * nGridPerSide * nGridPerSide + bufferCenter * nGridPerSide + bufferCenter] = 0.0;

  double meanStrainDecrement = 1.0 / nSite;
  double factor[MaxDimension];

  // strain buffer calculated by Fourier transform
  kiss_fft_cpx *inbuf = new kiss_fft_cpx[nSite];
  kiss_fft_cpx *outbuf = new kiss_fft_cpx[nSite];
 
  for (int strain=0; strain<nStrain; strain++){  
   
    for (int site=0; site<nSite; site++){
      // fill in the inbuf of FFT
      unstack(site, site_coor);

      i = site_coor[0];   // in the x direction
      i = (i > bufferCenter) ? i-nGridPerSide : i;  // PBC conditions
      qx = 2.0 * M_PI * i/L;
      qx2 = qx * qx;
      
      j = site_coor[1];   // in the y direction
      j = (j > bufferCenter) ? j-nGridPerSide : j;
      qy = 2.0 * M_PI * j/L;
      qy2 = qy * qy;
     
      if (d==2){
        q2 = qx2 + qy2;
        q4 = q2 * q2;
        if (strain==1)  // xy to xy
          inbuf[site].r = -4.0*qx2*qy2 / q4;
        else
          inbuf[site].r = -1.0*(qx2-qy2)*(qx2-qy2) / q4;
      }

      else{
        k = site_coor[2];   // in the z direction
        k = (k > bufferCenter) ? k-nGridPerSide : k;
        qz = 2.0 * M_PI * k/L;
        qz2 = qz *qz;
      
        q2 = qx2 + qy2 + qz2;
        q4 = q2 * q2;

        if (strain==0)  // xx to xx
          inbuf[site].r = (3.0*qx2*(q2-qx2) - q4) / q4;            
        else if(strain==1)  // xy to xy
          inbuf[site].r = (-4.0*qx2*qy2 - q2*qz2) / q4;
        else if(strain==2)  // xz to xz 
          inbuf[site].r = (-4.0*qx2*qz2 - q2*qy2) / q4;
        else if(strain==3)  // yy to yy
          inbuf[site].r = (3.0*qy2*(q2-qy2) - q4) / q4;
        else    // strain == 4, yz to yz
          inbuf[site].r = (-4.0*qy2*qz2 - q2*qx2) / q4; 
      }
      
      inbuf[site].i = 0.0;
    }
    
    inbuf[0].r = 0;   
    
    if (d==2){
      int temp[2] = {nGridPerSide, nGridPerSide};
      kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
      kiss_fftnd(st, inbuf, outbuf);
      free(st);
    }

    else{
      int temp[3] = {nGridPerSide, nGridPerSide, nGridPerSide};
      kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
      kiss_fftnd(st, inbuf, outbuf);
      free(st);
    }

    factor[strain] = (-1.0 + meanStrainDecrement) / outbuf[0].r;
    
    // Fills in dEBuffer
    // the interaction kernel at the rearranging site is -1.0 b/c
    // the rearranging site loses all of the redistributed strain
    
    for (int site=0; site<nSite; site++){
      unstack(site, site_coor);
      for (int dd=0; dd<d; dd++){
        i = site_coor[dd] - bufferCenter;
        i = (i<0)? i + nGridPerSide : i;
        site_coor[dd] = i;
      }
      stack(site_coor, site_periodic);
   
      dEBuffer[strain][site_periodic].x[strain] = outbuf[site].r * factor[strain] - meanStrainDecrement;
    }
  }

  // fill in the dEBuffer of insymmetry terms, like xy to xx and xx to xy
  
  for (int strain=0; strain<nStrain; strain++){
    for (int dim=0; dim<nStrain-1; dim++){  
      int buf_id = strain * (nStrain-1) + dim;
      for (int site=0; site<nSite; site++){
      
        unstack(site, site_coor);
        i = site_coor[0];   // in the x direction
        i = (i > bufferCenter) ? i-nGridPerSide : i;  // PBC conditions
        qx = 2.0 * M_PI * i/L;
        qx2 = qx * qx;
        j = site_coor[1];   // in the y direction
        j = (j > bufferCenter) ? j-nGridPerSide : j;
        qy = 2.0 * M_PI * j/L;
        qy2 = qy * qy;

        if (d==2){
          q2 = qx2 +qy2;
          q4 = q2 * q2;
          inbuf[site].r = -2.0*qx*qy*(qx2-qy2) / q4;
        }

        else {
          k = site_coor[2];   // in the z direction
          k = (k > bufferCenter) ? k-nGridPerSide : k;
          qz = 2.0 * M_PI * k/L;
          qz2 = qz *qz;
          q2 = qx2 + qy2 + qz2;
          q4 = q2 * q2;

          if (buf_id == 0)      // xx to xy
            inbuf[site].r = (-0.75*qx*qy*(q2+qx2)) / q4;
          else if (buf_id == 1) // xx to xz
            inbuf[site].r = (-0.75*qx*qz*(q2+qx2)) / q4;
          else if (buf_id == 2) // xx to yy
            inbuf[site].r = (0.5*q4 + 1.5*qx2*qy2) / q4;
          else if (buf_id == 3) // xx to yz
            inbuf[site].r = (1.5*qy*qz*qx2) / q4;

          else if (buf_id == 4) // xy to xx
            inbuf[site].r = (2.0*qx*qy*(q2-2.0*qx2)) / q4;
          else if (buf_id == 5) // xy to xz
            inbuf[site].r = (qy*qz*(q2-4.0*qx2)) / q4;
          else if (buf_id == 6) // xy to yy
            inbuf[site].r = (2.0*qx*qy*(q2-2.0*qy2)) / q4;
          else if (buf_id == 7) // xy to yz
            inbuf[site].r = (qx*qz*(q2-4.0*qy2)) / q4;

          else if (buf_id == 8)  // xz to xx
            inbuf[site].r = (2.0*qx*qz*(q2-2.0*qx2)) / q4;
          else if (buf_id == 9)  // xz to xy
            inbuf[site].r = (qy*qz*(q2-4.0*qx2))/ q4;
          else if (buf_id == 10) // xz to yy
            inbuf[site].r = (-4.0*qx*qz*qy2) / q4;
          else if (buf_id == 11) // xz to yz
            inbuf[site].r = (qx*qy*(q2-4.0*qz2)) / q4;

          else if (buf_id == 12) // yy to xx
            inbuf[site].r = (0.5*q4 + 1.5*qx2*qy2) / q4;
          else if (buf_id == 13) // yy to xy
            inbuf[site].r = (-0.75*qx*qy*(q2+qy2)) / q4;
          else if (buf_id == 14) // yy to xz
            inbuf[site].r = (1.5*qy2*qx*qz) / q4;
          else if (buf_id == 15) // yy to yz
            inbuf[site].r = (-0.75*qy*qz*(q2+qy2)) / q4;

          else if (buf_id == 16) // yz to xx
            inbuf[site].r = (-4.0*qy*qz*qx2) / q4;
          else if (buf_id == 17) // yz to xy
            inbuf[site].r = (qx*qz*(q2-4.0*qy2)) / q4;
          else if (buf_id == 18) // yz to xz
            inbuf[site].r = (qx*qy*(q2-4.0*qz2)) / q4;
          else if (buf_id == 19) // yz to yy
            inbuf[site].r = (2.0*qy*qz*(q2-2.0*qy2)) / q4;
        }

        inbuf[site].i = 0.0;
      }
      inbuf[0].r = 0;
      
      if (d==2){
        int temp[2] = {nGridPerSide, nGridPerSide};
        kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
        kiss_fftnd(st, inbuf, outbuf);
        free(st);
        for (int site=0; site<nSite; site++){
          unstack(site, site_coor);
          for (int dd=0; dd<d; dd++){
            i = site_coor[dd] - bufferCenter;
            i = (i<0)? i + nGridPerSide : i;
            site_coor[dd] = i;
          }
          stack(site_coor, site_periodic);
    
          if (strain == 0)
            dEBuffer[strain][site_periodic].x[dim+1] = outbuf[site].r * factor[strain];
        
          else
            dEBuffer[strain][site_periodic].x[dim] = outbuf[site].r * factor[strain];
        }
      }
       
      else{
        int temp[3] = {nGridPerSide, nGridPerSide, nGridPerSide};
        kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
        kiss_fftnd(st, inbuf, outbuf);
        free(st);
         
        for (int site=0; site<nSite; site++){
          unstack(site, site_coor);
          for (int dd=0; dd<d; dd++){
            i = site_coor[dd] - bufferCenter;
            i = (i<0)? i + nGridPerSide : i;
            site_coor[dd] = i;
          }
          stack(site_coor, site_periodic);
    
          if (strain == 0){
            dEBuffer[strain][site_periodic].x[dim+1] = outbuf[site].r*factor[strain];
          }
        
          else {
            if (dim < strain)
              dEBuffer[strain][site_periodic].x[dim] = outbuf[site].r*factor[strain];
            else 
              dEBuffer[strain][site_periodic].x[dim+1] = outbuf[site].r*factor[strain];
          }       
        }
      }
    }    
  }
 
  double dE_sum_1 = 0.0;   // diagonal
  double dE_sum_2 = 0.0;   // off-diagonal

  for (int na = 0; na<nStrain; na++){
    for (int nb = 0; nb<nStrain; nb++){
      for (int nc = 0; nc<nSite; nc++){
        if (nb == na)
          dE_sum_1 += dEBuffer[na][nc].x[nb];
        else
          dE_sum_2 += dEBuffer[na][nc].x[nb];
      }
    }
  }

  double dS_sum = 0.0;
  for (int na = 0; na<nSite; na++){
    dS_sum += dSBuffer[na];
  }

    
  printf ("%lf\n", dE_sum_1);   // output for chekcing whether diagonal sum is -5 (-2 for 2d),
  printf ("%lf\n", dE_sum_2);   // and off diagonal sum is 0
  printf ("%lf\n", dS_sum);     // output to check whether angular term sum is 0 (average softness change is 0)
  fflush(stdout);

  delete[] outbuf;
  delete[] inbuf;

}
