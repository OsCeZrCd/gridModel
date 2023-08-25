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
  double Lx = nxGrid * lGrid;
  double Ly = nyGrid * lGrid;
  double Lz = nzGrid * lGrid;


  // Initializes change of strain (dE) and change of softness (dS)
  // buffers as well as bufferCenter
  bufferCenterX = (nxGrid-1) / 2;
  bufferCenterY = (nyGrid-1) / 2;
  bufferCenterZ = (nzGrid-1) / 2;


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
    site_dx[0] = (site_coor[0] - bufferCenterX) * lGrid;
    site_dx[1] = (site_coor[1] - bufferCenterY) * lGrid;
    site_dx[2] = (site_coor[2] - bufferCenterZ) * lGrid;
    dSBuffer[site] = dsFromRearranger(site_dx);      // generate softness change
  }

  // set the ds for the center grid equals to 0
  if (d==2)
    dSBuffer[bufferCenterY * nxGrid + bufferCenterX] = 0.0;
  else
    dSBuffer[bufferCenterZ * nxGrid * nyGrid + bufferCenterY * nxGrid + bufferCenterX] = 0.0;

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
      i = (i > bufferCenterX) ? i-nxGrid : i;  // PBC conditions
      qx = 2.0 * M_PI * i/Lx;
      qx2 = qx * qx;
      
      j = site_coor[1];   // in the y direction
      j = (j > bufferCenterY) ? j-nyGrid : j;
      qy = 2.0 * M_PI * j/Ly;
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
        k = (k > bufferCenterZ) ? k-nzGrid : k;
        qz = 2.0 * M_PI * k/Lz;
        qz2 = qz *qz;
      
        q2 = qx2 + qy2 + qz2;
        q4 = q2 * q2;

        if (strain==0)  // xx to xx
          inbuf[site].r = 2.0*((q2-qx2)*qx2 + qx2*qz2)/q4 - 1;            
        else if(strain==1)  // xy to xy
          inbuf[site].r = (-4.0*qx2*qy2 - q2*qz2) / q4;
        else if(strain==2)  // xz to xz 
          inbuf[site].r = (-4.0*qx2*qz2 - q2*qy2) / q4;
        else if(strain==3)  // yy to yy
          inbuf[site].r = 2.0*((q2-qy2)*qy2 + qy2*qz2)/q4 - 1;
        else    // strain == 4, yz to yz
          inbuf[site].r = (-4.0*qy2*qz2 - q2*qx2) / q4; 
      }
      
      inbuf[site].i = 0.0;
    }
    
    inbuf[0].r = 0;   
    
    if (d==2){
      int temp[2] = {nxGrid, nyGrid};
      kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
      kiss_fftnd(st, inbuf, outbuf);
      free(st);
    }

    else{
      int temp[3] = {nxGrid, nyGrid, nzGrid};
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
      
      i = site_coor[0] + bufferCenterX;
      site_coor[0] = (i>=nxGrid)? i - nxGrid : i;
      i = site_coor[1] + bufferCenterY;
      site_coor[1] = (i>=nyGrid)? i - nyGrid : i;
      i = site_coor[2] + bufferCenterZ;
      site_coor[2] = (i>=nzGrid)? i - nzGrid : i;

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
        i = (i > bufferCenterX) ? i-nxGrid : i;  // PBC conditions
        qx = 2.0 * M_PI * i/Lx;
        qx2 = qx * qx;
        j = site_coor[1];   // in the y direction
        j = (j > bufferCenterY) ? j-nyGrid : j;
        qy = 2.0 * M_PI * j/Ly;
        qy2 = qy * qy;

        if (d==2){
          q2 = qx2 +qy2;
          q4 = q2 * q2;
          inbuf[site].r = -2.0*qx*qy*(qx2-qy2) / q4;
        }

        else {
          k = site_coor[2];   // in the z direction
          k = (k > bufferCenterZ) ? k-nzGrid : k;
          qz = 2.0 * M_PI * k/Lz;
          qz2 = qz *qz;
          q2 = qx2 + qy2 + qz2;
          q4 = q2 * q2;

          if (buf_id == 0)      // xx to xy
            inbuf[site].r = qx*qy*(q2 - 2*qx2 + 2*qz2) / q4;
          else if (buf_id == 1) // xx to xz
            inbuf[site].r = 2*qx*qz*(qz2-qx2) / q4;
          else if (buf_id == 2) // xx to yy
            inbuf[site].r = 2*qy2*(qz2-qx2) / q4;
          else if (buf_id == 3) // xx to yz
            inbuf[site].r = qy*qz*(-q2 - 2*qx2 + 2*qz2) / q4;

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
            inbuf[site].r = 2.0*qx2*(qz2-qy2) / q4;
          else if (buf_id == 13) // yy to xy
            inbuf[site].r = qx*qy*(q2 - 2*qy2 + 2*qz2) / q4;
          else if (buf_id == 14) // yy to xz
            inbuf[site].r = qx*qz*(-q2 - 2*qy2 + 2*qz2) / q4;
          else if (buf_id == 15) // yy to yz
            inbuf[site].r = 2*qy*qz*(qz2 - qy2) / q4;

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
        int temp[2] = {nxGrid, nyGrid};
        kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
        kiss_fftnd(st, inbuf, outbuf);
        free(st);
      }
       
      else{
        int temp[3] = {nxGrid, nyGrid, nzGrid};
        kiss_fftnd_cfg st = kiss_fftnd_alloc(temp, d, true, nullptr, nullptr);
        printf ("In: %d %lf\n", strain, inbuf[0].r);
        kiss_fftnd(st, inbuf, outbuf);
        printf ("Out: %d %lf\n", strain, outbuf[0].r);
        free(st);
         
        for (int site=0; site<nSite; site++){
          unstack(site, site_coor);

          i = site_coor[0] + bufferCenterX;
          site_coor[0] = (i>=nxGrid)? i - nxGrid : i;
          i = site_coor[1] + bufferCenterY;
          site_coor[1] = (i>=nyGrid)? i - nyGrid : i;
          i = site_coor[2] + bufferCenterZ;
          site_coor[2] = (i>=nzGrid)? i - nzGrid : i;

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
 
  if (d==3)
    printf ("%lf %lf %lf %lf %lf\n", factor[0], factor[1], factor[2], factor[3], factor[4]);

  double dE_sum_1 = 0.0;   // diagonal
  double dE_sum_2 = 0.0;   // off-diagonal

  for (int na = 0; na<nStrain; na++){
    for (int nb = 0; nb<nStrain; nb++){
      for (int nc = 0; nc<nSite; nc++){
        if (nb == na && nb == 1)
          dE_sum_1 += dEBuffer[na][nc].x[nb];
        else{
          //dEBuffer[na][nc].x[nb] = 0;
          dE_sum_2 += dEBuffer[na][nc].x[nb];
        }
      }
    }
  }

  double dS_sum = 0.0;
  for (int na = 0; na<nSite; na++){
    dS_sum += dSBuffer[na];
  }

    
  printf ("%lf\n", dE_sum_1);   // output for chekcing whether diagonal sum is -5 (-2 for 2d),
  printf ("%lf\n", dE_sum_2);   // and off diagonal sum is 0
  printf ("%lf\n", dS_sum);
  fflush(stdout);

  if (d==2){
    unstack(0, center);
    center[0] = bufferCenterX;
    center[1] = bufferCenterY;
    int center_id;
    stack(center, center_id);
    for (int na=0; na<nStrain; na++){
      for (int nb=0; nb<nStrain; nb++){
        printf ("%lf ", dEBuffer[na][center_id].x[nb]);
        fflush(stdout);
      }
      printf ("\n");
      fflush(stdout);
    }
  }

  else{
    unstack(0, center);
    center[0] = bufferCenterX;
    center[1] = bufferCenterY;
    center[2] = bufferCenterZ;
    int center_id;
    stack(center, center_id);
    for (int na=0; na<nStrain; na++){
      for (int nb=0; nb<nStrain; nb++){
        printf ("%lf ", dEBuffer[na][center_id].x[nb]);
        fflush(stdout);
      }
      printf ("\n");
      fflush(stdout);
    }
  }

  delete[] outbuf;
  delete[] inbuf;

}
