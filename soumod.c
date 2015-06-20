// soumod.c reads slip distribution with a certain discretisation step from given file sourceold, 
// resamples it to the given discretisation step disc and adjust the values that way 
// that a k²-distribution of slip is obtained. Hence, the low frequency content of the original
// slip distribution is kept and the high frequency content is varied randomly.
//
// if no slip distribution is given, the programm creates a random slip distribution with the parameters given
//
// The focal mechanisms of the different fault patches are varied. Variation of strike and dip are chosen, so that 
// the fault plane assumes the shape of a fractal surface with aspect ratio 0.1 (in respect to the shorter fault edge)
// Rake is varied by +-10 degree
//
//////////////////////////////////////////////// The rupture velocity is varied.
//
// Katrin Kieling, august 2010 
//
// last modified: Jun 25 2013
// 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "RECIPES/nrutil.c"
#include "RECIPES/ran1.c"
#include "RECIPES/dft2d.c"
#include "readidem3.c"
#include "readidem2.c"
#include "readidem.c"
#include "rupfront.c"
//  #include "readxyplane.c"
int
soumod (char *sourceold, char *modtype, double disc, double Me, double rupvel,
		char *taperswitch, int seed, double **density, int num_model_depth, double lathypo, double lonhypo, double dephypo,
		double seglat, double seglon, double segdepth, double *dip1, double moment)
{
  FILE *fin3, *fout, *fout2;
  char outname[] = "rupture.out";
  char *kfile;
  char data[256];
  double reflat, reflon, refdepth, strike, rake;
  double odiscx, odiscy, maxox=0, maxoy=0;
  double dk_x, dk_y, kc, kcy, kcy_int; 
  double *lon, *lat, *z, *ox, *oy, *REslipo, *rake_v, *strike_v, *dip, *lato, *lono, *depo;
  double *nx, *ny, *REslip, *IMslip, *k_x, *k_y, *REold, *IMold, *FToldRE,
    *FToldIM, *dipold, *rakeold, *strikeold, *latold, *lonold, *depthold;
  double   *REslipnew, *IMslipnew;
  double *REsurf, *IMsurf, *FTsurfRE, *FTsurfIM, *dstrike, *ddip, *drake,
    *FTrand2RE, *FTrand2IM;
  double height, maxsurf, minsurf, maxoslip=0, maxslip;
  double   *FTrandRE, *FTrandIM, *specREold, *specIMold, *specREnew,*specIMnew , *taper, *taper2;
//double *vel;
  double sma, smb, smc, bgc, bga, alphax, alphay;
  double refdepth1, ydis;
  double minslip, diffspec, velmean;
  double slipsum=0.0;
  double oxdim = 0.0;
  double oydim = 0.0;
  float aspect = 0.02;
  float roughness = 2.0;
  float roughslip =2.3;
//  float fraction = 0.85;	//fraction MUST be < 1, for fraction>1 rupture velocity may become negative
  int old = 0;
  int smooth = 0;
  int i, j, k, l;
  int ns, nns, nsx, nsy;
  long seed2;
  double pi = 4.0 * atan (1.0);
  float REARTH = 6.371 * pow (10, 3);	// radius of earth in km

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (strcmp (sourceold, "new") != 0)
	  { seed2 = -fabs (seed);
        old=1;
        if (strcmp (modtype, "idem") == 0)
	      {printf ("idem\n");
	       readidem (sourceold, &ns, &reflat, &reflon, &refdepth,
	    		     &rake, &strike, &dip, &oxdim, &oydim, &nsx, &nsy,
	    		     &odiscx, &odiscy, &ox, &oy, &REslipo,
	    		     &rake_v, &strike_v, &lato, &lono, &depo);
	      }
        else if (strcmp (modtype, "idem2") == 0)
	      { printf ("idem2\n");
	        readidem2 (sourceold, &ns, &reflat, &reflon, &refdepth, &rake,
		             &strike, &dip, &oxdim, &oydim, &nsx, &nsy, &odiscx,
		             &odiscy, &ox, &oy, &REslipo, &rake_v, &strike_v,&lato, &lono, &depo);
	      }
	    else if (strcmp (modtype, "idem3") == 0)
	      { printf ("idem3\n");
	        readidem3 (sourceold, &ns, &reflat, &reflon, &refdepth, &rake,
	        		 &strike, &dip, &oxdim, &oydim, &nsx, &nsy, &odiscx,
	        		 &odiscy, &ox, &oy, &REslipo, &rake_v, &strike_v,&lato, &lono, &depo);
	      }

//      else if (strcmp (modtype, "xyplane") == 0)
//	{
//	  printf ("xyplane\n");
//	  readxyplane (sourceold, data, &ns, &reflat, &reflon, &refdepth,
//		       &rake, &strike, &dip, &oxdim, &oydim, &nsx, &nsy,
//		       &odiscx, &odiscy, &ox, &oy, &REslipo, &rake_v, &lato, &lono, &depo);
//	}
        else
	      { printf ("unknown slip distribution format");
	      }
        printf ("source read\n");
        *dip1=0;
        for (i = 1; i <= ns; i++)
	      { *dip1=*dip1+dip[i];
//	        printf("REslipo=%lf, ox=%lf, oy=%lf\n", REslipo[i], ox[i], oy[i]);
	        if (REslipo[i] > maxoslip)
	          { maxoslip = REslipo[i];
	          }
	        if (ox[i] > maxox)
	          { maxox = ox[i];
	          }
	        if (oy[i] > maxoy)
	          { maxoy = oy[i];
	          }
	      }
        *dip1=*dip1/ns;
        printf(" dip1=%lf\n",*dip1);
 //       maxoy = maxoy + odiscy;
        printf ("disc=%lf, maxox=%lf, maxoy=%lf\n", disc, maxox, maxoy);
//resampling of slipmap
        nns = (int) (ceil (oxdim / disc) * ceil (oydim / disc));	// number of patches
//      printf (" nns+1=%d\n", nnsdd+1);
        nsx = (int) ceil (oxdim / disc);
        nsy = (int) ceil (oydim / disc);
        REold = dvector (1, nns+1);
        IMold = dvector (1, nns+1);
        latold= dvector (1,nns+1);
        lonold= dvector (1,nns+1);
        depthold= dvector (1,nns+1);
        dipold = dvector (1, nns+1);
        rakeold = dvector (1, nns+1);
        strikeold = dvector (1, nns+1);
        nx = dvector (1, nns+1);
        ny = dvector (1, nns+1);

        i = 1;
        for (k = 1; k <= nsy; k++)
	      { for (j = 1; j <= nsx; j++)
	          { nx[i] = (j - 1) * disc;	//relative position along strike
	            ny[i] = (k - 1) * disc;	//relative position along dip
	            i++;
	          }
	      }
        printf ("disc=%lf, maxox=%lf, maxoy=%lf\n", disc, maxox, maxoy);
        if (disc < odiscx && disc < odiscy)	// refine slip distribution
          { i = 1;
            l = 1;
	        while (l <= ns && i <= nns)
	          {	//printf("nx[i]=%lf, ox[l]-0.1*odiscx=%lf, ny[i]=%lf, oy[l]-0.1*odiscy=%lf\n",nx[i],ox[l]-0.1*odiscx,ny[i],oy[l]-0.1*odiscy);
	            if (ox[l] > maxox - 0.2 * odiscx)
		          {		//printf("AHA! maxox=ox\n");
		            if (oy[l] > maxoy - 0.3 * odiscy)
		            { REold[i] = REslipo[l];
		              latold[i] = lato[l];
		              lonold[i] = lono[l];
		              depthold[i] = depo[l];
		              IMold[i] = 0;
		              dipold[i] = dip[l];
		              rakeold[i] = rake_v[l];
		              strikeold[i] = strike_v[l];
		              i++;
		              l = 1;
		       //       printf("strike=%lf, REold=%lf, nx=%lf, ny=%lf\n", strikeold[i-1], REold[i-1], nx[i-1], ny[i-1]);
		            }
		            else if (ny[i] < oy[l] + odiscy && ny[i] >= oy[l] - 0.2 * odiscy)
		            { REold[i] = REslipo[l];
		              latold[i] = lato[l];
		              lonold[i] = lono[l];
		              depthold[i] = depo[l];
		              IMold[i] = 0;
		              dipold[i] = dip[l];
		              rakeold[i] = rake_v[l];
		              strikeold[i] = strike_v[l];
		              i++;
		        	  l = 1;
		    //    	  printf("strike=%lf, REold=%lf, nx=%lf, ny=%lf\n", strikeold[i-1], REold[i-1], nx[i-1], ny[i-1]);
		            }
		            else
		            { l++;
		            }
		    }
	      else if (oy[l] > maxoy - 0.1 * odiscy && nx[i] < ox[l] + odiscx && nx[i] >= ox[l] - 0.1 * odiscx)
		    { REold[i] = REslipo[l];
	          latold[i] = lato[l];
	          lonold[i] = lono[l];
	          depthold[i] = depo[l];
	          IMold[i] = 0;
	          dipold[i] = dip[l];
	          rakeold[i] = rake_v[l];
	          strikeold[i] = strike_v[l];
	          i++;
	          l = 1;
	     //     printf("strike=%lf, REold=%lf, nx=%lf, ny=%lf\n", strikeold[i-1], REold[i-1], nx[i-1], ny[i-1]);
		    }
	      else if (nx[i] < ox[l] + odiscx && ny[i] < oy[l] + odiscy && nx[i] >= ox[l] - 0.1 * odiscx && ny[i] >= oy[l] - 0.1 * odiscy)
		    { REold[i] = REslipo[l];
	          latold[i] = lato[l];
	          lonold[i] = lono[l];
	          depthold[i] = depo[l];
	          IMold[i] = 0;
	          dipold[i] = dip[l];
	          rakeold[i] = rake_v[l];
	          strikeold[i] = strike_v[l];
	          i++;
	          l = 1;
	        //  printf("strike=%lf, REold=%lf, nx=%lf, ny=%lf\n", strikeold[i-1], REold[i-1], nx[i-1], ny[i-1]);

		    }
	      else
		    {		//printf("nx=%lf, ox=%lf, ox+odisc=%lf, ny=%lf, oy=%lf, oy+odisc=%lf\n", nx[i], ox[l], ox[l]+odiscx, ny[i], oy[l], oy[l]+odiscy);
	          if (l == ns - 1)
	        	{ REold[i] = REold[i - 1];
	        	  latold[i] =latold[i-1];
	        	  lonold[i] = lonold[i-1];
	        	  depthold[i]= depthold[i-1];
	        	  IMold[i] = 0;
	        	  dipold[i] = dipold[i - 1];
	        	  rakeold[i] = rakeold[i - 1];
	        	  strikeold[i] = strikeold[i - 1];
	        	  i++;
	        	  l = 0;
	        	//  printf("strike=%lf, REold=%lf, nx=%lf, ny=%lf\n", strikeold[i-1], REold[i-1], nx[i-1], ny[i-1]);
	        	}
	          l++;
		    }

	    }
	 }
      else if (disc > odiscx && disc > odiscy)	// smooth slip distribution
	{ printf("smooth=1: Slipverteilung glätten\n");
	  smooth = 1;
	  i = 1;
	  l = 1;
	  for (i = 1; i <= nns; i++)
	    {
	      REold[i] = 0;
	      IMold[i] = 0;
	      k = 0;
	      while (l < ns)
		{
		  if (ox[l] < nx[i] + disc && oy[l] < ny[i] + disc
		      && ox[l] >= nx[i] && oy[l] >= ny[i])
		    {
		      REold[i] = REold[i] + REslipo[l];
		      latold[i]=latold[i] + lato[l];
		      lonold[i]=lonold[i] + lono[l];
		      depthold[i]=depthold[i] + depo[l];
		      dipold[i] = dipold[i] + dip[l];
		      rakeold[i] = rakeold[i] + rake_v[l];
		      strikeold[i] = strikeold[i] + strike_v[l];
		      k++;
		    }
		  l++;
		}
	      if (REold[i] > 0)
		{
		  REold[i] = REold[i] / k;
		  latold[i]=latold[i]/k;
		  lonold[i]=lonold[i]/k;
		  depthold[i]=depthold[i]/k;
		  dipold[i] = dipold[i] / k;
		  rakeold[i] = rakeold[i] / k;
		  strikeold[i] = strikeold[i] / k;
		}
	      l = 1;
	    }
	  //// calculated position for each fault patch
	  lat = dvector (1, nns+1);
	  lon = dvector (1, nns+1);
	  z = dvector (1, nns+1);
	  i = 1;
	  l = 1;
	  refdepth1 = refdepth;
	  ydis = 0.0;
	  for (k = 1; k <= nsx; k++)
	    {
	      for (j = 1; j <= nsy; j++)
	    	{
		      if (i > 1 && ny[i] > ny[i - 1])
		        { refdepth1 = z[i - 1];
		          printf ("refdepth=%f,dipold=%f,i=%d,nns=%d\n",refdepth1, dipold[i], i, nns+1);
		          ydis = ydis + disc * cos (dipold[i - 1] * pi / 180);
		        }
		      z[i] = refdepth1 + disc * sin (dipold[i] * pi / 180);	// absolute depth
		      sma = sqrt (pow (nx[i], 2) + pow (ydis + disc * cos (dipold[i] * pi / 180), 2)) / REARTH;
		      smb = 0.5 * pi - reflat * pi / 180;
		      bgc =strikeold[i] * pi / 180 +atan ((ydis + disc * cos (dipold[i] * pi / 180)) / nx[i]);
		      smc =acos (cos (sma) * cos (smb) + sin (sma) * sin (smb) * cos (bgc));
		      bga = asin (sin (sma) * sin (bgc) / sin (smc));
		      lat[i] = 90 - smc * 180 / pi;	//absolute latitude
		      lon[i] = fmod (reflon + bga * 180 / pi, 360);	//absolute longitude
		      if (nx[i] == 0 && ny[i] == 0)	// set outer border to zero
		        { lat[i] = reflat;
		          lon[i] = reflon;
		        }
		      i++;
	    	}
	    }
	  lat[1] = reflat;
	  lon[1] = reflon;
	}
      else
	   {
	     printf("choose discretisation either higher or lower to make clear between smoothing or refining! \n");
	     exit(0);
	   }
      printf ("smooth=%d \n", smooth);
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  else if (ns == 11 || ns == 12)	// case 2) create new random slip distribution
    {
      old = 0;
      fin3 = fopen (sourceold, "r");
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &reflat);	// reference latitude
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &reflon);	// reference longitude 
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &refdepth);	// refenrence depth in km
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &oxdim);	// length in km
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &oydim);	// width in km
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &strike);	// strike in degrees
      fgets (data, 256, fin3);
      sscanf (data, "%lf", dip1);	// dip in degrees
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &rake);	// rake in degrees
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &disc);	// discretisation step in km
      fgets (data, 256, fin3);
      sscanf (data, "%ld", &seed2);	// seed = start number of random walk
      fgets (data, 256, fin3);
      sscanf (data, "%lf", &velmean);	// mean rupture velocity
      fgets (data, 256, fin3);
      if (ns == 12)
	{
	  fgets (data, 256, fin3);
	  sscanf (data, "%lf", &Me);	// energy magnitude
	}
      else
	{
	  Me = 4.38 + 1.49 * log10 (oxdim);
	}
      fclose(fin3);
      nsx = ceil (oxdim / disc);
      nsy = ceil (oydim / disc);
      nns = (int) (ceil (oxdim / disc) * ceil (oydim / disc));	// number of patches
      ns = nsx * nsy;
      i = 1;
      l = 1;
      dipold = dvector (1, nns+1);
      rakeold = dvector (1, nns+1);
      strikeold = dvector (1, nns+1);
      while (i <= nns)
	{
	  dipold[i] = *dip1;
	  rakeold[i] = rake;
	  strikeold[i] = strike;
	  i++;
	}
      seed2 = -fabs (seed2);
    }
  else
    {
      printf ("number of input parameters too small or too large!\n");
      printf
	("a) For newslip distribution give reflat, refon, refdepth, length, width, strike, dip, surf, discretisation, the mean rupture velocity, and a random seed number. You may add an explicit value for the energy magnitude \n");
      printf
	("b) For resampling a given slip distribution give filename of old distribution, the modeltype (ff=USGS finite fault model, srcmod=models provided on www.seismo.ethz.ch/srcmod/Homepage.html, idem for input file in the same format as output file), the new discretisation step, the energy magnitude, and the mean rupture velocity\n");
    }
  if (smooth == 0)
    {
      nx = dvector (1, nns+1);
      ny = dvector (1, nns+1);
      i = 1;
      for (k = 1; k <= nsy; k++)
	{
	  for (j = 1; j <= nsx; j++)
	    {
	      nx[i] = (j - 1) * disc;	//relative postiion along strike
	      ny[i] = (k - 1) * disc;	//relative postiion along dip
	      i++;
	    }
	}
      //// calculated position for each fault patch
      lat = dvector (1, nns+1);
      lon = dvector (1, nns+1);
      z = dvector (1, nns+1);
      i = 1;
      l = 1;
      refdepth1 = refdepth;
      ydis = 0.0;
      for (k = 1; k <= nsx; k++)
	{
	  for (j = 1; j <= nsy; j++)
	    {
	      if (i > 1 && ny[i] > ny[i - 1])
		{
		  refdepth1 = z[i - 1];
		  ydis = ydis + disc * cos (dipold[i - 1] * pi / 180);
		}
	      z[i] = refdepth1 + disc * sin (dipold[i] * pi / 180);	// absolute depth
	      sma =
		sqrt (pow (nx[i], 2) +
		      pow (ydis + disc * cos (dipold[i] * pi / 180),
			   2)) / REARTH;
	      smb = 0.5 * pi - reflat * pi / 180;
	      bgc =
		strikeold[i] * pi / 180 +
		atan ((ydis + disc * cos (dipold[i] * pi / 180)) / nx[i]);
	      smc =
		acos (cos (sma) * cos (smb) +
		      sin (sma) * sin (smb) * cos (bgc));
	      bga = asin (sin (sma) * sin (bgc) / sin (smc));
	      lat[i] = 90 - smc * 180 / pi;	//absolute latitude
	      lon[i] = fmod (reflon + bga * 180 / pi, 360);	//absolute longitude
	      if (nx[i] == 0 && ny[i] == 0)	// set outer border to zero
		{
		  lat[i] = reflat;
		  lon[i] = reflon;
		}
	      i++;
	    }
	}
      lat[1] = reflat;
      lon[1] = reflon;

      // generate wavenumber vector
 //     kny_x = 1. / disc / 2.;	// Nyquist wave numbers=maximal wave number
 //     kny_y = 1. / disc / 2.;
      dk_x = 1. / oxdim;	// discretisation in wavenumber space
      dk_y = 1. / oydim;
      k_x = dvector (1, nns+1);
      k_y = dvector (1, nns+1);
      for (i = 0; i < nns; i++)
	{
	  l = floor (i / nsx);
	  if (l > nsy / 2)
	    {
	      l = floor (i / nsx) - nsy;
	    }
	  j = i % nsx;
	  if (j > nsx / 2)
	    {
	      j = i % nsx - nsx;
	    }
	  k_x[i + 1] = j * dk_x;
	  k_y[i + 1] = l * dk_y;
	}
      // 2dft of old distribution if given
      if (old == 1)
	{
	  FToldRE = dvector (1, nns+1);
	  FToldIM = dvector (1, nns+1);
	  dft2d (nsy, nsx, 0, REold, IMold, FToldRE, FToldIM);
	}

///////////////////////////////////////////////////////    
///// taper 
///// create taper (similar to tukey window but with only 1/4 of the cos period, alpha = fraction of length used for decay, alpha=1 = rectangular window, alpha=0 = similar Hann window)
///// if refdepth<1 there is no taper at the upper border
////////////////////////////////////////////////////////////////
      taper=dvector(1,nns);
      taper2=dvector(1,nns);
      alphax=0.51;
      alphay=0.51;
    /*  printf("alphax=%lf, alphay=%lf\n",alphax, alphay);
      if (nsx*(1-alphax)<9.0)
    	{ alphax=1.0-1.0/nsx*9.0;
    	}
      if (nsy*(1-alphay)<9.0)
    	{ alphay=1.0-1.0/nsy*9.0;
    	}  
    */  printf("alphax=%lf, alphay=%lf\n",alphax, alphay);

       for (i=1;i<=nns;i++)
    	     {  if (fabs(nx[i]-(nsx-1)/2*disc)>alphax*nsx/2*disc && fabs(ny[i]-(nsy-1)/2*disc)>alphay*nsy/2*disc)
                            {
                              taper[i]=cos(pi/2*(fabs(fabs(nx[i]-(nsx+1)/2*disc)-alphax*nsx/2*disc))/((1-alphax)*nsx/2*disc))*cos(pi/2*(fabs(fabs(ny[i]-(nsy+1)/2*disc)-alphay*nsy/2*disc))/((1-alphay)*nsy/2*disc));
                            }
                    else if (fabs(nx[i]-(nsx-1)/2*disc)>alphax*nsx/2*disc && fabs(ny[i]-(nsy-1)/2*disc)<=alphay*nsy/2*disc)
                            { taper[i]=cos(pi/2*(fabs(fabs(nx[i]-(nsx+1)/2*disc)-alphax*nsx/2*disc))/((1-alphax)*nsx/2*disc));
                            }
                    else if (fabs(nx[i]-(nsx-1)/2*disc)<=alphax*nsx/2*disc && fabs(ny[i]-(nsy-1)/2*disc)>alphay*nsy/2*disc)
                            { taper[i]=cos(pi/2*(fabs(fabs(ny[i]-(nsy+1)/2*disc)-alphay*nsy/2*disc))/((1-alphay)*nsy/2*disc));
                            }
                    else if (fabs(nx[i]-(nsx-1)/2*disc)<=alphax*nsx/2*disc && fabs(ny[i]-(nsy-1)/2*disc)<=alphay*nsy/2*disc)
                            {
                              taper[i]=1;
                            }
    		if (taper[i]<0)
    			{taper[i]=0;
    			}
                   }
            
       if (refdepth<2.0)
            { for (i=1;i<=nns;i++)
                 {  if ((ny[i]-(nsy-1)/2*disc)>alphay*nsy/2*disc)
                      { taper2[i]=cos(pi/2*(fabs(fabs(ny[i]-(nsy+1)/2*disc)-alphay*nsy/2*disc))/((1-alphay)*nsy/2*disc));
                      }
                    else 
                      { taper2[i]=1;
                      }
    		    if (taper2[i]<0)
    		      {taper2[i]=0;
    		      }
                 }
            }
    
//////////////////////////////////////////////////////
    //  create random field of size nns, apply taper and make 2D-DFT of it
      REslip   = dvector(1,nns+1);
      IMslip   = dvector(1,nns+1);
      REslipnew   = dvector(1,nns+1);
      IMslipnew   = dvector(1,nns+1);
      REsurf   = dvector(1,nns+1);
      IMsurf   = dvector(1,nns+1);
      FTsurfRE = dvector(1,nns+1);
      FTsurfIM = dvector(1,nns+1);
      FTrandRE = dvector(1,nns+1);
      FTrandIM = dvector(1,nns+1);
      FTrand2RE = dvector(1,nns+1);
      FTrand2IM = dvector(1,nns+1);
      specREold   = dvector(1,nns+1);
      specIMold   = dvector(1,nns+1);
      specREnew   = dvector(1,nns+1);
      specIMnew   = dvector(1,nns+1);
      for (i=1;i<=nns;i++)
            {// printf("seed=%ld\n",seed2);
              REslip[i]=ran1(&seed2);
              IMslip[i]=0;
            }
      dft2d(nsy,nsx,0, REslip, IMslip, FTrandRE, FTrandIM);

      // create a second random field for the variation of strike and dip
      for (i = 1; i <= nns; i++)
    	{
	  REslip[i] = ran1 (&seed2);
	  IMslip[i] = 0;
    	}
      dft2d (nsy, nsx, 0, REslip, IMslip, FTrand2RE, FTrand2IM);

    //  calculate spectrum with random phases and amplitude decay 1/(k)^(2) corresponding to a Hurst exponent of 1 for k higher than the corner wavenumber kc  
      kc=pow(10,(1.82-0.5*(Me))); 		// corner wavenumber according to Causse et al. (2010), 
      kcy=3.0*MAX(odiscx/oxdim, odiscy/oydim);
      kcy_int=kcy*0.8;
    // combine spectra of old and random distribution if old slip is given
      if(old==1)
    	{ // calculate difference between old and random spectrum to bring them to the same level
          l=0;
    	  diffspec=0;
    	  for(i=1;i<=nns;i++)
             { if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kcy+kcy_int && sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) > kcy-kcy_int)
    		 { diffspec=diffspec+log10(sqrt(FTrandRE[i]*FTrandRE[i]+FTrandIM[i]*FTrandIM[i])
                                     /( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)))
                                    -log10(sqrt(FToldRE[i]*FToldRE[i]+FToldIM[i]*FToldIM[i]));
    		   l++;
    		 }
             }
    	  diffspec=diffspec/l;
    	  printf("diffspec=%lf\n",diffspec);
          
    	  for(i=1;i<=nns;i++)
                      // for k< kcy-kcy_int keep only old spectrum
    		    { if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kcy-kcy_int)
                        { specREold[i]=FToldRE[i]*pow(10,diffspec);
    		          specIMold[i]=FToldIM[i]*pow(10,diffspec);
                          specREnew[i]=0.0;
    		          specIMnew[i]=0.0;
                        }
                      // for k in the combination interval decrease contribution of old spectrum with a sin-function
                      // and increase contribution of random spectrum with a sin function
                      // amplitude of the random spectrum is multiplied by 1/k² if k > kc
                    else if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kcy+kcy_int 
                             && sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) > kcy-kcy_int 
                             && sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) > kc)
                        { specREold[i]=FToldRE[i]*pow(10,diffspec)*1/2*(1-sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                          specIMold[i]=FToldIM[i]*pow(10,diffspec)*1/2*(1-sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                          specREnew[i]=FTrandRE[i]/pow(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)),roughslip)
                                       *1/2*(1+sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                          specIMnew[i]=FTrandIM[i]/pow(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)),roughslip)
                                       *1/2*(1+sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                        }
                      // for k in the combination interval decrease contribution of old spectrum with a sin-function
                      // and increase contribution of random spectrum with a sin function
                      // amplitude of random spectrum is multiplied by 1/kc² for k<kc
                    else if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kcy+kcy_int 
                             && sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) > kcy-kcy_int 
                             && sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kc)
                        { specREold[i]=FToldRE[i]*pow(10,diffspec)*1/2*(1-sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                          specIMold[i]=FToldIM[i]*pow(10,diffspec)*1/2*(1-sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                          specREnew[i]=FTrandRE[i]/(pow(fabs(kc), roughslip))*1/2*(1+sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                          specIMnew[i]=FTrandIM[i]/(pow(fabs(kc), roughslip))*1/2*(1+sin(pi*(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]),2))-kcy)/(2*kcy_int)));
                        }
                       // for k>kcy+kcy_int keep only random spectrum
                       // amplitude of the random spectrum is multiplied by 1/k² if k > kc
                     else
                        {  specREold[i]=0.0;
                           specIMold[i]=0.0;
                           specREnew[i]=FTrandRE[i]/pow(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)),roughslip);
                           specIMnew[i]=FTrandIM[i]/pow(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)),roughslip);
                        }
                     // do k-square sclaing for second random field
                     if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) > kc)
                        { FTsurfRE[i] =FTrand2RE[i] / pow(sqrt(pow (fabs (k_x[i]), 2.0) + pow (fabs (k_y[i]), 2.0)),roughness);
                	  FTsurfIM[i] =FTrand2IM[i] / pow(sqrt(pow (fabs (k_x[i]), 2.0) + pow (fabs (k_y[i]), 2.0)),roughness);
                        }
                     else 
                        { FTsurfRE[i] =FTrand2RE[i] / (pow (fabs (kc), roughness));
                          FTsurfIM[i] = FTrand2IM[i] / (pow (fabs (kc), roughness));
                        }
//                printf(" i=%d, FTsurfRE[i]=%lf, FTrand2RE=%lf \n",i, FTsurfRE[i],FTrand2RE[i]);
             }
          dft2d(nsy,nsx,1, specREold, specIMold, REslip, IMslip);
        }
      // for new slip distribution multiply with 1/kc² for k<kc and with k² for k>kc
      else
        { for(i=1;i<=nns;i++)
            { if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) < kc)
               { specREnew[i]=FTrandRE[i]/(pow(fabs(kc), roughslip));
                 specIMnew[i]=FTrandIM[i]/(pow(fabs(kc), roughslip));
                 FTsurfRE[i] =FTrand2RE[i] / (pow (fabs (kc), roughness));
                 FTsurfIM[i] = FTrand2IM[i] / (pow (fabs (kc), roughness));
               }
              else
               { specREnew[i]=FTrandRE[i]/pow(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)),roughslip);
                 specIMnew[i]=FTrandIM[i]/pow(sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)),roughslip);
                 FTsurfRE[i] =FTrand2RE[i] /pow(sqrt (pow (fabs (k_x[i]), 2.0) + pow (fabs (k_y[i]), 2.0)),roughness);
                 FTsurfIM[i] =FTrand2IM[i] /pow(sqrt (pow (fabs (k_x[i]), 2.0) + pow (fabs (k_y[i]), 2.0)),roughness);
                 
               }
             }
        }

      // inverse fourier transform to obtain final slip distribution
      dft2d (nsy, nsx, 1, specREnew, specIMnew, REslipnew, IMslipnew);
      dft2d (nsy, nsx, 1, FTsurfRE, FTsurfIM, REsurf, IMsurf);
    // set negativ values of slip to 0 by substracting the minimum value 
      minslip=REslip[1];
      for(i=1;i<=nns;i++)
        { if (minslip>REslipnew[i])
          {  minslip=REslipnew[i];
          }
        }
      if(old==1)
      // combine the two slip distributions
       {for(i=1;i<=nns;i++)
         { if (refdepth<2.0)
            { REslip[i]=REslip[i]*taper2[i]+(REslipnew[i]-minslip)*taper[i];
            }
           else
            { REslip[i]=REslip[i]+(REslipnew[i]-minslip)*taper[i];
            }
         }
       }
      else
       {for(i=1;i<=nns;i++)
           { REslip[i]=(REslipnew[i]-minslip)*taper[i];
           }
       }
      minslip=REslip[1];
      minsurf=REsurf[1];
      for(i=1;i<=nns;i++)
        { if (minslip>REslip[i])
            { minslip=REslip[i];
            }
	  if (REsurf[i] < minsurf)
	    { minsurf = REsurf[i];
//              printf("1. REsurf[i]=%lf\n",REsurf[i]);
	    }
	}
       for(i=1;i<=nns;i++)
        {  REslip[i]=REslip[i]-minslip;
           slipsum=slipsum+REslip[i];
        }
///// //////////////calkulate spectrum of final distribution for plotting
      dft2d(nsy,nsx,0, REslip, IMslip, specREnew, specIMnew);
//// print output to file
      kfile= "kxprofile.dat";
      fout=fopen(kfile,"w");
      kfile= "kyprofile.dat";
      fout2=fopen(kfile,"w");
      fprintf(fout,"#  kx slip_amplitude\n");
      fprintf(fout2,"#  ky slip_amplitude\n");
      printf("nns=%d",nns);
      for (k=1;k<=nns;k++)
        { if ( k_y[k]==0.0 )
            { fprintf(fout,"%8.4lf\t %8.4lf\n",k_x[k],sqrt(specREnew[k]*specREnew[k]+specIMnew[k]*specIMnew[k]));
            }
          if ( k_x[k]==0.0 )
            { fprintf(fout2,"%8.4lf\t %8.4lf\n",k_y[k],sqrt(specREnew[k]*specREnew[k]+specIMnew[k]*specIMnew[k]));
            }
	}
      fclose(fout2);
      fclose(fout);
      // adjust height of fractal surface to aspect ratio
      maxsurf = REsurf[1] - minsurf;
      for (i = 1; i <= nns; i++)
	      { REsurf[i] = REsurf[i] - minsurf;
	        if (REsurf[i] > maxsurf)
	          { maxsurf = REsurf[i];
	          }
	      }
      if (oxdim > oydim)
	    { height = oydim * aspect;
	    }
      else
	    { height = oxdim * aspect;
	    }
   
      fout = fopen ("surface", "w");
      fprintf (fout, "#  lat  lon x[km] y[km] heigth\n");
      for (i = 1; i <= nns; i++)
	    { REsurf[i] = REsurf[i] / maxsurf * height;
//	      fprintf (fout, "%8.4lf\t %8.4lf\t%lf\t%lf\t%lf\n",
//				   lat[i], lon[i], nx[i], ny[i], REsurf[i]);
	    }
      fclose(fout);
      // calculate deviations of strike and dip as local angles along the fractal surface
      printf ("disc=%lf \n", disc);
      dstrike = dvector (1, nns+1);
      ddip = dvector (1, nns+1);
      drake = dvector (1, nns+1);
      for (i = 1; i <= nns; i++)
	    { if (nx[i] == 0)
	        { dstrike[i] = strikeold[i] +atan ((REsurf[i + 1] - REsurf[i]) / disc) * 180.0 / pi;
	        }
	      else if (nx[i] == (nsx - 1) * disc)
	        { dstrike[i] = strikeold[i] + atan ((REsurf[i] - REsurf[i - 1]) / disc) * 180.0 / pi;
	        }
	      else
	        { dstrike[i] = strikeold[i] + atan ((REsurf[i + 1] - REsurf[i - 1]) / 2 / disc) * 180.0 / pi;
	        }
	  if (ny[i] == 0)
	    {
	      ddip[i] =
		dipold[i] +
		atan ((REsurf[i + (nsx)] - REsurf[i]) / disc) * 180.0 / pi;
	    }
	  else if (ny[i] == (nsy - 1) * disc)
	    {
	      ddip[i] =
		dipold[i] +
		atan ((REsurf[i] - REsurf[i - nsx]) / disc) * 180.0 / pi;
	    }
	  else
	    {
	      ddip[i] =
		dipold[i] +
		atan ((REsurf[i + nsx] -
		       REsurf[i - nsx]) / 2 / disc) * 180.0 / pi;
	    }
	  drake[i] = rakeold[i] + REsurf[i] / height * 20.0 - 10.0;
	}
      // create a third random field for the variation of rake
      for (i = 1; i <= nns; i++)
    	{
	  REsurf[i] = ran1 (&seed2);
	  IMsurf[i] = 0;
    	}
      dft2d (nsy, nsx, 0, REsurf, IMsurf, FTrand2RE, FTrand2IM);
      // do k-square sclaing for second random field
      for (i = 1; i <= nns; i++)
    	{ if (sqrt( pow(fabs(k_x[i]), 2.0)+ pow(fabs(k_y[i]), 2.0)) > kc)
            { FTsurfRE[i] =FTrand2RE[i] / (pow (fabs (k_x[i]), roughness) + pow (fabs (k_y[i]), roughness));
              FTsurfIM[i] = FTrand2IM[i] / (pow (fabs (k_x[i]), roughness) + pow (fabs (k_y[i]), roughness));
            }
          else 
            { FTsurfRE[i] = FTrand2RE[i] / (pow (fabs (kc), roughness));
              FTsurfIM[i] = FTrand2IM[i] / (pow (fabs (kc), roughness));
            }
        }
      dft2d (nsy, nsx, 1, FTsurfRE, FTsurfIM, REsurf, IMsurf);
      minsurf=REsurf[1];
      for(i=1;i<=nns;i++)
	{ if (REsurf[i] < minsurf)
	    { minsurf = REsurf[i];
	    }
	}
      // adjust height of fractal surface to aspect ratio
      maxsurf = REsurf[1] - minsurf;
      for (i = 1; i <= nns; i++)
	      { REsurf[i] = REsurf[i] - minsurf;
	        if (REsurf[i] > maxsurf)
	          { maxsurf = REsurf[i];
	          }
	      }
      for (i = 1; i <= nns; i++)
	    { drake[i] = rakeold[i] + REsurf[i] / maxsurf  * 20.0 - 10.0;
	    }

      //// determine a possible hypocentre: find patch with maximum slip, determine the nearest point at which slip = 1/2max
      maxslip = (REslip[1] - (minslip)) * taper[1] + 0.001;
      for (i = 1; i <= nns; i++)
	{
	  if (maxslip < (REslip[i] - (minslip)) * taper[i] + 0.001)
	    {
	      maxslip = (REslip[i] - (minslip)) * taper[i] + 0.001;
//	      maxlat = nx[i];
//	      maxlon = ny[i];
	    }
	}
/*      distance = oxdim;
      for (i = 1; i <= nns; i++)
	{
	  dist =
	    sqrt ((nx[i] - maxlat) * (nx[i] - maxlat) +
		  (ny[i] - maxlon) * (ny[i] - maxlon));
	  //         printf("lathypo=%lf lonhypo=%lf, dist=%lf, distance=%lf, 1/2*maxslip= %lf, slip=%lf, i=%d \n", lat1, lon1, dist, distance, 0.5*maxslip, ((REslip[i]-(minslip))*taper[i]+0.001),i);
	  if ((0.5 * maxslip > ((REslip[i] - (minslip)) * taper[i] + 0.001))
	      && (distance > dist) && taper[i] > 0.95)
	    {
	      lat1 = lat[i];
	      lon1 = lon[i];
	      xhypo = nx[i];
	      yhypo = ny[i];
	      distance = dist;
	      //printf("lathypo=%lf lonhypo=%lf, xhypo=%lf yhypo=%lf i=%d \n", lat1, lon1, xhypo, yhypo,i);
	    }
	}
      printf ("lathypo=%lf lonhypo=%lf, xhypo=%lf yhypo=%lf \n", lat1, lon1,
	      xhypo, yhypo);
   */

//////////////////////////////////////////////////////////////////////////////////////
   ///// calculate variation for rupture velocity
   ///////////////////////////////////////////////////////////////////////////////////
//   vel = dvector (1, nns+1);
//   for (i = 1; i <= nns; i++)
//	{
//	  vel[i] =
//	    ((REslip[i] - (minslip)) * taper[i] +
//	     0.001) / maxslip * velmean * fraction * 2 + velmean -
//	    velmean * fraction;
//	}
    }
//// print output to file
  strcat(sourceold, "_genso");
  fout = fopen (outname, "w");
  fout2 = fopen (sourceold, "w");
  fprintf (fout, "#  lat  lon x[km] y[km] vel[km/s]\n");
  fprintf (fout2,
	   "#  lat  lon  refdepth[km] slip[m] rake[grad] dstrike[grad]  ddip[grad] x[km] y[km] lat_old lon_old depthold \n");
  printf ("nns=%d\n", nns);
  if (smooth == 1) // output for smoothed slip distribution
          {
	    for (k = 1; k <= nns; k++)
		  { fprintf (fout2,
				   "%8.4lf\t %8.4lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				   lat[k], lon[k], z[k], REold[k], rakeold[k], strikeold[k],
				   dipold[k], nx[k], ny[k]);
		  }
            printf("Daten geschrieben\n");
	  }
  else if (smooth == 0) // output for refined slip distribution
	  { for (k = 1; k <= nns; k++)
			  { //printf("k=%d \n",k);
			       fprintf (fout2,
						   "%8.4lf\t %8.4lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n",
						   lat[k], lon[k], z[k],
						   REslip[k]/slipsum*moment+0.00001, drake[k],
						   dstrike[k], ddip[k], nx[k], ny[k], latold[k], lonold[k], depthold[k]);
//				    fprintf (fout, "%8.4lf\t %8.4lf\t %lf\t %lf\t %lf\t \n", lat[k],
//						   lon[k], vel[k], nx[k], ny[k]);

			  }

	  }
  fclose (fout2);
  fclose (fout);
  if (smooth == 0)
    {
      if (old == 1)
	{
	  free_dvector (FToldRE, 1, nns+1);
	  free_dvector (FToldIM, 1, nns+1);
	  free_dvector (REold, 1, nns+1);
	  free_dvector (IMold, 1, nns+1);
	  free_dvector (latold,1,nns+1);
	  free_dvector (lonold,1,nns+1);
	  free_dvector (depthold,1,nns+1);
      free_dvector (ox,1,ns);
      free_dvector (oy,1,ns);
      free_dvector (REslipo,1,ns);
      free_dvector (lato,1,ns);
      free_dvector (lono,1,ns);
      free_dvector (depo,1,ns);
      free_dvector (rake_v,1,ns);
      free_dvector (strike_v,1,ns);
      free_dvector (dip,1,ns);

	}
      free_dvector (nx, 1, nns+1);
      free_dvector (ny, 1, nns+1);
      free_dvector (k_x, 1, nns+1);
      free_dvector (k_y, 1, nns+1);
      free_dvector (lat, 1, nns+1);
      free_dvector (lon, 1, nns+1);
      free_dvector (z, 1, nns+1);
      free_dvector (taper, 1, nns+1);
      free_dvector (dstrike, 1, nns+1);
      free_dvector (ddip, 1, nns+1);
      free_dvector (drake, 1, nns+1);
      free_dvector (dipold, 1, nns+1);
      free_dvector (rakeold, 1, nns+1);
      free_dvector (strikeold, 1, nns+1);
      free_dvector (REslip, 1, nns+1);
      free_dvector (IMslip, 1, nns+1);
      free_dvector (FTrand2RE, 1, nns+1);
      free_dvector (FTrand2IM, 1, nns+1);
      free_dvector (FTrandRE, 1, nns+1);
      free_dvector (FTrandIM, 1, nns+1);
      free_dvector (specREnew, 1, nns+1);
      free_dvector (specIMnew, 1, nns+1);
      free_dvector (specREold, 1, nns+1);
      free_dvector (specIMold, 1, nns+1);
      free_dvector (REsurf, 1, nns+1);
      free_dvector (IMsurf, 1, nns+1);
      free_dvector (FTsurfRE, 1, nns+1);
      free_dvector (FTsurfIM, 1, nns+1);
    }
  else if (smooth == 1)
    {
      free_dvector (nx, 1, nns+1);
      free_dvector (ny, 1, nns+1);
      free_dvector (lat, 1, nns+1);
      free_dvector (lon, 1, nns+1);
      free_dvector (z, 1, nns+1);
      free_dvector (REold, 1, nns+1);
      free_dvector (dipold, 1, nns+1);
      free_dvector (rakeold, 1, nns+1);
      free_dvector (latold,1,nns+1);
      free_dvector (lonold,1,nns+1);
      free_dvector (depthold,1,nns+1);
    }
  printf ("... soumod end!\n");
  return nns;
}
