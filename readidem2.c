#include <math.h>
void readidem2( char *sourceold, int *pns, double *reflat, double *reflon, double *refdepth, double *rake,
		double *strike, double **dip, double *oxdim, double *oydim, int *nsx, int *nsy, double *odiscx, double *odiscy,
		double **ox, double **oy, double **REslipo, double **rake_v, double **strike_v, double **lato, double **lono,
		double **depo)
{ FILE *fin2;
  int i, k,l;
  char data[256];
  double lat1, lon1, lat2, lon2,depth1,dip1;
  double pi = 4.0*atan(1.0);
  double zref;
  fin2 = fopen(sourceold,"r");          // open old slip distribution
  i=0;
  while (fgets(data,256,fin2)!=NULL)
    { if (data[0]!='#')
	{ i++;
	}
    }
  *pns=i;
  rewind(fin2);
  fgets(data,256,fin2);
  *reflat=99999.0;
  while (*reflat==99999.0)
    {if (data[0]!='#')
	{ sscanf(data,"%lf %lf %lf %*f %lf %lf %lf\n",reflat,reflon,refdepth,rake,strike,&dip1);
	}
      else { fgets(data,256,fin2);
      }
    }
  rewind(fin2);
  *oxdim=0.0;
  *oydim=0.0;
  i=0;
  while (fgets(data,256,fin2)!=NULL)
    {
      if (data[0]!='#')
	{ sscanf(data,"%lf %lf %lf %*f %*f %lf %lf\n",&lat1,&lon1,&depth1,strike,&dip1);
	  if (depth1==*refdepth)
	    { lat2=lat1;
	      lon2=lon1;
	      i++;
	    }
	  else
		{*oxdim=180.0/pi*111.12*acos(sin(*reflat/180*pi)*sin(lat2/180*pi)+cos(*reflat/180*pi)*cos(lat2/180*pi)*cos((*reflon-lon2)/180*pi));
		 *oydim=(depth1-*refdepth)/sin(dip1/180*pi);
	  }
	}
    }
  rewind(fin2);
  *nsx=i;
  *nsy=(int) *pns/(*nsx);
  *odiscx=(*oxdim/(*nsx-1));                 //old discretisation along strike
  *oxdim=*oxdim+(*odiscx);
  *odiscy=  (*oydim/(*nsy-1)); //old discretisation along dip
  printf("oxdim=%lf, oydim=%lf, nsx=%d, nsy=%d, pns=%d\n",*oxdim,*oydim, *nsx ,*nsy,*pns);

  
  // caclculate positiopns in local cartesian coordinate system
  *ox        = dvector(1,*pns);
  *oy        = dvector(1,*pns);
  *REslipo   = dvector(1,*pns);
  *lato	     = dvector(1,*pns);
  *lono      = dvector(1,*pns);
  *depo      = dvector(1,*pns);
  *rake_v    = dvector(1,*pns);
  *strike_v    = dvector(1,*pns);
  *dip       = dvector(1,*pns);
  i=1;
  k=-1;
  l=0;
  zref=0;

  while (fgets(data,256,fin2)!=NULL)
    { if (data[0]!='#')
	  {sscanf(data,"%lf %lf %lf %lf %lf %lf %lf\n",&lat1,&lon1,&depth1,&(*REslipo)[i],&(*rake_v)[i],&(*strike_v)[i], &(*dip)[i]);
	   if (depth1>zref)
	    { zref=depth1;

	      k++;
	    }
  	  else if (depth1<=zref)
	    { zref=depth1;

	      k=0;
	      l++;
        }
      (*ox)[i]=*odiscx*l;
	  (*oy)[i]=*odiscy*k;
	  (*lato)[i]=lat1;
	  (*lono)[i]=lon1;
	  (*depo)[i]=depth1;
	//  printf("depo[%i]=%lf, ox[%d]=%lf, oy=%lf\n",i, (*depo)[i],i, (*ox)[i],i, *oy[i]);
	  i++;
	  }
    }
  printf("oxdim=%lf, oydim=%lf, odiscx=%lf, odiscy=%lf\n",*oxdim,*oydim, *odiscx, *odiscy);
  fclose(fin2);
}


