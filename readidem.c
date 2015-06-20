#include <math.h>
void readidem( char *sourceold, int *ns, double *reflat, double *reflon, double *refdepth,
		       double *rake, double *strike, double **dip, double *oxdim, double *oydim, int *nsx, int *nsy,
		       double *odiscx, double *odiscy, double **ox, double **oy, double **REslipo,
		       double **rake_v, double **strike_v, double **lato, double **lono, double **depo)
{ FILE *fin2;
  int i,k;
  char data[256];
  double lat1, lon1, lat2, lon2,depth1,dip1;
  double pi = 4.0*atan(1.0);
  double xref, yref, zref;
  fin2 = fopen(sourceold,"r");          // open old slip distribution
  i=0;
  while (fgets(data,256,fin2)!=NULL)
    { if (data[0]!='#')
	{ i++;
	}
    }
  *ns=i;
  rewind(fin2);
  fgets(data,256,fin2);
  *reflat=99999.0;
  while (*reflat==99999.0)
    {if (data[0]!='#')
	   { sscanf(data,"%lf %lf %lf %*f %lf %lf %lf\n",reflat,reflon,refdepth,rake,strike,&dip1);
	   }
     else
       { fgets(data,256,fin2);
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

	    }
	}
    }
  rewind(fin2);
  *nsx=i;
  *nsy=(int) *ns/(*nsx);
  *oydim=(depth1-*refdepth)/sin(dip1/180*pi);
  *odiscy=*oydim/(*nsy-1);                   //old discretisation along dip
  *oydim=*oydim+*odiscy;
  *odiscx=(*oxdim/(*nsx-1));                 //old discretisation along strike
  *oxdim=*oxdim+*odiscx;
  printf("oxdim=%lf, oydim=%lf, odiscx=%lf, odiscy=%lf, nsx=%d, nsy=%d\n",*oxdim,*oydim, *odiscx, *odiscy, *nsx, *nsy);

  
  // caclculate positions in local cartesian coordinate system
  *ox        = dvector(1,*ns);
  *oy        = dvector(1,*ns);
  *REslipo   = dvector(1,*ns);
  *lato      = dvector(1,*ns);
  *lono      = dvector(1,*ns);
  *depo      = dvector(1,*ns);
  *rake_v    = dvector(1,*ns);
  *strike_v  = dvector(1,*ns);
  *dip       = dvector(1,*ns);
  i=1;
  k=0;
  zref=0.0;
  xref=*reflat;
  yref=*reflon;
  while (fgets(data,256,fin2)!=NULL)
    { if (data[0]!='#')
	{ sscanf(data,"%lf %lf %lf %lf %lf %lf %lf\n",&lat1,&lon1,&depth1,&(*REslipo)[i],&(*rake_v)[i], &(*strike_v)[i],&(*dip)[i]);
	  //printf("lat1[%d]=%lf,lon1[%d]=%lf, depth=%lf, slip=%lf\n",i,lat1,i,lon1,depth1, (*REslipo)[i]);
	  if (depth1>zref)
	    { zref=depth1;
	      k++;
	      xref=lat1;
	      yref=lon1;
	    }
	  (*oy)[i]=(k-1)*(*odiscy);
          (*ox)[i]=180.0/pi*111.12*acos(sin(xref/180.0*pi)*sin(lat1/180.0*pi)+cos(xref/180.0*pi)*cos(lat1/180.0*pi)*cos((yref-lon1)/180.0*pi));
	  (*lato)[i]=lat1;
	  (*lono)[i]=lon1;
	  (*depo)[i]=depth1;
	  i++;
	}
    }
  fclose(fin2);
}

