#include <math.h>
void readxyplane( char *sourceold, char *data, int *ns, double *reflat, double *reflon, double *refdepth, double *rake, double *strike, double **dip, double *oxdim, double *oydim, int *nsx, int *nsy, double *odiscx, double *odiscy, double **ox, double **oy, double **REslipo, double **rake_v, double **lato, double **lono, double **depo)
 {
  FILE *fin2;
  int i;
  double minx,miny,maxy;
  double hypolat, hypolon, hypodepth, dip1;
  double sma, smb, smc, bga, bgc,pi, REARTH;
  pi      = 4.0*atan(1.0);
  REARTH  = 6.371*pow(10,3); // radius of earth in km

  fin2 = fopen(sourceold,"r");          // open old slip distribution
  i=0;
  while (fgets(data,256,fin2)!=NULL)
    { if (data[0]!='#')
	{ i++;
	}
    }
  *ns=i;
  rewind(fin2);

  *ox        = dvector(1,*ns);
  *oy        = dvector(1,*ns);
  *REslipo   = dvector(1,*ns);
  *lato      = dvector(1,*ns);
  *lono      = dvector(1,*ns);
  *depo      = dvector(1,*ns);
  *rake_v    = dvector(1,*ns);
  *dip       = dvector(1,*ns);
  minx=0.0;
  miny=0.0;
  maxy=0.0;
  fgets(data,256,fin2);
  sscanf(data,"%*s %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf",&hypolat,&hypolon,&hypodepth,strike,&dip1,rake);
  //printf("reflat %lf reflon %lf refdepth %lf strike %lf dip %lf rake %lf\n",hypolat,hypolon,hypodepth,*strike,dip1,*rake);

  i=1;
  while (fgets(data,256,fin2)!=NULL)
    {if (data[0]!='#')
	{ sscanf(data,"%lf %lf %lf %lf \n",&(*ox)[i],&(*oy)[i],&(*rake_v)[i],&(*REslipo)[i]);
	  if ((*ox)[i]<minx)
	    {minx=(*ox)[i];
	    }
	  if ((*oy)[i]<miny)
	    {miny=(*oy)[i];
	    }
	  if((*oy)[i]>maxy)
	    {maxy=(*oy)[i];
	    }
	  (*dip[i])=dip1;
	  i++;
	}
     }
  rewind(fin2); 
    
  if ((*ox)[1]!=(*ox)[2])
    { *odiscx=sqrt((*ox)[1]*(*ox)[1])-sqrt((*ox)[2]*(*ox)[2]);
      i=1;
      while ((*oy)[i]==(*oy)[1])
	{i++;
	}
      *nsx=i-1;
      *odiscy=sqrt((*oy)[1]*(*oy)[1])-sqrt((*oy)[i+1]*(*oy)[i+1]);
      *nsy=*ns/(*nsx);
    }
  else if ((*oy)[1]==(*oy)[2])
    { *odiscy=sqrt((*oy)[1]*(*oy)[1])-sqrt((*oy)[2]*(*oy)[2]);
      i=1;
      while ((*ox)[i]==(*ox)[1])
	{ i++;
	}
      *nsy=i-1;
      *odiscx=sqrt((*ox)[1]*(*ox)[1])-sqrt((*ox)[i+1]*(*ox)[i+1]);
      *nsx=*ns/(*nsy);
    }

  i=1;
  while(i<=*ns)
    {(*ox)[i]=(*ox)[i]-(minx);
      (*oy)[i]=((*oy)[i]-(maxy))*(-1);
      i++;
    }
  *oxdim=((*nsx)-1)*(*odiscx);
  *oydim=((*nsy)-1)*(*odiscy);

  *refdepth=hypodepth+miny*sin((90-dip1)*pi/180);  // absolute depth
  sma=sqrt(pow((minx),2)+pow((miny)*cos((90-dip1)*pi/180),2))/REARTH;	
  smb=0.5*pi-hypolat*pi/180;
  bgc=fmod(*strike+180.0,360.0)*pi/180+atan(miny*cos((90-dip1)*pi/180)/(minx));
  smc=acos(cos(sma)*cos(smb)+sin(sma)*sin(smb)*cos(bgc));
  bga=asin(sin(sma)*sin(bgc)/sin(smc));
  *reflat=90-smc*180/pi;			//absolute latitude
  *reflon=fmod(hypolon+bga*180/pi,360);   //absolute longitude

  fclose(fin2);
 }

