/*
 * readidem3.c
 *
 *  Created on: Jul 4, 2013
 *      Author: katrin
 */

void readidem3( char *sourceold, int *pns, double *reflat, double *reflon, double *refdepth, double *rake,
		double *strike, double **dip, double *oxdim, double *oydim, int *nsx, int *nsy, double *odiscx, double *odiscy,
		double **ox, double **oy, double **REslipo, double **rake_v, double **strike_v, double **lato, double **lono,
		double **depo)
{ FILE *fin2;
  int i, l, nsxy;
  char data[256];
  double lat1, lon1,depth1,dip1, ox1, oy1;
  printf("readidem3sou\n");
  fin2 = fopen(sourceold,"r");          // open old slip distribution

  printf("sourceold=%s\n",sourceold);
  i=0;
  if (sourceold == NULL) perror ("Error opening file");
  while (fgets(data,256,fin2)!=NULL)
    { //printf("i=%d\n",i);
      if (data[0]!='#')
	    { i++;
	    }
    }
  printf("i=%d\n",i);
  *pns=i;
  rewind(fin2);
  fgets(data,256,fin2);
  *reflat=99999.0;
  while (*reflat==99999.0)
    { printf("reflat=%lf\n", *reflat);
      if (data[0]!='#')
	    { sscanf(data,"%lf %lf %lf %*f %lf %*f %lf\n",reflat,reflon,refdepth,rake,&dip1);
	    }
      else
    	{ fgets(data,256,fin2);
        }
    }
  printf("reflat=%lf\n", *reflat);
  rewind(fin2);
  *oxdim=0.0;
  *oydim=0.0;
  i=0;
  l=0;
  while (fgets(data,256,fin2)!=NULL)
    {
      if (data[0]!='#')
	    { sscanf(data,"%lf %lf %lf %*f %*f %*f %*f %lf %lf  \n",&lat1,&lon1,&depth1,&ox1,&oy1);
	      if (*oxdim<ox1)
	        { *oxdim=ox1;
	          i++;
	        }
	      if (*oydim<oy1)
	    	{ *oydim=oy1;
	          l++;
	        }
	    }
    }
  rewind(fin2);
  *nsx=i;
  *nsy=l;
  if(*nsx*(*nsy)!=*pns)
	  { printf("nsx*nsy!=pns\n");
	    exit(1);
	  }

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

  while (fgets(data,256,fin2)!=NULL)
    { if (data[0]!='#')
	  {sscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&(*lato)[i],&(*lono)[i],&(*depo)[i],&(*REslipo)[i],
			  &(*rake_v)[i],&(*strike_v)[i], &(*dip)[i], &(*ox)[i], &(*oy)[i]);

	//  printf("depo[%i]=%lf, ox[%d]=%lf, oy=%lf\n",i, (*depo)[i],i, (*ox)[i],i, *oy[i]);
	  i++;
	  }
    }
  nsxy=i-1;
  *odiscx=((*ox)[i-1]-(*ox)[1])/(*nsx-1);                 //old discretisation along strike
  *odiscy=((*oy)[i-1]-(*oy)[1])/(*nsy-1); //old discretisation along dip
  for (i=1; i<=nsxy; i++)
    {(*ox)[i]=(*ox)[i]-(*odiscx);
     (*oy)[i]=(*oy)[i]-(*odiscy);
    }
  printf("oy[i-1]=%lf, oy[1]=%lf, nsy=%d\n",(*oy)[i-1],(*oy)[1],*nsy-1);
  printf("oxdim=%lf, oydim=%lf, odiscx=%lf, odiscy=%lf, i=%d\n",*oxdim,*oydim, *odiscx, *odiscy,i);
  fclose(fin2);
}

