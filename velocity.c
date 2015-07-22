/*
 * velocity.c
 * reads the velocity structure of Greens Functions created by QSCGRN
 *  Created on: Jun 17, 2013
 *      Author: katrin
 *
 *  required input: greensfile  -- File in which Greens function is stored
 *
 *  output: vectors s_depth, vp, vs and density containing
 *    		the used source depth and its according model vp vs and density
 *
 */
#include <stdlib.h>

void
velocity (char *greensfile, double **s_depth, double **vp, double **vs, double **density, int *num_model_depth,
		double *vprec, double *vsrec, double *rhorec, int *nt, double *window)
{ FILE *fin;
  char data[256];
  char dummy[256]="";
  char *pch;
  int i=0;

  fin = fopen(greensfile, "r");
  while (fgets (data, 256, fin) != NULL)
  	  { sscanf(data,"%*s %s %*s %*s %*s\n",dummy);
  	    if (strcmp (dummy, "depth[m]") == 0)
  	    	{ fgets (data, 256, fin);
  	    	  pch=strchr(data,'D');
  	    	  while (pch!=NULL)
  	    	    { *pch='E';
  	    		  pch=strchr(pch+1,'D');
  	    		}
  	    	  sscanf(data,"%*e %le %le %le \n", vprec, vsrec, rhorec);
  	    	}
  	    if (strcmp (dummy, "nt") == 0)
  	    	{ fgets (data, 256, fin);
  	    	  pch=strchr(data,'D');
  	    	  while (pch!=NULL)
  	    	  	{ *pch='E';
  	    	  	  pch=strchr(pch+1,'D');
  	    	  	}
  	    	  sscanf(data,"%d %lf %*s\n", nt, window);
  	    	}
		sscanf(data,"%*s %*s %*s %*s %s\n",dummy);
  	    if (strcmp (dummy, "depth:") == 0)
  	    	{ i++;
  	    	}
  	  }
  if (i==0)
	  { printf("Reading of velocity structure was not succesful!\n");
	  	exit(1);
	  }
  rewind(fin);
  *num_model_depth=i;
  *s_depth=dvector(1,i);
  *vp=dvector(1,i);
  *vs=dvector(1,i);
  *density=dvector(1,i);
  i=0;
  while (fgets (data, 256, fin) != NULL)
	  { sscanf(data,"%*s %*s %*s %*s %s\n",dummy);
	    if (strcmp (dummy, "depth:") == 0)
	    	{ i++;
	    	  fgets (data, 256, fin);
	    	  fgets (data, 256, fin);
	    	  pch=strchr(data,'D');
	    	  while (pch!=NULL)
	    	    { *pch='E';
	    	       pch=strchr(pch+1,'D');
	    	    }
	    	  sscanf(data,"%le %le %le %le \n",&(*s_depth)[i], &(*vp)[i], &(*vs)[i], &(*density)[i]);
	    	  printf("s_depth=%lf, vp=%lf, vs=%lf, density=%lf, k=%d \n", (*s_depth)[i], (*vp)[i], (*vs)[i], (*density)[i],i);
	    	}
	  }
  *num_model_depth=i;
  fclose (fin);
  }
