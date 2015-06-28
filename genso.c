/*
 * genso.c
 *
 * generates complete source description for a kinematic source model
 * reads input from file "infile" (user defined) (input file should be equivalent to the QSCMP input file)
 * outputs one file XXX_genso for each source segment input which includes refined slip distribution,
 * 			patch rupture times, patch rise times
 *
 * output: rise time is defined as 1/(2*pi*fc) , where fc is the corner frequency of the source time function.
 *                   The rise time is the time from rupture initiation to the peak moment rate.
 *
 *  Created on: Jun 6, 2013
 *  last modified: Jun 25, 2013
 *      Author: Katrin Kieling
 *
 *
 */
#include "soumod.c"
#include "velocity.c"
#include "risetime.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main (void)
  { FILE *fin, *fout;
    char infile[80];
    char greensfile[80]="";
    char data[256];
    char modtype[10], taperswitch[4];
    char sourcefile[80];
    char dummy[256];
    char *pch;
    char *fine_file[80];
    int  num_seg=0, num_model_depth=0, *num_patch, nt=0;
    int  sum_patch=0;
    int  seed,i,k,l;
    int kstart, kstop;
    double Me;
    double lathypo, lonhypo, dephypo;
    double rupvel;
    double *moment, *new_dstep, *seglat, *seglon, *segdepth;
    double *s_depth,*vp,*vs,*density;
    double vprec, vsrec, rhorec;
    double *lat, *lon, *depth, *slip, *rake, *strike, *dip, *time1, *nx, *ny,
            *rise;
    double *optrise, *opttime1;
    double *rupnx, *rupny, *rupz, *time2;
    double duration=0.0;
    double *dip1;
    double dt, window, fhighcut;

    //// read input file
    printf ("name of input file?\n");
    scanf ("%s", infile);
    fin = fopen (infile, "r");
    while (fgets (data, 256, fin) != NULL)
    	{ sscanf(data,"%s\n",dummy);
    	  if ( strncmp (dummy, "#",1 ) != 0)
    		  {strncat(greensfile,&(dummy[1]), strlen(dummy)-2);
    		   fgets (data, 256, fin);
    		   sscanf(data,"%s\n", dummy);
    		   strncat(greensfile, &(dummy[1]), strlen(dummy)-2);
    		   printf("greensfile %s \n",greensfile);
    		   break;
    		  }
    	}
    /*read velocity structure from greensfile*/
    velocity(greensfile, &s_depth, &vp, &vs,&density, &num_model_depth,
    		&vprec, &vsrec, &rhorec, &nt, &window);
    dt=window/nt;
    /* continue reading input file for source segments */
    i=0;
    while (fgets (data, 256, fin) != NULL)
    	{ sscanf(data,"%*2s %s %*s\n",dummy);
    	  if (strcmp (dummy, "n_seg") == 0)
    		  { fgets (data, 256, fin);
    		    fgets (data, 256, fin);

    			sscanf(data,"%d %lf %lf %lf %lf %lf %*s\n", &num_seg,&Me, &fhighcut, &lathypo, &lonhypo, &dephypo);
    			moment = dvector (1, num_seg);
    			new_dstep=dvector (1, num_seg);
    			seglat=dvector (1, num_seg);
    			seglon=dvector (1, num_seg);
    			segdepth=dvector (1, num_seg);
    			num_patch=ivector (1, num_seg);
    			dip1=dvector(1,num_seg);
    		    printf("num_seg=%d, Me=%lf \n",num_seg,Me);
    		  }
    	  if (strcmp (dummy, "file_name") == 0)
    		  { fgets (data, 256, fin);
    		    while (i<num_seg)
    		    	{ fgets (data, 256, fin);

    		    	  pch=strchr(data,'d');
    		    	  while (pch!=NULL)
    		    	   	{ *pch='E';
    		    	   	   pch=strchr(pch+1,'d');
    		    	   	}
    		    	  sscanf(data,"%*d %le %lf %lf %lf \n",&moment[i],&seglat[i], &seglon[i], &segdepth[i]);
    		          printf("segment=%d\n",i);
    		          fgets (data, 256, fin);
    		          fgets (data, 256, fin);
    		          sscanf(data,"%lf %d %s %s\n", &rupvel, &seed, taperswitch, modtype);
    		          fgets (data, 256, fin);
    		          sscanf(data,"%s %lf\n", &sourcefile[0], &new_dstep[i]);
     		          printf("sourcefile = %s, new_dstep=%lf seglat=%lf, seglon=%lf, segdepth=%lf, moment[is]=%le\n",
     		        		  &sourcefile[0], new_dstep[i], seglat[i],seglon[i],segdepth[i], moment[i]);

     		          /////////////////////////////////////////////////////
     		          //// refine slip distribution and calculate rupture times
     		          /////////////////////////////////////////////////////
    		          num_patch[i]= soumod(&sourcefile[0], modtype, new_dstep[i], Me, rupvel,
    		        		  taperswitch, seed, &density, num_model_depth, lathypo, lonhypo, dephypo,
    		        		  seglat[i], seglon[i], segdepth[i], &dip1[i], moment[i]);
    		          sum_patch=sum_patch+num_patch[i];
    		          printf("sum_patch=%d\n ", sum_patch);
    		          fine_file[i] = malloc(strlen(sourcefile)+1);
    		          if ( fine_file[i]== NULL)
    		        	{
    		              printf("stop: could not allocate fine_file\n");
    		              exit(1);
    		            }
    		          strcpy(fine_file[i], sourcefile);
    		          printf("fine_file=%s, i=%d\n", fine_file[i],i);
    		          i++;
    		    	}
    		  }
    	}
    // allocate memory for  all patches:

    fclose (fin);
    lat = dvector (1, sum_patch);
    lon = dvector (1, sum_patch);
    depth = dvector (1, sum_patch);
    slip = dvector (1, sum_patch);
    rake = dvector (1, sum_patch);
    strike = dvector (1, sum_patch);
    dip = dvector (1, sum_patch);
    opttime1 = dvector (1, sum_patch);
    nx = dvector (1, sum_patch);
    ny = dvector (1, sum_patch);
    optrise = dvector (1, sum_patch);
    time1 = dvector (1, sum_patch);
    rise = dvector (1, sum_patch);
        // read all patches
    k=1;
    for (i=0;i<num_seg; i++)
        { printf("fine_file=%s\n", fine_file[i]);
          fin = fopen (fine_file[i], "r");
          l=1;
          kstart=k;
          time2 = dvector (1, num_patch[i]);
          rupnx = dvector (1, num_patch[i]);
          rupny = dvector (1, num_patch[i]);
          rupz = dvector (1, num_patch[i]);
          while (fgets (data, 256, fin) != NULL)
        	{ if (data[0]!='#')
        		{ sscanf(data,"%lf %lf %lf %lf %lf %lf %lf %lf %lf \n",&lat[k],&lon[k],
        				&depth[k],&slip[k],&rake[k],&strike[k],&dip[k], &nx[k], &ny[k]);
        		  rupnx[l]=nx[k];
        		  rupny[l]=ny[k];
        		  rupz[l]=depth[k];

        		  k++;
        		  l++;
        		}
        	}
          printf("num_patch=%d, l=%d\n",num_patch[i],l);
          fclose(fin);
          kstop=k;
          //////////////////////////////////////////////////////////////////////////////////////
          /// Calculate rupture time
          /////////////////////////////////////////////////////////////////////////////////////
          printf(" dip1=%lf, lathypo=%lf, lonhypo=%lf, seglat=%lf, seglon=%lf\n",dip1[i], lathypo, lonhypo, seglat[i], seglon[i]);
          rupfront(&s_depth, &vs, num_model_depth, rupnx,rupny, rupz, num_patch[i], new_dstep[i],
                   lathypo, lonhypo, dephypo, seglat[i], seglon[i], segdepth[i],
                   dip1[i], rupvel,time2);
          l=1;
          for (k=kstart;k<kstop; k++)
        	  { //printf("k=%d, l=%d, time2=%lf\n", k,l,time2[k]);
                    time1[k]=time2[l];
        		l++;
			    if (time1[k]>duration)
			      { duration=time1[k];
			      }
        	  }
              free_dvector (time2,1,num_patch[i]);
    	}

//
    risetime(num_seg, sum_patch, moment, depth, slip, rake, strike, dip,
    		time1, num_model_depth, s_depth, vp, vs, density, Me,
    		vprec, vsrec, rhorec, duration, dt, fhighcut, rise,rupvel);
 
    fout = fopen ("source.out", "w");
    fprintf (fout,"%d\n", sum_patch-1);
    for (k=1;k<sum_patch; k++)
        { fprintf (fout,
				 "%8.4lf\t %8.4lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
				 lat[k], lon[k], depth[k], slip[k], strike[k],
				 dip[k],rake[k], time1[k],rise[k], nx[k], ny[k]);
        }
    fclose(fout);
      
    printf("Daten geschrieben\n");

    free_dvector (moment, 1, num_seg);
    free_dvector (new_dstep,1,num_seg);
    free_dvector (seglat,1,num_seg);
    free_dvector (seglon,1,num_seg);
    free_dvector (segdepth,1,num_seg);
    free_ivector (num_patch,1, num_seg);
    free_dvector (s_depth,1,num_model_depth);
    free_dvector (vs,1,num_model_depth);
    free_dvector (vp,1,num_model_depth);
    free_dvector (density,1,num_model_depth);
    free_dvector (rise,1,sum_patch);
    free_dvector (time1,1,sum_patch);
    free_dvector (lat,1,sum_patch);
    free_dvector (lon,1,sum_patch);
    free_dvector (depth,1,sum_patch);
    free_dvector (slip,1,sum_patch);
    free_dvector (rake,1,sum_patch);
    free_dvector (strike,1,sum_patch);
    free_dvector (dip,1,sum_patch);
    free_dvector (ny,1,sum_patch);
    free_dvector (nx,1,sum_patch);
    free_dvector (optrise,1,sum_patch);
    free_dvector (opttime1,1,sum_patch);

	return 0;
  }
