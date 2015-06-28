/*
 * risetime.c
 *  calculates patch rise times
 *  Created on: Jun 25, 2013
 *      Author: Katrin Kieling
 *
 * output: rise time is defined as 1/(2*pi*fc) , where fc is the corner frequency of the source time function.
 *                   The rise time is the time from rupture initiation to the peak moment rate.
 */
#include <complex.h>
void
risetime (int num_seg, int sum_patch, double *moment, double *depth, double *slip, double *rake, double *strike, double *dip,
		double *time1, int num_model_depth, double *s_depth, double *vp, double *vs, double *rho, double menergy,
		double vprec, double vsrec, double rhorec, double duration, double dt, double fhighcut, double *rise, double rupvel)
  { double momsum=0.0, dmsum;
    double sm[3][3];
    double *vps, *vss, *rhos;
    double *fcp;
    double st, ra, di;
    double es0, es, tau, scale;
    double risemean=0.0;
    double dfst;
    double maxfreq=20.0;
    int nt=512;
    double fst[nt];
    double t0;
    double complex sfct[nt];
    int i,k, jfst ;
    int ntst;
    int numite=0;
    FILE *fout;
    char *outfile, strrupvel[80];
    
    double  DEG2RAD=1.745329252e-02;
    double  PI= 4.0 * atan (1.0);
/////////////////////////////////////////////////////////
    /// Make correction for disturbed moment tensor: total moment should be correct
////////////////////////////////////////////////////////
    vps = dvector (1, sum_patch);
    vss= dvector (1, sum_patch);
    rhos = dvector (1, sum_patch);
    fcp = dvector (1, sum_patch);
    /////////////////should be checked////////////
    tau=duration/10.0;
    printf("nt=%d\n", nt);
   
    for (i=0; i<num_seg; i++)
    	{momsum=momsum+moment[i];
         printf("momsum=%le moment=%le\n", momsum, moment[i]);
    	}
    for(i=0; i<3; i++)
    	{ for(k=0; k<3; k++)
    		{ sm[i][k]=0.0;
    		}
    	}
    for (i=1; i<=sum_patch; i++)
    	{ st=strike[i]*DEG2RAD;
          di=dip[i]*DEG2RAD;
          ra=rake[i]*DEG2RAD;
          vss[i]=0.0;
          k=1;
          while(vss[i]==0.0)
   	    { 
              if(depth[i]*1000.0<=s_depth[k])
                { vss[i]= vs[k];
                  vps[i]=vp[k];
                  rhos[i]= rho[k];
                  k=1;
                }
              else if (depth[i]*1000.0<s_depth[num_model_depth])
                { k++;
                }
              else
                { vss[i]= vs[num_model_depth];
                  vps[i]= vp[num_model_depth];
                  rhos[i]=rho[num_model_depth];
                  k=1;
                }

            }
          fcp[i]=1.0/(2*PI*tau)*100;
          fcp[i]=fcp[i]*(vss[i]/vsrec)/sqrt(slip[i]);
          slip[i]=rhos[i]*vss[i]*vss[i]*slip[i];

          sm[0][0]=sm[0][0]+slip[i]*(-sin(di)*cos(ra)*sin(2.0*st)
        		  -sin(2.0*di)*sin(ra)*(sin(st)*sin(st)));
          sm[1][1]=sm[1][1]+slip[i]*( sin(di)*cos(ra)*sin(2.0*st)
                  -sin(2.0*di)*sin(ra)*(cos(st)*cos(st)));
          sm[0][1]=sm[0][1]+slip[i]*( sin(di)*cos(ra)*cos(2.0*st)
        		  +0.50*sin(2.0*di)*sin(ra)*sin(2.0*st));
          sm[1][2]=sm[1][2]+slip[i]*(-cos(di)*cos(ra)*sin(st)
                  +cos(2.0*di)*sin(ra)*cos(st));
          sm[2][0]=sm[2][0]+slip[i]*(-cos(di)*cos(ra)*cos(st)
        		  -cos(2.0*di)*sin(ra)*sin(st));
    	}
    sm[2][2]=-(sm[0][0]+sm[1][1]);
    dmsum=sqrt(0.50*(sm[0][0]*sm[0][0]+sm[1][1]*sm[1][1]
                        +sm[2][2]*sm[2][2])+sm[0][1]*sm[0][1]
                        +sm[1][2]*sm[1][2]+sm[2][0]*sm[2][0]);
    printf("%lf, %lf, %lf, %lf, %lf, %lf\n",sm[0][0],sm[1][1],sm[2][2],sm[0][1],sm[1][2],sm[2][0]);
    printf("momsum=%lf, dmsum=%lf\n", momsum, dmsum);
    for (i=1; i<=sum_patch; i++)
       	{slip[i]=slip[i]*momsum/dmsum;
     	}
        ////////////////////////////////////////////////////////
       /// radiated seismic energy, corner frequency and rise time
        ////////////////////////////////////////////////////////
    es0=pow(10.0,(1.50*menergy+4.40));

    printf("tau=%lf\n", tau);

    ntst=nt*2;
    dfst=0.01;
    //dfst=20.0/(ntst/2);
    printf("ntst=%d, dfst=%lf\n", ntst, dfst);
    for(jfst=0; jfst<nt;jfst++)
      {//equally spaced bins on logarithmic axes
        fst[jfst]=pow(10,(((log10(maxfreq)-log10(1*dfst))*jfst/nt)+log10(1*dfst)));
      }
    //calculation of source function and differenciation
    //    (multiply 2PIf in frequency domain)

    scale=1000.0;
    while((scale>1.01 || scale<0.99) && numite < 150)
      { numite++;
        for(jfst=0;jfst<ntst/2;jfst++)
          { sfct[jfst]=0.0;
        	for (i=1; i<=sum_patch; i++)
        	  { 
                    t0=time1[i];
                    sfct[jfst]=sfct[jfst]+(slip[i])*cexp((-2*PI*fst[jfst]*t0-2.0*atan2(2*PI*fst[jfst]/fcp[i],1))*I)
                                 /(1.0+4*PI*PI*fst[jfst]*fst[jfst]/(fcp[i]*fcp[i]))*sqrt((1.0+2.0/3.0*pow(vss[i]/vps[i],5))
                                                 /(10.0*PI*rhos[i]*pow(vss[i],5)));

                   }
        	sfct[jfst]=sfct[jfst]*(2*PI*fst[jfst]*I);
          }
  /*      if(fhighcut>0.0)
          { jh=1+(int)(fhighcut/dfst);
            printf("fhighcut=%lf\n", fhighcut);
            for(jfst=jh; jfst<ntst/2; jfst++)
              { sfct[jfst]=sfct[jfst]*(pow((fhighcut/fst[jfst]),2)+0.0*I);
              }
          }
          */
        ///////     energy integration^M
        es=0.0;
        printf("bla\n");
        for(jfst=1; jfst<ntst/2; jfst++)
        	{ es=es+pow(cabs(sfct[jfst]),2)*(fst[jfst]-fst[jfst-1]);
        	}
        es=es*2.0;
        ////     compare with es0 (which was calculated from the energy magnitude)
        ////     and update patch corner frequency

        scale=pow((es0/es),(1.0/3.0));
        tau=tau/scale;
        for (i=1; i<=sum_patch; i++)
        	{fcp[i]=scale*fcp[i];
        	}
        printf("scale=%lf, es0=%lf, es=%lf, ME=%lf\n",scale, es0, es,menergy);
      }
    if (numite>=150)
      { printf("Too many iterations\n");
        exit(1);
      }
    for (i=1; i<=sum_patch; i++ )
      {
        rise[i]=1.0/(fcp[i]);
        risemean=risemean+rise[i];
      }
    risemean=risemean/(sum_patch-1);
    printf("risemean= %lf\n", risemean);
/////////////////////////////////////////////
////  source function output for visualisation
    sprintf(strrupvel,"%lf",rupvel);
    outfile= strcat(strrupvel, "_f_sfct.dat");
    fout = fopen (outfile, "w");
    fprintf (fout, "f[Hz]         Spec\n");
    for(jfst=0; jfst<ntst/2; jfst++)
    	{ sfct[jfst]=0.0;
    	  for (i=1; i<=sum_patch; i++)
    		  { 
                  sfct[jfst]=sfct[jfst]+(slip[i])*cexp((-2*PI*fst[jfst]*time1[i]-2.0*atan2(2*PI*fst[jfst]/fcp[i],1))*I)
                                 /(1.0+4*PI*PI*fst[jfst]*fst[jfst]/(fcp[i]*fcp[i]));
    		  }
    	  fprintf (fout,"%lf \t %lf\n",fst[jfst],cabs(sfct[jfst]));
    	}
  /*  if(fhighcut>0.0)
      { jh=1+(int)(fhighcut/dfst);
        printf("fhighcut=%lf\n", fhighcut);
        for(jfst=jh; jfst<ntst/2; jfst++)
          { sfct[jfst]=sfct[jfst]*(pow((fhighcut/fst[jfst]),2)+0.0*I);
          }
      }
*/
    fclose(fout);
    
    free_dvector (vps,1,sum_patch);
    free_dvector (vss,1,sum_patch);
    free_dvector (rhos,1,sum_patch);
    free_dvector (fcp,1,sum_patch);
  
  }


