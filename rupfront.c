/*
 * rupfront.c
 *
 * calculate rupture times for each patch of a seismic source
 *
 *  Created on: Jun 6, 2013
 *  last modified: Jun 25, 2013
 *      Author: katrin
 *
 *  required input: s_depth, vs -- parameters of earth model velocity (which S-velocity vs is assigned to which depth s_depth
 *                  lat, lon, z -- positions of source patches
 */
#include "disazi.c"
int MINI(int a, int b)
	{ return (a < b) ? a : b;
	}
int MAXI(int a, int b)
	{ return (a > b) ? a : b;
	}

void
rupfront (double **s_depth, double **vs, int num_model_depth, double *nx, double *ny, double *z, int nns, double dstep,
		  double lathypo, double lonhypo, double dephypo, double seglat, double seglon, double segdepth,
		  double dip, double rupvel, double *time1)
	{ 
          float REARTH=6.371e+3;
          double  PI= 4.0 * atan (1.0); 
          double xdis=0, ydis=0;
          double dissq, diffstart;
          double diseps;
          double dis,sx,sy,s,radius,vsp, dage, det;
          double xs, ys, dxs[9], dys[9];
          double dx0, dy0;//, rx, ry;
          float fac=1.0;
          int startx, starty, stopx, stopy;
          double time0;
          double **vrup, **set, **age;

          double agemax=0.0, agemin;
          double depthsum=0, s_depth_max, vrupmean=0.0;
          double rupvel2=0.6;

          int i,k,ns;
          int ixhypo, iyhypo, ix, iy,inx, iny,ixy;
          int ixmin, iymin, ixmax, iymax, ix0, iy0, ix0old, iy0old;
          int nxset, nyset, ixs[9], iys[9];

         //printf("number of sourcepoints: %d\n", nns);

	  disazi(REARTH*1000.0,seglat,seglon, lathypo,lonhypo,&xdis,&ydis);
	  dissq=(xdis*xdis)+(ydis*ydis);
	  dstep=dstep*1000.0;
      diseps=1.0e-06*dstep;
      vrupmean=0.0;
      depthsum=0;
      printf("s_depth[%d]=%lf\n",num_model_depth-2, (*s_depth)[num_model_depth-2]);
      printf("s_depth[%d]=%lf\n",num_model_depth-1, (*s_depth)[num_model_depth-1]);
      printf("s_depth[%d]=%lf\n",num_model_depth, (*s_depth)[num_model_depth]);
      s_depth_max=MIN((*s_depth)[num_model_depth],60000.0);
      for (k=1;k<num_model_depth;k++) //calculate mean rupture velocity for the segment
        { if ((*s_depth)[k]<s_depth_max)
            { if(k==1)
           		{ vrupmean=vrupmean+(*vs)[k]*(*s_depth)[k];
           		  depthsum=(*s_depth)[k];
           		}
           	  else
           		{ vrupmean=vrupmean+(*vs)[k]*((*s_depth)[k]-(*s_depth)[k-1]);
           		  depthsum=depthsum+(*s_depth)[k]-(*s_depth)[k-1];
           		}
           	}
        }
      vrupmean=vrupmean/depthsum*rupvel;
///////determine hypocenter location xdis, ydis
/////// adjust location if hypocenter lies outside the fault (e.g. for multiple fault segments)
      printf(" dip=%lf\n",dip);
          if (dip>0.0)
            {
              ydis=((dephypo-segdepth)*1000.0)/sin(dip*PI/180.0);
              printf(" dip=%lf\n",dip);
              if (pow(dissq-((dephypo-segdepth)*1000.0/tan(dip*PI/180.0)),2)>0.0)
                { xdis=sqrt(dissq-pow(((dephypo-segdepth)*1000.0/tan(dip*PI/180.0)),2));
                }
              else
                { xdis=0.0;
                }

            }
          else
	    { printf("*dip < 0.0 \n");
              exit(0);
            }

//          printf("lat lon x y %lf %lf %lf %lf \n",lathypo,lonhypo,xdis,ydis);
          printf("xdis=%lf, ydis=%lf, dstep=%lf\n", xdis, ydis, dstep);
          printf("nx=%lf, ny=%lf,\n", nx[nns], ny[nns]);
          if(xdis>nx[nns]*1000.0)
            { diffstart=xdis-nx[nns]*1000.0;
              time0=diffstart/vrupmean;
              xdis=nx[nns]*1000.0;
            }
          else if(xdis<dstep)
        	{ diffstart=dstep-xdis;
              time0=diffstart/vrupmean;
              xdis=dstep;
        	}
          else
            { time0=0;
              xdis=(int)(xdis/dstep)*dstep;
            }

          if(ydis>ny[nns]*1000.0)
        	{ diffstart=ydis-ny[nns]*1000.0;
              time0=time0+diffstart/vrupmean;
              ydis=ny[nns]*1000.0;
        	}
          else if(ydis<dstep)
        	{ diffstart=dstep-ydis;
             time0=time0+diffstart/vrupmean;
             ydis=dstep;
        	}
          else
            {  ydis=(int)(ydis/dstep)*dstep;
            }
 //         printf("xdis=%lf, ydis=%lf, nx=%lf, ny=%lf, dstep=%lf\n", xdis, ydis, nx[nns], ny[nns], dstep);
           printf("hypocenter %lf %lf %lf\n",xdis, ydis, time0);

           ixhypo=(int)(xdis/dstep);
           iyhypo=(int)(ydis/dstep);
           xs=0.5*nx[nns];
           ys=0.5*ny[nns];

           i=2;
           printf("nx(i)=%lf nx(i-1)=%lf\n",nx[i], nx[i-1]);
           while(nx[i]>nx[i-1])
        	   { i++;
        	   //  printf("nx(i)=%lf nx(i-1)=%lf\n",nx[i], nx[i-1]);
        	   }
           inx=i-1;
           iny=nns/inx;
           printf("iny %d, inx %d , nns %d\n",inx, iny, nns);
/////     fix the velocity field
/////     velocity is set to rupvel times the local S-velocity
          vrup = (double**)calloc(inx,sizeof(double*));
          if ( vrup== NULL)
            { printf("stop: could not allocate vrup\n");
              exit(1);
            }
          for(i=0;i<inx;i++)
            { vrup[i] = (double*)calloc(iny,sizeof(double));
              if ( vrup[i]== NULL)
                { printf("stop: could not allocate vrup\n");
                  exit(1);
                }
            }
          set = (double**)calloc(inx,sizeof(double*));
          if ( vrup== NULL)
            { printf("stop: could not allocate set\n");
              exit(1);
            }
          for(i=0;i<inx;i++)
            { set[i] = (double*)calloc(iny,sizeof(double));
              if ( set[i]== NULL)
                { printf("stop: could not allocate vrup\n");
                  exit(1);
                }
            }
          age = (double**)calloc(inx,sizeof(double*));
          if ( vrup== NULL)
            { printf("stop: could not allocate age\n");
              exit(1);
            }
          for(i=0;i<inx;i++)
            { age[i] = (double*)calloc(iny,sizeof(double));
              if ( age[i]== NULL)
                { printf("stop: could not allocate age\n");
                  exit(1);
                }
            }
           for (i=1;i<=nns;i++)
           	 { ix=(int) round(nx[i]*1000.0/dstep);
           	   iy=(int) round(ny[i]*1000.0/dstep);
           	   k=1;
           	   //printf("iy %d, ix %d  %d\n",iy, ix,i);
           	   while(vrup[ix][iy]==0.0)
           		   { if(z[i]*1000.0<=(*s_depth)[k])
           		       { vrup[ix][iy]= rupvel*(*vs)[k];
           		         k=1;
           		       }
           	             else if (z[i]*1000.0<(*s_depth)[num_model_depth])
           	               { k++;
           		       }
           	             else
           	               { vrup[ix][iy]= rupvel*(*vs)[num_model_depth];
           	        	 k=1;
           	               }
                           }
                   if (z[i]<=5.0) 
                     { vrup[ix][iy]= vrup[ix][iy]*rupvel2;
                     }
                   else if (z[i]>5.0 && z[i]<8.0)
                     { vrup[ix][iy]= vrup[ix][iy]*((1-rupvel2)/3.*z[i]+1-8./3.*(1-rupvel2));
                       fac=((1-rupvel2)/3.*z[i]+1-8./3.*(1-rupvel2));
                     }
           	 }

           for(ix=ixhypo-1; ix<ixhypo+1; ix++)
        	 { for(iy=iyhypo-1; iy<iyhypo+1; iy++)
                 { dis=sqrt(pow(((ix-ixhypo)*dstep),2)+pow(((iy-iyhypo)*dstep),2));
                   age[ix-1][iy-1]=2.0*dis/(vrup[ix-1][iy-1]+vrup[ixhypo-1][iyhypo-1]);
                   set[ix-1][iy-1]=1;
                   agemax=MAX(agemax,age[ix-1][iy-1]);
                 }
        	 }
           ixmin=ixhypo-1;
           iymin=iyhypo-1;
           ixmax=ixhypo+1;
           iymax=iyhypo+1;
           printf("min und max: %d / %d, %d / %d\n", ixmin, iymin, ixmax, iymax);
/////////////////////////////////////////////////7           100
           agemin=2.0*agemax;
           ix0=0;
           iy0=0;

           for(ix=ixmin; ix<ixmax; ix++)
        	   { for(iy=iymin; iy<iymax; iy++)
        		   { if(set[ix-1][iy-1]==1)
        			   { nxset=set[ix-1-1][iy-1]+set[ix+1-1][iy-1];
                         nyset=set[ix-1][iy-1-1]+set[ix-1][iy+1-1];
                         if (ix==inx || ix==1)
                        	 { nxset=nxset+1;
                        	 }
                         if (iy==iny || iy==1)
                        	 { nyset=nyset+1;
                        	 }
                         if (nxset>0 && nyset>0 && nxset+nyset<4 && agemin>=age[ix-1][iy-1])
                           { agemin=age[ix-1][iy-1];
                             ix0=ix;
                             iy0=iy;
                           }
                       }
        		   }
        	   }
           while (ix0>0 || iy0>0)
             { //printf("agemin=%lf, ix=%d, iy=%d\n", agemin, ix, iy);
          /////////       the oldest edge patch is found: (ix0,iy0)
                 ns=0;
                 if(set[ix0-1-1][iy0-1]==1)
                   { ns=ns+1;
                     ixs[ns]=ix0-1;
                     iys[ns]=iy0;
                     dxs[ns]=-dstep;
                     dys[ns]=0.0;
                     if(set[ix0+1-1][iy0-1]==0)
                       { sx=(age[ix0-1][iy0-1]-age[ix0-1-1][iy0-1])/dstep;
                       }
                     else
                       { ns=ns+1;
                         ixs[ns]=ix0+1;
                         iys[ns]=iy0;
                         dxs[ns]=dstep;
                         dys[ns]=0.0;
                         sx=0.5*(age[ix0+1-1][iy0-1]-age[ix0-1-1][iy0-1])/dstep;
                       }
                   }
                 else
                   { ns=ns+1;
                     ixs[ns]=ix0+1;
                     iys[ns]=iy0;
                     dxs[ns]=dstep;
                     dys[ns]=0.0;
                     sx=(age[ix0+1-1][iy0-1]-age[ix0-1][iy0-1])/dstep;
                   }

                 if(set[ix0-1][iy0-1-1]==1)
                   { ns=ns+1;
                     ixs[ns]=ix0;
                     iys[ns]=iy0-1;
                     dxs[ns]=0.0;
                     dys[ns]=-dstep;
                   if(set[ix0-1][iy0-1+1]==0)
                     { sy=(age[ix0-1][iy0-1]-age[ix0-1][iy0-1-1])/dstep;
                     }
                   else
                     { ns=ns+1;
                       ixs[ns]=ix0;
                       iys[ns]=iy0+1;
                       dxs[ns]=0.0;
                       dys[ns]=dstep;
                       sy=0.5*(age[ix0-1][iy0+1-1]-age[ix0-1][iy0-1-1])/dstep;
                     }
                   }
                 else
                   { ns=ns+1;
                     ixs[ns]=ix0;
                     iys[ns]=iy0+1;
                     dxs[ns]=0.0;
                     dys[ns]=dstep;
                     sy=(age[ix0-1][iy0+1-1]-age[ix0-1][iy0-1])/dstep;
                   }

               s=sqrt(sx*sx+sy*sy);
               sx=sx/s;
               sy=sy/s;
   //            printf("s, sx, sy: %lf, %lf, %lf\n",s, sx, sy);
               radius=0.0;
               for( i=1; i<ns; i++)
                 { vsp=0.5*(vrup[ix0-1][iy0-1]+vrup[ixs[i]-1][iys[i]-1]);
                   dage=age[ixs[i]-1][iys[i]-1]-age[ix0-1][iy0-1];
                   det=2.0*(vsp*dage-dxs[i]*sx-dys[i]*sy);
                   if(abs(det)>diseps)
                     { radius=MAX(radius,
                        (dxs[i]*dxs[i]+dys[i]*dys[i]-(vsp*dage)*(vsp*dage))/det);
                     }
                 }
               xs=radius*sx;
               ys=radius*sy;
               startx=MAXI(1,ix0-1);
               starty=MAXI(1,iy0-1);
               stopx=MINI(ix0+1,inx);
               stopy=MINI(iy0+1,iny);

               for(ix=startx; ix<=stopx; ix++)
                 { dx0=(ix-ix0)*dstep;
                   for(iy=starty; iy<=stopy; iy++)
                     { dy0=(iy-iy0)*dstep;
                       if(set[ix-1][iy-1]==0)
                         { age[ix-1][iy-1]=age[ix0-1][iy0-1]+2.0*(sqrt((xs+dx0)*(xs+dx0)+(ys+dy0)*(ys+dy0))-radius)
                                       /(vrup[ix-1][iy-1]+vrup[ix0-1][iy0-1]);
                   //        printf("age=%lf,ix=%d,iy=%d, vel=%lf, vel(ix0,iy0)=%lf\n", age[ix-1][iy-1], ix, iy, vrup[ix-1][iy-1], vrup[ix0-1][iy0-1]);
                           agemax=MAX(agemax,age[ix-1][iy-1]);
                           set[ix-1][iy-1]=1;
                         }
                     }
                 }

               ixmin=MINI(ixmin,MAXI(1,ix0-1));
               iymin=MINI(iymin,MAXI(1,iy0-1));
               ixmax=MAXI(ixmax,MINI(inx-1,ix0+1));
               iymax=MAXI(iymax,MINI(iny-1,iy0+1));
  //         printf("min und max: %d / %d, %d / %d\n", ixmin, iymin, ixmax, iymax);
               
               agemin=2.0*agemax;
               ix0old=ix0;
               iy0old=iy0;
               ix0=0;
               iy0=0;

               for(ix=ixmin; ix<=ixmax; ix++)
                   { for(iy=iymin; iy<=iymax; iy++)
                           { if(set[ix-1][iy-1]==1)
                              { // printf("agemin=%lf\n",agemin);
                                if (ix==inx)
                                    { nxset=set[ix-1-1][iy-1]+1;
                                    }
                                else if (ix==1)
                                    { nxset=set[ix+1-1][iy-1]+1;
                                    }
                                else
                                    {nxset=set[ix-1-1][iy-1]+set[ix+1-1][iy-1];
                                    }
                                if (iy==iny)
                                    { nyset=set[ix-1][iy-1-1]+1;
                                    }
                                else if (iy==1)
                                    { nyset=set[ix-1][iy+1-1]+1;
                                    }
                                else
                                    { nyset=set[ix-1][iy-1-1]+set[ix-1][iy+1-1];
                                    }
                                if (ix==inx || ix==1)
                               	 { nxset=nxset+1;
                               	 }
                                if (iy==iny || iy==1)
                               	 { nyset=nyset+1;
                               	 }
                                if (nxset>0 && nyset>0 && nxset+nyset<4 && agemin>=age[ix-1][iy-1] && (ix!=ix0old || iy!=iy0old))
                                  { agemin=age[ix-1][iy-1];
                                    ix0=ix;
                                    iy0=iy;

                                  }
                              }
                           }
                   }

             }
//

           ///// add age to start time (time0) to determine absolute rupture time (time1)
           ixy=1;
           for(iy=1;iy<=iny;iy++)
             { //ry=(iy-1)*dstep;
               for(ix=1; ix<=inx; ix++)
                 {//rx=(ix-1)*dstep;
                   time1[ixy]=age[ix-1][iy-1]+time0;
      //             printf("%lf\t %lf\t %lf\t %lf\t %lf\n", rx, ry, nx[ixy], ny[ixy], time1[ixy]);
                   ixy=ixy+1;

                 }
             }
      }
