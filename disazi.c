/*
 * disazi.c
 *
 * calculates horizontal distances between the reference point lateq/loneq and the station point latst/lonst
 *
 *  Created on: Jun 17, 2013
 *      Author: katrin
 *      translated from Fortran code (Rongjiang Wang)
 *
 *
 */

double MIN(double a, double b)
	{ return (a < b) ? a : b;
	}
double MAX(double a, double b)
	{ return (a > b) ? a : b;
	}

void
disazi( float rearth, double lateq, double loneq, double latst, double lonst, double *xnorth, double *yeast)
	{

        int     iangle;
        double  latb,lonb,latc,lonc,angleb,anglec;
        double  dis,a,b,c,s,aa;
        double  PI= 4.0 * atan (1.0);
        double  PI2 = 2*PI;
        double  DEGTORAD=1.745329252e-02;

//       A is north pole
        latb=lateq*DEGTORAD;
        lonb=loneq*DEGTORAD;
        latc=latst*DEGTORAD;
        lonc=lonst*DEGTORAD;


        if(lonb<0.0)
        	{lonb=lonb+PI2;
        	}

        if(lonc<0.0)
        	{lonc=lonc+PI2;
        	}
        b=0.5*PI-latb;
        c=0.5*PI-latc;
        if(lonc>lonb)
        	{ aa=lonc-lonb;
        	  if(aa<=PI)
                { iangle=1;
                }
              else
            	{aa=PI2-aa;
                 iangle=-1;
            	}
        	}
        else
            { aa=lonb-lonc;
              if(aa<=PI)
                {iangle=-1;
                }
              else
                {aa=PI2-aa;
                 iangle=1;
                }
            }
        s=cos(b)*cos(c)+sin(b)*sin(c)*cos(aa);
        a=acos(copysign(MIN(fabs(s),1.0),s));

        dis=a*rearth;
        printf("s=%lf, a=%lf ,term=%lf\n", s, a, fabs(s));
        if(a*b*c==0.0)
          {angleb=0.0;
           anglec=0.0;
          }
        else
          {s=0.5*(a+b+c);
           a=MIN(a,s);
           b=MIN(b,s);
           c=MIN(c,s);
           anglec=2.0*asin(MIN(1.0,sqrt(sin(s-a)*sin(s-b)/(sin(a)*sin(b)))));
           angleb=2.0*asin(MIN(1.0,sqrt(sin(s-a)*sin(s-c)/(sin(a)*sin(c)))));
           if(iangle==1)
             {angleb=PI2-angleb;
             }
           else
        	 {anglec=PI2-anglec;
        	 }
           }
//
//       cartesian coordinates of the station
//
        *xnorth=dis*cos(anglec);
        *yeast=dis*sin(anglec);

    }
