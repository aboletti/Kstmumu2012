#include <iostream>
#include "TMath.h"

using namespace std;

double Pwave (double fl, double ft, double P1, double P5, double ctl, double ctk, double phi)
{
  return 4*fl*ctk*ctk*(1-ctl*ctl) + ft*(1-ctk*ctk)*(1+ctl*ctl) +
    P1*ft*(1-ctk*ctk)*(1-ctl*ctl)*cos(2*phi) +
    4*P5*ctk*cos(phi)*sqrt(fl*ft*(1-ctk*ctk)*(1-ctl*ctl));
}

double Swave (double fs, double as, double A5, double ctl, double ctk, double phi)
{
  return 4./3.*((fs+as*ctk)*(1-ctl*ctl) + 
		A5*sqrt((1-ctk*ctk)*(1-ctl*ctl))*cos(phi));
}

void plot_limits (int bin, int effe)
{
  double flarr[9] = {0.641004  ,0.799186 ,0.619384 ,0.503676 ,0,0.392124   ,0,0.476826 ,0.377081  };
  double fsarr[9] = {0.00196504,0.0101492,0.0113141,0.0227136,0,8.33273e-05,0,0.0108452,0.00948832};
  double asarr[9] = {0.0883043 ,-0.276216,-0.255672,-0.32599 ,0,0.0176166  ,0,0.220459 ,-0.178423 };

  double fl = flarr[bin];
  double ft = 1-fl;
  double fs = fsarr[bin];
  double as = asarr[bin];

  int max_effe = 3;

  // double iP1 = -1 + 0.01*P1indx;
  // cout<<"- "<<iP1<<endl;
  for (double iP1 = -1.; iP1<1.; iP1+=0.01) {
    for (double iP5 = -1.5; iP5<0; iP5+=0.01) {
      // cout<<"* "<<iP5<<endl;
      bool out = false;
      double A5smax = 0.89 * sqrt(3*fs*(1-fs)*ft*(1+iP1));
      double iA5 = effe * A5smax / max_effe;
      // for (double iA5 = -1*A5smax; iA5<=A5smax; iA5+=(A5smax>0?A5smax/2:1)) {
      // 	out = false;
	// cout<<"+ "<<iA5<<endl;
	for (double ctk = -1.; ctk<=1; ctk+=0.02) {
	  for (double ctl =0.; ctl<=1; ctl+=0.02) {
	    for (double phi =0; phi<TMath::Pi(); phi+=0.02) if ( (1-fs)*Pwave(fl,ft,iP1,iP5,ctl,ctk,phi)+Swave(fs,as,iA5,ctl,ctk,phi) < 0) {
		out = true;
		break;
	      }
	    if (out) break;
	  }
	  if (out) break;
	}
      // 	if (out) break;
      // }
      if (!out) {
	cout<<iP5-0.01<<" "<<iP1<<endl;
	break;
      }
    }
  }
  
  return;
}
