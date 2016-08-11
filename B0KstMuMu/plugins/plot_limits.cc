#include <iostream>
#include "TMath.h"

using namespace std;

void plot_limits (int bin)
{
  double flarr[9] = {0.641004,0.799186,0.619384,0.503676,0,0.392124,0,0.476826,0.377081};

  double fl = flarr[bin];
  double ft = 1-fl;
  
  for (double iP1 = -1.; iP1<1.; iP1+=0.01) {
    for (double iP5 = -1.5; iP5<0; iP5+=0.01) {
      bool out = false;   
      for (double ctk = -1.; ctk<=1; ctk+=0.02) {
	for (double ctl =0.; ctl<=1; ctl+=0.02) {
	  for (double phi =0; phi<TMath::Pi(); phi+=0.02) if (4*fl*ctk*ctk*(1-ctl*ctl) + ft*(1-ctk*ctk)*(1+ctl*ctl) +
							      iP1*ft*(1-ctk*ctk)*(1-ctl*ctl)*cos(2*phi) +
							      4*iP5*ctk*cos(phi)*sqrt(fl*ft*(1-ctk*ctk)*(1-ctl*ctl)) < 0) {
	      out = true;
	      break;
	    }
	  if (out) break;
	}
	if (out) break;
      }
      if (!out) {
	cout<<iP5-0.01<<" "<<iP1<<endl;
	break;
      }
    }
  }
  return;
}
