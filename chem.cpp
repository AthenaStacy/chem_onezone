
#include <math.h> 
#include <iostream> 
#include <fstream>
#include <string>
#include "ChemGlobals.H"
#include "ChemExternal.H"

using namespace std;
using namespace ChemNS;
		   
void ChemDriver(const double&, const double&, double&, double&, double*, double*);

FILE *outfile;
	
int main()
{

 outfile=fopen("output.dat","a");

  cout << "Starting chemistry" <<  endl;
  
  REAL xhe = 0.08;
  REAL gamma = 5.0/3.0;
  
  REAL dt = 1.0e15;
  REAL rho = 1.0e-16;
  REAL temperature = 1.9e4;//temperature = 950.0;
  REAL rho_0 = rho;

  REAL ei, nh, ntot;
  REAL rpar[nrpar];
  
  // starting mass fractions;
  REAL yIn[NINT], ysav[NINT];
  
  for(int i=0;i<nrpar;i++) rpar[i] = 0.0;

#if DO_CHEMISTRY == PRIMORDIAL

  double yhp, yh, yhm, yh2p, yh2, ydp, yd, ydm, yhdp, yhd, yd2p, yd2, yhep, yhe, yhepp;

/*
    yhp = 3.64e-08; 
    yhm = 2.22e-14; 
    yh2p = 5.03e-16; 
    yh2 = 1.e-3; 
    ydp = 1.78e-12; 
    yd = DeutAbund; 
    ydm = 1.e-40; 
    yhdp = 1.e-40; 
    yhd = 1.e-7; 
    yd2p = 1.e-40; 
    yd2 = 1.e-40; 
    yhep = 1.e-20; 
    yhe = 0.079; 
    yhepp = 1.e-20; 
*/

/*
    yhp = 1.e-9;//4.7e-11; 
    yhm = 5.855e-17; 
    yh2p = 1.3e-12; 
    yh2 = 0.1;//.498; 
    ydp = 1.3e-15; 
    yd = 3.74e-8; 
    ydm = 1.e-20; 
    yhdp = 1.5e-16; 
    yhd = 2.6e-5; 
    yd2p = 1.e-20; 
    yd2 = 3.5e-11; 
    yhep = 1.e-20; 
    yhe = 0.079; 
    yhepp = 1.e-20; 
*/
/*
    yhp = 1.; 
    yhm = 2.89678e-19; 
    yh2p = 2.3291e-07; 
    yh2 =3.08887e-09; 
    ydp = 9.90603e-05; 
    yd = 1e-20; 
    ydm = 1.e-20; 
    yhdp = 1.e-20; 
    yhd = 8.53085e-05; 
    yd2p = 7.88495e-06; 
    yd2 = 0.0244519; 
    yhep = 0.00108467; 
    yhe = 4.47086e-07; 
    yhepp = 0.0778623; 
*/
/*
    yhp = 0.099467; 
    yhm = 2.51847e-08; 
    yh2p = 1.1045e-07; 
    yh2 = 1.04949e-06; 
    ydp = 2.62282e-06; 
    yd = 2.33771e-05; 
    ydm = 6.70195e-13; 
    yhdp = 1.06685e-11; 
    yhd = 1.97779e-11; 
    yd2p = 1.56767e-16; 
    yd2 = 4.69013e-16; 
    yhep = 1.83336e-05; 
    yhe = 0.078929; 
    yhepp = 1.25586e-14; 
*/
    yhp = 7.e-10; 
    yhm = 9e-16; 
    yh2p = 1.e-15; 
    yh2 = 8e-3; 
    ydp = 2.e-14; 
    yd = 2.5e-5; 
    ydm = 2.3e-20; 
    yhdp = 1.6e-19; 
    yhd = 8.4e-7; 
    yd2p = 1.e-20; 
    yd2 = 4.4e-12; 
    yhep = 1.25e-20; 
    yhe = 0.078929; 
    yhepp = 1e-20; 

/*
  yhp = 1.0e-4;
  yhm = 1.0e-15;
  yh2p = 1.0e-15;
  yh2 = 1.0e-6;
  ydp = 1.0e-8;
  yd = 4.0e-5;
  ydm = 1.0e-15;
  yhdp = 1.0e-15;
  yhd = 1.0e-10;
  yh2p = 1.0e-15;
  yd2 = 1.0e-15;
  yhep = 1.0e-10;
  yhe = 0.08;
  yhepp = 1.0e-15;
*/

  yh = 1.0 - yhp - 2.0 * yh2;
  
  yIn[iHP] = yhp ;
  yIn[iH] = yh;
  yIn[iHM] = yhm;
  yIn[iH2P] = yh2p;
  yIn[iH2] = yh2;
  yIn[iDP] = ydp;
  yIn[iD] = yd;
  yIn[iDM] = ydm;
  yIn[iHDP] = yhdp;
  yIn[iHD] = yhd;
  yIn[iD2P] = yd2p;
  yIn[iD2] = yd2;
  yIn[iHEP] = yhep;
  yIn[iHE] = yhe;
  yIn[iHEPP] = yhepp;
  
#endif

#if DO_CHEMISTRY == CONTEMPORARY

 double yh, yhm, yh2, yhp, yh2p, yh3p, yhe, yhep, yc, ycp, yo, ychx, yohx, yhcop, yco, ym, ymp;
  yhm = 9.5e-13;
  yh2 = 1.0e-6;
  yhp = 7.5e-5;
  yh2p = 2.92e-13;
  yh3p = 5.92e-14;
  yhe = 0.08;
  yhep = 1.46e-5;
  yc = 1.6e-4;
  ycp = 1.0e-7;
  yo = 3.2e-4;
  ychx = 1.3e-12;
  yohx = 4.3e-16;
  yhcop = 4.56e-19;
  yco = 6.0e-10;
  ym = 3.9e-10;
  ymp = 1.99e-7;
  


  yh = 1.0 - yhp - 2.0 * yh2;


  yIn[iH] = yh;
 yIn[iHM] = yhm ;
 yIn[iHP] = yhp;
 yIn[iH2] = yh2;
 yIn[iH2P] = yh2p ;
 yIn[iH3P] = yh3p;
 yIn[iHE] = yhe;
 yIn[iHEP] = yhep;
 yIn[iC] = yc;
 yIn[iCP] = ycp;
 yIn[iO] = yo;
 yIn[iCHX] = ychx;
 yIn[iOHX] = yohx;
 yIn[iCO] = yco;
 yIn[iHCOP] = yhcop;
 yIn[iM] = ym;
 yIn[iMP] = ymp;


#endif

 
  rpar[g0_par] = 0.36;
  rpar[tdust_par] = 10.0;
  rpar[NH_par] = 1.0e10;
  rpar[NH2_par] = 1.0e10;
  rpar[NCO_par] = 1.0e10;
  rpar[zeta_par] = 1.0e-17;
  rpar[divv_par] = 1.3506e-13;
  
  //rpar[divv1_par] = 1.e-10; 
  //rpar[dx_par] = 5.e13;  
 
 rpar[divv1_par] = 3.e-8; 
 rpar[dx_par] = 5.e13; 

  double mu; 
  mu = (1. - yh2*h2conv - yhe*heconv) + yh2*h2conv/2 + yhe*heconv/4;
  mu = 1./mu; 
  if(yh2 > 0.2) gamma = 1.4;

  nh = rho / mh / (1.0 + 4.0*xhe);
  ntot = nh * (xhe + 1.0);
  //ei = temperature / (gamma - 1.0) * ntot * kboltz / rho / 1.3 * (1.0 + 0.08 - 0.5);
  ei = temperature / (gamma - 1.0) * ntot * kboltz / rho / mu;
  cout << "ei = " << ei << endl;
  
  cout << "starting T = " << temperature << " mu = " << mu << endl;

  REAL ei0;
  REAL step;
  REAL tff = 0.000476;
  REAL t=1.0;
  REAL duration = 3.e7; //duration=3.0e50; 
  REAL t_start;
  REAL rho_dot;
  REAL rho_prev;
  double gamma_here;
  gamma_here = 1.6666667;
  ofstream dump;
  string fout = "evol.dat";
  dump.open((fout).c_str());
  
  step = 5.0e10;
  rho = 1.0e-24;
  ei0 = ei;

  
  cout << "beginning abundances:" << endl;
  for(int i=0;i<NINT;i++){ 
    ysav[i] = yIn[i];
    cout << i << " " << yIn[i] << endl;
  }


  // for equilibrium temperature as function of density
  /*
  ofstream dump1;
  string fout1 = "temp_eq_simp.dat";
  dump1.open((fout1).c_str());

  int ndens = 200;
  double rhol = 2.0e-27;
  double rhoh = 3.0e-18;
  duration = 1.0e17;
  

  for(int i=0;i<ndens;i++) {

    //for(int s=0;s<NINT;s++) yIn[s] = ysav[s];

    rho = rhol * pow(10.0, log10(rhoh/rhol) * double(i) / double(ndens-1)   );
    cout << " rho = " << rho  <<endl;
    
    ei = 3.0e9;
    step = 1.0;
    t = 1.0;
    while(t<duration){
      step = step * 1.1;
      ChemDriver(step,rho,ei,gamma_here,yIn,rpar);  
      t = t + step;
    }
    
    temperature = ei * (gamma_here - 1.0) /  kb * mh * 1.3 / (1.0 + 0.08 - yIn[iH2]);
    dump1 << rho << " " << temperature << " " << yIn[iCO] << " " << yIn[iC] << " " << yIn[iCP] << " " << yIn[iCHX] << " " << yIn[iHCOP] << " " << yIn[iOHX] << " " <<  yIn[iO] << " " << yIn[iH3P] << endl;
    
  }
  
  dump1.close(); 
  return 0;
  */

  //
  
  //rho = 1.0e-16;
  //rho = 2.7e-10;
  //rho = 7.e-14;
  //rho = 6.e-17;
  rho = 4.e-15;
  rho_0 = rho;
  rho_prev = rho_0;
  step = 1.0;
  t = 1.0;

  while(t < duration) // Interates until dt is advanced
    {
      //rho = pow(pow(rho_0,-0.5) - tff * t / 2.0 ,-2.0);  //keep denisty constant to mimic Orion2
      rho_dot = tff * pow(rho,1.5);
      //rho = rho_0;
      //step = step * 1.1;
      step = 0.1 * rho / rho_dot;

      //if(yIn[iH2] > 0.2) gamma = 1.4;
      //else 
      gamma = 5./3.;     
 
      // adiabatically heat gas here
      ei = ei * pow(rho/rho_prev,gamma-1);
      cout << "input ei = " << ei << endl;

      rho_prev = rho;
     
       mu = (1. - yIn[iH2]*h2conv - yIn[iHE]*heconv) + yIn[iH2]*h2conv/2 + yIn[iHE]*heconv/4;
       mu = 1./mu;
 
      ChemDriver(step,rho,ei,gamma,yIn,rpar);
      //temperature = ei * (gamma - 1.0) /  kboltz * mh * 1.3 / (1.0 + 0.08 - 0.5);
      temperature = ei * (gamma - 1.0) /  kboltz * mh * mu;

      double heat_fac=1.0;
      //if(rho > 1.e-11 && temperature < 2500.) heat_fac = 2500./temperature;
      //ei = ei * heat_fac;
      //temperature = temperature * heat_fac;

      nh = rho / mh / (1.0 + 4.0*xhe); double dum=0.;
      fprintf(outfile, "%15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g %15.11g\n", nh, yIn[iH2], yIn[iHD], yIn[iDP], yIn[iHEP], yIn[iHP], dum, dum, tff, dum, temperature, dum, dum, dum, dum, dum);

      cout << "t = " << t << " rho = " << rho << ", temp = " << temperature << " h2 =" << yIn[iH2] << endl;
      dump << rho/1.67e-24 << " " << t << " " << temperature << " " << yIn[iH2] /2.0 * 1.3 << " " << yIn[iHP] << " " << yIn[iH] <<  endl;
      t = t + step;
      
      if(rho > 1.0e-8) break;


    }
  dump.close(); 
  cout << "ending mass fractions = " << endl;
  for(int i=0;i<NINT;i++) cout << i << " " << yIn[i] << endl;
  cout << " " << endl;
  cout << "ending T = " << ei * (gamma - 1.0) /  kboltz * mh * mu << endl;
  cout << "ending time = " << t << endl;
  
  fclose(outfile);

  return 0;
		 
		     
}
