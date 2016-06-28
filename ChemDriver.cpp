
#include <iostream>
#include <iomanip>
using namespace std;

#include "ChemDriver.H"
using namespace ChemNS;


void fill_species(REAL* ymass, REAL* ytot);

void Chem_Burner(REAL dt, REAL* eosIn, REAL* xIn, REAL* xOut, double* rpar);

void GetNetworkRates(const REAL Temp, const REAL Den,REAL* ys, REAL* RatRaw, double* rpar);


REAL yin_original[NINT];
REAL rpar_original[nrpar];
REAL ei_original;
REAL dt_original;
REAL dens_original;

void ChemDriver(const double &dt_in, const double &den_in, double &ei_in, double &gamma_here, double* y_fromOrion, double* rpar)
{
    int i;
    double smallx = 1.0e-20;
    
    // Initialize arrays (xIN will be read in, actually)
    REAL xIn[NINT] = {0.0};
    REAL yIn[NINT] = {0.0};
    REAL yOut[NINT] = {0.0};
    REAL eosIn[6];  // Density, Temperature, EngDen, Pressure, Metal Fraction, gamma

    // save original values for use in AbortChem
    for(i=0;i<NINT;i++) yin_original[i] = y_fromOrion[i];
    for(i=0;i<nrpar;i++) rpar_original[i] = rpar[i];
    ei_original = ei_in;
    dt_original = dt_in;
    dens_original = den_in;
    //

    
    // Initialization parameters
    // and other parameters that might become read in
    // CTSS: deprecated, J21 is now in rpar array
    // and shielding factors should be determined consistently
    //bool chem_cc_case = 1;
    //REAL fssh2 = 1.0;
    //REAL fsshd = 1.0;
    //REAL J21 = 0.0;  // J21 radiation
    //REAL RadParams[3] = {J21,fssh2,fsshd};
    
    
    // Load from Orion2
    for(i=0;i<NINT;i++) xIn[i] = y_fromOrion[i];
    REAL dt = dt_in;
    eosIn[0] = den_in;
    eosIn[1] = 0.0; //Temp
    eosIn[2] = ei_in;
    eosIn[3] = 0.0; //P
    eosIn[4] = 0.0; // mass frac
    eosIn[5] = gamma_here;
    rpar[gamma_par] = gamma_here;
    rpar[rho_par] = den_in;
    rpar[e0_par] = ei_in;
    rpar[temp_par] = -1.0;

    // check initial abundances
    for(i=0;i<NINT;i++)
    {
      if(xIn[i] < -1.0){
	cout << "big negative abundance upon entering chemistry! " << i << " " << xIn[i] << endl;
	AbortChem();
      }
    }
    if(xIn[iH] > 1.0){
      if(xIn[iH] > 2.0){
	cout << "too big H abundance too high upon entering chemistry" << endl;
	AbortChem();
      }
      xIn[iH] = 1.0 - xIn[iHP] - 2.0 * xIn[iH2];
    }
    if(xIn[iH2] > 0.5){
      if(xIn[iH2] > 1.0){
	cout << "too H2 abundance too high upon entering chemistry" << endl;
	AbortChem();
      }
      xIn[iH2] = 0.5 * (1.0 - xIn[iH] - xIn[iHP]);
    }
    
      
    

 
    // yIn = abundance array of tracked species (size NINT)
    // ytot = abundance array of ALL species (including electrons, size NINT+1)
    REAL ytot[NSPECIES];
    for(i=0;i<NINT;i++)
    {
      yIn[i] = xIn[i] ; 
	ytot[i] = yIn[i];
        if(yIn[i]< smallx ) yIn[i] = smallx;
    }
 

    fill_species(yIn,ytot);
    Chem_EOS(2,eosIn,ytot,NSPECIES,1);
    rpar[temp_par] = eosIn[1];
    rpar_original[temp_par] = rpar[temp_par];
    
 
    // Call BURNER
    Chem_Burner(dt,eosIn,yIn,yOut,rpar);
    
 
    fill_species(yIn,ytot);
    Chem_EOS(1,eosIn,ytot,NSPECIES,2);
        
    for(i=0;i<NINT;i++) y_fromOrion[i] = yOut[i]; 
    //den_in = eosIn[0];
    ei_in = eosIn[2];    
    

    

// check final abundances
    for(i=0;i<NINT;i++)
    {
      if(yOut[i] < 0.0){
	cout << "negative abundance upon leaving chemistry! " << i << " " << xIn[i] << endl;
	AbortChem();
      }
      if(yOut[iH] > 1.5){
	cout << "H abundance too high upon leaving chemistry" << endl;
	cout << "current abundances:" << endl;
	for(i=0;i<NINT;i++){
	  cout << i << " " << yOut[i] << endl;
	}
	AbortChem();
      }
      if(yOut[iH2] > 1.0){
	cout << "H2 abundance too high upon leaving chemistry" << endl;
	AbortChem();
      }
    }
    
    
    if(eosIn[1] > 1.0e9){
      cout << "leaving chemistry with temp too high, > 1.0e9. Temp = " << eosIn[1] << endl;
      AbortChem();
    }
    

}




void AbortChem()
{ 
  int i;
  cout << "Aborting chemistry!" << endl;
  cout << "original rpar: " << endl;
  for(i=0;i<nrpar;i++) cout << i << " " << rpar_original[i] << endl;
  cout << "original abundances:" << endl;
  for(i=0;i<NINT;i++) cout << i << " " << yin_original[i] << endl;
  cout << "original ei = " << ei_original << endl;
  cout << "original dt = " << dt_original << endl;
  cout << "original density = " << dens_original << endl;
  cout << "goodbye..."<<endl; 
  exit(0);
}
