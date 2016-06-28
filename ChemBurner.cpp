
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
using namespace std;

#include "ChemBurner.H"

using namespace ChemNS;


void ChemIntegrate(REAL start,REAL firstDT,REAL tmin,REAL dt,REAL* eosIn, REAL* ymass, void (*f)(REAL,REAL*,REAL* ,REAL*, double* ),double* rpar,REAL* RatRaw,int &jcounts);

//ARS moved these declrations outside of Chem_Burner.  Unsure why, but declaring them inside Chem_Burner can
//lead to memory issues where dt is somehow overwritten to zero...  Stack overflow?

double ysave[NINT];
double ymass[NINT];
double ytot[NSPECIES];
double RatRaw[NREACTION];
double ddyddx[NINT];

void Chem_Burner(REAL dt, REAL* eosIn, REAL* yIn, REAL* yOut, double* rpar)
{

    // Integration/Debug stuff
    int i;
    double smallx = 1.0e-20, bigx = 0.5, bigx2 = 1.0;
    double temp_before;
    double ttot; ttot = dt; //ARS added "ttot" to save original value of dt
    REAL tmin  = dt*1.e-20, dtmax = 0.2*ttot; 


    //REAL ysave[NINT],ymass[NINT];
    for(i=0;i<NINT;i++)
      {
      ymass[i] = yIn[i];
      ysave[i] = yIn[i];
    }
    //REAL ytot[NSPECIES];
    for(int i=0;i<NSPECIES;i++) {
      ytot[i] = 0.0;
    }
  

    // Prepare the network rates (inputs: temp, density, array)
    //REAL RatRaw[NREACTION];

    GetNetworkRates(eosIn[1],eosIn[0],ymass, RatRaw, rpar);

 
    // Get the values of derivatives for the ODE solver
    //REAL ddyddx[NINT];
    GetNetwork(dt,ddyddx,RatRaw,ymass,rpar);

    // Given these, we can estimate the first timestep to try
    double firstDT; firstDT = dt;
    double tempDT = 0.0;
    int isSmall = -1;
    for(i=0;i<NINT;i++)
    {
	tempDT = fabs(max(ymass[i],1.0e-4)/ddyddx[i]);
        if(tempDT < firstDT) if(ymass[i] > smallx)
        {
            firstDT = tempDT;
            isSmall = i;
        }
    }
    if(firstDT > dtmax) firstDT = dtmax;
    // CTSS: should dt also be based on cooling time here?

    
    int icounts = 0, jcounts = 0;
    double advT = 0.0;
    while(dt>0) // Main integration loop
    {
        icounts = icounts + 1;
        
        // We have now set the time to integrate before updating photochem
        
        tmin = firstDT*1.0e-20;
        double endT = firstDT;
        
	rpar[temp_par] = eosIn[1];
        ChemIntegrate(0,firstDT,tmin,endT,eosIn,ymass,GetNetwork,rpar,RatRaw,jcounts);
	for(i=0;i<NINT;i++){
	  if(ymass[i] != ymass[i]){
	    cout << "ChemIntegrate produced nan!" << endl;
	    for(i=0;i<NINT;i++) cout << i << " " << ymass[i] << endl;
	    cout << "eos[0] = " << eosIn[0] << endl;
	    cout << "eos[1] = " << eosIn[1] << endl;
	    cout << "eos[2] = " << eosIn[2] << endl;
	    AbortChem();
	  }
	}
	
	
        // Prevent negative values
        for(i=0;i<NINT;i++) if(ymass[i] < smallx ) ymass[i] = smallx;

        //ARS also putting a check on the H2 and hP abundance 
        #if DO_CHEMISTRY == PRIMORDIAL
        if(ymass[iH2] > bigx) ymass[iH2] = bigx;
        if(ymass[iHP] > bigx2) ymass[iHP] = bigx2;        
        #endif

	//----------
	// CTSS: let's turn this off for now, or at least figure a better way to do it:
        // Determines how much ionization occured from collisions vs. photoionizations
        //double sdotRate=0.0;
        //double photoc[NSPECIES] = {0.0};
        //Chemistry_photo(ymass,RatRaw,photoc);
        //for(i=0;i<NINT;i++) sdotRate = sdotRate + (ymass[i]-ysave[i])*bion[i]*photoc[i];
        // Update energy density from this photochem / Update chemEOS from there
        //eosIn[2] = eosIn[2] + sdotRate*conv;
        //if(eosIn[2]<0) cout << "ERROR: energy density is negative!!!" << endl;
        //chemEOS(ymass,eosIn,aion,gamma,Do_cool);
	//----------
        
        // Cooling
	temp_before = eosIn[1];
	Cool_gas(eosIn,ymass,RatRaw,rpar,endT); //updates e (=eosIn[2])
	fill_species(ymass,ytot);
	Chem_EOS(2,eosIn,ytot,NSPECIES,6);
        rpar[temp_par] = eosIn[1];  //ARS updating temp value in rpar
	if(eosIn[1] < 0.1){
	  cout << "temp below 0.1 after cool species! not ok ... " << eosIn[1] << endl;
	  cout << "temp before cool_gas = " << temp_before << endl;
	  cout << "endT = " << endT << endl;
	  cout << "current abundances = " << endl;
	  for(int i=0;i<NINT;i++) cout << i << " " << ymass[i] << endl;
	  AbortChem();
	}
        

    
        // Re-evaluate rates and network using new state values and species abundances
        GetNetworkRates(eosIn[1],eosIn[0],ymass,RatRaw,rpar);
        GetNetwork(0.0,ddyddx,RatRaw,ymass,rpar); // first entry time, not actually used
        
        // We have now successfully advanced 'endT' worth of time
        // Update time left
        advT = advT + endT;
        dt = dt - endT;
        firstDT = dt;
        
        // Re-evaluate firstDT
        for(i=0;i<NINT;i++)
        {
	    tempDT = fabs(max(ymass[i],1.0e-4)/ddyddx[i]);
            if(tempDT < firstDT) if(ymass[i] > smallx)
            {
                firstDT = tempDT;
                isSmall = i;
            }
        }
        if(firstDT > dtmax) firstDT = dtmax;

        //ARS adding a different temp cooling constraint
        double Tdot = 0.01*fabs(temp_before - eosIn[1])/firstDT;
        double tkDT = eosIn[1] / Tdot;
        if(tkDT < firstDT) firstDT = tkDT;

        if(firstDT > dt) firstDT = dt; // should never happen, but why not check.
        
    }

 
    // All done!
    for(i=0;i<NINT;i++)
    {
      yOut[i] = ymass[i];
      if(yOut[i] < smallx) yOut[i] = smallx;
    }
 
  //ARS also putting a check on the H2 and HPabundance
  #if DO_CHEMISTRY == PRIMORDIAL
   if(yOut[iH2] > bigx) yOut[iH2] = bigx;
   if(yOut[iHP] > bigx2) yOut[iHP] = bigx2;
  #endif

// ARS adding a few print statements for debugging purposes
    if(jcounts > 5000)
      {
      printf("Took a little while!  jcounts = %d, nh = %lg, dt = %lg, advT = %lg, ttot = %lg\n", jcounts, eosIn[0]/1.67e-24, dt, advT, ttot);
      for(i=0;i<6;i++) printf("4 i = %d, eos = %20.15g\n", i, eosIn[i]);
      for(i=0;i<NSPECIES;i++) printf("4 i = %d, xIn = %lg, xOut = %20.15g\n", i, yIn[i], yOut[i]);
      }

    double random, rfac = 1.e-6;
    random = rand();
    if(random < rfac*RAND_MAX)
      {
      printf("spot check!  nh = %lg dt = %lg firstDT = %lg advT = %lg ttot = %lg\n", eosIn[0]/1.67e-24, dt, firstDT, advT, ttot);
      for(i=0;i<6;i++) printf("4 i = %d, eos = %20.15g\n", i, eosIn[i]);
      for(i=0;i<NSPECIES;i++) printf("4 i = %d, xIn = %lg, xOut = %lg\n", i, yIn[i], yOut[i]);
      for(i=0;i<NINT;i++) printf("4 i = %d, ddyddx = %lg\n", i, ddyddx[i]);
      }
//////////////////////////////////////////////////  
  
};







