/* Chemistry Routine
 Based on Fortran code from Will Gray, re-written for C++ by Aaron Lee, 2014
 Heavily modified by Chalence Safranek-Shrader, 2016
 */

#include "ChemGlobals.H"
#include "ChemExternal.H"

using namespace ChemNS;

extern REAL yin_original[NINT];
extern REAL rpar_original[nrpar];
extern REAL ei_original;
extern REAL dt_original;
extern REAL dens_original;

void AbortChem();


void Chem_EOS(int mode, double* eosIn, double* y, int nSpec, int loc)
{
    //double mh = 1.67e-24;
    //double kb = 1.38e-16;
    double yhetot, fac;
    
    double den,temp,ei,pres,gc;
    den=eosIn[0];
    temp=eosIn[1];
    ei=eosIn[2];
    pres=eosIn[3];
    gc=eosIn[5];
    
    // CTSS: P = ntot * kb * T = rho * kb * T / (mu * mh)
    //       P = E * (gamma-1), where E = volumetric energy density
    //       P = ei / rho * (gamma-1), where ei = specific energy density
    //       To a very good approximation, ntot = (1 + yhe + ye - yh2) * nh
    //       where nh is the number density of hydrogen nuclei given by
    //       nh = rho / (mh * (1.0 + 4*yhe)) 
    //       ...which basically describes what's being done below:

    yhetot = y[iHE] + y[iHEP];
    fac = (1.0 + yhetot + y[iELEC] - y[iH2]) * kboltz / ((1.0 + 4.0*yhetot) * mh * (gc-1.0));
    if(mode==1) // determines ei and pressure from rho and T
    {
      ei = fac * temp;
      pres = (gc-1.0) * den * ei;
      
      eosIn[0]=den;
      eosIn[1]=temp;
      eosIn[2]=ei;
      eosIn[3]=pres;
    }
    else if(mode==2) //determines T and P from rho and ei
    {
      temp = ei / fac;
      pres = (gc-1.0) * den * ei;
      
      eosIn[0]=den;
      eosIn[1]=temp;
      eosIn[2]=ei;
      eosIn[3]=pres;
    }
    else if(mode==3) // determines T and ei from rho and P
    {
      // seems like this mode is never used...
      cout << "EOS mode 3 not done, and never really needed..." << endl;
      AbortChem();
    }
    else // doesn't change anything... not good.
    {
      cout << "Bad option in Chem_EOS " << mode << endl;
      AbortChem();
        
    }
    
    if(eosIn[1] != eosIn[1] || eosIn[2] != eosIn[2]) {
      cout << "nan after EOS! mode = " << mode << endl;
      cout << "loc = " << loc << endl;
      cout << "eos[0] = " << eosIn[0] << endl;
      cout << "eos[1] = " << eosIn[1] << endl;
      cout << "eos[2] = " << eosIn[2] << endl;
      cout << "gc = " << gc << endl;
      AbortChem();
    }
    
    if(temp > 1.0e9){
      cout << "temp too high after EOS. temp = " << temp << endl;
      cout << "mode = " << mode << endl;
      cout << "loc = " << loc << endl;
      cout << "eos[0] = " << eosIn[0] << endl;
      cout << "eos[1] = " << eosIn[1] << endl;
      cout << "eos[2] = " << eosIn[2] << endl;
      cout << "yh2 = " << y[iH2] << endl;
      cout << "yie = " << y[iELEC] << endl;
      AbortChem();
    
    }


    
}


