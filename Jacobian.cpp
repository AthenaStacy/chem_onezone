#include <iostream>
#include <iomanip>
#include <stdlib.h>
using namespace std;

#include "ChemGlobals.H"
using namespace ChemNS;

#include <math.h>

//REAL max(REAL x, REAL y) { if(x>y) return x; else return y; }
//REAL min(REAL x, REAL y) { if(x<y) return x; else return y; }

void fill_species(REAL* ymass, REAL* ytot);


void Jacobian(REAL x,REAL* yt,REAL (*dfdy)[NINT], REAL* ratraw, double* rpar)
{
    
    int i,j;
    
    for(i=0;i<NINT;i++) yt[i] = min(1.0 , max(yt[i] , 1.0e-20) );
    for(i=0;i<NINT;i++) for(j=0;j<NINT;j++) dfdy[i][j] = 0.0;
    
    REAL ys[NSPECIES];
    fill_species(yt,ys);
    
    // Jacobian Values

#if DO_CHEMISTRY == PRIMORDIAL
    
    dfdy[iHP][iHP]=  -ratraw[iR003]*ys[iH]-ratraw[iR005]*ys[iHM]-ratraw[iR007]*ys[iH2]
    -ratraw[iR013]*ys[iELEC]-ratraw[iR016]*ys[iHM]-ratraw[iR027]*ys[iHE]
    -ratraw[iR034]*ys[iD]-ratraw[iR041]*ys[iHD]-ratraw[iR042]*ys[iD]
    -ratraw[iR060]*ys[iDM]-ratraw[iR067]*ys[iDM]-ratraw[iR092]*ys[iHD]
    -ratraw[iR093]*ys[iHD]-ratraw[iR097]*ys[iD2]-ratraw[iR098]*ys[iD2]
    -ratraw[iR099]*ys[iD2];
    
    
    dfdy[iHP][iH]=  -ratraw[iR003]*ys[iHP]+ratraw[iR004]*ys[iH2P]+ratraw[iR012]*ys[iELEC]
    +ratraw[iR026]*ys[iHEP]+ratraw[iR035]*ys[iDP]+ratraw[iR038]*ys[iHDP]
    +ratraw[iR087]*ys[iD2P];
    
    dfdy[iHP][iHM]=  -ratraw[iR005]*ys[iHP]-ratraw[iR016]*ys[iHP];
    
    dfdy[iHP][iH2P]=  ratraw[iR004]*ys[iH]+ratraw[iR082]*ys[iD]+ratraw[iR118];
    
    dfdy[iHP][iH2]=  -ratraw[iR007]*ys[iHP]+ratraw[iR024]*ys[iHEP]+ratraw[iR039]*ys[iDP];
    
    dfdy[iHP][iDP]=  ratraw[iR035]*ys[iH]+ratraw[iR039]*ys[iH2]+ratraw[iR095]*ys[iHD];
    
    dfdy[iHP][iD]=  -ratraw[iR034]*ys[iHP]-ratraw[iR042]*ys[iHP]+ratraw[iR082]*ys[iH2P]
    +ratraw[iR085]*ys[iHDP];
    
    dfdy[iHP][iDM]=  -ratraw[iR060]*ys[iHP]-ratraw[iR067]*ys[iHP];
    
    dfdy[iHP][iHDP]=  ratraw[iR038]*ys[iH]+ratraw[iR085]*ys[iD]+ratraw[iR120];
    
    dfdy[iHP][iHD]=  -ratraw[iR041]*ys[iHP]-ratraw[iR092]*ys[iHP]-ratraw[iR093]*ys[iHP]
    +ratraw[iR095]*ys[iDP]+ratraw[iR102]*ys[iHEP];
    
    dfdy[iHP][iD2P]=  ratraw[iR087]*ys[iH];
    
    dfdy[iHP][iD2]=  -ratraw[iR097]*ys[iHP]-ratraw[iR098]*ys[iHP]-ratraw[iR099]*ys[iHP];
    
    dfdy[iHP][iHEP]=  ratraw[iR024]*ys[iH2]+ratraw[iR026]*ys[iH]+ratraw[iR102]*ys[iHD];
    
    dfdy[iHP][iHE]=  -ratraw[iR027]*ys[iHP];
    
    dfdy[iHP][iHEPP]= 0.0e0;
    
    //dfdy[iHP][iELEC]=  ratraw[iR012]*ys[iH]-ratraw[iR013]*ys[iHP];
    
    
    
    dfdy[iH][iHP]=  -ratraw[iR003]*ys[iH]+ratraw[iR005]*ys[iHM]+ratraw[iR005]*ys[iHM]
    +ratraw[iR007]*ys[iH2]+ratraw[iR013]*ys[iELEC]+ratraw[iR027]*ys[iHE]
    +ratraw[iR034]*ys[iD]+ratraw[iR067]*ys[iDM]+ratraw[iR092]*ys[iHD]
    +ratraw[iR099]*ys[iD2];
    
    dfdy[iH][iH]=  -ratraw[iR001]*ys[iELEC]-ratraw[iR002]*ys[iHM]-ratraw[iR003]*ys[iHP]
    -ratraw[iR004]*ys[iH2P]-ratraw[iR009]*ys[iH2]+ratraw[iR009]*ys[iH2]
    +ratraw[iR009]*ys[iH2]+ratraw[iR009]*ys[iH2]-ratraw[iR012]*ys[iELEC]
    -ratraw[iR015]*ys[iHM]+ratraw[iR015]*ys[iHM]+ratraw[iR015]*ys[iHM]
    -ratraw[iR026]*ys[iHEP]-ratraw[iR030]*ys[iH]*ys[iH]-ratraw[iR030]*ys[iH]*ys[iH]
    -ratraw[iR030]*ys[iH]*ys[iH]-ratraw[iR030]*ys[iH]*ys[iH]
    -ratraw[iR030]*ys[iH]*ys[iH]-ratraw[iR030]*ys[iH]*ys[iH]
    -ratraw[iR030]*ys[iH]*ys[iH]-ratraw[iR030]*ys[iH]*ys[iH]
    -ratraw[iR030]*ys[iH]*ys[iH]+ratraw[iR030]*ys[iH]*ys[iH]
    +ratraw[iR030]*ys[iH]*ys[iH]+ratraw[iR030]*ys[iH]*ys[iH]
    -ratraw[iR031]*ys[iH]*ys[iH2]-ratraw[iR031]*ys[iH]*ys[iH2]
    -ratraw[iR031]*ys[iH]*ys[iH2]-ratraw[iR031]*ys[iH]*ys[iH2]
    -ratraw[iR032]*ys[iH]*ys[iHE]-ratraw[iR032]*ys[iH]*ys[iHE]
    -ratraw[iR032]*ys[iH]*ys[iHE]-ratraw[iR032]*ys[iH]*ys[iHE]
    -ratraw[iR035]*ys[iDP]-ratraw[iR036]*ys[iD]-ratraw[iR038]*ys[iHDP]
    -ratraw[iR040]*ys[iHD]-ratraw[iR043]*ys[iDP]-ratraw[iR050]*ys[iHDP]
    -ratraw[iR052]*ys[iDM]-ratraw[iR055]*ys[iDM]-ratraw[iR064]*ys[iDM]
    +ratraw[iR064]*ys[iDM]-ratraw[iR083]*ys[iHDP]-ratraw[iR087]*ys[iD2P]
    -ratraw[iR088]*ys[iD2P]-ratraw[iR089]*ys[iD2P]-ratraw[iR107]*ys[iD2]
    -ratraw[iR108]*ys[iHD]+ratraw[iR108]*ys[iHD]+ratraw[iR108]*ys[iHD]
    -ratraw[iR112]*ys[iD2]+ratraw[iR112]*ys[iD2];
    
    dfdy[iH][iHM]=  -ratraw[iR002]*ys[iH]+ratraw[iR005]*ys[iHP]+ratraw[iR005]*ys[iHP]
    +ratraw[iR014]*ys[iELEC]-ratraw[iR015]*ys[iH]+ratraw[iR015]*ys[iH]
    +ratraw[iR015]*ys[iH]+ratraw[iR021]*ys[iH2P]+ratraw[iR022]*ys[iH2P]
    +ratraw[iR022]*ys[iH2P]+ratraw[iR022]*ys[iH2P]+ratraw[iR028]*ys[iHEP]
    +ratraw[iR029]*ys[iHE]+ratraw[iR053]*ys[iD]+ratraw[iR066]*ys[iDP]
    +ratraw[iR071]*ys[iHDP]+ratraw[iR072]*ys[iHDP]+ratraw[iR072]*ys[iHDP]
    +ratraw[iR075]*ys[iD2P]+ratraw[iR076]*ys[iD2P]+ratraw[iR116];
    
    dfdy[iH][iH2P]=  -ratraw[iR004]*ys[iH]+ratraw[iR006]*ys[iELEC]+ratraw[iR006]*ys[iELEC]
    +ratraw[iR021]*ys[iHM]+ratraw[iR022]*ys[iHM]+ratraw[iR022]*ys[iHM]
    +ratraw[iR022]*ys[iHM]+ratraw[iR048]*ys[iD]+ratraw[iR070]*ys[iDM]
    +ratraw[iR070]*ys[iDM]+ratraw[iR118];
    
    dfdy[iH][iH2]=  ratraw[iR007]*ys[iHP]+ratraw[iR008]*ys[iELEC]+ratraw[iR008]*ys[iELEC]
    -ratraw[iR009]*ys[iH]+ratraw[iR009]*ys[iH]+ratraw[iR009]*ys[iH]
    +ratraw[iR009]*ys[iH]+ratraw[iR010]*ys[iH2]+ratraw[iR010]*ys[iH2]
    +ratraw[iR010]*ys[iH2]+ratraw[iR010]*ys[iH2]+ratraw[iR011]*ys[iHE]
    +ratraw[iR011]*ys[iHE]+ratraw[iR023]*ys[iELEC]+ratraw[iR024]*ys[iHEP]
    -ratraw[iR031]*ys[iH]*ys[iH]-ratraw[iR031]*ys[iH]*ys[iH]
    +ratraw[iR037]*ys[iD]+ratraw[iR091]*ys[iDP]+ratraw[iR109]*ys[iHD]
    +ratraw[iR122]+ratraw[iR122];
    
    dfdy[iH][iDP]=  -ratraw[iR035]*ys[iH]-ratraw[iR043]*ys[iH]+ratraw[iR066]*ys[iHM]
    +ratraw[iR091]*ys[iH2]+ratraw[iR096]*ys[iHD];
    
    dfdy[iH][iD]=  ratraw[iR034]*ys[iHP]-ratraw[iR036]*ys[iH]+ratraw[iR037]*ys[iH2]
    +ratraw[iR048]*ys[iH2P]+ratraw[iR053]*ys[iHM]+ratraw[iR084]*ys[iHDP]
    +ratraw[iR106]*ys[iHD];
    
    dfdy[iH][iDM]=  -ratraw[iR052]*ys[iH]-ratraw[iR055]*ys[iH]-ratraw[iR064]*ys[iH]
    +ratraw[iR064]*ys[iH]+ratraw[iR067]*ys[iHP]+ratraw[iR070]*ys[iH2P]
    +ratraw[iR070]*ys[iH2P]+ratraw[iR074]*ys[iHDP];
    
    dfdy[iH][iHDP]=  -ratraw[iR038]*ys[iH]+ratraw[iR044]*ys[iELEC]-ratraw[iR050]*ys[iH]
    +ratraw[iR071]*ys[iHM]+ratraw[iR072]*ys[iHM]+ratraw[iR072]*ys[iHM]
    +ratraw[iR074]*ys[iDM]-ratraw[iR083]*ys[iH]+ratraw[iR084]*ys[iD]+ratraw[iR119];
    
    dfdy[iH][iHD]=  -ratraw[iR040]*ys[iH]+ratraw[iR057]*ys[iELEC]+ratraw[iR092]*ys[iHP]
    +ratraw[iR096]*ys[iDP]+ratraw[iR103]*ys[iHEP]+ratraw[iR106]*ys[iD]
    -ratraw[iR108]*ys[iH]+ratraw[iR108]*ys[iH]+ratraw[iR108]*ys[iH]
    +ratraw[iR109]*ys[iH2]+ratraw[iR110]*ys[iHE]+ratraw[iR111]*ys[iELEC]+ratraw[iR123];
    
    dfdy[iH][iD2P]=  ratraw[iR075]*ys[iHM]+ratraw[iR076]*ys[iHM]-ratraw[iR087]*ys[iH]
    -ratraw[iR088]*ys[iH]-ratraw[iR089]*ys[iH];
    
    dfdy[iH][iD2]=  ratraw[iR099]*ys[iHP]-ratraw[iR107]*ys[iH]-ratraw[iR112]*ys[iH]+ratraw[iR112]*ys[iH];
    
    dfdy[iH][iHEP]=  ratraw[iR024]*ys[iH2]-ratraw[iR026]*ys[iH]+ratraw[iR028]*ys[iHM]+ratraw[iR103]*ys[iHD];
    
    dfdy[iH][iHE]=  ratraw[iR011]*ys[iH2]+ratraw[iR011]*ys[iH2]+ratraw[iR027]*ys[iHP]
    +ratraw[iR029]*ys[iHM]-ratraw[iR032]*ys[iH]*ys[iH]-ratraw[iR032]*ys[iH]*ys[iH]
    +ratraw[iR110]*ys[iHD];
    
    dfdy[iH][iHEPP]= 0.0e0;
    
    //dfdy[iH][iELEC]=  -ratraw[iR001]*ys[iH]+ratraw[iR006]*ys[iH2P]+ratraw[iR006]*ys[iH2P]
    //+ratraw[iR008]*ys[iH2]+ratraw[iR008]*ys[iH2]-ratraw[iR012]*ys[iH]
    //+ratraw[iR013]*ys[iHP]+ratraw[iR014]*ys[iHM]+ratraw[iR023]*ys[iH2]
    //+ratraw[iR044]*ys[iHDP]+ratraw[iR057]*ys[iHD]+ratraw[iR111]*ys[iHD];
    
    
    
    dfdy[iHM][iHP]=  -ratraw[iR005]*ys[iHM]-ratraw[iR016]*ys[iHM];
    
    dfdy[iHM][iH]=  ratraw[iR001]*ys[iELEC]-ratraw[iR002]*ys[iHM]-ratraw[iR015]*ys[iHM]
    +ratraw[iR052]*ys[iDM];
    
    dfdy[iHM][iHM]=  -ratraw[iR002]*ys[iH]-ratraw[iR005]*ys[iHP]-ratraw[iR014]*ys[iELEC]
    -ratraw[iR015]*ys[iH]-ratraw[iR016]*ys[iHP]-ratraw[iR021]*ys[iH2P]
    -ratraw[iR022]*ys[iH2P]-ratraw[iR028]*ys[iHEP]-ratraw[iR029]*ys[iHE]
    -ratraw[iR053]*ys[iD]-ratraw[iR054]*ys[iD]-ratraw[iR061]*ys[iDP]
    -ratraw[iR066]*ys[iDP]-ratraw[iR071]*ys[iHDP]-ratraw[iR072]*ys[iHDP]
    -ratraw[iR075]*ys[iD2P]-ratraw[iR076]*ys[iD2P]-ratraw[iR116];
    
    dfdy[iHM][iH2P]=  -ratraw[iR021]*ys[iHM]-ratraw[iR022]*ys[iHM];
    
    dfdy[iHM][iH2]=  ratraw[iR023]*ys[iELEC];
    
    dfdy[iHM][iDP]=  -ratraw[iR061]*ys[iHM]-ratraw[iR066]*ys[iHM];
    
    dfdy[iHM][iD]=  -ratraw[iR053]*ys[iHM]-ratraw[iR054]*ys[iHM];
    
    dfdy[iHM][iDM]=  ratraw[iR052]*ys[iH];
    
    dfdy[iHM][iHDP]=  -ratraw[iR071]*ys[iHM]-ratraw[iR072]*ys[iHM];
    
    dfdy[iHM][iHD]=  ratraw[iR058]*ys[iELEC];
    
    dfdy[iHM][iD2P]=  -ratraw[iR075]*ys[iHM]-ratraw[iR076]*ys[iHM];
    
    dfdy[iHM][iD2]=  0.0e0;
    
    dfdy[iHM][iHEP]=  -ratraw[iR028]*ys[iHM];
    
    dfdy[iHM][iHE]=  -ratraw[iR029]*ys[iHM];
    
    dfdy[iHM][iHEPP]=  0.0e0;
    
    //dfdy[iHM][iELEC]=  ratraw[iR001]*ys[iH]-ratraw[iR014]*ys[iHM]
    //+ratraw[iR023]*ys[iH2]+ratraw[iR058]*ys[iHD];
    
    
    
    dfdy[iH2P][iHP]=  ratraw[iR003]*ys[iH]+ratraw[iR007]*ys[iH2]+ratraw[iR016]*ys[iHM]
    +ratraw[iR093]*ys[iHD];
    
    dfdy[iH2P][iH]=  ratraw[iR003]*ys[iHP]-ratraw[iR004]*ys[iH2P]+ratraw[iR050]*ys[iHDP];
    
    dfdy[iH2P][iHM]=  ratraw[iR016]*ys[iHP]-ratraw[iR021]*ys[iH2P]-ratraw[iR022]*ys[iH2P];
    
    dfdy[iH2P][iH2P]=  -ratraw[iR004]*ys[iH]-ratraw[iR006]*ys[iELEC]-ratraw[iR021]*ys[iHM]
    -ratraw[iR022]*ys[iHM]-ratraw[iR048]*ys[iD]-ratraw[iR069]*ys[iDM]
    -ratraw[iR070]*ys[iDM]-ratraw[iR081]*ys[iD]-ratraw[iR082]*ys[iD]
    -ratraw[iR118];
    
    dfdy[iH2P][iH2]=  ratraw[iR007]*ys[iHP]+ratraw[iR025]*ys[iHEP]+ratraw[iR090]*ys[iDP];
    
    dfdy[iH2P][iDP]=  ratraw[iR090]*ys[iH2];
    
    dfdy[iH2P][iD]=  -ratraw[iR048]*ys[iH2P]-ratraw[iR081]*ys[iH2P]-ratraw[iR082]*ys[iH2P];
    
    dfdy[iH2P][iDM]=  -ratraw[iR069]*ys[iH2P]-ratraw[iR070]*ys[iH2P];
    
    dfdy[iH2P][iHDP]=  ratraw[iR050]*ys[iH];
    
    dfdy[iH2P][iHD]=  ratraw[iR093]*ys[iHP];
    
    dfdy[iH2P][iD2P]= 0.0e0;
    
    dfdy[iH2P][iD2]=  0.0e0;
    
    dfdy[iH2P][iHEP]=  ratraw[iR025]*ys[iH2];
    
    dfdy[iH2P][iHE]=  0.0e0;
    
    dfdy[iH2P][iHEPP]=  0.0e0;
    
    //dfdy[iH2P][iELEC]=  -ratraw[iR006]*ys[iH2P];
    
    
    
    dfdy[iH2][iHP]=  -ratraw[iR007]*ys[iH2]+ratraw[iR041]*ys[iHD];
    
    dfdy[iH2][iH]=  ratraw[iR002]*ys[iHM]+ratraw[iR004]*ys[iH2P]-ratraw[iR009]*ys[iH2]
    +ratraw[iR030]*ys[iH]*ys[iH]+ratraw[iR030]*ys[iH]*ys[iH]
    +ratraw[iR030]*ys[iH]*ys[iH]-ratraw[iR031]*ys[iH]*ys[iH2]
    -ratraw[iR031]*ys[iH]*ys[iH2]+ratraw[iR031]*ys[iH]*ys[iH2]
    +ratraw[iR031]*ys[iH]*ys[iH2]+ratraw[iR031]*ys[iH]*ys[iH2]
    +ratraw[iR031]*ys[iH]*ys[iH2]+ratraw[iR032]*ys[iH]*ys[iHE]
    +ratraw[iR032]*ys[iH]*ys[iHE]+ratraw[iR040]*ys[iHD]+ratraw[iR083]*ys[iHDP];
    
    dfdy[iH2][iHM]=  ratraw[iR002]*ys[iH]+ratraw[iR021]*ys[iH2P];
    
    dfdy[iH2][iH2P]=  ratraw[iR004]*ys[iH]+ratraw[iR021]*ys[iHM]
    +ratraw[iR069]*ys[iDM]+ratraw[iR081]*ys[iD];
    
    dfdy[iH2][iH2]=  -ratraw[iR007]*ys[iHP]-ratraw[iR008]*ys[iELEC]-ratraw[iR009]*ys[iH]
    -ratraw[iR010]*ys[iH2]-ratraw[iR010]*ys[iH2]-ratraw[iR010]*ys[iH2]
    -ratraw[iR010]*ys[iH2]+ratraw[iR010]*ys[iH2]+ratraw[iR010]*ys[iH2]
    -ratraw[iR011]*ys[iHE]-ratraw[iR023]*ys[iELEC]-ratraw[iR024]*ys[iHEP]
    -ratraw[iR025]*ys[iHEP]-ratraw[iR031]*ys[iH]*ys[iH]+ratraw[iR031]*ys[iH]*ys[iH]
    +ratraw[iR031]*ys[iH]*ys[iH]-ratraw[iR037]*ys[iD]-ratraw[iR039]*ys[iDP]
    -ratraw[iR090]*ys[iDP]-ratraw[iR091]*ys[iDP]-ratraw[iR109]*ys[iHD]
    +ratraw[iR109]*ys[iHD]-ratraw[iR113]*ys[iD2]+ratraw[iR113]*ys[iD2]-ratraw[iR122];
    
    dfdy[iH2][iDP]=  -ratraw[iR039]*ys[iH2]-ratraw[iR090]*ys[iH2]-ratraw[iR091]*ys[iH2];
    
    dfdy[iH2][iD]=  -ratraw[iR037]*ys[iH2]+ratraw[iR081]*ys[iH2P];
    
    dfdy[iH2][iDM]=  ratraw[iR069]*ys[iH2P];
    
    dfdy[iH2][iHDP]=  ratraw[iR083]*ys[iH];
    
    dfdy[iH2][iHD]=  ratraw[iR040]*ys[iH]+ratraw[iR041]*ys[iHP]-ratraw[iR109]*ys[iH2]+ratraw[iR109]*ys[iH2];
    
    dfdy[iH2][iD2P]= 0.0e0;
    
    dfdy[iH2][iD2]=  -ratraw[iR113]*ys[iH2]+ratraw[iR113]*ys[iH2];
    
    dfdy[iH2][iHEP]=  -ratraw[iR024]*ys[iH2]-ratraw[iR025]*ys[iH2];
    
    dfdy[iH2][iHE]=  -ratraw[iR011]*ys[iH2]+ratraw[iR032]*ys[iH]*ys[iH];
    
    dfdy[iH2][iHEPP]=  0.0e0;
    
    //dfdy[iH2][iELEC]=  -ratraw[iR008]*ys[iH2]-ratraw[iR023]*ys[iH2];
    
    
    
    dfdy[iDP][iHP]=  ratraw[iR034]*ys[iD]+ratraw[iR041]*ys[iHD]+ratraw[iR097]*ys[iD2];
    
    dfdy[iDP][iH]=  -ratraw[iR035]*ys[iDP]-ratraw[iR043]*ys[iDP]+ratraw[iR083]*ys[iHDP]+ratraw[iR089]*ys[iD2P];
    
    dfdy[iDP][iHM]=  -ratraw[iR061]*ys[iDP]-ratraw[iR066]*ys[iDP];
    
    dfdy[iDP][iH2P]=  ratraw[iR081]*ys[iD];
    
    dfdy[iDP][iH2]=  -ratraw[iR039]*ys[iDP]-ratraw[iR090]*ys[iDP]-ratraw[iR091]*ys[iDP];
    
    dfdy[iDP][iDP]=  -ratraw[iR033]*ys[iELEC]-ratraw[iR035]*ys[iH]-ratraw[iR039]*ys[iH2]
    -ratraw[iR043]*ys[iH]-ratraw[iR047]*ys[iHE]-ratraw[iR061]*ys[iHM]
    -ratraw[iR062]*ys[iDM]-ratraw[iR066]*ys[iHM]-ratraw[iR068]*ys[iDM]
    -ratraw[iR080]*ys[iD]-ratraw[iR090]*ys[iH2]-ratraw[iR091]*ys[iH2]
    -ratraw[iR094]*ys[iHD]-ratraw[iR095]*ys[iHD]-ratraw[iR096]*ys[iHD]
    -ratraw[iR100]*ys[iD2];
    
    dfdy[iDP][iD]=  ratraw[iR034]*ys[iHP]+ratraw[iR045]*ys[iELEC]+ratraw[iR046]*ys[iHEP]
    +ratraw[iR049]*ys[iHDP]-ratraw[iR080]*ys[iDP]+ratraw[iR081]*ys[iH2P]
    +ratraw[iR086]*ys[iD2P];
    
    dfdy[iDP][iDM]=  -ratraw[iR062]*ys[iDP]-ratraw[iR068]*ys[iDP];
    
    dfdy[iDP][iHDP]=  ratraw[iR049]*ys[iD]+ratraw[iR083]*ys[iH]+ratraw[iR119];
    
    dfdy[iDP][iHD]=  ratraw[iR041]*ys[iHP]-ratraw[iR094]*ys[iDP]-ratraw[iR095]*ys[iDP]
    -ratraw[iR096]*ys[iDP]+ratraw[iR103]*ys[iHEP];
    
    dfdy[iDP][iD2P]=  ratraw[iR086]*ys[iD]+ratraw[iR089]*ys[iH]+ratraw[iR121];
    
    dfdy[iDP][iD2]=  ratraw[iR097]*ys[iHP]-ratraw[iR100]*ys[iDP]+ratraw[iR105]*ys[iHEP];
    
    dfdy[iDP][iHEP]=  ratraw[iR046]*ys[iD]+ratraw[iR103]*ys[iHD]+ratraw[iR105]*ys[iD2];
    
    dfdy[iDP][iHE]=  -ratraw[iR047]*ys[iDP];
    
    dfdy[iDP][iHEPP]= 0.0e0;
    
   // dfdy[iDP][iELEC]=  -ratraw[iR033]*ys[iDP]+ratraw[iR045]*ys[iD];
    
    
    
    dfdy[iD][iHP]=  -ratraw[iR034]*ys[iD]-ratraw[iR042]*ys[iD]+ratraw[iR067]*ys[iDM]
    +ratraw[iR093]*ys[iHD]+ratraw[iR098]*ys[iD2];
    
    dfdy[iD][iH]=  ratraw[iR035]*ys[iDP]-ratraw[iR036]*ys[iD]+ratraw[iR040]*ys[iHD]
    +ratraw[iR050]*ys[iHDP]+ratraw[iR052]*ys[iDM]+ratraw[iR064]*ys[iDM]
    +ratraw[iR088]*ys[iD2P]+ratraw[iR107]*ys[iD2]+ratraw[iR108]*ys[iHD]
    +ratraw[iR112]*ys[iD2]+ratraw[iR112]*ys[iD2];
    
    dfdy[iD][iHM]=  -ratraw[iR053]*ys[iD]-ratraw[iR054]*ys[iD]+ratraw[iR066]*ys[iDP]
    +ratraw[iR072]*ys[iHDP]+ratraw[iR076]*ys[iD2P]+ratraw[iR076]*ys[iD2P];
    
    dfdy[iD][iH2P]=  -ratraw[iR048]*ys[iD]+ratraw[iR069]*ys[iDM]+ratraw[iR070]*ys[iDM]
    -ratraw[iR081]*ys[iD]-ratraw[iR082]*ys[iD];
    
    dfdy[iD][iH2]=  -ratraw[iR037]*ys[iD]+ratraw[iR090]*ys[iDP]+ratraw[iR109]*ys[iHD]
    +ratraw[iR113]*ys[iD2]+ratraw[iR113]*ys[iD2];
    
    dfdy[iD][iDP]=  ratraw[iR033]*ys[iELEC]+ratraw[iR035]*ys[iH]+ratraw[iR047]*ys[iHE]
    +ratraw[iR066]*ys[iHM]+ratraw[iR068]*ys[iDM]+ratraw[iR068]*ys[iDM]
    -ratraw[iR080]*ys[iD]+ratraw[iR090]*ys[iH2]+ratraw[iR094]*ys[iHD]
    +ratraw[iR100]*ys[iD2];
    
    dfdy[iD][iD]=  -ratraw[iR034]*ys[iHP]-ratraw[iR036]*ys[iH]-ratraw[iR037]*ys[iH2]
    -ratraw[iR042]*ys[iHP]-ratraw[iR045]*ys[iELEC]-ratraw[iR046]*ys[iHEP]
    -ratraw[iR048]*ys[iH2P]-ratraw[iR049]*ys[iHDP]-ratraw[iR051]*ys[iELEC]
    -ratraw[iR053]*ys[iHM]-ratraw[iR054]*ys[iHM]-ratraw[iR056]*ys[iDM]
    -ratraw[iR080]*ys[iDP]-ratraw[iR081]*ys[iH2P]-ratraw[iR082]*ys[iH2P]
    -ratraw[iR084]*ys[iHDP]-ratraw[iR085]*ys[iHDP]-ratraw[iR086]*ys[iD2P]
    -ratraw[iR106]*ys[iHD];
    
    dfdy[iD][iDM]=  ratraw[iR052]*ys[iH]-ratraw[iR056]*ys[iD]+ratraw[iR063]*ys[iELEC]
    +ratraw[iR064]*ys[iH]+ratraw[iR065]*ys[iHE]+ratraw[iR067]*ys[iHP]
    +ratraw[iR068]*ys[iDP]+ratraw[iR068]*ys[iDP]+ratraw[iR069]*ys[iH2P]
    +ratraw[iR070]*ys[iH2P]+ratraw[iR073]*ys[iHDP]+ratraw[iR074]*ys[iHDP]
    +ratraw[iR074]*ys[iHDP]+ratraw[iR077]*ys[iD2P]+ratraw[iR078]*ys[iD2P]
    +ratraw[iR078]*ys[iD2P]+ratraw[iR078]*ys[iD2P]+ratraw[iR079]*ys[iHEP]
    +ratraw[iR117];
    
    dfdy[iD][iHDP]=  ratraw[iR044]*ys[iELEC]-ratraw[iR049]*ys[iD]+ratraw[iR050]*ys[iH]
    +ratraw[iR072]*ys[iHM]+ratraw[iR073]*ys[iDM]+ratraw[iR074]*ys[iDM]
    +ratraw[iR074]*ys[iDM]-ratraw[iR084]*ys[iD]-ratraw[iR085]*ys[iD]+ratraw[iR120];
    
    dfdy[iD][iHD]=  ratraw[iR040]*ys[iH]+ratraw[iR058]*ys[iELEC]+ratraw[iR093]*ys[iHP]
    +ratraw[iR094]*ys[iDP]+ratraw[iR102]*ys[iHEP]-ratraw[iR106]*ys[iD]
    +ratraw[iR108]*ys[iH]+ratraw[iR109]*ys[iH2]+ratraw[iR110]*ys[iHE]
    +ratraw[iR111]*ys[iELEC]+ratraw[iR123];
    
    dfdy[iD][iD2P]=  ratraw[iR076]*ys[iHM]+ratraw[iR076]*ys[iHM]+ratraw[iR077]*ys[iDM]
    +ratraw[iR078]*ys[iDM]+ratraw[iR078]*ys[iDM]+ratraw[iR078]*ys[iDM]
    -ratraw[iR086]*ys[iD]+ratraw[iR088]*ys[iH]+ratraw[iR121];
    
    dfdy[iD][iD2]=  ratraw[iR059]*ys[iELEC]+ratraw[iR098]*ys[iHP]+ratraw[iR100]*ys[iDP]
    +ratraw[iR105]*ys[iHEP]+ratraw[iR107]*ys[iH]+ratraw[iR112]*ys[iH]
    +ratraw[iR112]*ys[iH]+ratraw[iR113]*ys[iH2]+ratraw[iR113]*ys[iH2]
    +ratraw[iR114]*ys[iHE]+ratraw[iR114]*ys[iHE]+ratraw[iR115]*ys[iELEC]
    +ratraw[iR115]*ys[iELEC]+ratraw[iR124]+ratraw[iR124];
    
    dfdy[iD][iHEP]=  -ratraw[iR046]*ys[iD]+ratraw[iR079]*ys[iDM]
    +ratraw[iR102]*ys[iHD]+ratraw[iR105]*ys[iD2];
    
    dfdy[iD][iHE]=  ratraw[iR047]*ys[iDP]+ratraw[iR065]*ys[iDM]+ratraw[iR110]*ys[iHD]
    +ratraw[iR114]*ys[iD2]+ratraw[iR114]*ys[iD2];
    
    dfdy[iD][iHEPP]= 0.0e0;
    
    //dfdy[iD][iELEC]=  ratraw[iR033]*ys[iDP]+ratraw[iR044]*ys[iHDP]-ratraw[iR045]*ys[iD]
    //-ratraw[iR051]*ys[iD]+ratraw[iR058]*ys[iHD]+ratraw[iR059]*ys[iD2]
    //+ratraw[iR063]*ys[iDM]+ratraw[iR111]*ys[iHD]+ratraw[iR115]*ys[iD2]
    //+ratraw[iR115]*ys[iD2];
    
    
    
    dfdy[iDM][iHP]=  -ratraw[iR060]*ys[iDM]-ratraw[iR067]*ys[iDM];
    
    dfdy[iDM][iH]=  -ratraw[iR052]*ys[iDM]-ratraw[iR055]*ys[iDM]-ratraw[iR064]*ys[iDM];
    
    dfdy[iDM][iHM]=  ratraw[iR053]*ys[iD];
    
    dfdy[iDM][iH2P]=  -ratraw[iR069]*ys[iDM]-ratraw[iR070]*ys[iDM];
    
    dfdy[iDM][iH2]=  0.0e0;
    
    dfdy[iDM][iDP]=  -ratraw[iR062]*ys[iDM]-ratraw[iR068]*ys[iDM];
    
    dfdy[iDM][iD]=  ratraw[iR051]*ys[iELEC]+ratraw[iR053]*ys[iHM]-ratraw[iR056]*ys[iDM];
    
    dfdy[iDM][iDM]=  -ratraw[iR052]*ys[iH]-ratraw[iR055]*ys[iH]-ratraw[iR056]*ys[iD]
    -ratraw[iR060]*ys[iHP]-ratraw[iR062]*ys[iDP]-ratraw[iR063]*ys[iELEC]
    -ratraw[iR064]*ys[iH]-ratraw[iR065]*ys[iHE]-ratraw[iR067]*ys[iHP]
    -ratraw[iR068]*ys[iDP]-ratraw[iR069]*ys[iH2P]-ratraw[iR070]*ys[iH2P]
    -ratraw[iR073]*ys[iHDP]-ratraw[iR074]*ys[iHDP]-ratraw[iR077]*ys[iD2P]
    -ratraw[iR078]*ys[iD2P]-ratraw[iR079]*ys[iHEP]-ratraw[iR117];
    
    dfdy[iDM][iHDP]=  -ratraw[iR073]*ys[iDM]-ratraw[iR074]*ys[iDM];
    
    dfdy[iDM][iHD]=  ratraw[iR057]*ys[iELEC];
    
    dfdy[iDM][iD2P]=  -ratraw[iR077]*ys[iDM]-ratraw[iR078]*ys[iDM];
    
    dfdy[iDM][iD2]=  ratraw[iR059]*ys[iELEC];
    
    dfdy[iDM][iHEP]=  -ratraw[iR079]*ys[iDM];
    
    dfdy[iDM][iHE]=  -ratraw[iR065]*ys[iDM];
    
    dfdy[iDM][iHEPP]= 0.0e0;
    
    //dfdy[iDM][iELEC]=  ratraw[iR051]*ys[iD]+ratraw[iR057]*ys[iHD]+ratraw[iR059]*ys[iD2]-ratraw[iR063]*ys[iDM];
    
    
    
    dfdy[iHDP][iHP]=  ratraw[iR042]*ys[iD]+ratraw[iR060]*ys[iDM]+ratraw[iR092]*ys[iHD]+ratraw[iR098]*ys[iD2];
    
    dfdy[iHDP][iH]=  -ratraw[iR038]*ys[iHDP]+ratraw[iR043]*ys[iDP]-ratraw[iR050]*ys[iHDP]
    -ratraw[iR083]*ys[iHDP]+ratraw[iR088]*ys[iD2P];
    
    dfdy[iHDP][iHM]=  ratraw[iR061]*ys[iDP]-ratraw[iR071]*ys[iHDP]-ratraw[iR072]*ys[iHDP];
    
    dfdy[iHDP][iH2P]=  ratraw[iR048]*ys[iD];
    
    dfdy[iHDP][iH2]=  ratraw[iR091]*ys[iDP];
    
    dfdy[iHDP][iDP]=  ratraw[iR043]*ys[iH]+ratraw[iR061]*ys[iHM]+ratraw[iR091]*ys[iH2]+ratraw[iR094]*ys[iHD];
    
    dfdy[iHDP][iD]=  ratraw[iR042]*ys[iHP]+ratraw[iR048]*ys[iH2P]-ratraw[iR049]*ys[iHDP]
    -ratraw[iR084]*ys[iHDP]-ratraw[iR085]*ys[iHDP];
    
    dfdy[iHDP][iDM]=  ratraw[iR060]*ys[iHP]-ratraw[iR073]*ys[iHDP]-ratraw[iR074]*ys[iHDP];
    
    dfdy[iHDP][iHDP]=  -ratraw[iR038]*ys[iH]-ratraw[iR044]*ys[iELEC]-ratraw[iR049]*ys[iD]
    -ratraw[iR050]*ys[iH]-ratraw[iR071]*ys[iHM]-ratraw[iR072]*ys[iHM]
    -ratraw[iR073]*ys[iDM]-ratraw[iR074]*ys[iDM]-ratraw[iR083]*ys[iH]
    -ratraw[iR084]*ys[iD]-ratraw[iR085]*ys[iD]-ratraw[iR119]-ratraw[iR120];
    
    dfdy[iHDP][iHD]=  ratraw[iR092]*ys[iHP]+ratraw[iR094]*ys[iDP]+ratraw[iR101]*ys[iHEP];
    
    dfdy[iHDP][iD2P]=  ratraw[iR088]*ys[iH];
    
    dfdy[iHDP][iD2]=  ratraw[iR098]*ys[iHP];
    
    dfdy[iHDP][iHEP]=  ratraw[iR101]*ys[iHD];
    
    dfdy[iHDP][iHE]=  0.0e0;
    
    dfdy[iHDP][iHEPP]=  0.0e0;
    
    //dfdy[iHDP][iELEC]=  -ratraw[iR044]*ys[iHDP];
    
    
    
    dfdy[iHD][iHP]=  -ratraw[iR041]*ys[iHD]-ratraw[iR092]*ys[iHD]-ratraw[iR093]*ys[iHD]+ratraw[iR097]*ys[iD2];
    
    dfdy[iHD][iH]=  ratraw[iR036]*ys[iD]+ratraw[iR038]*ys[iHDP]-ratraw[iR040]*ys[iHD]
    +ratraw[iR055]*ys[iDM]+ratraw[iR089]*ys[iD2P]+ratraw[iR107]*ys[iD2]
    -ratraw[iR108]*ys[iHD];
    
    dfdy[iHD][iHM]=  ratraw[iR054]*ys[iD]+ratraw[iR071]*ys[iHDP];
    
    dfdy[iHD][iH2P]=  ratraw[iR082]*ys[iD];
    
    dfdy[iHD][iH2]=  ratraw[iR037]*ys[iD]+ratraw[iR039]*ys[iDP]-ratraw[iR109]*ys[iHD];
    
    dfdy[iHD][iDP]=  ratraw[iR039]*ys[iH2]-ratraw[iR094]*ys[iHD]-ratraw[iR095]*ys[iHD]-ratraw[iR096]*ys[iHD];
    
    dfdy[iHD][iD]=  ratraw[iR036]*ys[iH]+ratraw[iR037]*ys[iH2]+ratraw[iR049]*ys[iHDP]
    +ratraw[iR054]*ys[iHM]+ratraw[iR082]*ys[iH2P]-ratraw[iR106]*ys[iHD];
    
    dfdy[iHD][iDM]=  ratraw[iR055]*ys[iH]+ratraw[iR073]*ys[iHDP];
    
    dfdy[iHD][iHDP]=  ratraw[iR038]*ys[iH]+ratraw[iR049]*ys[iD]+ratraw[iR071]*ys[iHM]+ratraw[iR073]*ys[iDM];
    
    dfdy[iHD][iHD]=  -ratraw[iR040]*ys[iH]-ratraw[iR041]*ys[iHP]-ratraw[iR057]*ys[iELEC]
    -ratraw[iR058]*ys[iELEC]-ratraw[iR092]*ys[iHP]-ratraw[iR093]*ys[iHP]
    -ratraw[iR094]*ys[iDP]-ratraw[iR095]*ys[iDP]-ratraw[iR096]*ys[iDP]
    -ratraw[iR101]*ys[iHEP]-ratraw[iR102]*ys[iHEP]-ratraw[iR103]*ys[iHEP]
    -ratraw[iR106]*ys[iD]-ratraw[iR108]*ys[iH]-ratraw[iR109]*ys[iH2]
    -ratraw[iR110]*ys[iHE]-ratraw[iR111]*ys[iELEC]-ratraw[iR123];
    
    dfdy[iHD][iD2P]=  ratraw[iR089]*ys[iH];
    
    dfdy[iHD][iD2]=  ratraw[iR097]*ys[iHP]+ratraw[iR107]*ys[iH];
    
    dfdy[iHD][iHEP]=  -ratraw[iR101]*ys[iHD]-ratraw[iR102]*ys[iHD]-ratraw[iR103]*ys[iHD];
    
    dfdy[iHD][iHE]=  -ratraw[iR110]*ys[iHD];
    
    dfdy[iHD][iHEPP]= 0.0e0;
    
    //dfdy[iHD][iELEC]=  -ratraw[iR057]*ys[iHD]-ratraw[iR058]*ys[iHD]-ratraw[iR111]*ys[iHD];
    
    
    
    dfdy[iD2P][iHP]=  ratraw[iR099]*ys[iD2];
    
    dfdy[iD2P][iH]=  -ratraw[iR087]*ys[iD2P]-ratraw[iR088]*ys[iD2P]-ratraw[iR089]*ys[iD2P];
    
    dfdy[iD2P][iHM]=  -ratraw[iR075]*ys[iD2P]-ratraw[iR076]*ys[iD2P];
    
    dfdy[iD2P][iH2P]= 0.0e0;
    
    dfdy[iD2P][iH2]=  0.0e0;
    
    dfdy[iD2P][iDP]=  ratraw[iR062]*ys[iDM]+ratraw[iR080]*ys[iD]+ratraw[iR096]*ys[iHD]+ratraw[iR100]*ys[iD2];
    
    dfdy[iD2P][iD]=  ratraw[iR080]*ys[iDP]+ratraw[iR084]*ys[iHDP]-ratraw[iR086]*ys[iD2P];
    
    dfdy[iD2P][iDM]=  ratraw[iR062]*ys[iDP]-ratraw[iR077]*ys[iD2P]-ratraw[iR078]*ys[iD2P];
    
    dfdy[iD2P][iHDP]=  ratraw[iR084]*ys[iD];
    
    dfdy[iD2P][iHD]=  ratraw[iR096]*ys[iDP];
    
    dfdy[iD2P][iD2P]=  -ratraw[iR075]*ys[iHM]-ratraw[iR076]*ys[iHM]-ratraw[iR077]*ys[iDM]
    -ratraw[iR078]*ys[iDM]-ratraw[iR086]*ys[iD]-ratraw[iR087]*ys[iH]
    -ratraw[iR088]*ys[iH]-ratraw[iR089]*ys[iH]-ratraw[iR121];
    
    dfdy[iD2P][iD2]=  ratraw[iR099]*ys[iHP]+ratraw[iR100]*ys[iDP]+ratraw[iR104]*ys[iHEP];
    
    dfdy[iD2P][iHEP]=  ratraw[iR104]*ys[iD2];
    
    dfdy[iD2P][iHE]=  0.0e0;
    
    dfdy[iD2P][iHEPP]=  0.0e0;
    
    //dfdy[iD2P][iELEC]=  0.0e0;
    
    
    
    dfdy[iD2][iHP]=  -ratraw[iR097]*ys[iD2]-ratraw[iR098]*ys[iD2]-ratraw[iR099]*ys[iD2];
    
    dfdy[iD2][iH]=  ratraw[iR087]*ys[iD2P]-ratraw[iR107]*ys[iD2]-ratraw[iR112]*ys[iD2];
    
    dfdy[iD2][iHM]=  ratraw[iR075]*ys[iD2P];
    
    dfdy[iD2][iH2P]= 0.0e0;
    
    dfdy[iD2][iH2]=  -ratraw[iR113]*ys[iD2];
    
    dfdy[iD2][iDP]=  ratraw[iR095]*ys[iHD]-ratraw[iR100]*ys[iD2];
    
    dfdy[iD2][iD]=  ratraw[iR056]*ys[iDM]+ratraw[iR085]*ys[iHDP]+ratraw[iR086]*ys[iD2P]+ratraw[iR106]*ys[iHD];
    
    dfdy[iD2][iDM]=  ratraw[iR056]*ys[iD]+ratraw[iR077]*ys[iD2P];
    
    dfdy[iD2][iHDP]=  ratraw[iR085]*ys[iD];
    
    dfdy[iD2][iHD]=  ratraw[iR095]*ys[iDP]+ratraw[iR106]*ys[iD];
    
    dfdy[iD2][iD2P]=  ratraw[iR075]*ys[iHM]+ratraw[iR077]*ys[iDM]+ratraw[iR086]*ys[iD]+ratraw[iR087]*ys[iH];
    
    dfdy[iD2][iD2]=  -ratraw[iR059]*ys[iELEC]-ratraw[iR097]*ys[iHP]-ratraw[iR098]*ys[iHP]
    -ratraw[iR099]*ys[iHP]-ratraw[iR100]*ys[iDP]-ratraw[iR104]*ys[iHEP]
    -ratraw[iR105]*ys[iHEP]-ratraw[iR107]*ys[iH]-ratraw[iR112]*ys[iH]
    -ratraw[iR113]*ys[iH2]-ratraw[iR114]*ys[iHE]-ratraw[iR115]*ys[iELEC]-ratraw[iR124];
    
    dfdy[iD2][iHEP]=  -ratraw[iR104]*ys[iD2]-ratraw[iR105]*ys[iD2];
    
    dfdy[iD2][iHE]=  -ratraw[iR114]*ys[iD2];
    
    dfdy[iD2][iHEPP]= 0.0e0;
    
    //dfdy[iD2][iELEC]=  -ratraw[iR059]*ys[iD2]-ratraw[iR115]*ys[iD2];
    
    
    
    dfdy[iHEP][iHP]=  ratraw[iR027]*ys[iHE];
    
    dfdy[iHEP][iH]=  -ratraw[iR026]*ys[iHEP];
    
    dfdy[iHEP][iHM]=  -ratraw[iR028]*ys[iHEP];
    
    dfdy[iHEP][iH2P]= 0.0e0;
    
    dfdy[iHEP][iH2]=  -ratraw[iR024]*ys[iHEP]-ratraw[iR025]*ys[iHEP];
    
    dfdy[iHEP][iDP]=  ratraw[iR047]*ys[iHE];
    
    dfdy[iHEP][iD]=  -ratraw[iR046]*ys[iHEP];
    
    dfdy[iHEP][iDM]=  -ratraw[iR079]*ys[iHEP];
    
    dfdy[iHEP][iHDP]= 0.0e0;
    
    dfdy[iHEP][iHD]=  -ratraw[iR101]*ys[iHEP]-ratraw[iR102]*ys[iHEP]-ratraw[iR103]*ys[iHEP];
    
    dfdy[iHEP][iD2P]=  0.0e0;
    
    dfdy[iHEP][iD2]=  -ratraw[iR104]*ys[iHEP]-ratraw[iR105]*ys[iHEP];
    
    dfdy[iHEP][iHEP]=  -ratraw[iR018]*ys[iELEC]-ratraw[iR019]*ys[iELEC]-ratraw[iR024]*ys[iH2]
    -ratraw[iR025]*ys[iH2]-ratraw[iR026]*ys[iH]-ratraw[iR028]*ys[iHM]
    -ratraw[iR046]*ys[iD]-ratraw[iR079]*ys[iDM]-ratraw[iR101]*ys[iHD]
    -ratraw[iR102]*ys[iHD]-ratraw[iR103]*ys[iHD]-ratraw[iR104]*ys[iD2]
    -ratraw[iR105]*ys[iD2];
    
    dfdy[iHEP][iHE]=  ratraw[iR017]*ys[iELEC]+ratraw[iR027]*ys[iHP]+ratraw[iR047]*ys[iDP];
    
    dfdy[iHEP][iHEPP]=  ratraw[iR020]*ys[iELEC];
    
    //dfdy[iHEP][iELEC]=  ratraw[iR017]*ys[iHE]-ratraw[iR018]*ys[iHEP]-ratraw[iR019]*ys[iHEP]
    //+ratraw[iR020]*ys[iHEPP];
    
    
    
    dfdy[iHE][iHP]=  -ratraw[iR027]*ys[iHE];
    
    dfdy[iHE][iH]=  ratraw[iR026]*ys[iHEP]-ratraw[iR032]*ys[iH]*ys[iHE]-ratraw[iR032]*ys[iH]*ys[iHE]
    +ratraw[iR032]*ys[iH]*ys[iHE]+ratraw[iR032]*ys[iH]*ys[iHE];
    
    dfdy[iHE][iHM]=  ratraw[iR028]*ys[iHEP]-ratraw[iR029]*ys[iHE]+ratraw[iR029]*ys[iHE];
    
    dfdy[iHE][iH2P]= 0.0e0;
    
    dfdy[iHE][iH2]=  -ratraw[iR011]*ys[iHE]+ratraw[iR011]*ys[iHE]+ratraw[iR024]*ys[iHEP]+ratraw[iR025]*ys[iHEP];
    
    dfdy[iHE][iDP]=  -ratraw[iR047]*ys[iHE];
    
    dfdy[iHE][iD]=  ratraw[iR046]*ys[iHEP];
    
    dfdy[iHE][iDM]=  -ratraw[iR065]*ys[iHE]+ratraw[iR065]*ys[iHE]+ratraw[iR079]*ys[iHEP];
    
    dfdy[iHE][iHDP]= 0.0e0;
    
    dfdy[iHE][iHD]=  ratraw[iR101]*ys[iHEP]+ratraw[iR102]*ys[iHEP]+ratraw[iR103]*ys[iHEP]
    -ratraw[iR110]*ys[iHE]+ratraw[iR110]*ys[iHE];
    
    dfdy[iHE][iD2P]= 0.0e0;
    
    dfdy[iHE][iD2]=  ratraw[iR104]*ys[iHEP]+ratraw[iR105]*ys[iHEP]-ratraw[iR114]*ys[iHE]+ratraw[iR114]*ys[iHE];
    
    dfdy[iHE][iHEP]=  ratraw[iR019]*ys[iELEC]+ratraw[iR024]*ys[iH2]+ratraw[iR025]*ys[iH2]
    +ratraw[iR026]*ys[iH]+ratraw[iR028]*ys[iHM]+ratraw[iR046]*ys[iD]
    +ratraw[iR079]*ys[iDM]+ratraw[iR101]*ys[iHD]+ratraw[iR102]*ys[iHD]
    +ratraw[iR103]*ys[iHD]+ratraw[iR104]*ys[iD2]+ratraw[iR105]*ys[iD2];
    
    dfdy[iHE][iHE]=  -ratraw[iR011]*ys[iH2]+ratraw[iR011]*ys[iH2]-ratraw[iR017]*ys[iELEC]
    -ratraw[iR027]*ys[iHP]-ratraw[iR029]*ys[iHM]+ratraw[iR029]*ys[iHM]
    -ratraw[iR032]*ys[iH]*ys[iH]+ratraw[iR032]*ys[iH]*ys[iH]
    -ratraw[iR047]*ys[iDP]-ratraw[iR065]*ys[iDM]+ratraw[iR065]*ys[iDM]
    -ratraw[iR110]*ys[iHD]+ratraw[iR110]*ys[iHD]-ratraw[iR114]*ys[iD2]
    +ratraw[iR114]*ys[iD2];
    
    dfdy[iHE][iHEPP]= 0.0e0;
    
    //dfdy[iHE][iELEC]=  -ratraw[iR017]*ys[iHE]+ratraw[iR019]*ys[iHEP];
    
    
    
    dfdy[iHEPP][iHP]=  0.0e0;
    
    dfdy[iHEPP][iH]=  0.0e0;
    
    dfdy[iHEPP][iHM]=  0.0e0;
    
    dfdy[iHEPP][iH2P]=  0.0e0;
    
    dfdy[iHEPP][iH2]=  0.0e0;
    
    dfdy[iHEPP][iDP]=  0.0e0;
    
    dfdy[iHEPP][iD]=  0.0e0;
    
    dfdy[iHEPP][iDM]=  0.0e0;
    
    dfdy[iHEPP][iHDP]=  0.0e0;
    
    dfdy[iHEPP][iHD]=  0.0e0;
    
    dfdy[iHEPP][iD2P]=  0.0e0;
    
    dfdy[iHEPP][iD2]=  0.0e0;
    
    dfdy[iHEPP][iHEP]=  ratraw[iR018]*ys[iELEC];
    
    dfdy[iHEPP][iHE]=  0.0e0;
    
    dfdy[iHEPP][iHEPP]=  -ratraw[iR020]*ys[iELEC];
  
#endif

#if DO_CHEMISTRY == CONTEMPORARY
    
/*
  dydt[iH] = -ys[iH] * ys[iELEC] * ratraw[iR001]
      - ys[iHM] * ys[iH] * ratraw[iR002]
      - ys[iH] * ys[iHP] * ratraw[iR003]
      - ys[iH] * ys[iH2P] * ratraw[iR004]
      + 2.0 * ys[iHM] * ys[iHP] * ratraw[iR005]
      + 2.0 * ys[iH2P] * ys[iELEC] * ratraw[iR006]
      + ys[iH2] * ys[iHP] * ratraw[iR007]
      + 2.0 * ys[iH2] * ys[iELEC] * ratraw[iR008]
      + 2.0 * ys[iH2] * ys[iH] * ratraw[iR009]
      + 2.0 * ys[iH2] * ys[iH2] * ratraw[iR010]
      - ys[iH] * ys[iELEC] * ratraw[iR011]
      + ys[iHP] * ys[iELEC] * ratraw[iR012]
      + y[iHM] * ys[iELEC] * ratraw[iR013]
      + ys[iHM] * ys[iH] * ratraw[iR014]
      - ys[iHEP] * ys[iH] * ratraw[iR018]
      + ys[iHE] * ys[iHP] * ratraw[iR019]
      + 2.0 * ys[iH2] * ys[iHE] * ratraw[iR020]
      + ys[iH2] * ys[iHEP] * ratraw[iR022]
      + ys[iHEP] * ys[iHM] * ratraw[iR023]
      + ys[iH3P] * ys[iELEC] * ratraw[iR024]
      + 3.0 * ys[iH3P] * ys[iELEC] * ratraw[iR025]
      - 2.0 * ys[iH] * ys[iH] * ys[iH] * ratraw[iR027]
      - 2.0 * ys[iH] * ys[iH] * ys[iH2] * ratraw[iR028] 
      - 2.0 * ys[iH] * ys[iH] * ys[iHE] * ratraw[iR029] 
      - 2.0 * ys[iH] * ratraw[iR030]
      + ys[iCP] * ys[iH2] * ratraw[iR035]
      + ys[iO] * ys[iCHX] * ratraw[iR037]
      + ys[iC] * ys[iOHX] * ratraw[iR038]
      + ys[iHCOP] * ys[iELEC] * ratraw[iR040]
      + ys[iH2P] * ys[iH2] * ratraw[iR043]
      - ys[iH3P] * ys[iH] * ratraw[iR044]
      + ys[iCHX] * ratraw[iP001]
      + ys[iOHX] * ratraw[iP002]
      + ys[iHCOP] * ratraw[iP004]
      + ys[iHM] * ratraw[iP005]
      + ys[iH2P] * ratraw[iP006]
      + 2.0 * ys[iH2] * ratraw[iP007]
      + ys[iH3P] * ratraw[iP009]
      - ys[iH] * ratraw[iP012]
      + 2.0 * ys[iH2] * ratraw[iP016] ;
*/

    
    dfdy[iH][iH] = - ys[iELEC] * ratraw[iR001]
      - ys[iHM] * ratraw[iR002]
      - ys[iHP] * ratraw[iR003]
      - ys[iH2P] * ratraw[iR004]
      + 2.0 * ys[iH2] * ratraw[iR009]
      - ys[iELEC] * ratraw[iR011]
      + ys[iHM] *  ratraw[iR014]
      - ys[iHEP] * ratraw[iR018]
      - 6.0 * ys[iH] * ys[iH] * ratraw[iR027]
      - 4.0 * ys[iH] * ys[iH2] * ratraw[iR028] 
      - 4.0 * ys[iH] * ys[iHE] * ratraw[iR029] 
      - 2.0 * ratraw[iR030]
      - ys[iH3P] * ratraw[iR044]
      - ratraw[iP012];
  

  dfdy[iH][iHM] =  - ys[iH] * ratraw[iR002]
    + 2.0 * ys[iHP] * ratraw[iR005]
    + ys[iELEC] * ratraw[iR013]
    + ys[iH] * ratraw[iR014]
    + ys[iHEP] * ratraw[iR023]
    +  ratraw[iP005];
 

  dfdy[iH][iH2] =  ys[iHP] * ratraw[iR007]
      + 2.0 * ys[iELEC] * ratraw[iR008]
      + 2.0 * ys[iH] * ratraw[iR009]
      + 4.0 * ys[iH2] * ratraw[iR010]
      + 2.0 * ys[iHE] * ratraw[iR020]
      + ys[iHEP] * ratraw[iR022]
      - 2.0 * ys[iH] * ys[iH] * ratraw[iR028] 
      + ys[iCP] * ratraw[iR035]
      + ys[iH2P] * ratraw[iR043]
      + 2.0 * ratraw[iP007]
      + 2.0 * ratraw[iP016] ;

  dfdy[iH][iHP] = - ys[iH] * ratraw[iR003]
    + 2.0 * ys[iHM] *ratraw[iR005]
    + ys[iH2] * ratraw[iR007]
    + ys[iELEC] * ratraw[iR012]
    + ys[iHE] * ratraw[iR019];

  dfdy[iH][iH2P] =  - ys[iH] * ratraw[iR004]
      + 2.0 * ys[iELEC] * ratraw[iR006]
      + ys[iH2] * ratraw[iR043]
      + ratraw[iP006];

  dfdy[iH][iH3P] = ys[iELEC] * ratraw[iR024]
      + 3.0 * ys[iELEC] * ratraw[iR025]
      - ys[iH] * ratraw[iR044]
      + ratraw[iP009] ;

  dfdy[iH][iHE] = ys[iHP] * ratraw[iR019]
    + 2.0 * ys[iH2] * ratraw[iR020]
    - 2.0 * ys[iH] * ys[iH] * ratraw[iR029];

  dfdy[iH][iHEP] = - ys[iH] * ratraw[iR018]
    + ys[iH2] * ratraw[iR022]
    + ys[iHM] * ratraw[iR023] ;

  dfdy[iH][iC] =  ys[iOHX] * ratraw[iR038];

  dfdy[iH][iCP] =  ys[iH2] * ratraw[iR035];

  dfdy[iH][iO] =  ys[iCHX] * ratraw[iR037];

  dfdy[iH][iCHX] =  + ys[iO] * ratraw[iR037]
    + ratraw[iP001];

  dfdy[iH][iOHX] =  ys[iC] * ratraw[iR038]
    + ratraw[iP002];

 dfdy[iH][iHCOP] = ys[iELEC] * ratraw[iR040]
   + ratraw[iP004];

 dfdy[iH][iCO] = 0.0;

 dfdy[iH][iM] = 0.0;

 dfdy[iH][iMP] = 0.0;


 /*
    dydt[iHM] = ys[iH] * ys[iELEC] * ratraw[iR001]
      - ys[iHM] * ys[iH] * ratraw[iR002]
      - ys[iHM] * ys[iHP] * ratraw[iR005]
      - ys[iHM] * ys[iELEC] * ratraw[iR013] 
      - ys[iHM] * ys[iH] * ratraw[iR014]
      - ys[iHM] * ys[iHP] * ratraw[iR015]
      - ys[iHEP] * ys[iHM] * ratraw[iR023]
      - ys[iHM] * ratraw[iP005] ;
 */

 dfdy[iHM][iH] =  ys[iELEC] * ratraw[iR001]
   - ys[iHM] * ratraw[iR002]
   - ys[iHM] * ratraw[iR014];

  dfdy[iHM][iHM] = -  ys[iH] * ratraw[iR002]
      - ys[iHP] * ratraw[iR005]
      - ys[iELEC] * ratraw[iR013] 
      - ys[iH] * ratraw[iR014]
      - ys[iHP] * ratraw[iR015]
      - ys[iHEP] * ratraw[iR023]
      - ratraw[iP005];

  dfdy[iHM][iH2] = 0.0;

  dfdy[iHM][iHP] = - ys[iHM] * ratraw[iR005]
    - ys[iHM] * ratraw[iR015];

  dfdy[iHM][iH2P] = 0.0;

  dfdy[iHM][iH3P] = 0.0;

  dfdy[iHM][iHE] = 0.0;

  dfdy[iHM][iHEP] =  -ys[iHM] * ratraw[iR023];

  dfdy[iHM][iC] = 0.0;

  dfdy[iHM][iCP] = 0.0;

  dfdy[iHM][iO] = 0.0;

  dfdy[iHM][iCHX] = 0.0;

  dfdy[iHM][iOHX] = 0.0;

  dfdy[iHM][iHCOP] = 0.0;

  dfdy[iHM][iCO] = 0.0;

 dfdy[iHM][iM] = 0.0;

 dfdy[iHM][iMP] = 0.0;

 /*
    dydt[iH2] = ys[iHM] * ys[iH] * ratraw[iR002]
      + ys[iH] * ys[iH2P] * ratraw[iR004] 
      - ys[iH2] * ys[iHP] * ratraw[iR007]
      - ys[iH2] * ys[iELEC] * ratraw[iR008] 
      - ys[iH2] * ys[iH] * ratraw[iR009]
      - ys[iH2] * ys[iH2] * ratraw[iR010]
      - ys[iH2] * ys[iHE] * ratraw[iR020]
      - ys[iH2] * ys[iHEP] * ratraw[iR021]
      - ys[iH2] * ys[iHEP] * ratraw[iR022]
      + ys[iH3P] * ys[iELEC] * ratraw[iR024]
      - ys[iH2] * ys[iHP] * ratraw[iR026]
      + ys[iH] * ys[iH] * ys[iH] * ratraw[iR027] 
      + ys[iH] * ys[iH] * ys[iH2] * ratraw[iR028]
      + ys[iH] * ys[iH] * ys[iHE] * ratraw[iR029]
      + ys[iH] * ratraw[iR030] // dust grain formation
      + ys[iH3P] * ys[iC] * ratraw[iR031] 
      + ys[iH3P] * ys[iO] * ratraw[iR032]
      + ys[iH3P] * ys[iCO] * ratraw[iR033] 
      - ys[iCP] * ys[iH2] * ratraw[iR035] 
      + ys[iH3P] * ys[iM] * ratraw[iR042]
      - ys[iH2P] * ys[iH2] * ratraw[iR043]
      + ys[iH3P] * ys[iH] * ratraw[iR044] 
      - ys[iH2] * ratraw[iP007]
      + ys[iH3P] * ratraw[iP008]
      - ys[iH2] * ratraw[iP015]
      - ys[iH2] * ratraw[iP016] ;
 */


 dfdy[iH2][iH] = ys[iHM] *ratraw[iR002]
      + ys[iH2P] * ratraw[iR004] 
      - ys[iH2] * ratraw[iR009]
      + 3.0 * ys[iH] * ys[iH] * ratraw[iR027] 
      + 2.0 * ys[iH] * ys[iH2] * ratraw[iR028]
      + 2.0 * ys[iH] * ys[iHE] * ratraw[iR029]
      + ratraw[iR030] // dust grain formation
      + ys[iH3P] * ratraw[iR044];

 dfdy[iH2][iHM] = ys[iH] * ratraw[iR002];

  dfdy[iH2][iH2] =  - ys[iHP] * ratraw[iR007]
      - ys[iELEC] * ratraw[iR008] 
      - ys[iH] * ratraw[iR009]
      - 2.0 * ys[iH2] * ratraw[iR010]
      - ys[iHE] * ratraw[iR020]
      - ys[iHEP] * ratraw[iR021]
      - ys[iHEP] * ratraw[iR022]
      - ys[iHP] * ratraw[iR026]
      + ys[iH] * ys[iH] * ratraw[iR028]
      - ys[iCP] *  ratraw[iR035] 
      - ys[iH2P] *  ratraw[iR043]
      - ratraw[iP007]
      - ratraw[iP015]
      - ratraw[iP016];

  dfdy[iH2][iHP] = - ys[iH2] * ratraw[iR007]
    - ys[iH2] * ratraw[iR026];

  dfdy[iH2][iH2P] = ys[iH] * ratraw[iR004] 
    - ys[iH2] * ratraw[iR043];

  dfdy[iH2][iH3P] = ys[iELEC] * ratraw[iR024]
      + ys[iC] * ratraw[iR031] 
      + ys[iO] * ratraw[iR032]
      + ys[iCO] * ratraw[iR033] 
      + ys[iM] * ratraw[iR042]
      + ys[iH] * ratraw[iR044] 
      + ratraw[iP008];

  dfdy[iH2][iHE] =  - ys[iH2] * ratraw[iR020]
    + ys[iH] * ys[iH] * ratraw[iR029];

  dfdy[iH2][iHEP] = - ys[iH2] * ratraw[iR021]
    - ys[iH2] * ratraw[iR022];

  dfdy[iH2][iC] = ys[iH3P] * ratraw[iR031]; 

  dfdy[iH2][iCP] = - ys[iH2] * ratraw[iR035];

  dfdy[iH2][iO] = ys[iH3P] * ratraw[iR032];

  dfdy[iH2][iCHX] = 0.0;

  dfdy[iH2][iOHX] = 0.0;

  dfdy[iH2][iHCOP] = 0.0;

  dfdy[iH2][iCO] = ys[iH3P] * ratraw[iR033] ;

  dfdy[iH2][iM] = ys[iH3P] * ratraw[iR042];

  dfdy[iH2][iMP] = 0.0;
    

  /*
    dydt[iHP] = - ys[iH] * ys[iHP] * ratraw[iR003]
      + ys[iH] * ys[iH2P] * ratraw[iR004] 
      - ys[iHM] * ys[iHP] * ratraw[iR005]
      - ys[iH2] * ys[iHP] * ratraw[iR007]
      + ys[iH] * ys[iELEC] * ratraw[iR011] 
      - ys[iHP] * ys[iELEC] * ratraw[iR012]
      - ys[iHM] * ys[iHP] * ratraw[iR015]
      + ys[iHEP] * ys[iH] * ratraw[iR018]
      - ys[iHE] * ys[iHP] * ratraw[iR019]
      + ys[iH2] * ys[iHEP] * ratraw[iR022]
      - ys[iH2] * ys[iHP] * ratraw[iR026]
      + ys[iH2P] * ratraw[iP006]
      + ys[iH3P] * ratraw[iP008]
      + ys[iH] * ratraw[iP012];
  */


  dfdy[iHP][iH] = - ys[iHP] * ratraw[iR003]
      + ys[iH2P] * ratraw[iR004] 
      + ys[iELEC] * ratraw[iR011] 
      + ys[iHEP] *  ratraw[iR018]
      + ratraw[iP012];

  dfdy[iHP][iHM] = - ys[iHP] * ratraw[iR005]
    - ys[iHP] * ratraw[iR015];

  dfdy[iHP][iH2] = - ys[iHP] * ratraw[iR007]
      + ys[iHEP] * ratraw[iR022]
      - ys[iHP] * ratraw[iR026];

  dfdy[iHP][iHP] = - ys[iH] * ratraw[iR003]
      - ys[iHM] * ratraw[iR005]
      - ys[iH2] * ratraw[iR007]
      - ys[iELEC] * ratraw[iR012]
      - ys[iHM] * ratraw[iR015]
      - ys[iHE] * ratraw[iR019]
    - ys[iH2] * ratraw[iR026];

  dfdy[iHP][iH2P] =  ys[iH] * ratraw[iR004] 
    + ratraw[iP006];

  dfdy[iHP][iH3P] = ratraw[iP008];

  dfdy[iHP][iHE] = - ys[iHP] * ratraw[iR019];

  dfdy[iHP][iHEP] = ys[iH] * ratraw[iR018]
    + ys[iH2] * ratraw[iR022];

    dfdy[iHP][iC] = 0.0;

  dfdy[iHP][iCP] = 0.0;

  dfdy[iHP][iO] = 0.0;

  dfdy[iHP][iCHX] = 0.0;

  dfdy[iHP][iOHX] = 0.0;

  dfdy[iHP][iHCOP] = 0.0;

  dfdy[iHP][iCO] = 0.0;

  dfdy[iHP][iM] = 0.0;

 dfdy[iHP][iMP] = 0.0;

 /*
    dydt[iH2P] = ys[iH] * ys[iHP] * ratraw[iR003]
      - ys[iH] * ys[iH2P] * ratraw[iR004]
      - ys[iH2P] * ys[iELEC] * ratraw[iR006]
      + ys[iH2] * ys[iHP] * ratraw[iR007]
      + ys[iHM] * ys[iHP] * ratraw[iR015]
      +ys[iH2] * ys[iHEP] * ratraw[iR021]
      - ys[iH2P] * ys[iH2] * ratraw[iR043]
      + ys[iH3P] * ys[iH] * ratraw[iR044]
      - ys[iH2P] * ratraw[iP006]
      + ys[iH3P] * ratraw[iP009]
      + ys[iH2] * ratraw[iP015] ;
 */


 dfdy[iH2P][iH] = ys[iHP] * ratraw[iR003]
      - ys[iH2P] * ratraw[iR004]
      + ys[iH3P] * ratraw[iR044];

 dfdy[iH2P][iHM] =  ys[iHP] * ratraw[iR015];

  dfdy[iH2P][iH2] = ys[iHP] * ratraw[iR007]
      + ys[iHEP] * ratraw[iR021]
      - ys[iH2P] * ratraw[iR043]
      + ratraw[iP015];

  dfdy[iH2P][iHP] = ys[iH] *ratraw[iR003]
      + ys[iH2] * ratraw[iR007]
    + ys[iHM] * ratraw[iR015];

  dfdy[iH2P][iH2P] = - ys[iH] * ratraw[iR004]
      - ys[iELEC] * ratraw[iR006]
      - ys[iH2] * ratraw[iR043]
      - ratraw[iP006];

  dfdy[iH2P][iH3P] = ys[iH] * ratraw[iR044]
    + ratraw[iP009];

  dfdy[iH2P][iHE] = 0.0;

  dfdy[iH2P][iHEP] = ys[iH2] * ratraw[iR021];

  dfdy[iH2P][iC] = 0.0;

  dfdy[iH2P][iCP] = 0.0;

  dfdy[iH2P][iO] = 0.0;

  dfdy[iH2P][iCHX] = 0.0;

  dfdy[iH2P][iOHX] = 0.0;

  dfdy[iH2P][iHCOP] = 0.0;

  dfdy[iH2P][iCO] = 0.0;

  dfdy[iH2P][iM] = 0.0;

 dfdy[iH2P][iMP] = 0.0;
  
 /*  
    dydt[iH3P] = - ys[iH3P] * ys[iELEC] * ratraw[iR024]
      - ys[iH3P] * ys[iELEC] * ratraw[iR025]
      + ys[iH2] * ys[iHP] * ratraw[iR026] 
      - ys[iH3P] * ys[iC] * ratraw[iR031]
      - ys[iH3P] * ys[iO] * ratraw[iR032]
      - ys[iH3P] * ys[iCO] * ratraw[iR033]
      - ys[iH3P] * ys[iM] * ratraw[iR042]
      - ys[iH3P] * ratraw[iP008]
      - ys[iH3P] * ratraw[iP009]
      + ys[iH2P] * ys[iH2] * ratraw[iR043]
      - ys[iH3P] * ys[iH] * ratraw[iR044] ;
 */


 dfdy[iH3P][iH] = - ys[iH3P] * ratraw[iR044] ;

 dfdy[iH3P][iHM] = 0.0;

  dfdy[iH3P][iH2] = ys[iHP] * ratraw[iR026] 
    + ys[iH2P] * ratraw[iR043];

  dfdy[iH3P][iHP] =  ys[iH2] * ratraw[iR026];

  dfdy[iH3P][iH2P] = ys[iH2] * ratraw[iR043];

  dfdy[iH3P][iH3P] =  -ys[iELEC] * ratraw[iR024]
      -  ys[iELEC] * ratraw[iR025]
      - ys[iC] * ratraw[iR031]
      - ys[iO] * ratraw[iR032]
      - ys[iCO] * ratraw[iR033]
      - ys[iM] * ratraw[iR042]
      - ratraw[iP008]
      - ratraw[iP009]
      - ys[iH] * ratraw[iR044];

  dfdy[iH3P][iHE] = 0.0;

  dfdy[iH3P][iHEP] = 0.0;

  dfdy[iH3P][iC] = - ys[iH3P] * ratraw[iR031];

  dfdy[iH3P][iCP] = 0.0;

  dfdy[iH3P][iO] = - ys[iH3P] * ratraw[iR032];

  dfdy[iH3P][iCHX] =  0.0;

  dfdy[iH3P][iOHX] = 0.0;

  dfdy[iH3P][iHCOP] = 0.0;

  dfdy[iH3P][iCO] =  - ys[iH3P] * ratraw[iR033];

  dfdy[iH3P][iM] = - ys[iH3P] * ratraw[iR042];

  dfdy[iH3P][iMP] = 0.0;
    

  /*
    dydt[iHE] = -ys[iHE] * ys[iELEC] * ratraw[iR016]
      + ys[iHEP] * ys[iELEC] * ratraw[iR017]
      + ys[iHEP] * ys[iH] * ratraw[iR018] 
      - ys[iHE] * ys[iHP] * ratraw[iR019]
      + ys[iH2] * ys[iHEP] * ratraw[iR021]
      +ys[iH2] * ys[iHEP] * ratraw[iR022]
      + ys[iHEP] * ys[iHM] * ratraw[iR023]
      + ys[iHEP] * ys[iCO] * ratraw[iR034]
      - ys[iHE] * ratraw[iP013] ;
  */


  dfdy[iHE][iH] =  ys[iHEP] * ratraw[iR018];

  dfdy[iHE][iHM] = ys[iHEP] * ratraw[iR023];

  dfdy[iHE][iH2] = ys[iHEP] * ratraw[iR021]
    + ys[iHEP] * ratraw[iR022];

  dfdy[iHE][iHP] = - ys[iHE] * ratraw[iR019];

  dfdy[iHE][iH2P] = 0.0;

  dfdy[iHE][iH3P] = 0.0;

  dfdy[iHE][iHE] = - ys[iELEC] * ratraw[iR016]
      - ys[iHP] * ratraw[iR019]
      - ratraw[iP013] ;

  dfdy[iHE][iHEP] = ys[iELEC] * ratraw[iR017]
      + ys[iH] * ratraw[iR018] 
      + ys[iH2] * ratraw[iR021]
      + ys[iH2] * ratraw[iR022]
      + ys[iHM] * ratraw[iR023]
      + ys[iCO] * ratraw[iR034];

  dfdy[iHE][iC] = 0.0;

  dfdy[iHE][iCP] = 0.0;

  dfdy[iHE][iO] = 0.0;

  dfdy[iHE][iCHX] = 0.0;

  dfdy[iHE][iOHX] = 0.0;

  dfdy[iHE][iHCOP] = 0.0;

  dfdy[iHE][iCO] = ys[iHEP] * ratraw[iR034];

  dfdy[iHE][iM] = 0.0;

  dfdy[iHE][iMP] = 0.0;
    
  /*
    dydt[iHEP] = ys[iHE] * ys[iELEC] * ratraw[iR016] 
      - ys[iHEP] * ys[iELEC] * ratraw[iR017]
      - ys[iHEP] * ys[iH] * ratraw[iR018]
      + ys[iHE] * ys[iHP] * ratraw[iR019] 
      - ys[iH2] * ys[iHEP] * ratraw[iR021]
      - ys[iH2] * ys[iHEP] * ratraw[iR022] 
      - ys[iHEP] * ys[iHM] * ratraw[iR023] 
      - ys[iHEP] * ys[iCO] * ratraw[iR034]
      + ys[iHE] * ratraw[iP013] ;
  */


  dfdy[iHEP][iH] = - ys[iHEP] * ratraw[iR018];

  dfdy[iHEP][iHM] = - ys[iHEP] * ratraw[iR023] ;

  dfdy[iHEP][iH2] = - ys[iHEP] * ratraw[iR021]
    - ys[iHEP] * ratraw[iR022];

  dfdy[iHEP][iHP] =  ys[iHE] * ratraw[iR019] ;

  dfdy[iHEP][iH2P] = 0.0;

  dfdy[iHEP][iH3P] = 0.0;

  dfdy[iHEP][iHE] =  ys[iELEC] * ratraw[iR016] 
      + ys[iHP] * ratraw[iR019] 
      + ratraw[iP013];

  dfdy[iHEP][iHEP] = -  ys[iELEC] * ratraw[iR017]
      - ys[iH] * ratraw[iR018]
      - ys[iH2] * ratraw[iR021]
      - ys[iH2] * ratraw[iR022] 
      - ys[iHM] * ratraw[iR023] 
    - ys[iCO] * ratraw[iR034];

  dfdy[iHEP][iC] = 0.0;

  dfdy[iHEP][iCP] = 0.0;

  dfdy[iHEP][iO] = 0.0;

  dfdy[iHEP][iCHX] = 0.0;

  dfdy[iHEP][iOHX] = 0.0;

  dfdy[iHEP][iHCOP] = 0.0;

  dfdy[iHEP][iCO] = - ys[iHEP] * ratraw[iR034];

  dfdy[iHEP][iM] = 0.0;

  dfdy[iHEP][iMP] = 0.0;

  /*    
    dydt[iC] = - ys[iH3P] * ys[iC] * ratraw[iR031]
      - ys[iC] * ys[iOHX] * ratraw[iR038]
      + ys[iCP] * ys[iELEC] * ratraw[iR039]
      + ys[iCHX] * ratraw[iP001]
      - ys[iC] * ratraw[iP010]
      + ys[iCO] * ratraw[iP011]
      - ys[iC] * ratraw[iP014] ;
  */

  dfdy[iC][iH] = 0.0;

  dfdy[iC][iHM] = 0.0;

  dfdy[iC][iH2] = 0.0;

  dfdy[iC][iHP] = 0.0;

  dfdy[iC][iH2P] = 0.0;

  dfdy[iC][iH3P] = - ys[iC] * ratraw[iR031];

  dfdy[iC][iHE] = 0.0;

  dfdy[iC][iHEP] = 0.0;

  dfdy[iC][iC] = - ys[iH3P] *  ratraw[iR031]
      - ys[iOHX] * ratraw[iR038]
      - ratraw[iP010]
      - ratraw[iP014];

  dfdy[iC][iCP] = ys[iELEC] * ratraw[iR039];

  dfdy[iC][iO] = 0.0;

  dfdy[iC][iCHX] = ratraw[iP001];

  dfdy[iC][iOHX] = - ys[iC] * ratraw[iR038];

  dfdy[iC][iHCOP] = 0.0;

  dfdy[iC][iCO] = ratraw[iP011];

  dfdy[iC][iM] = 0.0;

  dfdy[iC][iMP] = 0.0;
    

  /*
    dydt[iCP] =  ys[iHEP] * ys[iCO] * ratraw[iR034] 
      - ys[iCP] * ys[iH2] * ratraw[iR035] 
      - ys[iCP] * ys[iOHX] * ratraw[iR036] 
      - ys[iCP] * ys[iELEC] * ratraw[iR039]
      + ys[iC] * ratraw[iP010]
      + ys[iC] * ratraw[iP014] ;
  */

  dfdy[iCP][iH] = 0.0;

  dfdy[iCP][iHM] = 0.0;

  dfdy[iCP][iH2] = - ys[iCP] * ratraw[iR035] ;

  dfdy[iCP][iHP] = 0.0;

  dfdy[iCP][iH2P] = 0.0;

  dfdy[iCP][iH3P] = 0.0;

  dfdy[iCP][iHE] = 0.0;

  dfdy[iCP][iHEP] = ys[iCO] * ratraw[iR034] ;

  dfdy[iCP][iC] =  ratraw[iP010]
      + ratraw[iP014] ;

  dfdy[iCP][iCP] = - ys[iH2] * ratraw[iR035] 
      - ys[iOHX] * ratraw[iR036] 
    - ys[iELEC] * ratraw[iR039];

  dfdy[iCP][iO] = 0.0;

  dfdy[iCP][iCHX] = 0.0;

  dfdy[iCP][iOHX] = - ys[iCP] * ratraw[iR036] ;

  dfdy[iCP][iHCOP] = 0.0;

  dfdy[iCP][iCO] =  ys[iHEP] * ratraw[iR034] ;;

  dfdy[iCP][iM] = 0.0;

  dfdy[iCP][iMP] = 0.0;

  /*
    dydt[iO] = - ys[iH3P] * ys[iO] * ratraw[iR032] 
      + ys[iHEP] * ys[iCO] * ratraw[iR034]
      - ys[iO] * ys[iCHX] * ratraw[iR037]
      + ys[iOHX] * ratraw[iP002]
      + ys[iCO] * ratraw[iP011] ;
  */


  dfdy[iO][iH] = 0.0;

  dfdy[iO][iHM] = 0.0;

  dfdy[iO][iH2] = 0.0;

  dfdy[iO][iHP] = 0.0;

  dfdy[iO][iH2P] = 0.0;

  dfdy[iO][iH3P] = - ys[iO] * ratraw[iR032] ;

  dfdy[iO][iHE] = 0.0;

  dfdy[iO][iHEP] = ys[iCO] * ratraw[iR034];

  dfdy[iO][iC] = 0.0;

  dfdy[iO][iCP] = 0.0;

  dfdy[iO][iO] = - ys[iH3P] *  ratraw[iR032] 
    - ys[iCHX] * ratraw[iR037];

  dfdy[iO][iCHX] = - ys[iO] * ratraw[iR037];

  dfdy[iO][iOHX] = ratraw[iP002];

  dfdy[iO][iHCOP] = 0.0;

  dfdy[iO][iCO] = ys[iHEP] * ratraw[iR034]
    + ratraw[iP011];
  
  dfdy[iO][iM] = 0.0;
  
  dfdy[iO][iMP] = 0.0;
    

  /*
    dydt[iCHX] = ys[iH3P] * ys[iC] * ratraw[iR031] 
      + ys[iCP] * ys[iH2] * ratraw[iR035] 
      - ys[iO] * ys[iCHX] * ratraw[iR037] 
      - ys[iCHX] * ratraw[iP001] ;
  */


  dfdy[iCHX][iH] = 0.0;

  dfdy[iCHX][iHM] = 0.0;

  dfdy[iCHX][iH2] = ys[iCP] * ratraw[iR035] ;

  dfdy[iCHX][iHP] = 0.0;

  dfdy[iCHX][iH2P] = 0.0;

  dfdy[iCHX][iH3P] = ys[iC] * ratraw[iR031] ;

  dfdy[iCHX][iHE] = 0.0;

  dfdy[iCHX][iHEP] = 0.0;

  dfdy[iCHX][iC] = ys[iH3P] * ratraw[iR031] ;

  dfdy[iCHX][iCP] = ys[iH2] * ratraw[iR035] ;

  dfdy[iCHX][iO] = - ys[iCHX] * ratraw[iR037] ;

  dfdy[iCHX][iCHX] = - ys[iO] * ratraw[iR037] 
      - ratraw[iP001];

  dfdy[iCHX][iOHX] = 0.0;

 dfdy[iCHX][iHCOP] = 0.0;

 dfdy[iCHX][iCO] = 0.0;

 dfdy[iCHX][iM] = 0.0;

 dfdy[iCHX][iMP] = 0.0;


 /*
    dydt[iOHX] = ys[iH3P] * ys[iO] * ratraw[iR032]
      - ys[iCP] * ys[iOHX] * ratraw[iR036]
      - ys[iC] * ys[iOHX] * ratraw[iR038] 
      - ys[iOHX] * ratraw[iP002] ; 
 */


  dfdy[iOHX][iH] = 0.0;

  dfdy[iOHX][iHM] = 0.0;

  dfdy[iOHX][iH2] = 0.0;

  dfdy[iOHX][iHP] = 0.0;

  dfdy[iOHX][iH2P] = 0.0;

  dfdy[iOHX][iH3P] = ys[iO] * ratraw[iR032];

  dfdy[iOHX][iHE] = 0.0;

  dfdy[iOHX][iHEP] = 0.0;

  dfdy[iOHX][iC] = - ys[iOHX] * ratraw[iR038];

  dfdy[iOHX][iCP] = - ys[iOHX] * ratraw[iR036];

  dfdy[iOHX][iO] = ys[iH3P] * ratraw[iR032];

  dfdy[iOHX][iCHX] = 0.0;

  dfdy[iOHX][iOHX] = - ys[iCP] * ratraw[iR036]
      - ys[iC] *  ratraw[iR038] 
      - ratraw[iP002];

 dfdy[iOHX][iHCOP] = 0.0;

 dfdy[iOHX][iCO] = 0.0;

 dfdy[iOHX][iM] = 0.0;

 dfdy[iOHX][iMP] = 0.0;
    

 /*
    dydt[iHCOP] = ys[iH3P] * ys[iCO] * ratraw[iR033]
      + ys[iCP] * ys[iOHX] * ratraw[iR036]
      - ys[iHCOP] * ys[iELEC] * ratraw[iR040]
      - ys[iHCOP] * ratraw[iP004] ;
 */


 dfdy[iHCOP][iH] = 0.0 ;

  dfdy[iHCOP][iHM] = 0.0;

  dfdy[iHCOP][iH2] = 0.0;

  dfdy[iHCOP][iHP] = 0.0;

  dfdy[iHCOP][iH2P] = 0.0;

  dfdy[iHCOP][iH3P] = ys[iCO] * ratraw[iR033];

  dfdy[iHCOP][iHE] = 0.0;

  dfdy[iHCOP][iHEP] = 0.0;

  dfdy[iHCOP][iC] = 0.0;

  dfdy[iHCOP][iCP] = ys[iOHX] * ratraw[iR036];

  dfdy[iHCOP][iO] = 0.0;

  dfdy[iHCOP][iCHX] = 0.0;

  dfdy[iHCOP][iOHX] = ys[iCP] * ratraw[iR036];

 dfdy[iHCOP][iHCOP] = - ys[iELEC] * ratraw[iR040]
      - ratraw[iP004];

 dfdy[iHCOP][iCO] = ys[iH3P] * ratraw[iR033];

 dfdy[iHCOP][iM] = 0.0;

 dfdy[iHCOP][iMP] = 0.0;
    

 /*
    dydt[iCO] = - ys[iH3P] * ys[iCO] * ratraw[iR033] 
      - ys[iHEP] * ys[iCO] * ratraw[iR034]
      + ys[iO] * ys[iCHX] * ratraw[iR037]
      + ys[iC] * ys[iOHX] * ratraw[iR038]
      + ys[iHCOP] * ys[iELEC] * ratraw[iR040]
      + ys[iHCOP] * ratraw[iP004]
      - ys[iCO] * ratraw[iP011] ;
 */


 dfdy[iCO][iH] = 0.0;

  dfdy[iCO][iHM] = 0.0;

  dfdy[iCO][iH2] = 0.0;

  dfdy[iCO][iHP] = 0.0;

  dfdy[iCO][iH2P] = 0.0;

  dfdy[iCO][iH3P] = - ys[iCO] * ratraw[iR033] ;

  dfdy[iCO][iHE] = 0.0;

  dfdy[iCO][iHEP] = - ys[iCO] * ratraw[iR034];

  dfdy[iCO][iC] = ys[iOHX] * ratraw[iR038];

  dfdy[iCO][iCP] = 0.0;

  dfdy[iCO][iO] = ys[iCHX] * ratraw[iR037];

  dfdy[iCO][iCHX] =  ys[iO] * ratraw[iR037];

  dfdy[iCO][iOHX] =  ys[iC] * ratraw[iR038];

  dfdy[iCO][iHCOP] =  ys[iELEC] * ratraw[iR040]
    + ratraw[iP004];

  dfdy[iCO][iCO] =- ys[iH3P] * ratraw[iR033] 
      - ys[iHEP] * ratraw[iR034]
      - ratraw[iP011] ;

  dfdy[iCO][iM] = 0.0;

  dfdy[iCO][iMP] = 0.0;
    
  /*
    dydt[iM] =  ys[iMP] * ys[iELEC] * ratraw[iR041]
      - ys[iH3P] * ys[iM] * ratraw[iR042] 
      - ys[iM] * ratraw[iP003] ;
  */


  dfdy[iM][iH] = 0.0;

  dfdy[iM][iHM] = 0.0;

  dfdy[iM][iH2] = 0.0;

  dfdy[iM][iHP] = 0.0;

  dfdy[iM][iH2P] = 0.0;

  dfdy[iM][iH3P] = - ys[iM] * ratraw[iR042] ;

  dfdy[iM][iHE] = 0.0;

  dfdy[iM][iHEP] = 0.0;

  dfdy[iM][iC] = 0.0;

  dfdy[iM][iCP] = 0.0;

  dfdy[iM][iO] = 0.0;

  dfdy[iM][iCHX] = 0.0;

  dfdy[iM][iOHX] = 0.0;

  dfdy[iM][iHCOP] = 0.0;

  dfdy[iM][iCO] = 0.0;

  dfdy[iM][iM] =  - ys[iH3P] * ratraw[iR042] 
      - ratraw[iP003];

  dfdy[iM][iMP] = ys[iELEC] * ratraw[iR041];
    

  /*
    dydt[iMP] = - ys[iMP] * ys[iELEC] * ratraw[iR041]
      + ys[iH3P] * ys[iM] * ratraw[iR042]
      + ys[iM] * ratraw[iP003];
  */


  dfdy[iMP][iH] = 0.0;

  dfdy[iMP][iHM] = 0.0;

  dfdy[iMP][iH2] = 0.0;

  dfdy[iMP][iHP] = 0.0 ;

  dfdy[iMP][iH2P] = 0.0;

  dfdy[iMP][iH3P] = ys[iM] * ratraw[iR042];

  dfdy[iMP][iHE] = 0.0;

  dfdy[iMP][iHEP] = 0.0;

  dfdy[iMP][iC] =0.0 ;

  dfdy[iMP][iCP] =0.0 ;

  dfdy[iMP][iO] = 0.0;

  dfdy[iMP][iCHX] = 0.0;

  dfdy[iMP][iOHX] = 0.0;

  dfdy[iMP][iHCOP] = 0.0;

  dfdy[iMP][iCO] = 0.0;

  dfdy[iMP][iM] =  ys[iH3P] * ratraw[iR042]
      + ratraw[iP003];
 
  dfdy[iMP][iMP] = - ys[iELEC] * ratraw[iR041];

  
#endif

    //cout << "Jac start" << endl;
    //for(i=0;i<NSPECIES;i++)for(j=0;j<NSPECIES;j++) cout << setprecision(15) <<dfdy[i][j] << endl;
    //cout << "Jac end" << endl;
    
};
