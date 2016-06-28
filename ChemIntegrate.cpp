#include "ChemIntegrate.H"

using namespace ChemNS;



void ChemIntegrate(REAL start,REAL firstDT,REAL tmin,REAL dt,REAL* eosIn, REAL* ymass, void (*derivs)(REAL,REAL*,REAL*,REAL*, double* ), double* rpar, REAL* RatRaw, int &jcounts)
{
    int i;
    
    //REAL dxsav = 1;
    //int kmax = 1;
    REAL odescal = 1.0e-6;
    int stepMax = 1000000, stepCur = 0;
    
    
    int nok = 0;
    int nbad = 0;
    
    // Store boundary conditions
    REAL y[NINT];
    REAL yscal[NINT];
    for(i=0;i<NINT;i++) y[i] = ymass[i];
    
    REAL dydx[NINT];
    
    REAL x = start;
    REAL h = firstDT, hdid = 0, hnext=0;
    
    for(stepCur=1;stepCur < stepMax ; stepCur++)
    {
     
        
        // Yet another check that everything is positive
        for(i=0;i<NINT;i++) y[i] = max(y[i],1.0e-20);
        
        // Update the derivatives
        //cout << "About to enter loop derivs" << endl;
        derivs(h,dydx,RatRaw,y,rpar);
        
        // Used in error checking 
        for(i=0;i<NINT;i++) yscal[i] = max(odescal,fabs(y[i]));
        
        // Check to make sure we're not overshooting the timestep
        if( (x+h-dt)*(x+h-start) > 0.0 ) h = dt-x;
        
        hdid = 0; hnext = 0;
        
        // Step away! (y,x,hdid,hnext,jcounts updated)
        StepMe(y,dydx,x,h,hdid,hnext,yscal,derivs,jcounts,RatRaw,rpar);
        
            
        // Tally up
        if( hdid == h ) nok = nok+1;
        else nbad = nbad + 1;
        
        
        // Checks
        if(fabs(h) < tmin)
        {
            // print an error message
            cout << "ERROR: with fabs(h) < tmin" << endl;
        }
        
        if(stepCur==stepMax-1)
        {
            // Exceed max number of steps...
            cout << "ERROR: Exceeded max number of steps" << endl;
        }
        
        // Are we done?
        if( stepCur==stepMax-1 || (x-dt)*(dt-start) >= 0)
        {
            for(i=0;i<NINT;i++) ymass[i] = y[i];
            break;
        }
        
        
        // If not, get ready to start over
        h = hnext;
        
        
    } // end of stepping loop
    
    //cout << "jcount = " << jcounts << " nok = " << nok << " and nbad = " << nbad << endl;
    return;
    
};
