
#include "stdlib.h"
//#include "pluto.h"

#define real double
#define REAL double

#define PRIMORDIAL 3
#define CONTEMPORARY 2
#define DO_CHEMISTRY PRIMORDIAL
//#define DO_CHEMISTRY CONTEMPORARY

//ARS defining a few basic constants
#define mh  1.67e-24
#define eV 1.60219e-12
#define kboltz 1.3806e-16

#if DO_CHEMISTRY == PRIMORDIAL

#define MRHE2H 0.24
#define MRH    0.76
#define DeutAbund 2.6e-5
#define abhe ((1. - MRH)/4./MRH)
#define h2conv (2.0/1.3158)
#define heconv (4.0/1.3158)

// Total number of tracked species:
#define NSPECIES 16
// Total number of *integrated* species:
#define NINT 15
#define NREACTION 120


#endif // #if NETWORK == PRIMORDIAL


#if DO_CHEMISTRY == CONTEMPORARY

#define NSPECIES      18
#define NINT          17
#define NCREACTION   200
#define NPREACTION   30
#define NREACTION    NCREACTION+NPREACTION

#define YMETALS     1.0e-7
#define YOXYGEN     3.2e-4
#define YCARBON     1.4e-4
#define YHELIUM     0.08


#endif


void ChemDriver(double, double&, double&, double*, double*);

namespace ChemNS {
    
  // non tracked species must come last
    enum Species
    { 
#if DO_CHEMISTRY == PRIMORDIAL
        iHP   = 0,
        iH    = 1,
        iHM   = 2,
        iH2P  = 3,
        iH2   = 4,
        iDP   = 5,
        iD    = 6,
        iDM   = 7,
        iHDP  = 8,
        iHD   = 9,
        iD2P  = 10,
        iD2   = 11,
        iHEP  = 12,
        iHE   = 13,
        iHEPP = 14,
        iELEC = 15
#endif


#if DO_CHEMISTRY == CONTEMPORARY
	iH   = 0,
        iHM    = 1,
        iH2   = 2,
	iHP   = 3,
        iH2P  = 4,
        iH3P   = 5,
        iHE   = 6,
        iHEP    = 7,
        iC   = 8,
        iCP  = 9,
        iO   = 10,
        iCHX  = 11,
        iOHX   = 12,
        iHCOP  = 13,
	iCO  = 14,
        iM   = 15,
        iMP = 16,
	iELEC = 17
#endif

    };
    
    enum Reactions{
        
        iR001 = 0, iR002, iR003, iR004, iR005, iR006, iR007, iR008, iR009, iR010,
        iR011, iR012, iR013, iR014, iR015, iR016, iR017, iR018, iR019, iR020,
        iR021, iR022, iR023, iR024, iR025, iR026, iR027, iR028, iR029, iR030,
        iR031, iR032, iR033, iR034, iR035, iR036, iR037, iR038, iR039, iR040,
        iR041, iR042, iR043, iR044, iR045, iR046, iR047, iR048, iR049, iR050,
        iR051, iR052, iR053, iR054, iR055, iR056, iR057, iR058, iR059, iR060,
        iR061, iR062, iR063, iR064, iR065, iR066, iR067, iR068, iR069, iR070,
        iR071, iR072, iR073, iR074, iR075, iR076, iR077, iR078, iR079, iR080,
        iR081, iR082, iR083, iR084, iR085, iR086, iR087, iR088, iR089, iR090,
        iR091, iR092, iR093, iR094, iR095, iR096, iR097, iR098, iR099, iR100,
        iR101, iR102, iR103, iR104, iR105, iR106, iR107, iR108, iR109, iR110,
        iR111, iR112, iR113, iR114, iR115, iR116, iR117, iR118, iR119, iR120,
        iR121, iR122, iR123, iR124, iR125, iR126, iR127, iR128, iR129, iR130, 
	iR131, iR132, iR133, iR134, iR135, iR136, iR137, iR138, iR139, iR140, 
	iR141, iR142, iR143, iR144, iR145, iR146, iR147, iR148, iR149, iR150, 
	iR151, iR152, iR153, iR154, iR155, iR156, iR157, iR158, iR159, iR160, 
	iR161, iR162, iR163, iR164, iR165, iR166, iR167, iR168, iR169, iR170, 
	iR171, iR172, iR173, iR174, iR175, iR176, iR177, iR178, iR179, iR180, 
	iR181, iR182, iR183, iR184, iR185, iR186, iR187, iR188, iR189, iR190, 
	iR191, iR192, iR193, iR194, iR195, iR196, iR197, iR198, iR199, iR200, 
    };

#if DO_CHEMISTRY == CONTEMPORARY
  enum PReactions{
    
    iP001 = NCREACTION+1, iP002, iP003, iP004, iP005, iP006, iP007, iP008, iP009, iP010,
    iP011, iP012, iP013, iP014, iP015, iP016, iP017, iP018, iP019, iP020,

  };
#endif
    
}



