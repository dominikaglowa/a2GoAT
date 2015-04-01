#ifndef __GParticleReconstruction_h__
#define __GParticleReconstruction_h__

#include <TCutG.h>


#include "GDataChecks.h"

#define	PDG_ROOTINO 0

class	GParticleReconstruction : public GDataChecks
{
public:
    enum ReconstructType
    {
        ReconstructNone         = 0,
        ReconstructCutProton    = 1,
        ReconstructCutPion      = 2,
        ReconstructCutElectron  = 4,
        ReconstructTimeOfFlight = 8,
        ReconstructClusterSize  = 16,
        ReconstructAllPhotons   = 32,
        ReconstructAllProtons   = 64
    };

private:
	std::string config;

    ReconstructType   typeCB;
    ReconstructType   typeTAPS;

    TCutG* 		cutProtonCB;
    TCutG* 		cutPionCB;
    TCutG*		cutElectronCB;
    TCutG* 		cutProtonTAPS;
    TCutG* 		cutPionTAPS;
    TCutG*		cutElectronTAPS;

    TCutG* 		cutTimeOfFlightCB;
    TCutG* 		cutClusterSizeCB;
    TCutG* 		cutTimeOfFlightTAPS;
    TCutG* 		cutClusterSizeTAPS;

    Double_t	chargedThetaMin;
    Double_t	chargedThetaMax;

    Int_t* 		identified;
    Int_t* 		charge;
    Int_t* 		hadron;

    Bool_t 		chargeIgnorePID;
    Bool_t 		chargeIgnoreMWPC0;
    Bool_t 		chargeIgnoreMWPC1;
    Bool_t 		chargeIgnoreVETO;

    Double_t    timeCutCB[2];
    Double_t    timeCutTAPS[2];

    Bool_t      doScalerCorrection;
    Bool_t      doTrigger;
    Double_t    energySum;
    Int_t       multiplicity;

    TCutG*	OpenCutFile(Char_t* filename, Char_t* cutname);
    Bool_t  Trigger();

protected:
            Bool_t  ProcessEventWithoutFilling();
    virtual void    ProcessEvent();
    virtual Bool_t  Start();

public:
    GParticleReconstruction();
    virtual ~GParticleReconstruction();

    Bool_t	Init();
    void    SetCBTimeCut(const Double_t min, const Double_t max)    {timeCutCB[0]=min; timeCutCB[1]=max;}
    void    SetTAPSTimeCut(const Double_t min, const Double_t max)  {timeCutTAPS[0]=min; timeCutTAPS[1]=max;}
    void    SetScalerCorrection(const Bool_t value)                 {doScalerCorrection = value;}
    void    SetTrigger(const Double_t esum, const Int_t mult)       {doTrigger = kTRUE; energySum = esum; multiplicity = mult;}
    void    SetThetaRangeChargedParticles(const Double_t min, const Double_t max)  {chargedThetaMin = min; chargedThetaMax = max;}
};

#endif
