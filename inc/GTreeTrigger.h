#ifndef __GTreeTrigger_h__
#define __GTreeTrigger_h__


#include "Rtypes.h"
#include "GTree.h"


#define GTreeTrigger_MAX 256


class  GTreeTrigger : public GTree
{
private:
    Double_t 	energySum;
    Int_t 		multiplicity;
    Int_t 		nTriggerPattern;
    Int_t 		triggerPattern[GTreeTrigger_MAX];
    Int_t 		nErrors;
    Int_t 		errorModuleID[GTreeTrigger_MAX];
    Int_t 		errorModuleIndex[GTreeTrigger_MAX];
    Int_t 		errorCode[GTreeTrigger_MAX];
    Bool_t      helicity;
    Long64_t    MC_evt_id;
    Long64_t    MC_rnd_id;
    Bool_t      hasHelicity;
    Bool_t      hasMCID;

   protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeTrigger(GTreeManager *Manager);
    virtual ~GTreeTrigger();


    virtual void        Clear()                       {nTriggerPattern = 0; nErrors = 0;}
            Int_t 		GetMultiplicity()       const {return multiplicity;}
            Double_t	GetEnergySum()          const {return energySum;}
            Int_t		GetNTriggerPattern()    const {return nTriggerPattern;}
    const   Int_t*		GetTriggerPattern()     const {return triggerPattern;}
            Int_t		GetNErrors()            const {return nErrors;}
    const   Int_t*		GetErrorModuleID()      const {return errorModuleID;}
    const   Int_t*		GetErrorModuleIndex()   const {return errorModuleIndex;}
    const   Int_t*		GetErrorCode()          const {return errorCode;}
            Bool_t 	    GetHelicity()    	    const {return helicity;}
            Bool_t 	    HasHelicity()    	    const {return hasHelicity;}

};

#endif
