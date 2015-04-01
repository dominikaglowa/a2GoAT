#ifndef __GTreeEventParameters_h__
#define __GTreeEventParameters_h__


#include "Rtypes.h"
#include "GTree.h"


#define GTreeEventParameters_MAX 128


class  GTreeEventParameters : public GTree
{
private:
    Int_t     eventNumber;
    Int_t     nReconstructed;

   protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeEventParameters(GTreeManager *Manager);
    virtual ~GTreeEventParameters();


    virtual void        Clear()                         {nReconstructed = 0;}
            Int_t       GetEventNumber()                const {return eventNumber;}
            Int_t       GetNReconstructed()    	        const {return nReconstructed;}
            void        SetEventNumber(const Int_t number)    {eventNumber = number;}
            void        SetNReconstructed(const Int_t number) {nReconstructed = number;}

};

#endif
