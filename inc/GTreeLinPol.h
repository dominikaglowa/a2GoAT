#ifndef __GTreeLinPol_h__
#define __GTreeLinPol_h__


#include "Rtypes.h"
#include "GTree.h"


#define GTreeLinPol_MAX 4096

class  GTreeLinPol : public GTree
{
private:
    Int_t	plane;
    Double_t	edge;
    Double_t	edgeSetting;    
    Double_t	polarizationTable[352];
    Double_t	enhancementTable[352];

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeLinPol(GTreeManager *Manager);
    virtual ~GTreeLinPol();

    virtual void        Clear()     {}
            Int_t 	    GetPolarizationPlane() 	             const	{return plane;}
            Double_t 	GetEdge() 	                         const	{return edge;}
            Double_t 	GetEdgeSetting()                     const	{return edgeSetting;}
     const  Double_t*	GetPolarizationTable()               const	{return polarizationTable;}
     const  Double_t	GetPolarizationDegree(Int_t channel) const	{return polarizationTable[channel];}
     const  Double_t*	GetEnhancementTable() 	             const	{return enhancementTable;}
     const  Double_t	GetEnhancementTable(Int_t channel)   const	{return enhancementTable[channel];}

};


#endif
