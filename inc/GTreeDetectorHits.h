#ifndef __GTreeDetectorHits_h__
#define __GTreeDetectorHits_h__


#include "Rtypes.h"
#include "GTree.h"


class  GTreeDetectorHits    : public GTree
{
private:
    Int_t		nNaIHits;
    Int_t		NaIHits[720];
    Int_t		NaICluster[720];
    Int_t		nPIDHits;
    Int_t		PIDHits[24];
    Int_t		nMWPCHits;
    Int_t		MWPCHits[860];
    Int_t		nBaF2Hits;
    Int_t		BaF2Hits[438];
    Int_t		BaF2Cluster[438];
    Int_t		nVetoHits;
    Int_t		VetoHits[438];

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeDetectorHits(GTreeManager *Manager);
    virtual ~GTreeDetectorHits();

    virtual void            Clear()                         {nNaIHits=0; nPIDHits=0; nMWPCHits=0; nBaF2Hits=0; nVetoHits=0;}

            Int_t		GetNNaIHits()              	const	{return nNaIHits;}
    const	Int_t*		GetNaIHits()           		const	{return NaIHits;}
            Int_t		GetNaIHits(const Int_t index)	const	{return NaIHits[index];}
    const	Int_t*		GetNaICluster()           		const	{return NaICluster;}
            Int_t		GetNaICluster(const Int_t index)	const	{return NaICluster[index];}

            Int_t		GetNPIDHits()      			const	{return nPIDHits;}
    const	Int_t*		GetPIDHits()               	const	{return PIDHits;}
            Int_t		GetPIDHits(const Int_t index)	const	{return PIDHits[index];}

            Int_t		GetNMWPCHits()       			const	{return nMWPCHits;}
    const	Int_t*		GetMWPCHits()                	const	{return MWPCHits;}
            Int_t		GetMWPCHits(const Int_t index)	const	{return MWPCHits[index];}

            Int_t		GetNBaF2Hits()                   const	{return nBaF2Hits;}
    const	Int_t*		GetBaF2Hits()                    const	{return BaF2Hits;}
            Int_t		GetBaF2Hits(const Int_t index)	const	{return BaF2Hits[index];}
    const	Int_t*		GetBaF2Cluster()                    const	{return BaF2Cluster;}
            Int_t		GetBaF2Cluster(const Int_t index)	const	{return BaF2Cluster[index];}

            Int_t		GetNVetoHits()                 const	{return nVetoHits;}
    const	Int_t*		GetVetoHits()                  const	{return VetoHits;}
            Int_t		GetVetoHits(const Int_t index)	const	{return VetoHits[index];}
};


#endif
