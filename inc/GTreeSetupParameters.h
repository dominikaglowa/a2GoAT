#ifndef __GTreeSetupParameters_h__
#define __GTreeSetupParameters_h__


#include "Rtypes.h"
#include "GTree.h"


class  GTreeSetupParameters : public GTree
{
private:
    Int_t		nTagger;
    Double_t	TaggerGlobalOffset;
    Double_t	TaggerTDCLoThr[352];
    Double_t	TaggerTDCHiThr[352];
    Double_t	TaggerTDCOffset[352];
    Double_t	TaggerElectronEnergy[352];
    Double_t	TaggerPhotonEnergy[352];
    Double_t    TaggerEnergyWidth[352];

    Int_t		nNaI;
    Double_t	NaIGlobalOffset;
    Double_t	NaIGlobalScale;
    Int_t       NaIMaxClusters;
    Double_t	NaIClusterThr;
    Double_t	NaIADCLoThr[720];
    Double_t	NaIADCHiThr[720];
    Double_t	NaIADCGain[720];
    Double_t	NaITDCLoThr[720];
    Double_t	NaITDCHiThr[720];
    Double_t	NaITDCOffset[720];

    Int_t		nPID;
    Double_t	PIDGlobalOffset;
    Double_t	PIDADCLoThr[24];
    Double_t	PIDADCHiThr[24];
    Double_t	PIDADCPedestal[24];
    Double_t	PIDADCGain[24];
    Double_t	PIDTDCLoThr[24];
    Double_t	PIDTDCHiThr[24];
    Double_t	PIDTDCOffset[24];
    Double_t	PIDPhi[24];

    Int_t		nBaF2;
    Double_t	BaF2GlobalOffset;
    Double_t	BaF2GlobalScale;
    Double_t	BaF2Distance;
    Int_t       BaF2MaxClusters;
    Double_t	BaF2ClusterThr;
    Double_t	BaF2ADCLoThr[438];
    Double_t	BaF2ADCHiThr[438];
    Double_t	BaF2ADCPedestal[438];
    Double_t	BaF2ADCGain[438];
    Double_t	BaF2TDCLoThr[438];
    Double_t	BaF2TDCHiThr[438];
    Double_t	BaF2TDCOffset[438];
    Double_t	BaF2TDCGain[438];

    Int_t		nVeto;
    Double_t	VetoGlobalOffset;
    Double_t	VetoDistance;
    Double_t	VetoADCLoThr[438];
    Double_t	VetoADCHiThr[438];
    Double_t	VetoADCPedestal[438];
    Double_t	VetoADCGain[438];
    Double_t	VetoTDCLoThr[438];
    Double_t	VetoTDCHiThr[438];
    Double_t	VetoTDCOffset[438];

    Bool_t      hasOverlaps;

protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();

public:
    GTreeSetupParameters(GTreeManager *Manager);
    virtual ~GTreeSetupParameters();

    virtual void    Clear()             {}

            Int_t           GetNTagger()                               const	{return nTagger;}
            Double_t        GetTaggerGlobalOffset()                    const	{return TaggerGlobalOffset;}
    const	Double_t*       GetTaggerTDCLoThr()                        const	{return TaggerTDCLoThr;}
            Double_t        GetTaggerTDCLoThr(const Int_t index)       const	{return TaggerTDCLoThr[index];}
    const	Double_t*       GetTaggerTDCHiThr()                        const	{return TaggerTDCHiThr;}
            Double_t        GetTaggerTDCHiThr(const Int_t index)       const	{return TaggerTDCHiThr[index];}
    const	Double_t*       GetTaggerTDCOffset()                       const	{return TaggerTDCOffset;}
            Double_t        GetTaggerTDCOffset(const Int_t index)      const	{return TaggerTDCOffset[index];}
    const	Double_t*       GetTaggerElectronEnergy()                  const	{return TaggerPhotonEnergy;}
            Double_t        GetTaggerElectronEnergy(const Int_t index) const	{return TaggerPhotonEnergy[index];}
    const	Double_t*       GetTaggerPhotonEnergy()                    const	{return TaggerPhotonEnergy;}
            Double_t        GetTaggerPhotonEnergy(const Int_t index)   const	{return TaggerPhotonEnergy[index];}
    const	Double_t*       GetTaggerEnergyWidth()                     const	{return TaggerEnergyWidth;}
            Double_t        GetTaggerEnergyWidth(const Int_t index)    const	{return TaggerEnergyWidth[index];}

            Int_t           GetNNaI()                                  const	{return nNaI;}
            Double_t        GetNaIGlobalOffset()                       const	{return NaIGlobalOffset;}
            Double_t        GetNaIGlobalScale()                        const	{return NaIGlobalScale;}
            Int_t           GetNaIMaxClusters()                        const	{return NaIMaxClusters;}
            Double_t        GetNaIClusterThr()                         const	{return NaIClusterThr;}
    const	Double_t*       GetNaIADCLoThr()                           const	{return NaIADCLoThr;}
            Double_t        GetNaIADCLoThr(const Int_t index)          const	{return NaIADCLoThr[index];}
    const	Double_t*       GetNaIADCHiThr()                           const	{return NaIADCHiThr;}
            Double_t        GetNaIADCHiThr(const Int_t index)          const	{return NaIADCHiThr[index];}
    const	Double_t*       GetNaIADCGain()                            const	{return NaIADCGain;}
            Double_t        GetNaIADCGain(const Int_t index)           const	{return NaIADCGain[index];}
    const	Double_t*       GetNaITDCLoThr()                           const	{return NaITDCLoThr;}
            Double_t        GetNaITDCLoThr(const Int_t index)          const	{return NaITDCLoThr[index];}
    const	Double_t*       GetNaITDCHiThr()                           const	{return NaITDCHiThr;}
            Double_t        GetNaITDCHiThr(const Int_t index)          const	{return NaITDCHiThr[index];}
    const	Double_t*       GetNaITDCOffset()                          const	{return NaITDCOffset;}
            Double_t        GetNaITDCOffset(const Int_t index)         const	{return NaITDCOffset[index];}

            Int_t           GetNPID()                                  const	{return nPID;}
            Double_t        GetPIDGlobalOffset()                       const	{return PIDGlobalOffset;}
    const	Double_t*       GetPIDADCLoThr()                           const	{return PIDADCLoThr;}
            Double_t        GetPIDADCLoThr(const Int_t index)          const	{return PIDADCLoThr[index];}
    const	Double_t*       GetPIDADCHiThr()                           const	{return PIDADCHiThr;}
            Double_t        GetPIDADCHiThr(const Int_t index)          const	{return PIDADCHiThr[index];}
    const	Double_t*       GetPIDADCPedestal()                        const	{return PIDADCPedestal;}
            Double_t        GetPIDADCPedestal(const Int_t index)       const	{return PIDADCPedestal[index];}
    const	Double_t*       GetPIDADCGain()                            const	{return PIDADCGain;}
            Double_t        GetPIDADCGain(const Int_t index)           const	{return PIDADCGain[index];}
    const	Double_t*       GetPIDTDCLoThr()                           const	{return PIDTDCLoThr;}
            Double_t        GetPIDTDCLoThr(const Int_t index)          const	{return PIDTDCLoThr[index];}
    const	Double_t*       GetPIDTDCHiThr()                           const	{return PIDTDCHiThr;}
            Double_t        GetPIDTDCHiThr(const Int_t index)          const	{return PIDTDCHiThr[index];}
    const	Double_t*       GetPIDTDCOffset()                          const	{return PIDTDCOffset;}
            Double_t        GetPIDTDCOffset(const Int_t index)         const	{return PIDTDCOffset[index];}
    const	Double_t*       GetPIDPhi()                                const	{return PIDPhi;}
            Double_t        GetPIDPhi(const Int_t index)               const	{return PIDPhi[index];}

            Int_t           GetNBaF2()                                 const	{return nBaF2;}
            Double_t        GetBaF2GlobalOffset()                      const	{return BaF2GlobalOffset;}
            Double_t        GetBaF2GlobalScale()                       const	{return BaF2GlobalScale;}
            Double_t        GetBaF2Distance()                          const	{return BaF2Distance;}
            Int_t           GetBaF2MaxClusters()                       const	{return BaF2MaxClusters;}
            Double_t        GetBaF2ClusterThr()                        const	{return BaF2ClusterThr;}
    const	Double_t*       GetBaF2ADCLoThr()                          const	{return BaF2ADCLoThr;}
            Double_t        GetBaF2ADCLoThr(const Int_t index)         const	{return BaF2ADCLoThr[index];}
    const	Double_t*       GetBaF2ADCHiThr()                          const	{return BaF2ADCHiThr;}
            Double_t        GetBaF2ADCHiThr(const Int_t index)         const	{return BaF2ADCHiThr[index];}
    const	Double_t*       GetBaF2ADCPedestal()                       const	{return BaF2ADCPedestal;}
            Double_t        GetBaF2ADCPedestal(const Int_t index)      const	{return BaF2ADCPedestal[index];}
    const	Double_t*       GetBaF2ADCGain()                           const	{return BaF2ADCGain;}
            Double_t        GetBaF2ADCGain(const Int_t index)          const	{return BaF2ADCGain[index];}
    const	Double_t*       GetBaF2TDCLoThr()                          const	{return BaF2TDCLoThr;}
            Double_t        GetBaF2TDCLoThr(const Int_t index)         const	{return BaF2TDCLoThr[index];}
    const	Double_t*       GetBaF2TDCHiThr()                          const	{return BaF2TDCHiThr;}
            Double_t        GetBaF2TDCHiThr(const Int_t index)         const	{return BaF2TDCHiThr[index];}
    const	Double_t*       GetBaF2TDCOffset()                         const	{return BaF2TDCOffset;}
            Double_t        GetBaF2TDCOffset(const Int_t index)        const	{return BaF2TDCOffset[index];}
    const	Double_t*       GetBaF2TDCGain()                           const	{return BaF2TDCGain;}
            Double_t        GetBaF2TDCGain(const Int_t index)          const	{return BaF2TDCGain[index];}

            Int_t           GetNVeto()                                 const	{return nVeto;}
            Double_t        GetVetoGlobalOffset()                      const	{return VetoGlobalOffset;}
            Double_t        GetVetoDistance()                          const	{return VetoDistance;}
    const	Double_t*       GetVetoADCLoThr()                          const	{return VetoADCLoThr;}
            Double_t        GetVetoADCLoThr(const Int_t index)         const	{return VetoADCLoThr[index];}
    const	Double_t*       GetVetoADCHiThr()                          const	{return VetoADCHiThr;}
            Double_t        GetVetoADCHiThr(const Int_t index)         const	{return VetoADCHiThr[index];}
    const	Double_t*       GetVetoADCPedestal()                       const	{return VetoADCPedestal;}
            Double_t        GetVetoADCPedestal(const Int_t index)      const	{return VetoADCPedestal[index];}
    const	Double_t*       GetVetoADCGain()                           const	{return VetoADCGain;}
            Double_t        GetVetoADCGain(const Int_t index)          const	{return VetoADCGain[index];}
    const	Double_t*       GetVetoTDCLoThr()                          const	{return VetoTDCLoThr;}
            Double_t        GetVetoTDCLoThr(const Int_t index)         const	{return VetoTDCLoThr[index];}
    const	Double_t*       GetVetoTDCHiThr()                          const	{return VetoTDCHiThr;}
            Double_t        GetVetoTDCHiThr(const Int_t index)         const	{return VetoTDCHiThr[index];}
    const	Double_t*       GetVetoTDCOffset()                         const	{return VetoTDCOffset;}
            Double_t        GetVetoTDCOffset(const Int_t index)        const	{return VetoTDCOffset[index];}

            Bool_t          HasOverlaps()                              const    {return hasOverlaps;}
};

#endif
