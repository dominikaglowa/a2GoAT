#include "GTreeSetupParameters.h"
#include "GTreeManager.h"

GTreeSetupParameters::GTreeSetupParameters(GTreeManager *Manager)    :
    GTree(Manager, TString("setupParameters"), kFALSE, kTRUE),
    nTagger(0),
    TaggerGlobalOffset(0),
    nNaI(0),
    NaIGlobalOffset(0),
    NaIGlobalScale(0),
    NaIMaxClusters(0),
    NaIClusterThr(0),
    nPID(0),
    PIDGlobalOffset(0),
    nBaF2(0),
    BaF2GlobalOffset(0),
    BaF2GlobalScale(0),
    BaF2Distance(0),
    BaF2MaxClusters(0),
    BaF2ClusterThr(0),
    nVeto(0),
    VetoGlobalOffset(0),
    VetoDistance(0),
    hasOverlaps(0)
{
    for(Int_t i=0; i<352; i++)
    {
        TaggerTDCLoThr[i] = 0;
        TaggerTDCHiThr[i] = 0;
        TaggerTDCOffset[i] = 0;
        TaggerElectronEnergy[i] = 0;
        TaggerPhotonEnergy[i]   = 0;
        TaggerEnergyWidth[i]    = 0;
    }
    for(Int_t i=0; i<720; i++)
    {
        NaIADCLoThr[i] = 0;
        NaIADCHiThr[i] = 0;
        NaIADCGain[i] = 0;
        NaITDCLoThr[i] = 0;
        NaITDCHiThr[i] = 0;
        NaITDCOffset[i] = 0;
    }
    for(Int_t i=0; i<24; i++)
    {
        PIDADCLoThr[i] = 0;
        PIDADCHiThr[i] = 0;
        PIDADCPedestal[i] = 0;
        PIDADCGain[i] = 0;
        PIDTDCLoThr[i] = 0;
        PIDTDCHiThr[i] = 0;
        PIDTDCOffset[i] = 0;
        PIDPhi[i] = 0;
    }
    for(Int_t i=0; i<438; i++)
    {
        BaF2ADCLoThr[i] = 0;
        BaF2ADCHiThr[i] = 0;
        BaF2ADCPedestal[i] = 0;
        BaF2ADCGain[i] = 0;
        BaF2TDCLoThr[i] = 0;
        BaF2TDCHiThr[i] = 0;
        BaF2TDCOffset[i] = 0;
        BaF2TDCGain[i] = 0;
    }
    for(Int_t i=0; i<438; i++)
    {
        VetoADCLoThr[i] = 0;
        VetoADCHiThr[i] = 0;
        VetoADCPedestal[i] = 0;
        VetoADCGain[i] = 0;
        VetoTDCLoThr[i] = 0;
        VetoTDCHiThr[i] = 0;
        VetoTDCOffset[i] = 0;
    }
}

GTreeSetupParameters::~GTreeSetupParameters()
{

}

void    GTreeSetupParameters::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nTagger", &nTagger);
    inputTree->SetBranchAddress("TaggerGlobalOffset", &TaggerGlobalOffset);
    inputTree->SetBranchAddress("TaggerTDCLoThr", TaggerTDCLoThr);
    inputTree->SetBranchAddress("TaggerTDCHiThr", TaggerTDCHiThr);
    inputTree->SetBranchAddress("TaggerTDCOffset", TaggerTDCOffset);
    inputTree->SetBranchAddress("TaggerElectronEnergy", TaggerElectronEnergy);
    inputTree->SetBranchAddress("TaggerPhotonEnergy", TaggerPhotonEnergy);
    if(inputTree->GetBranch("TaggerEnergyWidth"))
    {
        inputTree->SetBranchAddress("TaggerEnergyWidth", TaggerEnergyWidth);
        hasOverlaps = true;
    }

    inputTree->SetBranchAddress("nNaI", &nNaI);
    inputTree->SetBranchAddress("NaIGlobalOffset", &NaIGlobalOffset);
    inputTree->SetBranchAddress("NaIGlobalScale", &NaIGlobalScale);
    inputTree->SetBranchAddress("NaIMaxClusters", &NaIMaxClusters);
    inputTree->SetBranchAddress("NaIClusterThr", &NaIClusterThr);
    inputTree->SetBranchAddress("NaIADCLoThr", NaIADCLoThr);
    inputTree->SetBranchAddress("NaIADCHiThr", NaIADCHiThr);
    inputTree->SetBranchAddress("NaIADCGain", NaIADCGain);
    inputTree->SetBranchAddress("NaITDCLoThr", NaITDCLoThr);
    inputTree->SetBranchAddress("NaITDCHiThr", NaITDCHiThr);
    inputTree->SetBranchAddress("NaITDCOffset", NaITDCOffset);

    inputTree->SetBranchAddress("nPID", &nPID);
    inputTree->SetBranchAddress("PIDGlobalOffset", &PIDGlobalOffset);
    inputTree->SetBranchAddress("PIDADCLoThr", PIDADCLoThr);
    inputTree->SetBranchAddress("PIDADCHiThr", PIDADCHiThr);
    inputTree->SetBranchAddress("PIDADCPedestal", PIDADCPedestal);
    inputTree->SetBranchAddress("PIDADCGain", PIDADCGain);
    inputTree->SetBranchAddress("PIDTDCLoThr", PIDTDCLoThr);
    inputTree->SetBranchAddress("PIDTDCHiThr", PIDTDCHiThr);
    inputTree->SetBranchAddress("PIDTDCOffset", PIDTDCOffset);
    inputTree->SetBranchAddress("PIDPhi", PIDPhi);

    inputTree->SetBranchAddress("nBaF2", &nBaF2);
    inputTree->SetBranchAddress("BaF2GlobalOffset", &BaF2GlobalOffset);
    inputTree->SetBranchAddress("BaF2GlobalScale", &BaF2GlobalScale);
    inputTree->SetBranchAddress("BaF2Distance", &BaF2Distance);
    inputTree->SetBranchAddress("BaF2MaxClusters", &BaF2MaxClusters);
    inputTree->SetBranchAddress("BaF2ClusterThr", &BaF2ClusterThr);
    inputTree->SetBranchAddress("BaF2ADCLoThr", BaF2ADCLoThr);
    inputTree->SetBranchAddress("BaF2ADCHiThr", BaF2ADCHiThr);
    inputTree->SetBranchAddress("BaF2ADCPedestal", BaF2ADCPedestal);
    inputTree->SetBranchAddress("BaF2ADCGain", BaF2ADCGain);
    inputTree->SetBranchAddress("BaF2TDCLoThr", BaF2TDCLoThr);
    inputTree->SetBranchAddress("BaF2TDCHiThr", BaF2TDCHiThr);
    inputTree->SetBranchAddress("BaF2TDCOffset", BaF2TDCOffset);
    inputTree->SetBranchAddress("BaF2TDCGain", BaF2TDCGain);

    inputTree->SetBranchAddress("nVeto", &nVeto);
    inputTree->SetBranchAddress("VetoGlobalOffset", &VetoGlobalOffset);
    inputTree->SetBranchAddress("VetoDistance", &VetoDistance);
    inputTree->SetBranchAddress("VetoADCLoThr", VetoADCLoThr);
    inputTree->SetBranchAddress("VetoADCHiThr", VetoADCHiThr);
    inputTree->SetBranchAddress("VetoADCPedestal", VetoADCPedestal);
    inputTree->SetBranchAddress("VetoADCGain", VetoADCGain);
    inputTree->SetBranchAddress("VetoTDCLoThr", VetoTDCLoThr);
    inputTree->SetBranchAddress("VetoTDCHiThr", VetoTDCHiThr);
    inputTree->SetBranchAddress("VetoTDCOffset", VetoTDCOffset);
}

void    GTreeSetupParameters::SetBranches()
{
    outputTree->Branch("nTagger", &nTagger, "nTagger/I");
    outputTree->Branch("TaggerGlobalOffset", &TaggerGlobalOffset, "TaggerGlobalOffset/D");
    outputTree->Branch("TaggerTDCLoThr", TaggerTDCLoThr, "TaggerTDCLoThr[nTagger]/D");
    outputTree->Branch("TaggerTDCHiThr", TaggerTDCHiThr, "TaggerTDCHiThr[nTagger]/D");
    outputTree->Branch("TaggerTDCOffset", TaggerTDCOffset, "TaggerTDCOffset[nTagger]/D");
    outputTree->Branch("TaggerElectronEnergy", TaggerElectronEnergy, "TaggerElectronEnergy[nTagger]/D");
    outputTree->Branch("TaggerPhotonEnergy", TaggerPhotonEnergy, "TaggerPhotonEnergy[nTagger]/D");
    if(hasOverlaps) outputTree->Branch("TaggerEnergyWidth", TaggerEnergyWidth, "TaggerEnergyWidth[nTagger]/D");

    outputTree->Branch("nNaI", &nNaI, "nNaI/I");
    outputTree->Branch("NaIGlobalOffset", &NaIGlobalOffset, "NaIGlobalOffset/D");
    outputTree->Branch("NaIGlobalScale", &NaIGlobalScale, "NaIGlobalScale/D");
    outputTree->Branch("NaIMaxClusters", &NaIMaxClusters, "NaIMaxClusters/I");
    outputTree->Branch("NaIClusterThr", &NaIClusterThr, "NaIClusterThr/D");
    outputTree->Branch("NaIADCLoThr", NaIADCLoThr, "NaIADCLoThr[nNaI]/D");
    outputTree->Branch("NaIADCHiThr", NaIADCHiThr, "NaIADCHiThr[nNaI]/D");
    outputTree->Branch("NaIADCGain", NaIADCGain, "NaIADCGain[nNaI]/D");
    outputTree->Branch("NaITDCLoThr", NaITDCLoThr, "NaITDCLoThr[nNaI]/D");
    outputTree->Branch("NaITDCHiThr", NaITDCHiThr, "NaITDCHiThr[nNaI]/D");
    outputTree->Branch("NaITDCOffset", NaITDCOffset, "NaITDCOffset[nNaI]/D");

    outputTree->Branch("nPID", &nPID, "nPID/I");
    outputTree->Branch("PIDGlobalOffset", &PIDGlobalOffset, "PIDGlobalOffset/D");
    outputTree->Branch("PIDADCLoThr", PIDADCLoThr, "PIDADCLoThr[nPID]/D");
    outputTree->Branch("PIDADCHiThr", PIDADCHiThr, "PIDADCHiThr[nPID]/D");
    outputTree->Branch("PIDADCPedestal", PIDADCPedestal, "PIDADCPedestal[nPID]/D");
    outputTree->Branch("PIDADCGain", PIDADCGain, "PIDADCGain[nPID]/D");
    outputTree->Branch("PIDTDCLoThr", PIDTDCLoThr, "PIDTDCLoThr[nPID]/D");
    outputTree->Branch("PIDTDCHiThr", PIDTDCHiThr, "PIDTDCHiThr[nPID]/D");
    outputTree->Branch("PIDTDCOffset", PIDTDCOffset, "PIDTDCOffset[nPID]/D");
    outputTree->Branch("PIDPhi", PIDPhi, "PIDPhi[nPID]/D");

    outputTree->Branch("nBaF2", &nBaF2, "nBaF2/I");
    outputTree->Branch("BaF2GlobalOffset", &BaF2GlobalOffset, "BaF2GlobalOffset/D");
    outputTree->Branch("BaF2GlobalScale", &BaF2GlobalScale, "BaF2GlobalScale/D");
    outputTree->Branch("BaF2Distance", &BaF2Distance, "BaF2Distance/D");
    outputTree->Branch("BaF2MaxClusters", &BaF2MaxClusters, "BaF2MaxClusters/I");
    outputTree->Branch("BaF2ClusterThr", &BaF2ClusterThr, "BaF2ClusterThr/D");
    outputTree->Branch("BaF2ADCLoThr", BaF2ADCLoThr, "BaF2ADCLoThr[nBaF2]/D");
    outputTree->Branch("BaF2ADCHiThr", BaF2ADCHiThr, "BaF2ADCHiThr[nBaF2]/D");
    outputTree->Branch("BaF2ADCPedestal", BaF2ADCPedestal, "BaF2ADCPedestal[nBaF2]/D");
    outputTree->Branch("BaF2ADCGain", BaF2ADCGain, "BaF2ADCGain[nBaF2]/D");
    outputTree->Branch("BaF2TDCLoThr", BaF2TDCLoThr, "BaF2TDCLoThr[nBaF2]/D");
    outputTree->Branch("BaF2TDCHiThr", BaF2TDCHiThr, "BaF2TDCHiThr[nBaF2]/D");
    outputTree->Branch("BaF2TDCOffset", BaF2TDCOffset, "BaF2TDCOffset[nBaF2]/D");
    outputTree->Branch("BaF2TDCGain", BaF2TDCGain, "BaF2TDCGain[nBaF2]/D");

    outputTree->Branch("nVeto", &nVeto, "nVeto/I");
    outputTree->Branch("VetoGlobalOffset", &VetoGlobalOffset, "VetoGlobalOffset/D");
    outputTree->Branch("VetoDistance", &VetoDistance, "VetoDistance/D");
    outputTree->Branch("VetoADCLoThr", VetoADCLoThr, "VetoADCLoThr[nVeto]/D");
    outputTree->Branch("VetoADCHiThr", VetoADCHiThr, "VetoADCHiThr[nVeto]/D");
    outputTree->Branch("VetoADCPedestal", VetoADCPedestal, "VetoADCPedestal[nVeto]/D");
    outputTree->Branch("VetoADCGain", VetoADCGain, "VetoADCGain[nVeto]/D");
    outputTree->Branch("VetoTDCLoThr", VetoTDCLoThr, "VetoTDCLoThr[nVeto]/D");
    outputTree->Branch("VetoTDCHiThr", VetoTDCHiThr, "VetoTDCHiThr[nVeto]/D");
    outputTree->Branch("VetoTDCOffset", VetoTDCOffset, "VetoTDCOffset[nVeto]/D");
}
