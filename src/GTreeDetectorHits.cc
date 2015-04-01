#include "GTreeDetectorHits.h"



GTreeDetectorHits::GTreeDetectorHits(GTreeManager *Manager)    :
    GTree(Manager, TString("detectorHits")),
    nNaIHits(0),
    nPIDHits(0),
    nMWPCHits(0),
    nBaF2Hits(0),
    nVetoHits(0)
{
}

GTreeDetectorHits::~GTreeDetectorHits()
{
}

void    GTreeDetectorHits::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nNaIHits", &nNaIHits);
    inputTree->SetBranchAddress("NaIHits", NaIHits);
    inputTree->SetBranchAddress("NaICluster", NaICluster);
    inputTree->SetBranchAddress("nPIDHits", &nPIDHits);
    inputTree->SetBranchAddress("PIDHits", PIDHits);
    inputTree->SetBranchAddress("nMWPCHits", &nMWPCHits);
    inputTree->SetBranchAddress("MWPCHits", MWPCHits);
    inputTree->SetBranchAddress("nBaF2Hits", &nBaF2Hits);
    inputTree->SetBranchAddress("BaF2Hits", BaF2Hits);
    inputTree->SetBranchAddress("BaF2Cluster", BaF2Cluster);
    inputTree->SetBranchAddress("nVetoHits", &nVetoHits);
    inputTree->SetBranchAddress("VetoHits", VetoHits);
}

void    GTreeDetectorHits::SetBranches()
{
    outputTree->Branch("nNaIHits", &nNaIHits, "nNaIHits/I");
    outputTree->Branch("NaIHits", NaIHits, "NaIHits[nNaIHits]/I");
    outputTree->Branch("NaICluster", NaICluster, "NaICluster[nNaIHits]/I");
    outputTree->Branch("nPIDHits", &nPIDHits, "nPIDHits/I");
    outputTree->Branch("PIDHits", PIDHits, "PIDHits[nPIDHits]/I");
    outputTree->Branch("nMWPCHits", &nMWPCHits, "nMWPCHits/I");
    outputTree->Branch("MWPCHits", MWPCHits, "MWPCHits[nMWPCHits]/I");
    outputTree->Branch("nBaF2Hits", &nBaF2Hits, "nBaF2Hits/I");
    outputTree->Branch("BaF2Hits", BaF2Hits, "BaF2Hits[nBaF2Hits]/I");
    outputTree->Branch("BaF2Cluster", BaF2Cluster, "BaF2Cluster[nBaF2Hits]/I");
    outputTree->Branch("nVetoHits", &nVetoHits, "nVetoHits/I");
    outputTree->Branch("VetoHits", VetoHits, "VetoHits[nVetoHits]/I");
}


