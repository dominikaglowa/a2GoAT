#include "GTreeParticle.h"
#include "GTreeManager.h"

#include <TCanvas.h>

using namespace std;


GTreeParticle::GTreeParticle(GTreeManager *Manager, const TString& _Name)    :
    GTreeTrack(Manager,_Name),
    nParticles(0)
{
    for(Int_t i=0; i<GTreeTrack_MAX; i++)
    {
        mass[i]       = 0;
        trackIndex[i] = -1;
    }
}

GTreeParticle::~GTreeParticle()
{
}

void    GTreeParticle::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nParticles",&nParticles);
    inputTree->SetBranchAddress("clusterEnergy",  clusterEnergy);
    inputTree->SetBranchAddress("theta", theta);
    inputTree->SetBranchAddress("phi",  phi);
    inputTree->SetBranchAddress("mass", mass);
    inputTree->SetBranchAddress("time", time);
    inputTree->SetBranchAddress("clusterSize", clusterSize);
    inputTree->SetBranchAddress("centralCrystal", centralCrystal);
    inputTree->SetBranchAddress("centralVeto", centralVeto);
    inputTree->SetBranchAddress("detectors", detectors);
    inputTree->SetBranchAddress("vetoEnergy", vetoEnergy);
    inputTree->SetBranchAddress("MWPC0Energy", MWPC0Energy);
    inputTree->SetBranchAddress("MWPC1Energy", MWPC1Energy);
    inputTree->SetBranchAddress("trackIndex", trackIndex);
}

void    GTreeParticle::SetBranches()
{
    outputTree->Branch("nParticles",&nParticles, "nParticles/I");
    outputTree->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
    outputTree->Branch("theta", theta, "theta[nParticles]/D");
    outputTree->Branch("phi", phi, "phi[nParticles]/D");
    outputTree->Branch("mass", mass, "mass[nParticles]/D");
    outputTree->Branch("time", time, "time[nParticles]/D");
    outputTree->Branch("clusterSize", clusterSize, "clusterSize[nParticles]/I");
    outputTree->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
    outputTree->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
    outputTree->Branch("detectors", detectors, "detectors[nParticles]/I");
    outputTree->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
    outputTree->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
    outputTree->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
    outputTree->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");

}


Bool_t	GTreeParticle::Write()
{
    if(!manager->outputFile)          return kFALSE;
    if(!outputTree)                   return kFALSE;
    if(!IsOpenForOutput())          return kFALSE;

    TCanvas c1("c1");
    c1.cd();
    outputTree->Draw("nParticles>>htmp(1,0,1)");
    TH1F* hist = (TH1F*)gDirectory->Get("htmp");
    //std::cout << hist->GetEntries()<< "\t" << hist->GetBinContent(hist->GetBin(1)) << std::endl;
    if(hist->GetEntries() == hist->GetBinContent(hist->GetBin(1)))
    {
        std::cout << "tree " << GetName() << " has not been written to disk. All Events have 0 " << GetName() << "." << std::endl;
        saveToFile  = false;
        return kTRUE;
    }

    return GTree::Write();
}

void    GTreeParticle::AddParticle(const Double_t _clusterEnergy, const Double_t _theta, const Double_t _phi, const Double_t _mass, const Double_t _time, const Int_t _clusterSize, const Int_t _centralCrystal, const Int_t _centralVeto, const Int_t _detectors, const Double_t _vetoEnergy, const Double_t _MWPC0Energy, const Double_t _MWPC1Energy, const Int_t _trackIndex)
{
    clusterEnergy[nParticles]  = _clusterEnergy;
    theta[nParticles]          = _theta;
    phi[nParticles]            = _phi;
    mass[nParticles]           = _mass;
    time[nParticles]           = _time;
    clusterSize[nParticles]    = _clusterSize;
    centralCrystal[nParticles] = _centralCrystal;
    centralVeto[nParticles]    = _centralVeto;
    detectors[nParticles]      = _detectors;
    vetoEnergy[nParticles]     = _vetoEnergy;
    MWPC0Energy[nParticles]    = _MWPC0Energy;
    MWPC1Energy[nParticles]    = _MWPC1Energy;
    trackIndex[nParticles]     = _trackIndex;
    nParticles++;
    manager->countReconstructed++;
}


void    GTreeParticle::RemoveAllParticles()
{
    manager->countReconstructed -= nParticles;
    Clear();
}


void    GTreeParticle::RemoveParticles(const Int_t nIndices, const Int_t* indices)
{
    Int_t*  sort    = new Int_t[nIndices];
    TMath::Sort(nIndices, indices, sort);
    for(Int_t i=0; i<nIndices; i++)
    {
        nParticles--;
        if(sort[i] != nParticles)
        {
            clusterEnergy[sort[i]]  = clusterEnergy[nParticles];
            theta[sort[i]]          = theta[nParticles];
            phi[sort[i]]            = phi[nParticles];
            mass[sort[i]]           = mass[nParticles];
            time[sort[i]]           = time[nParticles];
            clusterSize[sort[i]]    = clusterSize[nParticles];
            centralCrystal[sort[i]] = centralCrystal[nParticles];
            centralVeto[sort[i]]    = centralVeto[nParticles];
            detectors[sort[i]]   = detectors[nParticles];
            vetoEnergy[sort[i]]  = vetoEnergy[nParticles];
            MWPC0Energy[sort[i]] = MWPC0Energy[nParticles];
            MWPC1Energy[sort[i]] = MWPC1Energy[nParticles];
            trackIndex[sort[i]]  = trackIndex[nParticles];
        }
    }
    manager->countReconstructed -= nIndices;
}


void    GTreeParticle::PrintParticle(const Int_t i) const
{
    cout << "\tParticle " << i << ": " << endl;
    cout << "\tPx: " << Particle(i).Px() << "   Py: " << Particle(i).Py() << "   Pz: " << Particle(i).Pz() << "   E: " << Particle(i).E() << endl;
    cout << "\tDetectors: " << detectors[i] << "   clusterSize: " << clusterSize[i]  << endl;
    cout << "\tvetoEnergy: " << vetoEnergy[i] << "   MWPC0Energy: " << MWPC0Energy[i] << "   MWPC1Energy: " << MWPC1Energy[i]  << endl;
}

void    GTreeParticle::Print() const
{
    GTree::Print();
    cout << "nParticles: " << nParticles << endl;
    for(Int_t i=0; i<nParticles; i++)
        PrintParticle(i);
}
