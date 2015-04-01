#include "GTreeTrack.h"

GTreeTrack::GTreeTrack(GTreeManager *Manager, const TString& _Name)    :
    GTree(Manager,_Name),
    nTracks(0)
{
    for(Int_t i=0; i<GTreeTrack_MAX; i++)
    {
        clusterEnergy[i] = 0;
        theta[i] = 0;
        phi[i] = 0;
        time[i] = 0;
        clusterSize[i] = 0;
        centralCrystal[i] = -1;
        centralVeto[i] = -1;
        detectors[i] = 0;
        //Charged detector energies
        vetoEnergy[i] = 0;
        MWPC0Energy[i] = 0;
        MWPC1Energy[i] = 0;
    }
}

GTreeTrack::~GTreeTrack()
{

}

void    GTreeTrack::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nTracks",&nTracks);
    inputTree->SetBranchAddress("clusterEnergy",  clusterEnergy);
    inputTree->SetBranchAddress("theta", theta);
    inputTree->SetBranchAddress("phi",  phi);
    inputTree->SetBranchAddress("time", time);
    inputTree->SetBranchAddress("clusterSize", clusterSize);
    inputTree->SetBranchAddress("centralCrystal", centralCrystal);
    inputTree->SetBranchAddress("centralVeto", centralVeto);
    inputTree->SetBranchAddress("detectors", detectors);
    inputTree->SetBranchAddress("vetoEnergy", vetoEnergy);
    inputTree->SetBranchAddress("MWPC0Energy", MWPC0Energy);
    inputTree->SetBranchAddress("MWPC1Energy", MWPC1Energy);
    inputTree->SetBranchAddress("pseudoVertexX", pseudoVertexX);
    inputTree->SetBranchAddress("pseudoVertexY", pseudoVertexY);
    inputTree->SetBranchAddress("pseudoVertexZ", pseudoVertexZ);
}

void    GTreeTrack::SetBranches()
{
    outputTree->Branch("nTracks",&nTracks,"nTracks/I");
    outputTree->Branch("clusterEnergy",  clusterEnergy,  "clusterEnergy[nTracks]/D");
    outputTree->Branch("theta",  theta,  "theta[nTracks]/D");
    outputTree->Branch("phi",  phi,  "phi[nTracks]/D");
    outputTree->Branch("time", time, "time[nTracks]/D");
    outputTree->Branch("clusterSize", clusterSize, "clusterSize[nTracks]/I");
    outputTree->Branch("centralCrystal", centralCrystal, "centralCrystal[nTracks]/I");
    outputTree->Branch("centralVeto", centralVeto, "centralVeto[nTracks]/I");
    outputTree->Branch("detectors", detectors, "detectors[nTracks]/I");
    outputTree->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nTracks]/D");
    outputTree->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nTracks]/D");
    outputTree->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nTracks]/D");
    outputTree->Branch("pseudoVertexX", pseudoVertexX, "pseudoVertexX[nTracks]/D");
    outputTree->Branch("pseudoVertexY", pseudoVertexY, "pseudoVertexY[nTracks]/D");
    outputTree->Branch("pseudoVertexZ", pseudoVertexZ, "pseudoVertexZ[nTracks]/D");
}

void    GTreeTrack::Print(const Bool_t All) const
{
    std::cout << "GTreeTrack: nTracks->" << nTracks << std::endl;
}
