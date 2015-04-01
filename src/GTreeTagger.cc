#include "GTreeTagger.h"
#include "GTreeManager.h"

GTreeTagger::GTreeTagger(GTreeManager *Manager)    :
    GTree(Manager, TString("tagger")),
    nTagged(0),
    hasEnergy(0)
{
    for(Int_t i=0; i<GTreeTagger_MAX; i++)
    {
        taggedChannel[i] = 0;
        taggedTime[i]    = 0;
        taggedEnergy[i]  = 0;
    }
    for(Int_t i=0; i<352; i++) calibration[i] = 0;
}

GTreeTagger::~GTreeTagger()
{

}

void    GTreeTagger::SetBranchAdresses()
{
    inputTree->SetBranchAddress("nTagged", 	   &nTagged);
    inputTree->SetBranchAddress("taggedChannel", taggedChannel);
    inputTree->SetBranchAddress("taggedTime",    taggedTime);
    if(inputTree->GetBranch("taggedEnergy"))
    {
        inputTree->SetBranchAddress("taggedEnergy",  taggedEnergy);
        hasEnergy = true;
    }

}

void    GTreeTagger::SetBranches()
{
    outputTree->Branch("nTagged",       &nTagged,      "nTagged/I");
    outputTree->Branch("taggedChannel", taggedChannel, "taggedChannel[nTagged]/I");
    outputTree->Branch("taggedTime",    taggedTime,    "taggedTime[nTagged]/D");
    if(hasEnergy) outputTree->Branch("taggedEnergy",  taggedEnergy,  "taggedEnergy[nTagged]/D");
}

void    GTreeTagger::SetCalibration(const Int_t nChan, const Double_t *energy)
{
    for(Int_t i=0; i<nChan; i++) calibration[i] = energy[i];
}

TLorentzVector  GTreeTagger::GetVectorProtonTarget(const Int_t index)    const
{
    return TLorentzVector(0, 0, taggedEnergy[index], taggedEnergy[index] + (manager->pdgDB->GetParticle("proton")->Mass()*1000));
}
