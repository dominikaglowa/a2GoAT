#include "GTreeTrigger.h"



GTreeTrigger::GTreeTrigger(GTreeManager* Manager)    :
    GTree(Manager, TString("trigger")),
    energySum(0),
    multiplicity(0),
    nTriggerPattern(0),
    helicity(0),
    nErrors(0),
    MC_evt_id(-1),
    MC_rnd_id(-1),
    hasHelicity(0),
    hasMCID(0)
{
    for(Int_t i=0; i<GTreeTrigger_MAX; i++)
    {
        triggerPattern[i] = 0;
        errorModuleID[i] = 0;
        errorModuleIndex[i] = 0;
        errorCode[i] = 0;
    }
}

GTreeTrigger::~GTreeTrigger()
{

}

void    GTreeTrigger::SetBranchAdresses()
{
    inputTree->SetBranchAddress("energySum", 	      &energySum);
    inputTree->SetBranchAddress("multiplicity", 	  &multiplicity);
    inputTree->SetBranchAddress("nTriggerPattern",  &nTriggerPattern);
    inputTree->SetBranchAddress("triggerPattern",   triggerPattern);
    inputTree->SetBranchAddress("nErrors", 	      &nErrors);
    inputTree->SetBranchAddress("errorModuleID",    errorModuleID);
    inputTree->SetBranchAddress("errorModuleIndex", errorModuleIndex);
    inputTree->SetBranchAddress("errorCode",        errorCode);
    if(inputTree->GetBranch("helicity"))
    {
        inputTree->SetBranchAddress("helicity",  &helicity);
        hasHelicity = true;
    }
    if(inputTree->GetBranch("mc_evt_id"))
    {
        inputTree->SetBranchAddress("mc_evt_id", &MC_evt_id);
        inputTree->SetBranchAddress("mc_rnd_id", &MC_rnd_id);
        hasMCID = true;
    }
}

void    GTreeTrigger::SetBranches()
{
    outputTree->Branch("energySum",        &energySum,       "energySum/D");
    outputTree->Branch("multiplicity",     &multiplicity,    "multiplicity/I");
    outputTree->Branch("nTriggerPattern",  &nTriggerPattern, "nTriggerPattern/I");
    outputTree->Branch("triggerPattern",   triggerPattern,   "triggerPattern[nTriggerPattern]/I");
    outputTree->Branch("nErrors",          &nErrors,         "nErrors/I");
    outputTree->Branch("errorModuleID",    errorModuleID,    "errorModuleID[nErrors]/I");
    outputTree->Branch("errorModuleIndex", errorModuleIndex, "errorModuleIndex[nErrors]/I");
    outputTree->Branch("errorCode",        errorCode,        "errorCode[nErrors]/I");
    if(hasHelicity) outputTree->Branch("helicity",         &helicity,        "helicity/O");
    if(hasMCID)     outputTree->Branch("mc_evt_id",        &MC_evt_id,       "mc_evt_id/L");
    if(hasMCID)     outputTree->Branch("mc_rnd_id",        &MC_rnd_id,       "mc_rnd_id/L");
}


