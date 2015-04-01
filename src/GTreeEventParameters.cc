#include "GTreeEventParameters.h"



GTreeEventParameters::GTreeEventParameters(GTreeManager* Manager)    :
    GTree(Manager, TString("eventParameters")),
    eventNumber(0),
    nReconstructed(0)
{

}

GTreeEventParameters::~GTreeEventParameters()
{

}

void    GTreeEventParameters::SetBranchAdresses()
{
    inputTree->SetBranchAddress("eventNumber", 	&eventNumber);
    inputTree->SetBranchAddress("nReconstructed", 	&nReconstructed);

}

void    GTreeEventParameters::SetBranches()
{
    outputTree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    outputTree->Branch("nReconstructed", &nReconstructed, "nReconstructed/I");

  }


