#include "GTreeLinPol.h"

#include <TLeaf.h>


GTreeLinPol::GTreeLinPol(GTreeManager *Manager)    :
    GTree(Manager, TString("linPol"), kTRUE),
    plane(0),
    edge(0),
    edgeSetting(0)
{
    for(Int_t i=0; i<352; i++)
    {
        polarizationTable[i] = 0;
        enhancementTable[i] = 0;
    }
    
}

GTreeLinPol::~GTreeLinPol()
{

}

void    GTreeLinPol::SetBranchAdresses()
{
    inputTree->SetBranchAddress("plane", &plane);
    inputTree->SetBranchAddress("edge", &edge);
    inputTree->SetBranchAddress("edgeSetting", &edgeSetting);
    inputTree->SetBranchAddress("polarizationTable", polarizationTable);
    inputTree->SetBranchAddress("enhancementTable", enhancementTable);
    
}

void    GTreeLinPol::SetBranches()
{
    outputTree->Branch("plane", &plane, "plane/I");
    outputTree->Branch("edge", &edge, "edge/D");
    outputTree->Branch("edgeSetting", &edgeSetting, "edgeSetting/D");
    outputTree->Branch("polarizationTable", polarizationTable, "polarizationTable[352]/D");
    outputTree->Branch("enhancementTable", enhancementTable, "enhancementTable[352]/D");
}
