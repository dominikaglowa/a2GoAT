void RenameGoAT(TString sFile){

  if(sFile.EndsWith(".root")) RunFile(sFile);
  else{
    TString sDir = gSystem->ExpandPathName(sFile.Data());

    Int_t iLength;
    TSystemFile *sfFile;
    const char* cFile;

    gSystem->cd(sDir);
    
    TSystemDirectory *sdFile = new TSystemDirectory("files","./");
    TList *lFile = sdFile->GetListOfFiles();
    lFile->Sort();
    TIter itFile(lFile);
    
    while((sfFile=(TSystemFile*)itFile())){
      cFile = sfFile->GetName();
      sFile = cFile;
      iLength = strlen(cFile);
      if((iLength >= 5) && (strcmp(&cFile[iLength-5],".root") == 0)){
	RunFile(sFile);
      }
    }
  }
}

void RunFile(TString sFile){
  
  printf("Renaming %s\n",sFile.Data());
  TStopwatch timer;
  timer.Start();

  TFile fOld(sFile,"READ");

  // Acqu Trees
  TTree *taggerOld       = (TTree*)fOld.Get("treeTagger");
  TTree *tracksOld       = (TTree*)fOld.Get("treeRawEvent");
  TTree *detectorHitsOld = (TTree*)fOld.Get("treeDetectorHits");
  TTree *linPolOld       = (TTree*)fOld.Get("treeLinPol");
  TTree *triggerOld      = (TTree*)fOld.Get("treeTrigger");
  TTree *scalersOld       = (TTree*)fOld.Get("treeScaler");

  // GoAT Trees
  TTree *eventParametersOld = (TTree*)fOld.Get("treeEventParameters");
  TTree *rootinosOld        = (TTree*)fOld.Get("rootino");
  TTree *photonsOld         = (TTree*)fOld.Get("gamma");
  TTree *electronsOld       = (TTree*)fOld.Get("e-");
  TTree *chargedPionsOld    = (TTree*)fOld.Get("pi+");
  TTree *protonsOld         = (TTree*)fOld.Get("proton");
  TTree *neutronsOld        = (TTree*)fOld.Get("neutron");
  TTree *neutralPionsOld    = (TTree*)fOld.Get("pi0");
  TTree *etasOld            = (TTree*)fOld.Get("eta");
  TTree *etaPrimesOld       = (TTree*)fOld.Get("eta'");

  TClonesArray *particles   = new TClonesArray("TLorentzVector", 64);

  sFile.ReplaceAll(".root","_renamed.root");
  TFile fNew(sFile,"RECREATE");

  if(taggerOld){
    printf("         treeTagger to tagger\n");

    Int_t nTagged          = 0;
    Double_t *taggedEnergy = new Double_t[4096];
    Int_t *taggedChannel   = new Int_t[4096];
    Double_t *taggedTime   = new Double_t[4096];

    TTree *tagger = new TTree("tagger", "tagger");

    if(taggerOld->GetBranch("nTagged")){
      taggerOld->SetBranchAddress("nTagged", &nTagged);
      tagger->Branch("nTagged", &nTagged,"nTagged/I");
      if(taggerOld->GetBranch("tagged_ch")){
	taggerOld->SetBranchAddress("tagged_ch", taggedChannel);
	tagger->Branch("taggedChannel", taggedChannel, "taggedChannel[nTagged]/I");
      }
      if(taggerOld->GetBranch("tagged_t")){
	taggerOld->SetBranchAddress("tagged_t", taggedTime);
	tagger->Branch("taggedTime", taggedTime, "taggedTime[nTagged]/D");
      }
      if(taggerOld->GetBranch("photonbeam_E")){
	taggerOld->SetBranchAddress("photonbeam_E", taggedEnergy);
	tagger->Branch("taggedEnergy", taggedEnergy, "taggedEnergy[nTagged]/D");
      }
    }

    for(Int_t i=0; i<taggerOld->GetEntries(); i++){
      taggerOld->GetEvent(i);
      tagger->Fill();
    }
    tagger->Write();

    delete tagger;

    delete taggedEnergy;
    delete taggedChannel;
    delete taggedTime;
  }

  UInt_t nParticles       = 0;
  Double_t *clusterEnergy = new Double_t[128];
  Double_t *theta         = new Double_t[128];
  Double_t *phi           = new Double_t[128];
  Double_t *mass          = new Double_t[128];
  Double_t *time          = new Double_t[128];
  UChar_t *clusterSize    = new UChar_t[128];
  Int_t *centralCrystal   = new Int_t[128];
  Int_t *centralVeto      = new Int_t[128];
  UChar_t *apparatus      = new UChar_t[128];
  Double_t *vetoEnergy    = new Double_t[128];
  Double_t *MWPC0Energy   = new Double_t[128];
  Double_t *MWPC1Energy   = new Double_t[128];
  Double_t *pseudoVertexX = new Double_t[128];
  Double_t *pseudoVertexY = new Double_t[128];
  Double_t *pseudoVertexZ = new Double_t[128];
  Int_t *trackIndex       = new Int_t[128];

  Int_t nParticlesNew   = 0;
  Int_t *clusterSizeNew = new Int_t[128];
  Int_t *detectors   = new Int_t[128];

  if(tracksOld){
    printf("         treeRawEvent to tracks\n");

    Int_t nTracks           = 0;

    TTree *tracks = new TTree("tracks", "tracks");

    if(tracksOld->GetBranch("nParticles")){
      tracksOld->SetBranchAddress("nParticles", &nTracks);
      tracks->Branch("nTracks", &nTracks, "nTracks/I");
      if(tracksOld->GetBranch("Ek")){
	tracksOld->SetBranchAddress("Ek", clusterEnergy);
	tracks->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nTracks]/D");
      }
      if(tracksOld->GetBranch("Theta")){
	tracksOld->SetBranchAddress("Theta", theta);
	tracks->Branch("theta", theta, "theta[nTracks]/D");
      }
      if(tracksOld->GetBranch("Phi")){
	tracksOld->SetBranchAddress("Phi", phi);
	tracks->Branch("phi", phi, "phi[nTracks]/D");
      }
      if(tracksOld->GetBranch("time")){
	tracksOld->SetBranchAddress("time", time);
	tracks->Branch("time", time, "time[nTracks]/D");
      }
      if(tracksOld->GetBranch("clusterSize")){
	tracksOld->SetBranchAddress("clusterSize", clusterSize);
	tracks->Branch("clusterSize", clusterSizeNew, "clusterSize[nTracks]/I");
      }
      if(tracksOld->GetBranch("centralCrys")){
	tracksOld->SetBranchAddress("centralCrys", centralCrystal);
	tracks->Branch("centralCrystal", centralCrystal, "centralCrystal[nTracks]/I");
      }
      if(tracksOld->GetBranch("centralVeto")){
	tracksOld->SetBranchAddress("centralVeto", centralVeto);
	tracks->Branch("centralVeto", centralVeto, "centralVeto[nTracks]/I");
      }
      if(tracksOld->GetBranch("Apparatus")){
	tracksOld->SetBranchAddress("Apparatus", apparatus);
	tracks->Branch("detectors", detectors, "detectors[nTracks]/I");
      }
      if(tracksOld->GetBranch("d_E")){
	tracksOld->SetBranchAddress("d_E", vetoEnergy);
	tracks->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nTracks]/D");
      }
      if(tracksOld->GetBranch("WC0_E")){
	tracksOld->SetBranchAddress("WC0_E", MWPC0Energy);
	tracks->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nTracks]/D");
      }
      if(tracksOld->GetBranch("WC1_E")){
	tracksOld->SetBranchAddress("WC1_E", MWPC1Energy);
	tracks->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nTracks]/D");
      }
      if(tracksOld->GetBranch("WC_Vertex_X")){
	tracksOld->SetBranchAddress("WC_Vertex_X", pseudoVertexX);
	tracks->Branch("pseudoVertexX", pseudoVertexX, "pseudoVertexX[nTracks]/D");
      }
      if(tracksOld->GetBranch("WC_Vertex_Y")){
	tracksOld->SetBranchAddress("WC_Vertex_Y", pseudoVertexY);
	tracks->Branch("pseudoVertexY", pseudoVertexY, "pseudoVertexY[nTracks]/D");
      }
      if(tracksOld->GetBranch("WC_Vertex_Z")){
	tracksOld->SetBranchAddress("WC_Vertex_Z", pseudoVertexZ);
	tracks->Branch("pseudoVertexZ", pseudoVertexZ, "pseudoVertexZ[nTracks]/D");
      }
    }

    for(Int_t i=0; i<tracksOld->GetEntries(); i++){
      tracksOld->GetEvent(i);      
      for(Int_t j=0; j<nTracks; j++){
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
      }
      tracks->Fill();
    }
    tracks->Write();

    delete tracks;
  }

  if(detectorHitsOld){
    printf("         treeDetectorHits to detectorHits\n");

    Int_t nNaIHits          = 0;
    Int_t *NaIHits          = new Int_t[860];
    Int_t *NaICluster       = new Int_t[860];
    Int_t nPIDHits          = 0;
    Int_t *PIDHits          = new Int_t[860];
    Int_t nMWPCHits         = 0;
    Int_t *MWPCHits         = new Int_t[860];
    Int_t nBaF2Hits    = 0;
    Int_t *BaF2Hits    = new Int_t[860];
    Int_t *BaF2Cluster = new Int_t[860];
    Int_t nVetoHits         = 0;
    Int_t *VetoHits         = new Int_t[860];
	
    TTree *detectorHits = new TTree("detectorHits", "detectorHits");

    if(detectorHitsOld->GetBranch("nNaI_Hits")){
      detectorHitsOld->SetBranchAddress("nNaI_Hits", &nNaIHits);
      detectorHits->Branch("nNaIHits", &nNaIHits, "nNaIHits/I");
      if(detectorHitsOld->GetBranch("NaI_Hits")){
	detectorHitsOld->SetBranchAddress("NaI_Hits", NaIHits);
	detectorHits->Branch("NaIHits", NaIHits, "NaIHits[nNaIHits]/I");
      }
      if(detectorHitsOld->GetBranch("NaI_Cluster")){
	detectorHitsOld->SetBranchAddress("NaI_Cluster", NaICluster);
	detectorHits->Branch("NaICluster", NaICluster, "NaICluster[nNaIHits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nPID_Hits")){
      detectorHitsOld->SetBranchAddress("nPID_Hits", &nPIDHits);
      detectorHits->Branch("nPIDHits", &nPIDHits, "nPIDHits/I");
      if(detectorHitsOld->GetBranch("PID_Hits")){
	detectorHitsOld->SetBranchAddress("PID_Hits", PIDHits);
	detectorHits->Branch("PIDHits", PIDHits, "PIDHits[nPIDHits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nWC_Hits")){
      detectorHitsOld->SetBranchAddress("nWC_Hits", &nMWPCHits);
      detectorHits->Branch("nMWPCHits", &nMWPCHits, "nMWPCHits/I");
      if(detectorHitsOld->GetBranch("WC_Hits")){
	detectorHitsOld->SetBranchAddress("WC_Hits", MWPCHits);
	detectorHits->Branch("MWPCHits", MWPCHits, "MWPCHits[nMWPCHits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nBaF2_PbWO4_Hits")){
      detectorHitsOld->SetBranchAddress("nBaF2_PbWO4_Hits", &nBaF2Hits);
      detectorHits->Branch("nBaF2Hits", &nBaF2Hits, "nBaF2Hits/I");
      if(detectorHitsOld->GetBranch("BaF2_PbWO4_Hits")){
	detectorHitsOld->SetBranchAddress("BaF2_PbWO4_Hits", BaF2Hits);
	detectorHits->Branch("BaF2Hits", BaF2Hits, "BaF2Hits[nBaF2Hits]/I");
      }
      if(detectorHitsOld->GetBranch("BaF2_PbWO4_Cluster")){
	detectorHitsOld->SetBranchAddress("BaF2_PbWO4_Cluster", BaF2Cluster);
	detectorHits->Branch("BaF2Cluster", BaF2Cluster, "BaF2Cluster[nBaF2Hits]/I");
      }
    }
    if(detectorHitsOld->GetBranch("nVeto_Hits")){
      detectorHitsOld->SetBranchAddress("nVeto_Hits", &nVetoHits);
      detectorHits->Branch("nVetoHits", &nVetoHits, "nVetoHits/I");
      if(detectorHitsOld->GetBranch("Veto_Hits")){
	detectorHitsOld->SetBranchAddress("Veto_Hits", VetoHits);
	detectorHits->Branch("VetoHits", VetoHits, "VetoHits[nVetoHits]/I");
      }
    }

    for(Int_t i=0; i<detectorHitsOld->GetEntries(); i++){
      detectorHitsOld->GetEvent(i);
      detectorHits->Fill();
    }
    detectorHits->Write();

    delete detectorHits;

    delete NaIHits;
    delete NaICluster;
    delete PIDHits;
    delete MWPCHits;
    delete BaF2Hits;
    delete BaF2Cluster;
    delete VetoHits;
  }

  if(linPolOld){
    printf("         treeLinPol to linPol\n");

    Int_t plane                 = 0;
    Double_t edge               = 0;
    Double_t edgeSetting        = 0;
    Double_t *polarizationTable = new Double_t[352];
    Double_t *enhancementTable  = new Double_t[352];

    TTree *linPol = new TTree("linPol", "linPol");

    if(linPolOld->GetBranch("plane")){
      linPolOld->SetBranchAddress("plane", &plane);
      linPol->Branch("plane", &plane, "plane/I");
    }
    if(linPolOld->GetBranch("edge")){
      linPolOld->SetBranchAddress("edge", &edge);
      linPol->Branch("edge", &edge, "edge/D");
    }
    if(linPolOld->GetBranch("edgeSetting")){
      linPolOld->SetBranchAddress("edgeSetting", &edgeSetting);
      linPol->Branch("edgeSetting", &edgeSetting, "edgeSetting/D");
    }
    if(linPolOld->GetBranch("polTable")){
      linPolOld->SetBranchAddress("polTable", polarizationTable);
      linPol->Branch("polarizationTable", polarizationTable, "polarizationTable[352]/D");
    }
    if(linPolOld->GetBranch("enhTable")){
      linPolOld->SetBranchAddress("enhTable", enhancementTable);
      linPol->Branch("enhancementTable", enhancementTable, "enhancementTable[352]/D");
    }

    for(Int_t i=0; i<linPolOld->GetEntries(); i++){
      linPolOld->GetEvent(i);
      linPol->Fill();
    }
    linPol->Write();

    delete linPol;

    delete polarizationTable;
    delete enhancementTable;
  }

  if(triggerOld){
    printf("         treeTrigger to trigger\n");

    Double_t energySum      = 0;
    Int_t multiplicity      = 0;
    Bool_t helicity         = false;
    Int_t nErrors           = 0;
    Int_t nErrors           = 0;
    Int_t *errorModuleID    = new Int_t[256];
    Int_t *errorModuleIndex = new Int_t[256];
    Int_t *errorCode        = new Int_t[256];
    Int_t nTriggerPattern   = 0;
    Int_t *triggerPattern   = new Int_t[32];

    TTree *trigger = new TTree("trigger", "trigger");

    if(triggerOld->GetBranch("ESum")){
      triggerOld->SetBranchAddress("ESum", &energySum);
      trigger->Branch("energySum", &energySum, "energySum/D");
    }
    if(triggerOld->GetBranch("Mult")){
      triggerOld->SetBranchAddress("Mult", &multiplicity);
      trigger->Branch("multiplicity", &multiplicity, "multiplicity/I");
    }
    if(triggerOld->GetBranch("Helicity")){
      triggerOld->SetBranchAddress("Helicity", &helicity);
      trigger->Branch("helicity", &helicity, "helicity/O");
    }
    if(triggerOld->GetBranch("nError")){
      triggerOld->SetBranchAddress("nError", &nErrors);
      trigger->Branch("nErrors", &nErrors, "nErrors/I");
      if(triggerOld->GetBranch("ErrModID")){
	triggerOld->SetBranchAddress("ErrModID", errorModuleID);
	trigger->Branch("errorModuleID", errorModuleID, "errorModuleID[nErrors]/I");
      }
      if(triggerOld->GetBranch("ErrModIndex")){
	triggerOld->SetBranchAddress("ErrModIndex", errorModuleIndex);
	trigger->Branch("errorModuleIndex", errorModuleIndex, "errorModuleIndex[nErrors]/I");
      }
      if(triggerOld->GetBranch("ErrCode")){
	triggerOld->SetBranchAddress("ErrCode", errorCode);
	trigger->Branch("errorCode", errorCode, "errorCode[nErrors]/I");
      }
    }
    if(triggerOld->GetBranch("nTriggerPattern")){
      triggerOld->SetBranchAddress("nTriggerPattern", &nTriggerPattern);
      trigger->Branch("nTriggerPattern", &nTriggerPattern, "nTriggerPattern/I");
      if(triggerOld->GetBranch("TriggerPattern")){
	triggerOld->SetBranchAddress("TriggerPattern", triggerPattern);
	trigger->Branch("triggerPattern", triggerPattern, "triggerPattern[nTriggerPattern]/I");
      }
    }

    for(Int_t i=0; i<triggerOld->GetEntries(); i++){
      triggerOld->GetEvent(i);
      trigger->Fill();
    }
    trigger->Write();

    delete trigger;

    delete errorModuleID;
    delete errorModuleIndex;
    delete errorCode;
    delete triggerPattern;
  }

  Int_t eventNumber = 0;

  if(scalersOld){
    printf("         treeScaler to scalers\n");

    TString sScaler = scalersOld->GetBranch("Scaler")->GetTitle();
    sScaler.Remove(0,sScaler.First("[")+1);
    sScaler.Remove(sScaler.First("]"));
    const Int_t nScaler = sScaler.Atoi();
    Char_t str[256];
    sprintf(str, "scalers[%d]/i", nScaler);

    Int_t eventID = 0;
    UInt_t *scalerArray = new UInt_t[nScaler];

    TTree *scalers = new TTree("scalers", "scalers");
    
    if(scalersOld->GetBranch("eventNumber")){
      scalersOld->SetBranchAddress("eventNumber", &eventNumber);
      scalers->Branch("eventNumber", &eventNumber, "eventNumber/I");
    }
    if(scalersOld->GetBranch("eventID")){
      scalersOld->SetBranchAddress("eventID", &eventID);
      scalers->Branch("eventID", &eventID, "eventID/I");
    }
    if(scalersOld->GetBranch("Scaler")){
      scalersOld->SetBranchAddress("Scaler", scalerArray);
      scalers->Branch("scalers", scalerArray, str);
    }

    for(Int_t i=0; i<scalersOld->GetEntries(); i++){
      scalersOld->GetEvent(i);
      scalers->Fill();
    }
    scalers->Write();

    delete scalers;

    delete scalerArray;
  }

  if(eventParametersOld){
    printf("         treeEventParameters to eventParameters\n");

    UInt_t eventNumberU = 0;
    UChar_t nReconstructed = 0;
    Int_t nReconstructedNew = 0;

    TTree *eventParameters = new TTree("eventParameters", "eventParameters");

    if(eventParametersOld->GetBranch("EventNumber")){
      eventParametersOld->SetBranchAddress("EventNumber", &eventNumberU);
      eventParameters->Branch("eventNumber", &eventNumber, "eventNumber/I");
    }
    if(eventParametersOld->GetBranch("nReconstructed")){
      eventParametersOld->SetBranchAddress("nReconstructed", &nReconstructed);
      eventParameters->Branch("nReconstructed", &nReconstructedNew, "nReconstructed/I");
    }

    for(Int_t i=0; i<eventParametersOld->GetEntries(); i++){
      eventParametersOld->GetEvent(i);
      eventNumber = (Int_t)eventNumberU;
      nReconstructedNew = (Int_t)nReconstructed;
      eventParameters->Fill();
    }
    eventParameters->Write();

    delete eventParameters;
  }

  for(Int_t i=0; i<128; i++){
    centralCrystal[i] = -1;
    centralVeto[i] = -1;
    trackIndex[i] = -1;
  }

  TLorentzVector *particle;

  if(rootinosOld){
    printf("         rootino to rootinos\n");

    TTree *rootinos = new TTree("rootinos", "rootinos");

    if(rootinosOld->GetBranch("nParticles")){
      rootinosOld->SetBranchAddress("nParticles", &nParticles);
      rootinos->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(rootinosOld->GetBranch("particles.")){
	rootinosOld->SetBranchAddress("particles.", &particles);
	rootinos->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	rootinos->Branch("theta", theta, "theta[nParticles]/D");
	rootinos->Branch("phi", phi, "phi[nParticles]/D");
	rootinos->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(rootinosOld->GetBranch("time")){
	rootinosOld->SetBranchAddress("time", time);
	rootinos->Branch("time", time, "time[nParticles]/D");
      }
      if(rootinosOld->GetBranch("clusterSize")){
	rootinosOld->SetBranchAddress("clusterSize", clusterSize);
	rootinos->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      rootinos->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      rootinos->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(rootinosOld->GetBranch("Apparatus")){
	rootinosOld->SetBranchAddress("Apparatus", apparatus);
	rootinos->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(rootinosOld->GetBranch("d_E")){
	rootinosOld->SetBranchAddress("d_E", vetoEnergy);
	rootinos->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(rootinosOld->GetBranch("WC0_E")){
	rootinosOld->SetBranchAddress("WC0_E", MWPC0Energy);
	rootinos->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(rootinosOld->GetBranch("WC1_E")){
	rootinosOld->SetBranchAddress("WC1_E", MWPC1Energy);
	rootinos->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      rootinos->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<rootinosOld->GetEntries(); i++){
      rootinosOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
      }
      rootinos->Fill();
    }
    rootinos->Write();

    delete rootinos;
  }

  if(photonsOld){
    printf("         gamma to photons\n");

    TTree *photons = new TTree("photons", "photons");

    if(photonsOld->GetBranch("nParticles")){
      photonsOld->SetBranchAddress("nParticles", &nParticles);
      photons->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(photonsOld->GetBranch("particles.")){
	photonsOld->SetBranchAddress("particles.", &particles);
	photons->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	photons->Branch("theta", theta, "theta[nParticles]/D");
	photons->Branch("phi", phi, "phi[nParticles]/D");
	photons->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(photonsOld->GetBranch("time")){
	photonsOld->SetBranchAddress("time", time);
	photons->Branch("time", time, "time[nParticles]/D");
      }
      if(photonsOld->GetBranch("clusterSize")){
	photonsOld->SetBranchAddress("clusterSize", clusterSize);
	photons->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      photons->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      photons->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(photonsOld->GetBranch("Apparatus")){
	photonsOld->SetBranchAddress("Apparatus", apparatus);
	photons->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(photonsOld->GetBranch("d_E")){
	photonsOld->SetBranchAddress("d_E", vetoEnergy);
	photons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(photonsOld->GetBranch("WC0_E")){
	photonsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	photons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(photonsOld->GetBranch("WC1_E")){
	photonsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	photons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      photons->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<photonsOld->GetEntries(); i++){
      photonsOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
      }
      photons->Fill();
    }
    photons->Write();

    delete photons;
  }

  if(electronsOld){
    printf("         e- to electrons\n");

    TTree *electrons = new TTree("electrons", "electrons");

    if(electronsOld->GetBranch("nParticles")){
      electronsOld->SetBranchAddress("nParticles", &nParticles);
      electrons->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(electronsOld->GetBranch("particles.")){
	electronsOld->SetBranchAddress("particles.", &particles);
	electrons->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	electrons->Branch("theta", theta, "theta[nParticles]/D");
	electrons->Branch("phi", phi, "phi[nParticles]/D");
	electrons->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(electronsOld->GetBranch("time")){
	electronsOld->SetBranchAddress("time", time);
	electrons->Branch("time", time, "time[nParticles]/D");
      }
      if(electronsOld->GetBranch("clusterSize")){
	electronsOld->SetBranchAddress("clusterSize", clusterSize);
	electrons->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      electrons->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      electrons->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(electronsOld->GetBranch("Apparatus")){
	electronsOld->SetBranchAddress("Apparatus", apparatus);
	electrons->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(electronsOld->GetBranch("d_E")){
	electronsOld->SetBranchAddress("d_E", vetoEnergy);
	electrons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(electronsOld->GetBranch("WC0_E")){
	electronsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	electrons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(electronsOld->GetBranch("WC1_E")){
	electronsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	electrons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      electrons->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<electronsOld->GetEntries(); i++){
      electronsOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
      }
      electrons->Fill();
    }
    electrons->Write();

    delete electrons;
  }

  if(chargedPionsOld){
    printf("         pi+ to chargedPions\n");

    TTree *chargedPions = new TTree("chargedPions", "chargedPions");

    if(chargedPionsOld->GetBranch("nParticles")){
      chargedPionsOld->SetBranchAddress("nParticles", &nParticles);
      chargedPions->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(chargedPionsOld->GetBranch("particles.")){
	chargedPionsOld->SetBranchAddress("particles.", &particles);
	chargedPions->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	chargedPions->Branch("theta", theta, "theta[nParticles]/D");
	chargedPions->Branch("phi", phi, "phi[nParticles]/D");
	chargedPions->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(chargedPionsOld->GetBranch("time")){
	chargedPionsOld->SetBranchAddress("time", time);
	chargedPions->Branch("time", time, "time[nParticles]/D");
      }
      if(chargedPionsOld->GetBranch("clusterSize")){
	chargedPionsOld->SetBranchAddress("clusterSize", clusterSize);
	chargedPions->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      chargedPions->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      chargedPions->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(chargedPionsOld->GetBranch("Apparatus")){
	chargedPionsOld->SetBranchAddress("Apparatus", apparatus);
	chargedPions->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(chargedPionsOld->GetBranch("d_E")){
	chargedPionsOld->SetBranchAddress("d_E", vetoEnergy);
	chargedPions->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(chargedPionsOld->GetBranch("WC0_E")){
	chargedPionsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	chargedPions->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(chargedPionsOld->GetBranch("WC1_E")){
	chargedPionsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	chargedPions->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      chargedPions->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<chargedPionsOld->GetEntries(); i++){
      chargedPionsOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
      }
      chargedPions->Fill();
    }
    chargedPions->Write();

    delete chargedPions;
  }

  if(protonsOld){
    printf("         proton to protons\n");

    TTree *protons = new TTree("protons", "protons");

    if(protonsOld->GetBranch("nParticles")){
      protonsOld->SetBranchAddress("nParticles", &nParticles);
      protons->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(protonsOld->GetBranch("particles.")){
	protonsOld->SetBranchAddress("particles.", &particles);
	protons->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	protons->Branch("theta", theta, "theta[nParticles]/D");
	protons->Branch("phi", phi, "phi[nParticles]/D");
	protons->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(protonsOld->GetBranch("time")){
	protonsOld->SetBranchAddress("time", time);
	protons->Branch("time", time, "time[nParticles]/D");
      }
      if(protonsOld->GetBranch("clusterSize")){
	protonsOld->SetBranchAddress("clusterSize", clusterSize);
	protons->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      protons->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      protons->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(protonsOld->GetBranch("Apparatus")){
	protonsOld->SetBranchAddress("Apparatus", apparatus);
	protons->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(protonsOld->GetBranch("d_E")){
	protonsOld->SetBranchAddress("d_E", vetoEnergy);
	protons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(protonsOld->GetBranch("WC0_E")){
	protonsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	protons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(protonsOld->GetBranch("WC1_E")){
	protonsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	protons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      protons->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<protonsOld->GetEntries(); i++){
      protonsOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
      }
      protons->Fill();
    }
    protons->Write();

    delete protons;
  }

  if(neutronsOld){
    printf("         neutron to neutrons\n");

    TTree *neutrons = new TTree("neutrons", "neutrons");

    if(neutronsOld->GetBranch("nParticles")){
      neutronsOld->SetBranchAddress("nParticles", &nParticles);
      neutrons->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(neutronsOld->GetBranch("particles.")){
	neutronsOld->SetBranchAddress("particles.", &particles);
	neutrons->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	neutrons->Branch("theta", theta, "theta[nParticles]/D");
	neutrons->Branch("phi", phi, "phi[nParticles]/D");
	neutrons->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(neutronsOld->GetBranch("time")){
	neutronsOld->SetBranchAddress("time", time);
	neutrons->Branch("time", time, "time[nParticles]/D");
      }
      if(neutronsOld->GetBranch("clusterSize")){
	neutronsOld->SetBranchAddress("clusterSize", clusterSize);
	neutrons->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      neutrons->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      neutrons->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(neutronsOld->GetBranch("Apparatus")){
	neutronsOld->SetBranchAddress("Apparatus", apparatus);
	neutrons->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(neutronsOld->GetBranch("d_E")){
	neutronsOld->SetBranchAddress("d_E", vetoEnergy);
	neutrons->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(neutronsOld->GetBranch("WC0_E")){
	neutronsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	neutrons->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(neutronsOld->GetBranch("WC1_E")){
	neutronsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	neutrons->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      neutrons->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<neutronsOld->GetEntries(); i++){
      neutronsOld->GetEvent(i);
      neutrons->Fill();
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
      }
    }
    neutrons->Write();

    delete neutrons;
  }

  for(Int_t i=0; i<128; i++) trackIndex[i] = 0;

  UChar_t *nSubParticles    = new UChar_t[64];
  UChar_t *nSubRootinos     = new UChar_t[64];
  UChar_t *nSubPhotons      = new UChar_t[64];
  UChar_t *nSubChargedPions = new UChar_t[64];

  Int_t *nSubParticlesNew    = new Int_t[64];
  Int_t *nSubRootinosNew     = new Int_t[64];
  Int_t *nSubPhotonsNew      = new Int_t[64];
  Int_t *nSubChargedPionsNew = new Int_t[64];

  if(neutralPionsOld){
    printf("         pi0 to neutralPions\n");

    TTree *neutralPions = new TTree("neutralPions", "neutralPions");

    if(neutralPionsOld->GetBranch("nParticles")){
      neutralPionsOld->SetBranchAddress("nParticles", &nParticles);
      neutralPions->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(neutralPionsOld->GetBranch("particles.")){
	neutralPionsOld->SetBranchAddress("particles.", &particles);
	neutralPions->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	neutralPions->Branch("theta", theta, "theta[nParticles]/D");
	neutralPions->Branch("phi", phi, "phi[nParticles]/D");
	neutralPions->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("time")){
	neutralPionsOld->SetBranchAddress("time", time);
	neutralPions->Branch("time", time, "time[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("clusterSize")){
	neutralPionsOld->SetBranchAddress("clusterSize", clusterSize);
	neutralPions->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      neutralPions->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      neutralPions->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(neutralPionsOld->GetBranch("Apparatus")){
	neutralPionsOld->SetBranchAddress("Apparatus", apparatus);
	neutralPions->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("d_E")){
	neutralPionsOld->SetBranchAddress("d_E", vetoEnergy);
	neutralPions->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("WC0_E")){
	neutralPionsOld->SetBranchAddress("WC0_E", MWPC0Energy);
	neutralPions->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("WC1_E")){
	neutralPionsOld->SetBranchAddress("WC1_E", MWPC1Energy);
	neutralPions->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      if(neutralPionsOld->GetBranch("nSubParticles")){
	neutralPionsOld->SetBranchAddress("nSubParticles", nSubParticles);
	neutralPions->Branch("nSubParticles", nSubParticlesNew, "nSubParticles[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("nSubRootinos")){
	neutralPionsOld->SetBranchAddress("nSubRootinos", nSubRootinos);
	neutralPions->Branch("nSubRootinos", nSubRootinosNew, "nSubRootinos[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("nSubPhotons")){
	neutralPionsOld->SetBranchAddress("nSubPhotons", nSubPhotons);
	neutralPions->Branch("nSubPhotons", nSubPhotonsNew, "nSubPhotons[nParticles]/I");
      }
      if(neutralPionsOld->GetBranch("nSubChargedPi")){
	neutralPionsOld->SetBranchAddress("nSubChargedPi", nSubChargedPions);
	neutralPions->Branch("nSubChargedPions", nSubChargedPionsNew, "nSubChargedPions[nParticles]/I");
      }
      neutralPions->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<neutralPionsOld->GetEntries(); i++){
      neutralPionsOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
	else if(apparatus[j] == 3){
	  if(clusterEnergy[j] > 0) detectors[j] += 9;
	  if(vetoEnergy[j] > 0) detectors[j] += 34;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	nSubParticlesNew[j] = (Int_t)nSubParticles[j];
	nSubRootinosNew[j] = (Int_t)nSubRootinos[j];
	nSubPhotonsNew[j] = (Int_t)nSubPhotons[j];
	nSubChargedPionsNew[j] = (Int_t)nSubChargedPions[j];
      }
      neutralPions->Fill();
    }
    neutralPions->Write();

    delete neutralPions;
  }

  if(etasOld){
    printf("         eta to etas\n");

    TTree *etas = new TTree("etas", "etas");

    if(etasOld->GetBranch("nParticles")){
      etasOld->SetBranchAddress("nParticles", &nParticles);
      etas->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(etasOld->GetBranch("particles.")){
	etasOld->SetBranchAddress("particles.", &particles);
	etas->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	etas->Branch("theta", theta, "theta[nParticles]/D");
	etas->Branch("phi", phi, "phi[nParticles]/D");
	etas->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(etasOld->GetBranch("time")){
	etasOld->SetBranchAddress("time", time);
	etas->Branch("time", time, "time[nParticles]/D");
      }
      if(etasOld->GetBranch("clusterSize")){
	etasOld->SetBranchAddress("clusterSize", clusterSize);
	etas->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      etas->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      etas->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(etasOld->GetBranch("Apparatus")){
	etasOld->SetBranchAddress("Apparatus", apparatus);
	etas->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(etasOld->GetBranch("d_E")){
	etasOld->SetBranchAddress("d_E", vetoEnergy);
	etas->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(etasOld->GetBranch("WC0_E")){
	etasOld->SetBranchAddress("WC0_E", MWPC0Energy);
	etas->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(etasOld->GetBranch("WC1_E")){
	etasOld->SetBranchAddress("WC1_E", MWPC1Energy);
	etas->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      if(etasOld->GetBranch("nSubParticles")){
	etasOld->SetBranchAddress("nSubParticles", nSubParticles);
	etas->Branch("nSubParticles", nSubParticlesNew, "nSubParticles[nParticles]/I");
      }
      if(etasOld->GetBranch("nSubRootinos")){
	etasOld->SetBranchAddress("nSubRootinos", nSubRootinos);
	etas->Branch("nSubRootinos", nSubRootinosNew, "nSubRootinos[nParticles]/I");
      }
      if(etasOld->GetBranch("nSubPhotons")){
	etasOld->SetBranchAddress("nSubPhotons", nSubPhotons);
	etas->Branch("nSubPhotons", nSubPhotonsNew, "nSubPhotons[nParticles]/I");
      }
      if(etasOld->GetBranch("nSubChargedPi")){
	etasOld->SetBranchAddress("nSubChargedPi", nSubChargedPions);
	etas->Branch("nSubChargedPions", nSubChargedPionsNew, "nSubChargedPions[nParticles]/I");
      }
      etas->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<etasOld->GetEntries(); i++){
      etasOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
	else if(apparatus[j] == 3){
	  if(clusterEnergy[j] > 0) detectors[j] += 9;
	  if(vetoEnergy[j] > 0) detectors[j] += 34;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	nSubParticlesNew[j] = (Int_t)nSubParticles[j];
	nSubRootinosNew[j] = (Int_t)nSubRootinos[j];
	nSubPhotonsNew[j] = (Int_t)nSubPhotons[j];
	nSubChargedPionsNew[j] = (Int_t)nSubChargedPions[j];
      }
      etas->Fill();
    }
    etas->Write();

    delete etas;
  }

  if(etaPrimesOld){
    printf("         eta' to etaPrimes\n");

    TTree *etaPrimes = new TTree("etaPrimes", "etaPrimes");

    if(etaPrimesOld->GetBranch("nParticles")){
      etaPrimesOld->SetBranchAddress("nParticles", &nParticles);
      etaPrimes->Branch("nParticles", &nParticlesNew, "nParticles/I");
      if(etaPrimesOld->GetBranch("particles.")){
	etaPrimesOld->SetBranchAddress("particles.", &particles);
	etaPrimes->Branch("clusterEnergy", clusterEnergy, "clusterEnergy[nParticles]/D");
	etaPrimes->Branch("theta", theta, "theta[nParticles]/D");
	etaPrimes->Branch("phi", phi, "phi[nParticles]/D");
	etaPrimes->Branch("mass", mass, "mass[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("time")){
	etaPrimesOld->SetBranchAddress("time", time);
	etaPrimes->Branch("time", time, "time[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("clusterSize")){
	etaPrimesOld->SetBranchAddress("clusterSize", clusterSize);
	etaPrimes->Branch("clusterSize", clusterSizeNew, "clusterSize[nParticles]/I");
      }
      etaPrimes->Branch("centralCrystal", centralCrystal, "centralCrystal[nParticles]/I");
      etaPrimes->Branch("centralVeto", centralVeto, "centralVeto[nParticles]/I");
      if(etaPrimesOld->GetBranch("Apparatus")){
	etaPrimesOld->SetBranchAddress("Apparatus", apparatus);
	etaPrimes->Branch("detectors", detectors, "detectors[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("d_E")){
	etaPrimesOld->SetBranchAddress("d_E", vetoEnergy);
	etaPrimes->Branch("vetoEnergy", vetoEnergy, "vetoEnergy[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("WC0_E")){
	etaPrimesOld->SetBranchAddress("WC0_E", MWPC0Energy);
	etaPrimes->Branch("MWPC0Energy", MWPC0Energy, "MWPC0Energy[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("WC1_E")){
	etaPrimesOld->SetBranchAddress("WC1_E", MWPC1Energy);
	etaPrimes->Branch("MWPC1Energy", MWPC1Energy, "MWPC1Energy[nParticles]/D");
      }
      if(etaPrimesOld->GetBranch("nSubParticles")){
	etaPrimesOld->SetBranchAddress("nSubParticles", nSubParticles);
	etaPrimes->Branch("nSubParticles", nSubParticlesNew, "nSubParticles[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("nSubRootinos")){
	etaPrimesOld->SetBranchAddress("nSubRootinos", nSubRootinos);
	etaPrimes->Branch("nSubRootinos", nSubRootinosNew, "nSubRootinos[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("nSubPhotons")){
	etaPrimesOld->SetBranchAddress("nSubPhotons", nSubPhotons);
	etaPrimes->Branch("nSubPhotons", nSubPhotonsNew, "nSubPhotons[nParticles]/I");
      }
      if(etaPrimesOld->GetBranch("nSubChargedPi")){
	etaPrimesOld->SetBranchAddress("nSubChargedPi", nSubChargedPions);
	etaPrimes->Branch("nSubChargedPions", nSubChargedPionsNew, "nSubChargedPions[nParticles]/I");
      }
      etaPrimes->Branch("trackIndex", trackIndex, "trackIndex[nParticles]/I");
    }

    for(Int_t i=0; i<etaPrimesOld->GetEntries(); i++){
      etaPrimesOld->GetEvent(i);
      nParticlesNew = (Int_t)nParticles;
      for(Int_t j=0; j<nParticlesNew; j++){
	particle = (TLorentzVector*)particles->At(j);
	clusterEnergy[j] = (particle->E()-particle->M());
	theta[j] = (TMath::RadToDeg()*particle->Theta());
	phi[j] = (TMath::RadToDeg()*particle->Phi());
	mass[j] = particle->M();
	clusterSizeNew[j] = (Int_t)clusterSize[j];
	detectors[j] = 0;
	if(apparatus[j] == 1){
	  if(clusterEnergy[j] > 0) detectors[j] += 1;
	  if(vetoEnergy[j] > 0) detectors[j] += 2;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	else if(apparatus[j] == 2){
	  if(clusterEnergy[j] > 0) detectors[j] += 8;
	  if(vetoEnergy[j] > 0) detectors[j] += 32;
	}
	else if(apparatus[j] == 3){
	  if(clusterEnergy[j] > 0) detectors[j] += 9;
	  if(vetoEnergy[j] > 0) detectors[j] += 34;
	  if(MWPC0Energy[j] > 0 || MWPC1Energy[j] > 0) detectors[j] += 4;
	}
	nSubParticlesNew[j] = (Int_t)nSubParticles[j];
	nSubRootinosNew[j] = (Int_t)nSubRootinos[j];
	nSubPhotonsNew[j] = (Int_t)nSubPhotons[j];
	nSubChargedPionsNew[j] = (Int_t)nSubChargedPions[j];
      }
      etaPrimes->Fill();
    }
    etaPrimes->Write();

    delete etaPrimes;
  }

  delete clusterEnergy;
  delete theta;
  delete phi;
  delete mass;
  delete time;
  delete clusterSize;
  delete clusterSizeNew;
  delete centralCrystal;
  delete centralVeto;
  delete apparatus;
  delete detectors;
  delete vetoEnergy;
  delete MWPC0Energy;
  delete MWPC1Energy;
  delete trackIndex;

  delete nSubParticles;
  delete nSubRootinos;
  delete nSubPhotons;
  delete nSubChargedPions;

  delete nSubParticlesNew;
  delete nSubRootinosNew;
  delete nSubPhotonsNew;
  delete nSubChargedPionsNew;

  if(fOld.GetListOfKeys()->Contains("CheckCB")){
    fNew.mkdir("CheckCB");
    fNew.cd("CheckCB");

    TH2F *Check_CBdE_E = (TH2F*)fOld.Get("CheckCB/Check_CBdE_E");
    Check_CBdE_E->Write();
    TH2F *Check_CBPhiCorr = (TH2F*)fOld.Get("CheckCB/Check_CBPhiCorr");	
    Check_CBPhiCorr->Write();
    TH2F *Check_CBdE_E_1PID = (TH2F*)fOld.Get("CheckCB/Check_CBdE_E_1PID");	
    Check_CBdE_E_1PID->Write();
    TH2F *Check_CBPhiCorr_1PID = (TH2F*)fOld.Get("CheckCB/Check_CBPhiCorr_1PID");	
    Check_CBPhiCorr_1PID->Write();
    TH2F *Check_CBdE_E_pi0 = (TH2F*)fOld.Get("CheckCB/Check_CBdE_E_pi0");	
    Check_CBdE_E_pi0->Write();
    TH2F *Check_CBPhiCorr_pi0 = (TH2F*)fOld.Get("CheckCB/Check_CBPhiCorr_pi0");	
    Check_CBPhiCorr_pi0->Write();

    delete Check_CBdE_E;
    delete Check_CBPhiCorr;
    delete Check_CBdE_E_1PID;
    delete Check_CBPhiCorr_1PID;
    delete Check_CBdE_E_pi0;
    delete Check_CBPhiCorr_pi0;

    fNew.cd();
  }

  if(fOld.GetListOfKeys()->Contains("CheckTAPS")){
    fNew.mkdir("CheckTAPS");
    fNew.cd("CheckTAPS");

    TH2F *Check_TAPSdE_E = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSdE_E");
    Check_TAPSdE_E->Write();
    TH2F *Check_TAPSPhiCorr = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSPhiCorr");
    Check_TAPSPhiCorr->Write();
    TH2F *Check_TAPSdE_E_1Veto = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSdE_E_1PID");
    Check_TAPSdE_E_1Veto->SetName("Check_TAPSdE_E_1Veto");
    Check_TAPSdE_E_1Veto->Write();
    TH2F *Check_TAPSPhiCorr_1Veto = (TH2F*)fOld.Get("CheckTAPS/Check_TAPSPhiCorr_1Veto");
    Check_TAPSPhiCorr_1Veto->Write();

    delete Check_TAPSdE_E;
    delete Check_TAPSPhiCorr;
    delete Check_TAPSdE_E_1Veto;
    delete Check_TAPSPhiCorr_1Veto;

    fNew.cd();
  }

  if(fOld.GetListOfKeys()->Contains("CheckHitPatterns")){
    fNew.mkdir("CheckHitPatterns");
    fNew.cd("CheckHitPatterns");	

    TH2F *Check_CBHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_CBHits");
    Check_CBHits->Write();
    TH2F *Check_CBADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_CBADCHits");
    Check_CBADCHits->Write();
    TH2F *Check_CBTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_CBTDCHits");
    Check_CBTDCHits->Write();

    TH2F *Check_PIDHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_PIDHits");
    Check_PIDHits->Write();
    TH2F *Check_PIDADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_PIDADCHits");
    Check_PIDADCHits->Write();
    TH2F *Check_PIDTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_PIDTDCHits");
    Check_PIDTDCHits->Write();

    TH2F *Check_TAPSHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_TAPSHits");
    Check_TAPSHits->Write();
    TH2F *Check_TAPSADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_TAPSADCHits");
    Check_TAPSADCHits->Write();
    TH2F *Check_TAPSTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_TAPSTDCHits");
    Check_TAPSTDCHits->Write();

    TH2F *Check_VetoHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_VetoHits");
    Check_VetoHits->Write();
    TH2F *Check_VetoADCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_VetoADCHits");
    Check_VetoADCHits->Write();
    TH2F *Check_VetoTDCHits = (TH2F*)fOld.Get("CheckHitPatterns/Check_VetoTDCHits");
    Check_VetoTDCHits->Write();

    delete Check_CBHits;
    delete Check_CBADCHits;
    delete Check_CBTDCHits;

    delete Check_PIDHits;
    delete Check_PIDADCHits;
    delete Check_PIDTDCHits;

    delete Check_TAPSHits;
    delete Check_TAPSADCHits;
    delete Check_TAPSTDCHits;

    delete Check_VetoHits;
    delete Check_VetoADCHits;
    delete Check_VetoTDCHits;

    fNew.cd();
  }

  if(fOld.GetListOfKeys()->Contains("GoAT_File")){
    TNamed *GoAT_File = (TNamed*)fOld.Get("GoAT_File");
    GoAT_File->Write();
    
    delete GoAT_File;
  }

  if(fOld.GetListOfKeys()->Contains("CountScalerValid")){
    TH1I *CountScalerValid = (TH1I*)fOld.Get("CountScalerValid");
    CountScalerValid->Write();

    delete CountScalerValid;
  }

  delete taggerOld;
  delete tracksOld;
  delete detectorHitsOld;
  delete linPolOld;
  delete triggerOld;
  delete scalersOld;

  delete eventParametersOld;
  delete rootinosOld;
  delete photonsOld;
  delete electronsOld;
  delete chargedPionsOld;
  delete protonsOld;
  delete neutronsOld;
  delete neutralPionsOld;
  delete etasOld;
  delete etaPrimesOld;

  delete particles;

  fNew.Close();
  fOld.Close();

  timer.Stop();
  printf("Renaming took %.1f seconds.\n\n", timer.RealTime());

}
