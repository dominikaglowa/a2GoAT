void DrawCluster(TString sData, Int_t iMinEvn=0, Int_t iMinPart=0, Int_t iMinNaI=0, Int_t iMinBaF=0, TString sNaI="NaI.dat", TString sBaF="BaF2-PbWO4.dat"){

  gROOT->Clear();

  TGraph *gBox = new TGraph(4);
  gBox->SetName("gBox");
  gBox->SetPoint(0,  0,-180);
  gBox->SetPoint(1,  0, 180);
  gBox->SetPoint(2,180,-180);
  gBox->SetPoint(3,180, 180);
  gBox->SetMarkerStyle(1);
  gBox->SetMarkerColor(0);

  TGraph *gBox2 = new TGraph();
  gBox2->SetName("gBox2");
  gBox2->SetMarkerStyle(24);

  TGraph2D *gBox3 = new TGraph2D(8);
  gBox3->SetName("gBox3");
  gBox3->SetPoint(0,-60,-60,-50);
  gBox3->SetPoint(1,-60,-60,150);
  gBox3->SetPoint(2,-60, 60,-50);
  gBox3->SetPoint(3,-60, 60,150);
  gBox3->SetPoint(4, 60,-60,-50);
  gBox3->SetPoint(5, 60,-60,150);
  gBox3->SetPoint(6, 60, 60,-50);
  gBox3->SetPoint(7, 60, 60,150);
  gBox3->SetMarkerStyle(1);
  gBox3->SetMarkerColor(0);
  
  Int_t i,j;

  TTree *tNaI = new TTree("tNaI","NaI");
  tNaI->ReadFile(sNaI,"Ch:X:Y:Z");
  Int_t iNaI = tNaI->Draw("X:Y:Z","","goff");
  Double_t *dNaIX = tNaI->GetV1();
  Double_t *dNaIY = tNaI->GetV2();
  Double_t *dNaIZ = tNaI->GetV3();

  const Int_t cNaI = iNaI;

  Double_t dNaITh[cNaI];
  Double_t dNaIPh[cNaI];
  for(i=0; i<iNaI; i++){
    dNaITh[i] = (TMath::RadToDeg()*TMath::ACos(dNaIZ[i]/TMath::Sqrt(TMath::Power(dNaIX[i],2)+TMath::Power(dNaIY[i],2)+TMath::Power(dNaIZ[i],2))));
    dNaIPh[i] = (TMath::RadToDeg()*TMath::ATan2(dNaIY[i],dNaIX[i]));
  }

  TTree *tBaF = new TTree("tBaF","BaF");
  tBaF->ReadFile(sBaF,"Ch:X:Y:Z");
  Int_t iBaF = tBaF->Draw("X:Y:Z","","goff");
  Double_t *dBaFX = tBaF->GetV1();
  Double_t *dBaFY = tBaF->GetV2();
  Double_t *dBaFZ = tBaF->GetV3();

  const Int_t cBaF = iBaF;

  Double_t dBaFTh[cBaF];
  Double_t dBaFPh[cBaF];
  for(i=0; i<iBaF; i++){
    dBaFTh[i] = (TMath::RadToDeg()*TMath::ACos(dBaFZ[i]/TMath::Sqrt(TMath::Power(dBaFX[i],2)+TMath::Power(dBaFY[i],2)+TMath::Power(dBaFZ[i],2))));
    dBaFPh[i] = (TMath::RadToDeg()*TMath::ATan2(dBaFY[i],dBaFX[i]));
    gBox2->SetPoint(i,dBaFX[i],dBaFY[i]);
  }

  TFile fData(sData,"READ");
  TTree *tData = (TTree*)fData.Get("detectorHits");

  Int_t iNNaIHits;
  Int_t* iNaIHits;
  Int_t* iNaICluster;

  iNaIHits = new Int_t[720];
  iNaICluster = new Int_t[720];

  tData->SetBranchAddress("nNaIHits",&iNNaIHits);
  tData->SetBranchAddress("NaIHits",iNaIHits);
  tData->SetBranchAddress("NaICluster",iNaICluster);

  Int_t iNBaFHits;
  Int_t* iBaFHits;
  Int_t* iBaFCluster;

  iBaFHits = new Int_t[720];
  iBaFCluster = new Int_t[720];

  tData->SetBranchAddress("nBaF2Hits",&iNBaFHits);
  tData->SetBranchAddress("BaF2Hits",iBaFHits);
  tData->SetBranchAddress("BaF2Cluster",iBaFCluster);

  TTree *tEvent = (TTree*)fData.Get("tracks");

  Int_t iNParticles;
  UChar_t* iClusterSize;
  UChar_t* iApparatus;

  iClusterSize = new UChar_t[128];
  iApparatus = new UChar_t[128];

  tEvent->SetBranchAddress("nParticles",&iNParticles);
  tEvent->SetBranchAddress("clusterSize",iClusterSize);
  tEvent->SetBranchAddress("apparatus",iApparatus);

  Bool_t bPart = true;
  Bool_t bNaI = true;
  Bool_t bBaF = true;

  Int_t iEvnN;
  for(iEvnN=iMinEvn; iEvnN<(tEvent->GetEntries()); iEvnN++){
    tEvent->GetEvent(iEvnN);

    if(iMinPart) bPart = false;
    if(iMinNaI) bNaI = false;
    if(iMinBaF) bBaF = false;

    if(iNParticles >= iMinPart) bPart = true;
    for(j=0; j<iNParticles; j++){
      if(iClusterSize[j] >= iMinNaI && iApparatus[j] == 1) bNaI = true;
      if(iClusterSize[j] >= iMinBaF && iApparatus[j] == 2) bBaF = true;
    }
    if(bPart && bNaI && bBaF) break;
  }

  cout << endl << "Plotting event " << iEvnN << endl << endl;

  tData->GetEvent(iEvnN);

  j=-1;
  for(i=0; i<iNNaIHits; i++){
    if(iNaICluster[i] > j) j = iNaICluster[i];
  }
  const Int_t iNNaICluster = j+1;

  j=-1;
  for(i=0; i<iNBaFHits; i++){
    if(iBaFCluster[i] > j) j = iBaFCluster[i];
  }
  const Int_t iNBaFCluster = j+1;

  cout << iNNaICluster << " NaI Clusters, " << iNBaFCluster << " BaF Clusters" << endl;

  gROOT->cd();

  const Int_t iNData = iNNaICluster+iNBaFCluster+1;
  Int_t iDataN;

  TGraph *gData[iNData];
  for(i=0; i<iNData; i++){
    gData[i] = new TGraph();
    gData[i]->SetName(Form("gData_%d",i));
    gData[i]->SetMarkerStyle(20);
    gData[i]->SetMarkerColor(i+1);
  }

  TGraph *gData2[iNData];
  for(i=0; i<iNData; i++){
    gData2[i] = new TGraph();
    gData2[i]->SetName(Form("gData_%d",i));
    gData2[i]->SetMarkerStyle(20);
    gData2[i]->SetMarkerColor(i+1);
  }

  TGraph2D *gData3[iNData];
  for(i=0; i<iNData; i++){
    gData3[i] = new TGraph2D();
    gData3[i]->SetName(Form("gData2_%d",i));
    gData3[i]->SetMarkerStyle(20);
    gData3[i]->SetMarkerColor(i+1);
  }

  cout << endl << "NaI" << endl << endl;
  for(Int_t i=0; i<iNNaIHits; i++){
    cout << i << "\t" << iNaIHits[i] << "\t" << iNaICluster[i] << "\t" << dNaIX[iNaIHits[i]] << "\t" << dNaIY[iNaIHits[i]] << "\t" << dNaIZ[iNaIHits[i]] << endl;
    iDataN = iNaICluster[i]+1;
    gData[iDataN]->SetPoint(gData[iDataN]->GetN(),dNaITh[iNaIHits[i]],dNaIPh[iNaIHits[i]]);
    gData3[iDataN]->SetPoint(gData3[iDataN]->GetN(),dNaIX[iNaIHits[i]],dNaIY[iNaIHits[i]],dNaIZ[iNaIHits[i]]);
  }
  cout << endl << "BaF" << endl << endl;
  for(Int_t i=0; i<iNBaFHits; i++){
    cout << i << "\t" << iBaFHits[i] << "\t" << iBaFCluster[i] << "\t" << dBaFX[iBaFHits[i]] << "\t" << dBaFY[iBaFHits[i]] << "\t" << dBaFZ[iBaFHits[i]] << "\t" << iBaFCluster[i]+iNNaICluster+1 << endl;
    if(iBaFCluster[i] < 0) iDataN = 0;
    else iDataN = iBaFCluster[i]+iNNaICluster+1;
    gData[iDataN]->SetPoint(gData[iDataN]->GetN(),dBaFTh[iBaFHits[i]],dBaFPh[iBaFHits[i]]);
    gData2[iDataN]->SetPoint(gData2[iDataN]->GetN(),dBaFX[iBaFHits[i]],dBaFY[iBaFHits[i]]);
    gData3[iDataN]->SetPoint(gData3[iDataN]->GetN(),dBaFX[iBaFHits[i]],dBaFY[iBaFHits[i]],dBaFZ[iBaFHits[i]]);
    gBox2->SetPoint(iBaFHits[i],-100,0);
  }

  if(gROOT->FindObject("c1")) delete (TCanvas*)gROOT->FindObject("c1");
  if(gROOT->FindObject("c2")) delete (TCanvas*)gROOT->FindObject("c2");
  if(gROOT->FindObject("c3")) delete (TCanvas*)gROOT->FindObject("c3");

  TCanvas *c1 = new TCanvas("c1","Clusters in Theta and Phi",0,0,800,600);
  gBox->Draw("ap");
  gBox->SetTitle("Clusters in Theta and Phi");
  gBox->GetXaxis()->SetTitle("Theta (deg)");
  gBox->GetXaxis()->SetTitleOffset(1.2);
  gBox->GetYaxis()->SetTitle("Phi (deg)");
  gBox->GetYaxis()->SetTitleOffset(1.2);
  for(i=0; i<iNData; i++){
    if(gData[i]->GetN() == 1) gData[i]->SetPoint(1,-180,0);
    if(gData[i]->GetN() == 0) gData[i]->SetPoint(0,-180,0);
    gData[i]->Draw("same p");
  }  

  TCanvas *c2 = new TCanvas("c2","Clusters in TAPS",800,0,600,600);
  gBox2->Draw("ap");
  gBox2->SetTitle("Clusters in TAPS");
  gBox2->GetXaxis()->SetRangeUser(-70,70);
  gBox2->GetXaxis()->SetTitle("X (cm)");
  gBox2->GetXaxis()->SetTitleOffset(1.2);
  gBox2->GetYaxis()->SetRangeUser(-70,70);
  gBox2->GetYaxis()->SetTitle("Y (cm)");
  gBox2->GetYaxis()->SetTitleOffset(1.2);
  for(i=0; i<iNData; i++){
    if(gData2[i]->GetN() == 1) gData2[i]->SetPoint(1,-100,0);
    if(gData2[i]->GetN() == 0) gData2[i]->SetPoint(0,-100,0);
    gData2[i]->Draw("same p");
  }

  TCanvas *c3 = new TCanvas("c3","Clusters in X, Y, and Z",300,300,800,600);
  gBox3->Draw("p");
  gBox3->SetTitle("Clusters in X, Y, and Z");
  gBox3->GetXaxis()->SetTitle("X (cm)");
  gBox3->GetXaxis()->SetTitleOffset(1.5);
  gBox3->GetYaxis()->SetTitle("Y (cm)");
  gBox3->GetYaxis()->SetTitleOffset(1.5);
  gBox3->GetZaxis()->SetTitle("Z (cm)");
  gBox3->GetZaxis()->SetTitleOffset(1.5);
  for(i=0; i<iNData; i++){
    if(gData3[i]->GetN() == 1) gData3[i]->SetPoint(1,0,0,-1000);
    if(gData3[i]->GetN() == 0) gData3[i]->SetPoint(0,0,0,-1000);
    gData3[i]->Draw("same p");
  }
}
