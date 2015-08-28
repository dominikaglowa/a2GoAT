void TaggEffBgSub(TString sBeam, TString sBkg1, TString sBkg2="", Bool_t bFreeScal=false){
  TFile fBeam(sBeam.Data(),"READ");
  TFile fBkg1(sBkg1.Data(),"READ");  
  gROOT->cd();

  TH1D *hBeamAllHits = (TH1D*)fBeam.Get("TaggerAllHits");
  TH1D *hBeamSingles = (TH1D*)fBeam.Get("TaggerSingles");
  TH1D *hBeamDoubles = (TH1D*)fBeam.Get("TaggerDoubles");
  TH1D *hBeamAccScal = (TH1D*)fBeam.Get("TaggerAccScal");
  TH1D *hBeamLiveTime = (TH1D*)fBeam.Get("LiveTimeScal");

  Double_t dBeamClock = hBeamLiveTime->GetBinContent(hBeamLiveTime->GetXaxis()->FindBin("Clock"));
  Double_t dBeamInhib = hBeamLiveTime->GetBinContent(hBeamLiveTime->GetXaxis()->FindBin("Inhibited"));

  cout << "Clock=" << dBeamClock << "   Inhib=" << dBeamInhib << endl;

  TH1D *hBkg1AccScal = (TH1D*)fBkg1.Get("TaggerAccScal");
  TH1D *hBkg1LiveTime = (TH1D*)fBkg1.Get("LiveTimeScal");

  TH1D *hBackAccScal = (TH1D*)hBkg1AccScal->Clone("hBackAccScal");

  Double_t dBackClock = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Clock"));
  Double_t dBackInhib = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Inhibited"));
  
  if(sBkg2 != ""){
    TFile fBkg2(sBkg2.Data(),"READ");
    gROOT->cd();

    TH1D *hBkg2AccScal = (TH1D*)fBkg2.Get("TaggerAccScal");
    TH1D *hBkg2LiveTime = (TH1D*)fBkg2.Get("LiveTimeScal");

    hBackAccScal->Add(hBkg2AccScal);

    dBackClock += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Clock"));
    dBackInhib += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Inhibited"));
  }

  cout << "Bcg Clock=" << dBackClock << "   Gcg Inhib=" << dBackInhib << endl;
  TFile fOut("TaggEff_sn124.root","RECREATE");

  TH1D *hEffAllHits = (TH1D*)hBeamAllHits->Clone("hEffAllHits");
  TH1D *hEffAllHits2 = (TH1D*)hBeamAllHits->Clone("hEffAllHits2");
  TH1D *hEffSingles = (TH1D*)hBeamSingles->Clone("hEffSingles");
  TH1D *hEffDoubles = (TH1D*)hBeamDoubles->Clone("hEffDoubles");
  TH1D *hEffAccScal = (TH1D*)hBeamAccScal->Clone("hEffAccScal");
  TH1D *hEffAccScal2 = (TH1D*)hBeamAccScal->Clone("hEffAccScal2");
  TH1D *hBckgAccScal = (TH1D*)hBackAccScal->Clone("hBckgAccScal");
  hEffAccScal->Sumw2();
  hEffAccScal2->Sumw2();

  if(bFreeScal){
    hEffAccScal->Add(hBackAccScal,(-dBeamClock/dBackClock));
    hEffAccScal->Scale(dBeamInhib/dBeamClock);
  }
  else hEffAccScal->Add(hBackAccScal,(-dBeamInhib/dBackClock));

  hEffAllHits->Divide(hEffAccScal);
  hEffSingles->Divide(hEffAccScal);
  hEffDoubles->Divide(hEffAccScal);
  hEffAllHits2->Divide(hEffAccScal2);
		     
  fOut.Write();
  fOut.Close();
}
