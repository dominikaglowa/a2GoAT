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

  TH1D *hBkg1AccScal = (TH1D*)fBkg1.Get("TaggerAccScal");
  TH1D *hBkg1LiveTime = (TH1D*)fBkg1.Get("LiveTimeScal");

  TH1D *hBackAccScal = (TH1D*)hBkg1AccScal->Clone("hBackAccScal");

  Double_t dBackClock = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Clock"));
  Double_t dBackInhib = hBkg1LiveTime->GetBinContent(hBkg1LiveTime->GetXaxis()->FindBin("Inhibited"));
  
  if(sBkg2 != ""){
    TFile fBkg2(sBkg2,"READ");
    gROOT->cd();

    TH1D *hBkg2AccScal = (TH1D*)fBkg2.Get("TaggerAccScal");
    TH1D *hBkg2LiveTime = (TH1D*)fBkg2.Get("LiveTimeScal");

    hBackAccScal->Add(hBkg2AccScal);

    dBackClock += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Clock"));
    dBackInhib += hBkg2LiveTime->GetBinContent(hBkg2LiveTime->GetXaxis()->FindBin("Inhibited"));
  }

  TFile fOut("TaggEff_sn120_4c.root","RECREATE");

  TH1D *hEffAllHits = (TH1D*)hBeamAllHits->Clone("TaggEffAllHits");
  TH1D *hEffSingles = (TH1D*)hBeamSingles->Clone("TaggEffSingles");
  TH1D *hEffDoubles = (TH1D*)hBeamDoubles->Clone("TaggEffDoubles");
  TH1D *hEffAccScal = (TH1D*)hBeamAccScal->Clone("TaggEffAccScal");
  hEffAccScal->Sumw2();

  if(bFreeScal){
    hEffAccScal->Add(hBackAccScal,(-dBeamClock/dBackClock));
    hEffAccScal->Scale(dBeamInhib/dBeamClock);
  }
  else hEffAccScal->Add(hBackAccScal,(-dBeamInhib/dBackInhib));

  hEffAllHits->Divide(hEffAccScal);
  hEffSingles->Divide(hEffAccScal);
  hEffDoubles->Divide(hEffAccScal);
		     
  fOut.Write();
  fOut.Close();
}
