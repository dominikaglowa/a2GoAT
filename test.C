{
TFile *_file0 = TFile::Open("Physics_CB_331_renamed.root");
TH1F *h2dFile_pr=(TH1F*)_file0->Get("Test_histo_h1_prompt");
TH1F *h2dFile_rn=(TH1F*)_file0->Get("Test_histo_h1_random");
  double random_factor = 0.055;
h2dFile_pr->Add(h2dFile_rn,-random_factor);
h2dFile_pr->Draw();
TH1F *h2dFile=(TH1F*)_file0->Get("DeltaE_CM_D_2g");
h2dFile->SetLineColor(2);
h2dFile->Rebin(4);
h2dFile->Draw("SAME");

}
