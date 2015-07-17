{
  TChain chain("tagger");
   chain.Add("./goat_empty_vertex_corr/GoAT_CB_*");

   TH1F *hist = new TH1F("hist","taggedChannel",100,0,350);

   double taggedChannel[200]; 
   int ntagged;
   chain.SetBranchAddress("taggedChannel", &taggedChannel);
   chain.SetBranchAddress("nTagged", &ntagged);
   Int_t ntaggedChannel = chain.GetEntries();

   // for (Int_t i=0; i<ntaggedEnergy; i++){
   for (Int_t i=0; i<ntaggedChannel; i++){
     chain.GetEvent(i);
     cout << i << " " << taggedChannel[0] <<endl;
     for (int j=0; j<ntagged; j++) {
       hist->Fill(taggedChannel[j]);
     }
   }
   TFile f2("goat_total_empty_taggedchannel2.root","RECREATE");
   f2.cd();
   hist->Write();
   f2.Close();
}
