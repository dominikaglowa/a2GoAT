{
  TChain chain("tagger");
   chain.Add("./goat_empty_vertex_corr/GoAT_CB_*");

   TH1F *hist = new TH1F("hist","taggedEnergy",100,0,1000);

   double taggedEnergy[100]; 
   int ntagged;
   chain.SetBranchAddress("taggedEnergy", &taggedEnergy);
   chain.SetBranchAddress("nTagged", &ntagged);
   Int_t ntaggedEnergy = chain.GetEntries();

   // for (Int_t i=0; i<ntaggedEnergy; i++){
   for (Int_t i=0; i<ntaggedEnergy; i++){
     chain.GetEvent(i);
     cout << i << " " << taggedEnergy[0] <<endl;
     for (int j=0; j<ntagged; j++) {
       hist->Fill(taggedEnergy[j]);
     }
   }
   TFile f2("goat_total_empty_energychannel.root","RECREATE");
   f2.cd();
   hist->Write();
   f2.Close();
}
