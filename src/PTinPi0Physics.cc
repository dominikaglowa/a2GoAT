#include "PTinPi0Physics.h"

PTinPi0Physics::PTinPi0Physics()
{ 
    time 	= new GH1("time", 	"time", 	1400, -700, 700);
    time_cut 	= new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

    time_2g 	= new GH1("time_2g",	"time_2g", 	1400, -700, 700);
    time_2g_cut = new GH1("time_2g_cut","time_2g_cut", 	1400, -700, 700);

    IM 		= new GH1("IM", 	"IM", 		400,   0, 400);
    IM_2g 	= new GH1("IM_2g", 	"IM_2g", 	400,   0, 400);
  
    MM		= new GH1("MM", 	"MM", 	 	400,   107800, 108200);     
    MM_2g	= new GH1("MM_2g", 	"MM_2g", 	400,  107800, 108200);

    TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);

    //my added histograms

    MMom		= new GH1("MMom", 	"MMom;q(fm^{-1})", 	 	400, 0., 2.0);     
    MMom_2g	= new GH1("MMom_2g", 	"MMom_2g;q(fm^{-1})", 	400,   0, 2.0);
    MMomD		= new GH1("MMomD", 	"MMom;q(fm^{-1})", 	 	400, 0., 2.0);     
    MMomD_2g	= new GH1("MMomD_2g", 	"MMom_2g;q(fm^{-1})", 	400,   0, 2.0);

    DeltaE_CM_D_2g	= new GH1("DeltaE_CM_D_2g", 	"DeltaE_{CM} 2 #gamma ; MeV", 	400, -60., 60.); 
    Test_histo_g1  =  new GH1("Test_histo_g1", "DeltaE_{CM} 2 #gamma ; MeV",   100, -60., 60.);
    Test_histo_h1_prompt  =  new TH1F("Test_histo_h1_prompt", "Promp DeltaE_{CM} 2 #gamma ; MeV",   100, -60., 60.);
    Test_histo_h1_random  =  new TH1F("Test_histo_h1_random", "Random DeltaE_{CM} 2 #gamma ; MeV",   100, -60., 60.);
    
    // int bin_q=300;
    // int bin_e=23;
    // char Title[60];
    // char Title2[160];
    // for (int i=0; i<bin_q; i++) {
    //   for (int j=0; j<bin_e; j++) {
    // 	//	sprintf(Title,"pimissen_q_%d_%.3f_%.3f",j,i*0.005,(i*0.005+0.005)); //original q binning
    // 		sprintf(Title,"pimissen_q_%d_%.3f_%.3f",j,i*0.02,(i*0.02+0.02)); //rebin 4
    // 	//sprintf(Title,"pimissen_q_%d_%.3f_%.3f",j,i*0.05,(i*0.05+0.05)); //rebin 10
    // 	sprintf(Title2,"DeltaE_{CM} 2 #gamma for q_{bin}=%i and E_{bin}=%i; DeltaE_{CM}(MeV)",i,j);
    // 	DeltaE_Missmom_BeamE[i*bin_e+j] =  new TH1(Title,Title2,100,-60,60);
    //   }
    // }

    int bin_q=300;
    int bin_e=23;
    char Title[60];
    char Title2[160];
    for (int j=0; j<bin_e; j++) {
      
      // sprintf(Title,"pimissen_pr_DE_E_%d",j);
      // sprintf(Title2,"DeltaE_{CM} 2 #gamma for E_{bin}=%i; DeltaE_{CM}(MeV); q(fm^{-1})",j);
      // DeltaE_Missmom_BeamE[j] = new GH2(Title,Title2,100,-60.0,60.0,300,0.0,1.5,352);



      sprintf(Title,"pimissen_pr_DE_E_%d",j);
      sprintf(Title2,"Prompt DeltaE_{CM} 2 #gamma for E_{bin}=%i; DeltaE_{CM}(MeV); q(fm^{-1})",j);
      DeltaE_Missmom_BeamE_Prompt[j] = new TH2F(Title,Title2,100,-60,60,300,0.0,1.5);

      sprintf(Title,"pimissen_r_DE_E_%d",j);
      sprintf(Title2,"random DeltaE_{CM} 2 #gamma for E_{bin}=%i; DeltaE_{CM}(MeV); q(fm^{-1})",j);
      DeltaE_Missmom_BeamE_Random[j] = new TH2F(Title,Title2,100,-60,60,300,0.0,1.5);


      // for (int i=0; i<bin_q; i++) {
      // 	sprintf(Title,"pimissen_%i,pr_DE_E_%d",i,j);
      // 	sprintf(Title2,"Prompt DeltaE_{CM} 2 #gamma for E_{bin}=%i e mom_{bin}=%i; DeltaE_{CM}(MeV)",j,i);
      // 	Test_histo[bin_e] = new TH2F(Title,Title2,100,-60,60,300,0.0,1.5);



      // }
    }

int bin_theta=180;
int bin_en=23;
    char Title3[60];
    char Title4[160];
    for (int j=0; j<bin_en; j++) {

       sprintf(Title3,"pimissen_pr_DE_theta_E_%d",j);
      sprintf(Title4,"Prompt DeltaE_{CM} 2 #gamma for E_{bin}=%i; DeltaE_{CM}(MeV); theta(deg)",j);
      DeltaE_Thetapi0_BeamE_Prompt[j] = new TH2F(Title3,Title4,100,-60,60,180,0,180);

      sprintf(Title3,"pimissen_r_DE_theta_E_%d",j);
      sprintf(Title4,"random DeltaE_{CM} 2 #gamma for E_{bin}=%i; DeltaE_{CM}(MeV); theta(deg)",j);
      DeltaE_Thetapi0_BeamE_Random[j] = new TH2F(Title3,Title4,100,-60,60,180,0,180);

}

}

PTinPi0Physics::~PTinPi0Physics()
{
}

Bool_t	PTinPi0Physics::Init()
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;

	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PTinPi0Physics::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();

    return kTRUE;
}

void	PTinPi0Physics::ProcessEvent()
{
	// fill time diff (tagger - pi0), all pi0
    FillTime(*GetNeutralPions(),time);
    FillTimeCut(*GetNeutralPions(),time_cut);
	
	// fill missing mass, all pi0
    FillMissingMass(*GetNeutralPions(),MM);
	
	// fill invariant mass, all pi0
    FillMass(*GetNeutralPions(),IM);

    //my added histograms

	// fill missing momentum, all pi0
    FillMissingMomentum(*GetNeutralPions(),MMom);
	// fill missing momentum calculated using Dan routine, all pi0
    FillMissingMomentumD(*GetNeutralPions(),MMomD);
	//
	
    int i=0;
    // Some neutral decays
    // for (Int_t i = 0; i < GetNeutralPions()->GetNParticles(); i++)
    // {
        // Fill MM for 2 photon decay
      if ( (GetNeutralPions()->GetNParticles()==1) && (GetNeutralPions()->GetNSubParticles(i) == 2) && (GetNeutralPions()->GetNSubPhotons(i) == 2))
        {
		// fill time diff (tagger - pi0), this pi0
        FillTime(*GetNeutralPions(),i,time_2g);
        FillTimeCut(*GetNeutralPions(),i,time_2g_cut);

	int n_track_h = GetTracks()->GetNTracks();
	for (int j=0; j<  n_track_h ; j++) {
	  TLorentzVector v4_h =  GetTracks()->GetVector(j);
	  //cout << v4_h.Px() << endl;
	}
			
		// fill missing mass, this pi0
	FillMissingMass(*GetNeutralPions(),i,MM_2g);
            
		// fill invariant mass, this pi0
            FillMass(*GetNeutralPions(),i,IM_2g);

      //my added histograms

	    FillDeltaE_Missmom(*GetNeutralPions());
	    FillDeltaE_Thetapi0(*GetNeutralPions());
	    FillDeltaE(*GetNeutralPions(),i,DeltaE_CM_D_2g);
	    
		// fill missing momentum, this pi0
	    FillMissingMomentum(*GetNeutralPions(),i,MMom_2g);
		// fill missing momentum calculated using Dan routine, this pi0
	    FillMissingMomentumD(*GetNeutralPions(),i,MMomD_2g);
		//

        }

      //    }

}

void	PTinPi0Physics::ProcessScalerRead()
{
	// Fill Tagger Scalers
	FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	PTinPi0Physics::Write()
{
    // Write all GH1's and TObjects defined in this class
	GTreeManager::Write();
}


void PTinPi0Physics::FillDeltaE_Missmom(const GTreeMeson& tree)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        	FillDeltaE_Missmom(tree, i, j);
	}
    }
}

void PTinPi0Physics::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        	FillDeltaE_Missmom(tree, particle_index, i);
	}
}

void PTinPi0Physics::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index)
{ 

  //  if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
  // if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;
    
  Double_t time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
  TLorentzVector beam2(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));  

  double qmin=0.0; //min and max in fm^-1
  double qmax=1.5;
  int nqbin = 300;
  int qbin;
  int nEbin=23;
  int Ebin=-1;
  double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};
    
    // calc particle time diff
  //    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    Double_t q=CalcMissingMomentumD(tree,particle_index,tagger_index);
    Double_t E_beam=beam2.E();
    Double_t DeltaE=CalcDeltaED(tree,particle_index,tagger_index);

    for (int i=0; i<nEbin; i++) {
      if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
    }
    qbin = int(TMath::Floor(q*50));
    if ( Ebin != -1 && qbin >0 && qbin< nqbin ) {
      DeltaE_CM_D_2g->Fill(DeltaE,time);
      //    cout << q << " " << DeltaE << " " << time << " " << Ebin << endl;
      //     DeltaE_Missmom_BeamE_Prompt[Ebin]->Fill(DeltaE,q);
      // gHist[qbin*nEbin+Ebin]->Fill(DeltaE,time);
      //      cout << "DeltaE=" << DeltaE << "  q=" << q << "  time=" << time << " IsPrompt=" << GHistBGSub::IsPrompt(time) << " IsRandom=" << GHistBGSub::IsRandom(time) << endl; 

      //      DeltaE_Missmom_BeamE[Ebin]->Fill(DeltaE,q,time);
      if (GHistBGSub::IsPrompt(time)) {
	DeltaE_Missmom_BeamE_Prompt[Ebin]->Fill(DeltaE,q);
	Test_histo_h1_prompt->Fill(DeltaE);
      }
      if (GHistBGSub::IsRandom(time)) {
	DeltaE_Missmom_BeamE_Random[Ebin]->Fill(DeltaE,q);
	Test_histo_h1_random->Fill(DeltaE);
      }
    }


						
}

void PTinPi0Physics::FillDeltaE_Thetapi0(const GTreeMeson& tree)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  	FillDeltaE_Thetapi0(tree, i, j);
	}
    }
}


void PTinPi0Physics::FillDeltaE_Thetapi0(const GTreeMeson& tree, Int_t particle_index)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
	  	FillDeltaE_Thetapi0(tree, particle_index, i);
	}
}


void PTinPi0Physics::FillDeltaE_Thetapi0(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index)
{ 

  Double_t time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
  TLorentzVector beam2(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));  

  double thetamin=0.0; //min and max in fm^-1
  double thetamax=180;
  int nthetabin = 180;
  int thetabin;
  int nEbin=23;
  int Ebin=-1;
  double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};
  TLorentzVector particle	= tree.Particle(particle_index); 
  Double_t theta=particle.Theta()*TMath::RadToDeg();

  Double_t E_beam=beam2.E();
  Double_t DeltaE=CalcDeltaED(tree,particle_index,tagger_index);

    for (int i=0; i<nEbin; i++) {
      if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
    }

    if ( Ebin != -1) { 

      if (GHistBGSub::IsPrompt(time)) {
   	DeltaE_Thetapi0_BeamE_Prompt[Ebin]->Fill(DeltaE,theta);
    	//Test_histo_h1_prompt->Fill(DeltaE);
      }

      if (GHistBGSub::IsRandom(time)) {
    	DeltaE_Thetapi0_BeamE_Random[Ebin]->Fill(DeltaE,theta);
    	//Test_histo_h1_random->Fill(DeltaE);
      }
    }
}

void PTinPi0Physics::FillDeltaE_Thetapi0_corr1(const GTreeMeson& tree)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  	FillDeltaE_Thetapi0(tree, i, j);
	}
    }
}


void PTinPi0Physics::FillDeltaE_Thetapi0_corr1(const GTreeMeson& tree, Int_t particle_index)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
	  	FillDeltaE_Thetapi0(tree, particle_index, i);
	}
}


void PTinPi0Physics::FillDeltaE_Thetapi0_corr1(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index)
{ 

  Double_t time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
  TLorentzVector beam2(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));  

  double thetamin=0.0; //min and max in fm^-1
  double thetamax=180;
  int nthetabin = 180;
  int thetabin;
  int nEbin=23;
  int Ebin=-1;
  double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};
  TLorentzVector particle	= tree.Particle(particle_index); 
  Double_t theta=particle.Theta()*TMath::RadToDeg();
  
  TVector3 v_shift(0.0,0.0,0.44); 
  TVector3 v_part,v_part_new;
  v_part.SetMagThetaPhi(45.411,particle.Theta(),particle.Phi());
  v_part_new = v_part - v_shift;
  theta = v_part_new.Theta()*TMath::RadToDeg();

  Double_t E_beam=beam2.E();
  Double_t DeltaE=CalcDeltaED(tree,particle_index,tagger_index);

    for (int i=0; i<nEbin; i++) {
      if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
    }

    if ( Ebin != -1) { 

      if (GHistBGSub::IsPrompt(time)) {
   	DeltaE_Thetapi0_BeamE_Prompt[Ebin]->Fill(DeltaE,theta);
    	//Test_histo_h1_prompt->Fill(DeltaE);
      }

      if (GHistBGSub::IsRandom(time)) {
    	DeltaE_Thetapi0_BeamE_Random[Ebin]->Fill(DeltaE,theta);
    	//Test_histo_h1_random->Fill(DeltaE);
      }
    }
}

void PTinPi0Physics::FillDeltaE_Thetapi0_corr2(const GTreeMeson& tree)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
	  	FillDeltaE_Thetapi0(tree, i, j);
	}
    }
}


void PTinPi0Physics::FillDeltaE_Thetapi0_corr2(const GTreeMeson& tree, Int_t particle_index)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
	  	FillDeltaE_Thetapi0(tree, particle_index, i);
	}
}


void PTinPi0Physics::FillDeltaE_Thetapi0_corr2(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index)
{ 

  Double_t time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
  TLorentzVector beam2(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));  

  double thetamin=0.0; //min and max in fm^-1
  double thetamax=180;
  int nthetabin = 180;
  int thetabin;
  int nEbin=23;
  int Ebin=-1;
  double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};
  TLorentzVector particle	= tree.Particle(particle_index); 
  Double_t theta=particle.Theta()*TMath::RadToDeg();

  TVector3 v_shift(0.0,0.0,0.44); 
  TVector3 v_part,v_part_new;
  v_part.SetMagThetaPhi(45.411,particle.Theta(),particle.Phi());
  v_part_new = v_part - v_shift;
  theta = v_part_new.Theta()*TMath::RadToDeg();

  Double_t E_beam=beam2.E();
  Double_t DeltaE=CalcDeltaED_corr(tree,particle_index,tagger_index);

    for (int i=0; i<nEbin; i++) {
      if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
    }

    if ( Ebin != -1) { 

      if (GHistBGSub::IsPrompt(time)) {
   	DeltaE_Thetapi0_BeamE_Prompt[Ebin]->Fill(DeltaE,theta);
    	//Test_histo_h1_prompt->Fill(DeltaE);
      }

      if (GHistBGSub::IsRandom(time)) {
    	DeltaE_Thetapi0_BeamE_Random[Ebin]->Fill(DeltaE,theta);
    	//Test_histo_h1_random->Fill(DeltaE);
      }
    }
}



//

