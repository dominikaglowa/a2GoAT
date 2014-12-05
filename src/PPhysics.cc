#ifndef __CINT__

#include "PPhysics.h"

PPhysics::PPhysics() 
{
	TC_cut_min = 0;
	TC_cut_max = 352;
}

PPhysics::~PPhysics()
{
}

Bool_t	PPhysics::Init(const Char_t *configfile)
{
	return kTRUE;
}

void	PPhysics::Reconstruct()
{
}

// ----------------------------------------------------------------------------------------
// TH1 routines
// ----------------------------------------------------------------------------------------
void PPhysics::FillScalers(Int_t low_scaler_number, Int_t high_scaler_number, TH1* hist)
{
	Int_t nFillScalers = high_scaler_number - low_scaler_number + 1;

	if( nFillScalers < hist->GetNbinsX())
	{
	    cout << "Error: FillScalers - histogram has insufficient bins for range" << endl;
	    return;
	}
	
	// To properly accumulate, create a histogram for this scaler read
	// cloning input histogram means the axis will be equivalent
	TH1* hist_current_SR = (TH1D*) hist->Clone();
	hist_current_SR->Reset();

	// Loop over scaler range, don't pull anything higher than the real # scalers
	if (low_scaler_number  < 0)
	{
		cout << "FillScalers given scaler number outside range: " << low_scaler_number << endl;
		cout << "Setting lower limit to zero and continuing" << endl;
		low_scaler_number = 0;
	}
	if (high_scaler_number > scalers->GetNScaler())
	{
		cout << "FillScalers given scaler number outside range: " << high_scaler_number << endl;
		cout << "Setting upper limit to "<< high_scaler_number << " and continuing" << endl;	
		high_scaler_number = scalers->GetNScaler();
	}

	for (int i = low_scaler_number; i <= high_scaler_number; i++) 
	{
		Int_t bin = i - low_scaler_number;
		hist_current_SR->SetBinContent(bin,scalers->GetScaler(i));
	}

	// Add to accumulated
	hist->Add(hist_current_SR);
}	

void PPhysics::FillMissingMass(const GTreeParticle& tree, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
	for (Int_t j = 0; j < tagger->GetNTagged(); j++)
	{
        	FillMissingMass(tree, i, j, Hprompt, Hrandom);
	}
    }
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        	FillMissingMass(tree, particle_index, i, Hprompt, Hrandom);
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom)
{
	time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
	missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

	if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(missingp4.M());
	if (GHistBGSub::IsRandom(time)) Hrandom->Fill(missingp4.M());						
}

void PPhysics::FillTime(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			// Is tagger channel rejected by user?
			if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
			if(tagger->GetTagged_ch(j) > TC_cut_max) continue;

           		time = tagger->GetTagged_t(j) - tree.GetTime(i);
			Hist->Fill(time);
		}
	}
}

void PPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
	for (Int_t j = 0; j < tagger->GetNTagged(); j++)
	{
	// Is tagger channel rejected by user?
		if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
		if(tagger->GetTagged_ch(j) > TC_cut_max) continue;
	
		time = tagger->GetTagged_t(j) - tree.GetTime(particle_index);
	Hist->Fill(time);
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			// Is tagger channel rejected by user?
			if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
			if(tagger->GetTagged_ch(j) > TC_cut_max) continue;

            		time = tagger->GetTagged_t(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
		}
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
	for (Int_t j = 0; j < tagger->GetNTagged(); j++)
	{
		// Is tagger channel rejected by user?
		if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
		if(tagger->GetTagged_ch(j) > TC_cut_max) continue;

		time = tagger->GetTagged_t(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		Hist->Fill(tree.Particle(i).M());
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
	Hist->Fill(tree.Particle(particle_index).M());
}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, Hprompt, Hrandom, MM_min, MM_max);
    }

}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
//    cout << tagger->GetTagged_ch(tagger_index) << " " << TC_cut_min << " " << TC_cut_max << endl;
    if(tagger->GetTagged_ch(tagger_index) < TC_cut_min) return;
    if(tagger->GetTagged_ch(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
//    cout << "time " << time << endl; 
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
 //   cout << "MM " << missingp4.M() << endl;     
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;

   // Calc phi
   double phi = tree.Particle(particle_index).Phi() * TMath::RadToDeg();
//    cout << "phi " << phi << endl;     
   
   if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(phi); //cout << "prompt" << endl;}
   if (GHistBGSub::IsRandom(time)) Hrandom->Fill(phi);	//cout << "random" << endl;}
}

Double_t PPhysics::CalcCoplanarity(const GTreeParticle& tree1, Int_t particle_index1, const GTreeParticle& tree2, Int_t particle_index2)
{
   double phi1 = tree1.Particle(particle_index1).Phi() * TMath::RadToDeg();
   double phi2 = tree2.Particle(particle_index2).Phi() * TMath::RadToDeg();
   double phidiff = TMath::Abs(phi1 - phi2);

   return phidiff;
}

// ----------------------------------------------------------------------------------------
// GH1 routines
// ----------------------------------------------------------------------------------------

void PPhysics::FillMissingMass(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			FillMissingMass(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillMissingMass(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
    if(tagger->GetTagged_ch(tagger_index) < TC_cut_min) return;
    if(tagger->GetTagged_ch(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(missingp4.M(),time, tagger->GetTagged_ch(tagger_index));
   else gHist->Fill(missingp4.M(),time);

}
//
//stuff from the old file
//
void PPhysics::FillMissingMomentum(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			FillMissingMomentum(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillMissingMomentum(tree, particle_index, i, gHist, TaggerBinning);
	}
}
void PPhysics::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

	// Fill GH1
	gHist->Fill(missingp4.Rho()/197.3,time);					

}

void PPhysics::FillMissingMomentumD(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
		  FillMissingMomentumD(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillMissingMomentumD(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
  
  
    double theta_pi0 = particle.Theta(); // theta pi0
    double mpi0 = 134.9766;
    double costheta = TMath::Cos(theta_pi0);
    double beta = beam.E()/(beam.E()+target.M());
    double Egamma_c = beam.E()*TMath::Sqrt( (1-beta)/(1+beta) );
    double gamma_l = 1./(TMath::Sqrt(1-beta*beta));
    double minv = 2*beam.E()*target.M() + target.M2();
    double Epi_c = TMath::Sqrt(minv)/2. + (mpi0*mpi0 - target.M2())/(2.*TMath::Sqrt(minv));
    double Epi = (Epi_c + TMath::Sqrt( Epi_c*Epi_c - ( 1 - beta*beta*costheta*costheta )*( gamma_l*gamma_l*beta*beta*mpi0*mpi0*costheta*costheta + Epi_c*Epi_c ) ) ) / ( gamma_l*(1-beta*beta*costheta*costheta) );

    double qsq = (Egamma_c - Epi_c)*(Egamma_c - Epi_c) + 2.*beam.E()*(Epi - TMath::Sqrt(Epi*Epi - mpi0*mpi0)*costheta) - mpi0*mpi0;
    Double_t q = TMath::Sqrt(qsq)/197.3;

	// Fill GH1
	gHist->Fill(q,time);					
}

void PPhysics::FillDeltaE(const GTreeMeson& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			FillDeltaE(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillDeltaE(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    

    // Fill GH1
    gHist->Fill(CalcDeltaED(tree,particle_index,tagger_index),time);					

}

void PPhysics::FillDeltaE_Missmom(const GTreeMeson& tree, GH1** gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
		  FillDeltaE_Missmom(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, GH1** gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
	  FillDeltaE_Missmom(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, GH1** gHist, Bool_t TaggerBinning)
{
  double qmin=0.0; //min and max in fm^-1
  double qmax=1.5;
  int nqbin = 300;
  int qbin;
  int nEbin=23;
  int Ebin=-1;
  double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};
    
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    Double_t q=CalcMissingMomentumD(tree,particle_index,tagger_index);
    Double_t E_beam=CalcBeamE(tagger_index);
    Double_t DeltaE=CalcDeltaED(tree,particle_index,tagger_index);

    for (int i=0; i<nEbin; i++) {
      if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
    }
    qbin = int(TMath::Floor(q*200)); //change this when rebinning in q so qbin=1
    if ( Ebin != -1 && qbin >0 && qbin< nqbin ) {
      gHist[qbin*nEbin+Ebin]->Fill(DeltaE,time);
    }
}


//
//
Double_t PPhysics::CalcMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.M();
}

Double_t PPhysics::CalcMissingEnergy(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree,particle_index, tagger_index);

	return missingp4.T();
}

TLorentzVector PPhysics::CalcMissingP4(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
	missingp4 	= beam + target - particle;						

	return missingp4;
}

//
//from the old file
//
Double_t PPhysics::CalcDeltaED(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index)
{
  particle	= tree.Particle(particle_index);
  beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
 
  int n_sub_part;
  int n_sub_phot;
  Double_t  E1, E2;
  Double_t E_diff=10000; // Value out of range in case one does not have just 2 photons

  // number of meson is ever 1 on the tree
  n_sub_part = tree.GetNSubParticles(0);
  n_sub_phot = tree.GetNSubPhotons(0);
  TLorentzVector gamma[n_sub_phot];
  if ( n_sub_phot ==2) {
    for (int i=0; i< n_sub_phot; i++) {
      gamma[i] = tree.SubPhotons(0,i);
    }

    E1 = gamma[0].E(); // photon 1 energy
    E2 = gamma[1].E(); // photon 2 energy
    Double_t beta = (beam.E()/(beam.E() + target.M())); 
    Double_t lorentz_gamma =  1/(sqrt(1 - beta*beta));
    Double_t Xform = (E1 - E2)/(E1 + E2);
    Double_t psi = gamma[0].Vect().Angle(gamma[1].Vect());
    Double_t costheta1 = gamma[0].CosTheta();
    Double_t costheta2 = gamma[1].CosTheta();
    Double_t mpi0 = 134.9766;
    Double_t M = target.M();
    Double_t Egamma=beam.E();
    E_diff = lorentz_gamma*((sqrt(2*mpi0*mpi0/((1-Xform*Xform)*(1-cos(psi))))) 
			    - ( beta*(E1*costheta1 + E2*costheta2)) ) 
      - ( (2*Egamma*M + mpi0*mpi0)/(2*sqrt(2*Egamma*M + M*M)) );
  }
  return E_diff;
  
}

Double_t PPhysics::CalcBeamE(Int_t tagger_index){
  return tagger->GetPhotonBeam_E(tagger_index);
}

Double_t PPhysics::CalcMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
  particle	= tree.Particle(particle_index);
  beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
  
  
  double theta_pi0 = particle.Theta(); // theta pi0
  double mpi0 = 134.9766;
  double costheta = TMath::Cos(theta_pi0);
  double beta = beam.E()/(beam.E()+target.M());
  double Egamma_c = beam.E()*TMath::Sqrt( (1-beta)/(1+beta) );
  double gamma_l = 1./(TMath::Sqrt(1-beta*beta));
  double minv = 2*beam.E()*target.M() + target.M2();
  double Epi_c = TMath::Sqrt(minv)/2. + (mpi0*mpi0 - target.M2())/(2.*TMath::Sqrt(minv));
  double Epi = (Epi_c + TMath::Sqrt( Epi_c*Epi_c - ( 1 - beta*beta*costheta*costheta )*( gamma_l*gamma_l*beta*beta*mpi0*mpi0*costheta*costheta + Epi_c*Epi_c ) ) ) / ( gamma_l*(1-beta*beta*costheta*costheta) );

  double qsq = (Egamma_c - Epi_c)*(Egamma_c - Epi_c) + 2.*beam.E()*(Epi - TMath::Sqrt(Epi*Epi - mpi0*mpi0)*costheta) - mpi0*mpi0;

  Double_t q = TMath::Sqrt(qsq)/197.3;
    
  return q;
  
}
//

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, gHist, TaggerBinning);
    }

}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
    if(tagger->GetTagged_ch(tagger_index) < TC_cut_min) return;
    if(tagger->GetTagged_ch(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;

   // Calc phi and Fill GH1
   double phi = tree.Particle(particle_index).Phi() * TMath::RadToDeg();
   
   if(TaggerBinning)   gHist->Fill(phi,time,tagger->GetTagged_ch(tagger_index));
   else gHist->Fill(phi,time);

}

void PPhysics::FillTime(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			// Is tagger channel rejected by user?
			if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
			if(tagger->GetTagged_ch(j) > TC_cut_max) continue;

           		time = tagger->GetTagged_t(j) - tree.GetTime(i);
			gHist->Fill(time);
		}
	}
}

void PPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
	for (Int_t j = 0; j < tagger->GetNTagged(); j++)
	{
		// Is tagger channel rejected by user?
		if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
		if(tagger->GetTagged_ch(j) > TC_cut_max) continue;
	
		time = tagger->GetTagged_t(j) - tree.GetTime(particle_index);
		gHist->Fill(time);
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			// Is tagger channel rejected by user?
			if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
			if(tagger->GetTagged_ch(j) > TC_cut_max) continue;

            		time = tagger->GetTagged_t(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
		}
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
	for (Int_t j = 0; j < tagger->GetNTagged(); j++)
	{
		// Is tagger channel rejected by user?
		if(tagger->GetTagged_ch(j) < TC_cut_min) continue;
		if(tagger->GetTagged_ch(j) > TC_cut_max) continue;

		time = tagger->GetTagged_t(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		gHist->Fill(tree.Particle(i).M());
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
	gHist->Fill(tree.Particle(particle_index).M());
}

Bool_t 	PPhysics::Write()
{
	return kTRUE;
}

// Some common initialisation stuff
Bool_t 	PPhysics::InitBackgroundCuts()
{
	// Set background cuts
	Double_t p1, p2, r1, r2;
	string config = ReadConfig("Set-Prompt-Cut");
	if(strcmp(config.c_str(), "nokey") == 0) 
		cout << "No BG subtraction - At least 1 prompt and random cut required" << endl;
	else if(sscanf( config.c_str(), "%lf %lf\n", &p1, &p2) == 2)
	{
	   config = ReadConfig("Add-Random-Cut",0);
	   if(strcmp(config.c_str(), "nokey") == 0) 
	   	cout << "No BG subtraction - At least 1 random cut required" << endl;
	   else if(sscanf( config.c_str(), "%lf %lf\n", &r1, &r2) == 2)
	   {
		cout << "Init BG cuts:" << endl;
		cout << "prompt(" << p1 << "," << p2 << ") ";
		cout << "random(" << r1 << "," << r2 << ") " << endl;

		GHistBGSub::InitCuts(p1,p2,r1,r2);

		// Look for additional random windows
		Int_t instance = 1;
		do
		{
			config = ReadConfig("Add-Random-Cut",instance);
			if(sscanf( config.c_str(), "%lf %lf\n", &r1, &r2) == 2)
			{
				cout << "Adding random cuts: ";
				cout << "random(" << r1 << "," << r2 << ") " << endl;

				GHistBGSub::AddRandCut(r1,r2);
			}		
			instance++;
		} while (strcmp(config.c_str(), "nokey") != 0);
	   }
	   else {cout << "Random window not set correctly" << endl; return kFALSE;}
	}
	else {cout << "Prompt window not set correctly" << endl; return kFALSE;}

	cout << endl;
	return kTRUE;

}

Bool_t 	PPhysics::InitTargetMass()
{
	Double_t mass;
	string config = ReadConfig("Target-Mass");
	if(strcmp(config.c_str(), "nokey") == 0)
	{
		cout << "Target mass unknown!" << endl;
	}
	else if(sscanf( config.c_str(), "%lf\n", &mass) == 1)
	{
		cout << "Setting Target mass: " << mass << " MeV" << endl;
		SetTarget(mass);		
	}
	else 
	{
		cout << "Target Mass not set correctly" << endl; 
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}

Bool_t 	PPhysics::InitTaggerChannelCuts()
{
	Double_t tc1, tc2;
	string config = ReadConfig("Tagger-Channel-Cut");
	if(sscanf( config.c_str(), "%lf %lf\n", &tc1, &tc2) == 2)
	{
		if ((tc1 < 0) || (tc1 > 352))
		{
		   cout << "Invalid tagger channel cut: " << tc1 << endl;
		   return kFALSE;
		}
		else if ((tc2 < 0) || (tc2 > 352))
		{
		   cout << "Invalid tagger channel cut: " << tc2 << endl;
		   return kFALSE;
		}
		
		cout << "Setting cut on tagger channels: " << tc1 << " to " << tc2 << endl;
		SetTC_cut(tc1,tc2);
	}
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
		cout << "Tagger Channel cut not set correctly" << endl; 
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}

Bool_t 	PPhysics::InitTaggerScalers()
{
	Int_t sc1, sc2;
	string config = ReadConfig("Tagger-Scalers");
	if(sscanf( config.c_str(), "%d %d\n", &sc1, &sc2) == 2)
	{
		cout << "Setting Tagger scaler channels: " << sc1 << " to " << sc2 << endl;
		SetTC_scalers(sc1,sc2);
	}
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
		cout << "Tagger Channel scalers not set correctly" << endl; 
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}
#endif
