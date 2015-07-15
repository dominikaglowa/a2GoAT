#ifndef __CINT__

#include "PPhysics.h"
#include "PPhysicsTin.h"
PPhysicsTin::PPhysicsTin() 
{
	TC_cut_min = 200;
	TC_cut_max = 300;
}

PPhysicsTin::~PPhysicsTin()
{
}

Bool_t	PPhysicsTin::Init()
{
	return kTRUE;
}

void	PPhysicsTin::Reconstruct()
{
}

// ----------------------------------------------------------------------------------------
// TH1 routines
// ----------------------------------------------------------------------------------------
void PPhysicsTin::FillScalers(Int_t low_scaler_number, Int_t high_scaler_number, TH1* hist)
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

    // Loop over scaler range, don't pull anything higher than the real # GetScalers()
	if (low_scaler_number  < 0)
	{
		cout << "FillScalers given scaler number outside range: " << low_scaler_number << endl;
		cout << "Setting lower limit to zero and continuing" << endl;
		low_scaler_number = 0;
	}
    if (high_scaler_number > GetScalers()->GetNScalers())
	{
		cout << "FillScalers given scaler number outside range: " << high_scaler_number << endl;
		cout << "Setting upper limit to "<< high_scaler_number << " and continuing" << endl;	
        high_scaler_number = GetScalers()->GetNScalers();
	}

    for (Int_t i = low_scaler_number; i <= high_scaler_number; i++)
	{
		Int_t bin = i - low_scaler_number;
        hist_current_SR->SetBinContent(bin,GetScalers()->GetScaler(i));
	}

	// Add to accumulated
	hist->Add(hist_current_SR);
}	

void PPhysicsTin::FillMissingMass(const GTreeParticle& tree, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        	FillMissingMass(tree, i, j, Hprompt, Hrandom);
	}
    }
}

void PPhysicsTin::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        	FillMissingMass(tree, particle_index, i, Hprompt, Hrandom);
	}
}

void PPhysicsTin::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom)
{
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

	if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(missingp4.M());
	if (GHistBGSub::IsRandom(time)) Hrandom->Fill(missingp4.M());						
}

void PPhysicsTin::FillTime(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			Hist->Fill(time);
		}
	}
}

void PPhysicsTin::FillTime(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
    // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;
	
        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
	Hist->Fill(time);
	}
}

void PPhysicsTin::FillTimeCut(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                    time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
		}
	}
}

void PPhysicsTin::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
	}
}

void PPhysicsTin::FillMass(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        Hist->Fill(tree.GetMass(i));
	}
}

void PPhysicsTin::FillMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    Hist->Fill(tree.GetMass(particle_index));
}

void PPhysicsTin::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, Hprompt, Hrandom, MM_min, MM_max);
    }

}

void PPhysicsTin::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
  //    cout << GetTagger()->GetTaggedChannel(tagger_index) << " " << TC_cut_min << " " << TC_cut_max << endl;
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
//    cout << "time " << time << endl; 
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
 //   cout << "MM " << missingp4.M() << endl;     
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;
   
   if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(tree.GetPhi(particle_index)); //cout << "prompt" << endl;}
   if (GHistBGSub::IsRandom(time)) Hrandom->Fill(tree.GetPhi(particle_index));	//cout << "random" << endl;}
}

Double_t PPhysicsTin::CalcCoplanarity(const GTreeParticle& tree1, Int_t particle_index1, const GTreeParticle& tree2, Int_t particle_index2)
{
   Double_t phi1 = tree1.GetPhi(particle_index1);
   Double_t phi2 = tree2.GetPhi(particle_index2);
   Double_t phidiff = TMath::Abs(phi1 - phi2);

   return phidiff;
}

//my TH1 histograms array
// void PPhysicsTin::FillDeltaE_Missmom(const GTreeMeson& tree, TH2F** Hprompt, TH2F** Hrandom)
// {
//     for (Int_t i = 0; i < tree.GetNParticles(); i++)
//     {
//     for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
// 	{
//         	FillDeltaE_Missmom(tree, i, j, Hprompt, Hrandom);
// 	}
//     }
// }

// void PPhysicsTin::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, TH2F** Hprompt, TH2F** Hrandom)
// {
//     for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
// 	{
//         	FillDeltaE_Missmom(tree, particle_index, i, Hprompt, Hrandom);
// 	}
// }

// void PPhysicsTin::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, TH2F** Hprompt, TH2F** Hrandom)
// { 

//   if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
//   if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;
    
//   time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
//   beam = TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));  

// {
//   double qmin=0.0; //min and max in fm^-1
//   double qmax=1.5;
//   int nqbin = 300;
//   int qbin;
//   int nEbin=23;
//   int Ebin=-1;
//   double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};
    
//     // calc particle time diff
//     time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
//     Double_t q=CalcMissingMomentumD(tree,particle_index,tagger_index);
//     Double_t E_beam=beam.E();
//     Double_t DeltaE=CalcDeltaED(tree,particle_index,tagger_index);

//     for (int i=0; i<nEbin; i++) {
//       if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
//     }
//     qbin = int(TMath::Floor(q*50));
//     if ( Ebin != -1 && qbin >0 && qbin< nqbin ) {
//       // gHist[qbin*nEbin+Ebin]->Fill(DeltaE,time);
//       if (GHistBGSub::IsPrompt(time)) Hprompt[Ebin]->Fill(DeltaE,q);
//       //      if (GHistBGSub::IsRandom(time)) Hrandom[Ebin]->Fill(DeltaE,q);
//     }
// }

						
// }
// //

// ----------------------------------------------------------------------------------------
// GH1 routines
// ----------------------------------------------------------------------------------------

void PPhysicsTin::FillMissingMass(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	  {//cout<<tree.GetNParticles()<<endl;

        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMass(tree, i, j, gHist, TaggerBinning);
			//cout<<i<<endl;
		}
	}
	  GetTracks()->GetNTracks();
	  //cout<<GetTracks()->GetNTracks()<<endl;
}

void PPhysicsTin::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMass(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysicsTin::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
  //    cout << GetTagger()->GetTaggedChannel(tagger_index) << " " << TC_cut_min << " " << TC_cut_max << endl;
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
    //    cout << "Missing momentum=" << missingp4.M() << endl;

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(missingp4.M(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(missingp4.M(),time);

}

//my added histograms
void PPhysicsTin::FillMissingMomentum(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMomentum(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysicsTin::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
  for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMomentum(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysicsTin::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
	if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;
	
    // calc particle time diff
	time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);

    // calc missing p4
	missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

    // Fill GH1
	if(TaggerBinning)   gHist->Fill(missingp4.Rho()/197.3,time, GetTagger()->GetTaggedChannel(tagger_index));
	else gHist->Fill(missingp4.Rho()/197.3,time);					

}

Double_t PPhysicsTin::CalcMissingMomentum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.Rho()/197.3;
}

void PPhysicsTin::FillMissingMomentumD(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMomentum(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysicsTin::FillMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
  for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMomentum(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysicsTin::FillMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
	if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;
	
    // calc particle time diff
	time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);

	particle	= tree.Particle(particle_index);
	beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
  
  
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
	if(TaggerBinning)   gHist->Fill(q,time, GetTagger()->GetTaggedChannel(tagger_index));
	else gHist->Fill(q,time);
				

}

void PPhysicsTin::FillDeltaE(const GTreeMeson& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillDeltaE(tree, i, j, gHist,TaggerBinning);
		}
	}
}

void PPhysicsTin::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillDeltaE(tree, particle_index, i, gHist,TaggerBinning);
	}
}

void PPhysicsTin::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{

  // Is tagger channel rejected by user?
  if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
  if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;
    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // Fill GH1
    gHist->Fill(CalcDeltaED(tree,particle_index,tagger_index),time);					

}

// void PPhysicsTin::FillDeltaE_Missmom(const GTreeMeson& tree, GH2** gHist, Bool_t TaggerBinning)
// {
// 	for (Int_t i = 0; i < tree.GetNParticles(); i++)
// 	{
// 		for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
// 		{
// 		  FillDeltaE_Missmom(tree, i, j, gHist,TaggerBinning);
// 		}
// 	}
// }

// void PPhysicsTin::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, GH2** gHist, Bool_t TaggerBinning)
// {
//     for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
// 	{
// 	  FillDeltaE_Missmom(tree, particle_index, i, gHist, TaggerBinning);
// 	}
// }

// void PPhysicsTin::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, GH2** gHist, Bool_t TaggerBinning)
// {

//   if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
//   if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;
  
//   double qmin=0.0; //min and max in fm^-1
//   double qmax=1.5;
//   int nqbin = 300;
//   int qbin;
//   int nEbin=23;
//   int Ebin=-1;
//   double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};

//   beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
    
//     // calc particle time diff
//     time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
//     Double_t q=CalcMissingMomentumD(tree,particle_index,tagger_index);
//     Double_t E_beam=beam.E()               ;
//     Double_t DeltaE=CalcDeltaED(tree,particle_index,tagger_index);

//     for (int i=0; i<nEbin; i++) {
//       if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
//     }
//     qbin = int(TMath::Floor(q*50));
//     if ( Ebin != -1 && qbin >0 && qbin< nqbin ) {
//       gHist[Ebin]->Fill(DeltaE,q,time);
//     }
// }

Double_t PPhysicsTin::CalcDeltaED(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index) {
  
  beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
 
 
  int n_sub_part;
  
  Double_t  E1, E2;
  Double_t E_diff=10000; // Value out of range in case one does not have just 2 photons

  // number of meson is ever 1 on the tree
  n_sub_part = tree.GetNSubParticles(0);
  //  Int_t n_sub_phot_ind[n_sub_part];
  std::vector<Int_t> n_sub_phot_ind=tree.GetTrackIndexList(0);
  
  // TLorentzVector gamma[n_sub_part];
  //int i=0;
  //int j;
  //if ( n_sub_part ==2) {
  //for(std::vector<Int_t>::iterator j=n_sub_phot_ind.begin();j!=n_sub_phot_ind.end();++j) {
  //gamma[i] =   GetTracks()->GetVector(*j);
  //i++;
  // }
    // std::for_each( n_sub_phot_ind.begin(), n_sub_phot_ind.end(),[&](Int_t n){
    // 	gamma[i] =   GetTracks()->GetVector(n);
    // 	i++;
    //   });
  TLorentzVector gamma1 = GetTracks()->GetVector(n_sub_phot_ind[0]);
  TLorentzVector gamma2 = GetTracks()->GetVector(n_sub_phot_ind[1]);


    E1 = gamma1.E(); // photon 1 energy
    E2 = gamma2.E(); // photon 2 energy
    
    //cout << E1 << " " << gamma1.Px() << " " << E2 << " " << gamma2.Px() << endl;

    Double_t beta = (beam.E()/(beam.E() + target.M())); 
    Double_t lorentz_gamma =  1/(sqrt(1 - beta*beta));
    Double_t Xform = (E1 - E2)/(E1 + E2);
    Double_t psi = gamma1.Vect().Angle(gamma2.Vect());
    Double_t costheta1 = gamma1.CosTheta();
    Double_t costheta2 = gamma2.CosTheta();
    Double_t theta1 = gamma1.Theta()/TMath::Pi()*180;
    Double_t theta2 = gamma2.Theta()/TMath::Pi()*180;
    Double_t mpi0 = 134.9766;
    Double_t M = target.M();
    Double_t Egamma=beam.E();
    //cout <<"theta2"<< theta2 << "   theta1" << theta1 << endl;
    if (theta1 > 25 && theta1<155 && theta2>25 &&theta2<155)	    E_diff = lorentz_gamma*((sqrt(2*mpi0*mpi0/((1-Xform*Xform)*(1-cos(psi))))) 
			    - ( beta*(E1*costheta1 + E2*costheta2)) ) 
      - ( (2*Egamma*M + mpi0*mpi0)/(2*sqrt(2*Egamma*M + M*M)) );
      else E_diff = 80.; // return a value out of our range for histrograms if the value is not in CrystalBall, also assuming 5 degrees for full acceptance of the photon
  return E_diff;
  
}


Double_t PPhysicsTin::CalcDeltaED_corr(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index) {
  
  beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
 
 
  int n_sub_part;
  
  Double_t  E1, E2;
  Double_t E_diff=10000; // Value out of range in case one does not have just 2 photons

  // number of meson is ever 1 on the tree
  n_sub_part = tree.GetNSubParticles(0);
  //  Int_t n_sub_phot_ind[n_sub_part];
  std::vector<Int_t> n_sub_phot_ind=tree.GetTrackIndexList(0);
  
  // TLorentzVector gamma[n_sub_part];
  //int i=0;
  //int j;
  //if ( n_sub_part ==2) {
  //for(std::vector<Int_t>::iterator j=n_sub_phot_ind.begin();j!=n_sub_phot_ind.end();++j) {
  //gamma[i] =   GetTracks()->GetVector(*j);
  //i++;
  // }
    // std::for_each( n_sub_phot_ind.begin(), n_sub_phot_ind.end(),[&](Int_t n){
    // 	gamma[i] =   GetTracks()->GetVector(n);
    // 	i++;
    //   });
  TLorentzVector gamma1 = GetTracks()->GetVector(n_sub_phot_ind[0]);
  TLorentzVector gamma2 = GetTracks()->GetVector(n_sub_phot_ind[1]);


    E1 = gamma1.E(); // photon 1 energy
    E2 = gamma2.E(); // photon 2 energy
    
    TVector3 v3_gamma1 = gamma1.Vect();
    TVector3 v3_gamma2 = gamma2.Vect();
    TVector3 v_shift(0.,0.,0.44);
    v3_gamma1 = v3_gamma1 - v_shift;
    v3_gamma2 = v3_gamma2 - v_shift;
    
    gamma1.SetVect(v3_gamma1);
    gamma2.SetVect(v3_gamma2);
    // cout << E1 << " " << gamma1.Px() << " " << E2 << " " << gamma2.Px() << endl;

    Double_t beta = (beam.E()/(beam.E() + target.M())); 
    Double_t lorentz_gamma =  1/(sqrt(1 - beta*beta));
    Double_t Xform = (E1 - E2)/(E1 + E2);
    Double_t psi = gamma1.Vect().Angle(gamma2.Vect());
    Double_t costheta1 = gamma1.CosTheta();
    Double_t costheta2 = gamma2.CosTheta();
    Double_t theta1 = gamma1.Theta()/TMath::Pi()*180;
    Double_t theta2 = gamma2.Theta()/TMath::Pi()*180;
    Double_t mpi0 = 134.9766;
    Double_t M = target.M();
    Double_t Egamma=beam.E();
    if (theta1 > 25 && theta1<155 && theta2>25 &&theta2<155)	    E_diff = lorentz_gamma*((sqrt(2*mpi0*mpi0/((1-Xform*Xform)*(1-cos(psi))))) 
			    - ( beta*(E1*costheta1 + E2*costheta2)) ) 
      - ( (2*Egamma*M + mpi0*mpi0)/(2*sqrt(2*Egamma*M + M*M)) );
      else E_diff = 80.; // return a value out of our range for histrograms if the value is not in CrystalBall, also assuming 5 degrees for full acceptance of the photon
  return E_diff;
  
}

Double_t PPhysicsTin::CalcMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.Rho()/197.3;
}

//
Double_t PPhysicsTin::CalcMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.M();
}

Double_t PPhysicsTin::CalcMissingEnergy(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree,particle_index, tagger_index);

	return missingp4.T();
}

TLorentzVector PPhysicsTin::CalcMissingP4(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
	missingp4 	= beam + target - particle;						

	return missingp4;
}

void PPhysicsTin::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, gHist, TaggerBinning);
    }

}

void PPhysicsTin::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;
   
   if(TaggerBinning)   gHist->Fill(tree.GetPhi(particle_index),time,GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(tree.GetPhi(particle_index),time);

}

void PPhysicsTin::FillTime(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			gHist->Fill(time);
		}
	}
}

void PPhysicsTin::FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;
	
        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		gHist->Fill(time);
	}
}

void PPhysicsTin::FillTimeCut(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            // Is tagger channel rejected by user?
            if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
            if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

                    time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
		}
	}
}

void PPhysicsTin::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        // Is tagger channel rejected by user?
        if(GetTagger()->GetTaggedChannel(j) < TC_cut_min) continue;
        if(GetTagger()->GetTaggedChannel(j) > TC_cut_max) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
	}
}

void PPhysicsTin::FillMass(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        gHist->Fill(tree.GetMass(i));
	}
}

void PPhysicsTin::FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    gHist->Fill(tree.GetMass(particle_index));
}

Bool_t 	PPhysicsTin::Write()
{
	return kTRUE;
}

// Some common initialisation stuff
Bool_t 	PPhysicsTin::InitBackgroundCuts()
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

Bool_t 	PPhysicsTin::InitTargetMass()
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

Bool_t 	PPhysicsTin::InitTaggerChannelCuts()
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

Bool_t 	PPhysicsTin::InitTaggerScalers()
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
        cout << "Tagger Channel GetScalers() not set correctly" << endl;
		return kFALSE;
	}

	cout << endl;
	return kTRUE;

}
#endif
