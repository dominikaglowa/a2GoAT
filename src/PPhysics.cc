#ifndef __CINT__

#include "PPhysics.h"

PPhysics::PPhysics() 
{
	TC_cut_min = 0;
	TC_cut_max = 352;
    TC_scaler_min = 204;
    TC_scaler_max = 555;
    LT_scaler_clock = 0;
    LT_scaler_inhib = 1;
    scalerHists = new TObjArray();
    nScalerHists = 0;
}

PPhysics::~PPhysics()
{
}

Bool_t	PPhysics::Init()
{
    if(!InitDecodeDoubles()) return kFALSE;
    if(!InitTaggerChannelCuts()) return kFALSE;
    if(!InitTaggerScalers()) return kFALSE;
    if(!InitLiveTimeScalers()) return kFALSE;
    if(!InitDisplayScalers()) return kFALSE;
	return kTRUE;
}

void	PPhysics::Reconstruct()
{
}

// ----------------------------------------------------------------------------------------
// TH1 routines
// ----------------------------------------------------------------------------------------
void PPhysics::FillScalers(Int_t low_scaler_number, Int_t high_scaler_number, TH1* hist, Int_t first_bin)
{
    Int_t nFillScalers = high_scaler_number - low_scaler_number + first_bin;

    if( nFillScalers > hist->GetNbinsX())
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
        low_scaler_number = 0;
        cout << "Setting lower limit to " << low_scaler_number << " and continuing" << endl;
	}
    if (high_scaler_number >= GetScalers()->GetNScalers())
	{
		cout << "FillScalers given scaler number outside range: " << high_scaler_number << endl;
        high_scaler_number = ((GetScalers()->GetNScalers()) - 1);
        cout << "Setting upper limit to " << high_scaler_number << " and continuing" << endl;
	}

    for (Int_t i = low_scaler_number; i <= high_scaler_number; i++)
	{
        Int_t bin = (i - low_scaler_number + first_bin);
        hist_current_SR->SetBinContent(bin,GetScalers()->GetScaler(i));
	}

	// Add to accumulated
	hist->Add(hist_current_SR);
}

void PPhysics::AddScalerHist(const char* name, Int_t lo, Int_t hi)
{
    TH1D *h1;

    if(nScalerHists && strcmp(name,scalerHists->At(nScalerHists-1)->GetName())==0)
    {
        cout << "    Appending to: " << name << " with scalers " << lo << " to " << hi << endl;
        h1 = (TH1D*)scalerHists->At(nScalerHists-1);
        h1->SetBins((h1->GetNbinsX()+(hi-lo+1)),0,(h1->GetNbinsX()+(hi-lo+1)));
        nScalerSets.at(nScalerHists-1) += 1;
    }
    else
    {
        cout << "Adding histogram: " << name << " with scalers " << lo << " to " << hi << endl;
        h1 = new TH1D(name,name,(hi-lo+1),0,(hi-lo+1));
        scalerHists->Add(h1);
        nScalerSets.push_back(1);
        nScalerHists++;
    }

    scalerChanL.push_back(lo);
    scalerChanH.push_back(hi);
}

void PPhysics::AddScalerHist(const char* name, Int_t scal, const char* label)
{
    TH1D *h1;

    if(nScalerHists && strcmp(name,scalerHists->At(nScalerHists-1)->GetName())==0)
    {
        cout << "    Appending to: " << name << " with scaler " << scal << " named " << label << endl;
        h1 = (TH1D*)scalerHists->At(nScalerHists-1);
        h1->SetBins((h1->GetNbinsX()+1),0,(h1->GetNbinsX()+1));
        h1->GetXaxis()->SetBinLabel(h1->GetNbinsX(),label);
        nScalerSets.at(nScalerHists-1) += 1;
    }
    else
    {
        cout << "Adding histogram: " << name << " with scaler " << scal << " named " << label << endl;
        h1 = new TH1D(name,name,1,0,1);
        h1->GetXaxis()->SetBinLabel(1,label);
        scalerHists->Add(h1);
        nScalerSets.push_back(1);
        nScalerHists++;
    }

    scalerChanL.push_back(scal);
    scalerChanH.push_back(scal);
}

void PPhysics::GoosyTagger(TH1* hist)
{
    Int_t nBins = hist->GetNbinsX();

    TH1* temp = (TH1D*) hist->Clone("temp");

    Int_t bin;
    for (Int_t i = 1; i <= nBins; i++)
    {
        if(i <= 24)
        {
            if(TMath::Even(i)) bin = (i/2);
            else bin = (13 + (i/2));
        }
        else bin = i;

        hist->SetBinContent(i,temp->GetBinContent(bin));
    }
    delete temp;
}

void PPhysics::GoosyVuprom(TH1* hist)
{
    Int_t nBins = hist->GetNbinsX();

    TH1* temp = (TH1D*) hist->Clone("temp");

    Int_t bin;
    for (Int_t i = 1; i <= nBins; i++)
    {
        if(i <= 12) bin = (2*i);
        else if(i <= 24) bin = (1+(2*(i-13)));
        else bin = i;

        hist->SetBinContent(i,temp->GetBinContent(bin));
    }
    delete temp;
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
    {
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
        {
                FillMissingMass(tree, i, j, Hprompt, Hrandom);
        }
    }
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        	FillMissingMass(tree, particle_index, i, Hprompt, Hrandom);
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom)
{
    if(RejectTagged(tagger_index)) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);

    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

    // Fill TH1
    if (GHistBGSub::IsPrompt(time)) Hprompt->Fill(missingp4.M());
	if (GHistBGSub::IsRandom(time)) Hrandom->Fill(missingp4.M());						
}

void PPhysics::FillTime(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            if(RejectTagged(j)) continue;

            time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			Hist->Fill(time);
		}
	}
}

void PPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        if(RejectTagged(j)) continue;
	
        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
        Hist->Fill(time);
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            if(RejectTagged(j)) continue;

            time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
		}
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        if(RejectTagged(j)) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) Hist->Fill(time);
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, TH1* Hist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        Hist->Fill(tree.GetMass(i));
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, TH1* Hist)
{
    Hist->Fill(tree.GetMass(particle_index));
}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, Hprompt, Hrandom, MM_min, MM_max);
    }

}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, TH1* Hprompt, TH1* Hrandom, Double_t MM_min, Double_t MM_max)
{
    if(RejectTagged(tagger_index)) return;

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

Double_t PPhysics::CalcCoplanarity(const GTreeParticle& tree1, Int_t particle_index1, const GTreeParticle& tree2, Int_t particle_index2)
{
   Double_t phi1 = tree1.GetPhi(particle_index1);
   Double_t phi2 = tree2.GetPhi(particle_index2);
   Double_t phidiff = TMath::Abs(phi1 - phi2);

   return phidiff;
}

// ----------------------------------------------------------------------------------------
// GH1 routines
// ----------------------------------------------------------------------------------------

void PPhysics::FillMissingMass(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMass(tree, i, j, gHist, TaggerBinning);
			//cout<<i<<endl;
		}
	}
	  GetTracks()->GetNTracks();
	  //cout<<GetTracks()->GetNTracks()<<endl;
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMass(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{
    if(RejectTagged(tagger_index)) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
    //    cout << "Missing momentum=" << missingp4.M() << endl;

   // Fill GH1
   if(TaggerBinning)   gHist->Fill(missingp4.M(),time, GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(missingp4.M(),time);

}

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
    beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
	missingp4 	= beam + target - particle;						

	return missingp4;
}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
    {
        FillBeamAsymmetry(tree, particle_index, i, gHist, TaggerBinning);
    }

}

void PPhysics::FillBeamAsymmetry(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning, Double_t MM_min, Double_t MM_max)
{
    if(RejectTagged(tagger_index)) return;

    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);
    if((missingp4.M() < MM_min) || (missingp4.M() > MM_max)) return;
   
   if(TaggerBinning)   gHist->Fill(tree.GetPhi(particle_index),time,GetTagger()->GetTaggedChannel(tagger_index));
   else gHist->Fill(tree.GetPhi(particle_index),time);

}

void PPhysics::FillTime(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            if(RejectTagged(j)) continue;

            time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			gHist->Fill(time);
		}
	}
}

void PPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        if(RejectTagged(j)) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		gHist->Fill(time);
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
            if(RejectTagged(j)) continue;

            time = GetTagger()->GetTaggedTime(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
		}
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
	{
        if(RejectTagged(j)) continue;

        time = GetTagger()->GetTaggedTime(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        gHist->Fill(tree.GetMass(i));
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    gHist->Fill(tree.GetMass(particle_index));
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
        cout << "Target mass unknown!" << endl << endl;
	}
	else if(sscanf( config.c_str(), "%lf\n", &mass) == 1)
	{
        cout << "Setting target mass: " << mass << " MeV" << endl << endl;
		SetTarget(mass);		
	}
	else 
	{
        cout << "Target mass not set correctly" << endl << endl;
		return kFALSE;
	}

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
           cout << "Invalid Tagger channel cut: " << tc1 << endl << endl;
		   return kFALSE;
		}
		else if ((tc2 < 0) || (tc2 > 352))
		{
           cout << "Invalid Tagger channel cut: " << tc2 << endl << endl;
		   return kFALSE;
		}
		
        cout << "Setting cut on Tagger channels: " << tc1 << " to " << tc2 << endl << endl;
		SetTC_cut(tc1,tc2);
	}
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
        cout << "Tagger channel cut not set correctly" << endl << endl;
		return kFALSE;
	}

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
        AddScalerHist("TaggerAccScal",sc1,sc2);
        cout << endl;
    }
	else if(strcmp(config.c_str(), "nokey") != 0)
	{
        cout << "Tagger scaler channels not set correctly" << endl << endl;
		return kFALSE;
	}

	return kTRUE;

}

Bool_t 	PPhysics::InitLiveTimeScalers()
{
    Int_t sc1, sc2, sc3;
    string config = ReadConfig("Live-Time-Scalers");
    if(sscanf( config.c_str(), "%d %d %d\n", &sc1, &sc2, &sc3) == 2)
    {
        cout << "Setting live time scaler channels: clock is " << sc1 << " and inhibited is " << sc2 << endl;
        SetLT_scalers(sc1,sc2);
        AddScalerHist("LiveTimeScal",sc1,"Clock");
        AddScalerHist("LiveTimeScal",sc2,"Inhibited");
        cout << endl;
    }
    else if(sscanf( config.c_str(), "%d %d %d\n", &sc1, &sc2, &sc3) == 3)
    {
        cout << "Setting live time scaler channels: clock is " << sc1 << ", inhibited is " << sc2 << ", and tagger inhibited is " << sc3 << endl;
        SetLT_scalers(sc1,sc2,sc3);
        AddScalerHist("LiveTimeScal",sc1,"Clock");
        AddScalerHist("LiveTimeScal",sc2,"Inhibited");
        AddScalerHist("LiveTimeScal",sc3,"Tagger Inhib");
        cout << endl;
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Live time scaler channels not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t  PPhysics::InitDisplayScalers()
{
    string config;
    Int_t instance = 0;
    char name[256];
    Int_t lo, hi;

    do
    {
        config = ReadConfig("Display-Scalers",instance);

        if(sscanf( config.c_str(), "%s %d %d\n", name, &lo, &hi) == 3)
        {
            AddScalerHist(name,lo,hi);
            instance++;
        }
        else if(strcmp(config.c_str(), "nokey") != 0)
        {
            cout << "Display scalers not set correctly" << endl;
            return kFALSE;
        }
    } while (strcmp(config.c_str(), "nokey") != 0);

    if(instance) cout << endl;

    return kTRUE;
}



Bool_t 	PPhysics::InitDecodeDoubles()
{
    Int_t dd;
    string config = ReadConfig("Decode-Doubles");
    if(sscanf( config.c_str(), "%d\n", &dd) == 1)
    {
        cout << "Setting decoding of doubles: " << dd << endl << endl;
        SetDecodeDoubles(dd);
    }
    else if(strcmp(config.c_str(), "nokey") != 0)
    {
        cout << "Decoding of doubles not set correctly" << endl << endl;
        return kFALSE;
    }

    return kTRUE;

}

Bool_t 	PPhysics::RejectTagged(Int_t tagger_index)
{
    Bool_t reject = false;

    // Is tagger channel rejected by user?
    if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) reject = true;
    if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) reject = true;

    // Is tagger hit a decoded double?
    if(GetTagger()->GetTaggedDouble(tagger_index)) reject = true;

    return reject;
}

Bool_t 	PPhysics::RejectDouble(Int_t tagger_index)
{
    Bool_t reject = false;

    // Is tagger channel rejected by user?
    if(GetTagger()->GetDoubleRandom(tagger_index) < TC_cut_min) reject = true;
    if(GetTagger()->GetDoubleRandom(tagger_index) > TC_cut_max) reject = true;

    return reject;
}

void	PPhysics::ProcessScalerRead()
{
    // Fill Scaler Histograms
    Int_t i = 0;
    for(Int_t j=0; j<nScalerHists; j++)
    {
        Int_t firstBin = 1;
        for(Int_t k=0; k<nScalerSets.at(j); k++)
        {
            FillScalers(scalerChanL.at(i),scalerChanH.at(i),(TH1*)scalerHists->At(j),firstBin);
            firstBin += (scalerChanH.at(i)-scalerChanL.at(i)+1);
            i++;
        }
    }
}

//my added histograms
void PPhysics::FillMissingMomentum(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMomentum(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
  for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMomentum(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
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

Double_t PPhysics::CalcMissingMomentum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.Rho()/197.3;
}

void PPhysics::FillMissingMomentumD(const GTreeParticle& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillMissingMomentum(tree, i, j, gHist, TaggerBinning);
		}
	}
}

void PPhysics::FillMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
  for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillMissingMomentum(tree, particle_index, i, gHist, TaggerBinning);
	}
}

void PPhysics::FillMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
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

void PPhysics::FillDeltaE(const GTreeMeson& tree, GH1* gHist, Bool_t TaggerBinning)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < GetTagger()->GetNTagged(); j++)
		{
			FillDeltaE(tree, i, j, gHist,TaggerBinning);
		}
	}
}

void PPhysics::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, GH1* gHist, Bool_t TaggerBinning)
{
    for (Int_t i = 0; i < GetTagger()->GetNTagged(); i++)
	{
        FillDeltaE(tree, particle_index, i, gHist,TaggerBinning);
	}
}

void PPhysics::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist, Bool_t TaggerBinning)
{

  // Is tagger channel rejected by user?
  if(GetTagger()->GetTaggedChannel(tagger_index) < TC_cut_min) return;
  if(GetTagger()->GetTaggedChannel(tagger_index) > TC_cut_max) return;
    // calc particle time diff
    time = GetTagger()->GetTaggedTime(tagger_index) - tree.GetTime(particle_index);
    
    // Fill GH1
    gHist->Fill(CalcDeltaED(tree,particle_index,tagger_index),time);					

}

Double_t PPhysics::CalcDeltaED(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index) {
  
  beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
 
 
  int n_sub_part;
  
  Double_t  E1, E2;
  Double_t E_diff=10000; // Value out of range in case one does not have just 2 photons

  // number of meson is ever 1 on the tree
  n_sub_part = tree.GetNSubParticles(0);
  //  Int_t n_sub_phot_ind[n_sub_part];
  std::vector<Int_t> n_sub_phot_ind=tree.GetTrackIndexList(0);

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

Double_t PPhysics::CalcDeltaED_corr(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index) {
  
  beam 		= TLorentzVector(0.,0.,GetTagger()->GetTaggedEnergy(tagger_index),GetTagger()->GetTaggedEnergy(tagger_index));
 
 
  int n_sub_part;
  
  Double_t  E1, E2;
  Double_t E_diff=10000; // Value out of range in case one does not have just 2 photons

  // number of meson is ever 1 on the tree
  n_sub_part = tree.GetNSubParticles(0);
  //  Int_t n_sub_phot_ind[n_sub_part];
  std::vector<Int_t> n_sub_phot_ind=tree.GetTrackIndexList(0);
  
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

Double_t PPhysics::CalcMissingMomentumD(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.Rho()/197.3;
}

#endif
