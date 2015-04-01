#include "GMesonReconstruction.h"
#include "GParticleReconstruction.h"

using namespace std;


GMesonReconstruction::GMesonReconstruction()    :
    widthNeutralPion(0),
    widthEta(0),
    widthEtaPrime(0)
{
}

GMesonReconstruction::~GMesonReconstruction()
{

}

Bool_t GMesonReconstruction::Start()
{
    GetNeutralPions()->CloseForInput();
    GetEtas()->CloseForInput();
    GetEtaPrimes()->CloseForInput();

    if(!TraverseEntries(0, GetPhotons()->GetNEntries()))		return kFALSE;

    return kTRUE;
}

Bool_t	GMesonReconstruction::Init()
{

    string  config = ReadConfig("Do-Meson-Reconstruction");
    if (strcmp(config.c_str(), "nokey") == 0)
    {
        mesonThetaMin = 0.0;
        mesonThetaMax = 180.0;
    }
    else if(sscanf( config.c_str(), "%*d %lf %lf\n", &mesonThetaMin, &mesonThetaMax) == 2)
    {
        cout << "meson reconstruction is active over theta range [" <<
        mesonThetaMin << "," << mesonThetaMax <<"]" << endl;
    }
    else
    {
        mesonThetaMin = 0.0;
        mesonThetaMax = 180.0;
    }

    config = ReadConfig("Cut-IM-Width-Pi0");
    sscanf( config.c_str(), "%lf\n", &widthNeutralPion);
    if(widthNeutralPion) cout << "Pi0 IM width cut set to " << widthNeutralPion << " MeV" << endl;
    else
    {
        widthNeutralPion = DEFAULT_WIDTH_NEUTRAL_PION;
        cout << "Pi0 IM width cut set to default (" << widthNeutralPion << " MeV)" << endl;
    }

    config = ReadConfig("Cut-IM-Width-Eta");
    sscanf( config.c_str(), "%lf\n", &widthEta);
    if(widthNeutralPion) cout << "Eta IM width cut set to " << widthEta << " MeV" << endl;
    else
    {
        widthEta = DEFAULT_WIDTH_ETA;
        cout << "Pi0 IM width cut set to default (" << widthEta << " MeV)" << endl;
    }

    config = ReadConfig("Cut-IM-Width-Eta-Prime");
    sscanf( config.c_str(), "%lf\n", &widthEtaPrime);
    if(widthEtaPrime) cout << "Eta-Prime IM width cut set to " << widthEtaPrime << " MeV" << endl;
    else
    {
        widthEtaPrime = DEFAULT_WIDTH_ETA_PRIME;
        cout << "Eta-Prime IM width cut set to default (" << widthEtaPrime << " MeV)" << endl;
    }
    cout << endl;

    return kTRUE;
}

Bool_t  GMesonReconstruction::ProcessEventWithoutFilling()
{
    GetNeutralPions()->Clear();
    GetEtas()->Clear();
    GetEtaPrimes()->Clear();

    Int_t       maxSubs = GetRootinos()->GetNParticles() + GetPhotons()->GetNParticles() + GetChargedPions()->GetNParticles();

    Int_t		index1	  [maxSubs * maxSubs];
    Int_t		index2	  [maxSubs * maxSubs];
    Int_t 		tempID    [maxSubs * maxSubs];
    Double_t 	diff_meson[maxSubs * maxSubs];
    Int_t 		sort_index[maxSubs * maxSubs];
    Bool_t 		is_meson  [maxSubs];

    Int_t       countRootinos	= 0;
    Int_t       countPhotons	= 0;
    Int_t       countChargedPi	= 0;


    TLorentzVector  daughter_list[maxSubs];
    Int_t       	daughter_index[maxSubs];
    Int_t       	pdg_list[maxSubs];

    // one list including only GetRootinos() and GetPhotons()
    Int_t 			ndaughter   = 0;
    
	// another list including charged pions ((and electrons)) too (eta/eta')
    Int_t 			ndaughter_full   = 0;

//###################################
//************From Here**************
//###################################

    TLorentzVector	reaction_p4;
    TLorentzVector	reaction_p4_full;
    
    // LEVEL 1:
    // Test full reaction 4 momentum
    // We will ignoring protons and neutrons in all reconstructions
    // We will include only GetPhotons() and GetRootinos() in the pi0 reconstruction
    // This is to test the following complex decays
    // n' -> pi+  pi-  n
    // n' -> pi0  pi0  n
    // n  -> pi0  pi0  pi0
    // n  -> pi0  pi+  pi-
    // n  -> pi0 (pi+  pi-  g)  - omega meson intermediate state
    // n  -> pi+  pi-  g		- direct n decay
    // 							    (or rho_0 intermediate state)    

    for (Int_t i = 0; i < GetRootinos()->GetNParticles(); i++)
    {
        if (GetRootinos()->GetTheta(i) < mesonThetaMin) continue; // user rejected theta region
        if (GetRootinos()->GetTheta(i) >= mesonThetaMax) continue; // user rejected theta region
		
		// reject zero energy particles (no CB cluster involved)
        if (GetRootinos()->GetClusterEnergy(i) == 0) continue;
		
        is_meson[ndaughter] = kFALSE;

        daughter_list[ ndaughter] = GetRootinos()->Particle(i);
        daughter_index[ndaughter] = i;
        pdg_list[ndaughter] 	  = PDG_ROOTINO;

        reaction_p4 	 += GetRootinos()->Particle(i);
        reaction_p4_full += GetRootinos()->Particle(i);
        
		ndaughter++;   
        ndaughter_full++;
        
        countRootinos++;

    }
    for (Int_t i = 0; i < GetPhotons()->GetNParticles(); i++)
    {
        if (GetPhotons()->GetTheta(i) < mesonThetaMin) continue; // user rejected theta region
        if (GetPhotons()->GetTheta(i) >= mesonThetaMax) continue; // user rejected theta region

		// reject zero energy particles (no CB cluster involved)
        if (GetPhotons()->GetClusterEnergy(i) == 0) continue;

        is_meson[ndaughter] = kFALSE;

        daughter_list[ ndaughter] = GetPhotons()->Particle(i);
        daughter_index[ndaughter] = i;
        pdg_list[ndaughter] 	  = pdgDB->GetParticle("gamma")->PdgCode();;

        reaction_p4 	 += GetPhotons()->Particle(i);
        reaction_p4_full += GetPhotons()->Particle(i);
        
		ndaughter++;   
        ndaughter_full++;
        
        countPhotons++;
    }
    for (Int_t i = 0; i < GetChargedPions()->GetNParticles(); i++)
    {
        if (GetChargedPions()->GetTheta(i) < mesonThetaMin) continue; // user rejected theta region
        if (GetChargedPions()->GetTheta(i) >= mesonThetaMax) continue; // user rejected theta region

		// reject zero energy particles (no CB cluster involved)
		// note, not currently possible for charged pions 
		// included for completeness
        if (GetChargedPions()->GetClusterEnergy(i) == 0) continue;

        is_meson[ndaughter_full] = kFALSE;

        daughter_list[ ndaughter_full] = GetChargedPions()->Particle(i);
        daughter_index[ndaughter_full] = i;
        pdg_list[ndaughter_full] 	   = pdgDB->GetParticle("pi+")->PdgCode();;

        reaction_p4_full += GetChargedPions()->Particle(i);
        ndaughter_full++;
      
        countChargedPi++;

    }

    Double_t diff_pi0  = TMath::Abs( reaction_p4.M() 	  - (pdgDB->GetParticle("pi0" )->Mass()*1000) )/widthNeutralPion;
    Double_t diff_eta  = TMath::Abs( reaction_p4_full.M() - (pdgDB->GetParticle("eta" )->Mass()*1000) )/widthEta;
    Double_t diff_etap = TMath::Abs( reaction_p4_full.M() - (pdgDB->GetParticle("eta'")->Mass()*1000) )/widthEtaPrime;

    if ((diff_pi0 <= 1.0) && (diff_pi0 < diff_eta) && (diff_pi0 < diff_etap) && (ndaughter >= 2))
    {
		// Add pi0 (charged pion list not included)
        GetNeutralPions()->AddParticle(countRootinos, countPhotons, 0, daughter_index, daughter_list);
						
        // Remove GetRootinos() and GetPhotons() from original particle list
        GetRootinos()->RemoveParticles(countRootinos, daughter_index);
        GetPhotons()->RemoveParticles(countPhotons, &daughter_index[countRootinos]);

		// Complete reaction satisfied, return kTRUE
        return kTRUE;
    }
    else if ((diff_eta <= 1.0) && (diff_eta < diff_pi0) && (diff_eta < diff_etap) && (ndaughter_full >= 2))
    {
        // Add eta 
        GetEtas()->AddParticle(countRootinos, countPhotons, countChargedPi, daughter_index, daughter_list);

		// Remove daughters from original particle list
        GetRootinos()->RemoveParticles(countRootinos, daughter_index);
        GetPhotons()->RemoveParticles(countPhotons, &daughter_index[countRootinos]);
        GetChargedPions()->RemoveParticles(countChargedPi, &daughter_index[countRootinos+countPhotons]);
        
        // Complete reaction satisfied, return kTRUE
        return kTRUE;
    }
    else if ((diff_etap <= 1.0) && (diff_etap < diff_pi0) && (diff_etap < diff_eta) && (ndaughter_full >= 2))
    {
        // Add eta prime
        GetEtaPrimes()->AddParticle(countRootinos, countPhotons, countChargedPi, daughter_index, daughter_list);

		// Remove daughters from original particle list
        GetRootinos()->RemoveParticles(countRootinos, daughter_index);
        GetPhotons()->RemoveParticles(countPhotons, &daughter_index[countRootinos]);
        GetChargedPions()->RemoveParticles(countChargedPi, &daughter_index[countRootinos+countPhotons]);
        
        // Complete reaction satisfied, return kTRUE
        return kTRUE;
    }

    // LEVEL 2:
    // Well that didn't work, let's try to make some 2 particle checks
    // Loop over possible 2-particle combinations (skip i=j, ij = ji)
    // to check pi0 -> 2g, n -> 2g , n' -> 2g
    // Find all pairs within IM limits and sort by best Chi
    // Use only GetPhotons() and GetRootinos()
    // Don't double count in sorting!

    Int_t k = 0;
    for (Int_t i = 0; i < ndaughter; i++)
    {
        for (Int_t j = i+1; j < ndaughter; j++)
        {
            TLorentzVector p4 = daughter_list[i] + daughter_list[j];

            Double_t diff_pi0  = TMath::Abs( p4.M() - (pdgDB->GetParticle("pi0" )->Mass()*1000) )/widthNeutralPion;
            Double_t diff_eta  = TMath::Abs( p4.M() - (pdgDB->GetParticle("eta" )->Mass()*1000) )/widthEta;
            Double_t diff_etap = TMath::Abs( p4.M() - (pdgDB->GetParticle("eta'")->Mass()*1000) )/widthEtaPrime;

            if ((diff_pi0 <= 1.0) && (diff_pi0 < diff_eta) && (diff_pi0 < diff_etap))
            {
                diff_meson[k] 	= diff_pi0;
                tempID[k] 		= pdgDB->GetParticle("pi0")->PdgCode();
                index1[k]		= i;
                index2[k]		= j;
                k++;
            }
            else if ((diff_eta <= 1.0) && (diff_eta < diff_pi0) && (diff_eta < diff_etap))
            {
                diff_meson[k]	= diff_eta;
                tempID[k] 		= pdgDB->GetParticle("eta")->PdgCode();
                index1[k]		= i;
                index2[k]		= j;
                k++;
            }
            else if ((diff_etap <= 1.0) && (diff_etap < diff_pi0) && (diff_etap < diff_eta))
            {
                diff_meson[k]	= diff_etap;
                tempID[k] 		= pdgDB->GetParticle("eta'")->PdgCode();
                index1[k]		= i;
                index2[k]		= j;
                k++;
            }
        }
    }

    TMath::Sort(k, diff_meson, sort_index, kFALSE);

    Int_t		nIndex_rootino_delete    = 0;
    Int_t		nIndex_photon_delete    = 0;
 //   Int_t		nIndex_chargedPion_delete = 0;
    Int_t		index_rootino_delete[GTreeTrack_MAX];
    Int_t		index_photon_delete[GTreeTrack_MAX];
  //  Int_t		index_chargedPion_delete[GTreeParticle_MaxEntries];

    for (Int_t i = 0; i < k; i++)
    {
        //particle pair already involved in a meson reconstruction?
        if(is_meson[index1[sort_index[i]]] == kTRUE)  continue;
        if(is_meson[index2[sort_index[i]]] == kTRUE)  continue;

        // New meson identified!
        is_meson[index1[sort_index[i]]] = kTRUE;
        is_meson[index2[sort_index[i]]] = kTRUE;

        // Add to particle list
        if(tempID[sort_index[i]] == pdgDB->GetParticle("pi0")->PdgCode())
        {
            GetNeutralPions()->AddParticle(daughter_index[index1[sort_index[i]]], daughter_list[index1[sort_index[i]]], pdg_list[index1[sort_index[i]]], daughter_index[index2[sort_index[i]]], daughter_list[index2[sort_index[i]]], pdg_list[index2[sort_index[i]]]);
            if(index1[sort_index[i]] < countRootinos)
            {
                index_rootino_delete[nIndex_rootino_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_rootino_delete++;
            }
            else if(index1[sort_index[i]] < (countRootinos+countPhotons))
            {
                index_photon_delete[nIndex_photon_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_photon_delete++;
            }
 /*           else
            {
                index_chargedPion_delete[nIndex_chargedPion_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_chargedPion_delete++;
            }*/
            if(index2[sort_index[i]] < countRootinos)
            {
                index_rootino_delete[nIndex_rootino_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_rootino_delete++;
            }
            else if(index2[sort_index[i]] < (countRootinos+countPhotons))
            {
                index_photon_delete[nIndex_photon_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_photon_delete++;
            }
/*            else
            {
                index_chargedPion_delete[nIndex_chargedPion_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_chargedPion_delete++;
            }*/
        }
        else if(tempID[sort_index[i]] == pdgDB->GetParticle("eta")->PdgCode())
        {
            GetEtas()->AddParticle(daughter_index[index1[sort_index[i]]], daughter_list[index1[sort_index[i]]], pdg_list[index1[sort_index[i]]], daughter_index[index2[sort_index[i]]], daughter_list[index2[sort_index[i]]], pdg_list[index2[sort_index[i]]]);
            if(index1[sort_index[i]] < countRootinos)
            {
                index_rootino_delete[nIndex_rootino_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_rootino_delete++;
            }
            else if(index1[sort_index[i]] < (countRootinos+countPhotons))
            {
                index_photon_delete[nIndex_photon_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_photon_delete++;
            }
/*          else
            {
                index_chargedPion_delete[nIndex_chargedPion_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_chargedPion_delete++;
            }*/
            if(index2[sort_index[i]] < countRootinos)
            {
                index_rootino_delete[nIndex_rootino_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_rootino_delete++;
            }
            else if(index2[sort_index[i]] < (countRootinos+countPhotons))
            {
                index_photon_delete[nIndex_photon_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_photon_delete++;
            }
/*            else
            {
                index_chargedPion_delete[nIndex_chargedPion_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_chargedPion_delete++;
            }*/
        }
        else if(tempID[sort_index[i]] == pdgDB->GetParticle("eta'")->PdgCode())
        {
            GetEtaPrimes()->AddParticle(daughter_index[index1[sort_index[i]]], daughter_list[index1[sort_index[i]]], pdg_list[index1[sort_index[i]]], daughter_index[index2[sort_index[i]]], daughter_list[index2[sort_index[i]]], pdg_list[index2[sort_index[i]]]);
            if(index1[sort_index[i]] < countRootinos)
            {
                index_rootino_delete[nIndex_rootino_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_rootino_delete++;
            }
            else if(index1[sort_index[i]] < (countRootinos+countPhotons))
            {
                index_photon_delete[nIndex_photon_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_photon_delete++;
            }
 /*           else
            {
                index_chargedPion_delete[nIndex_chargedPion_delete] = daughter_index[index1[sort_index[i]]];
                nIndex_chargedPion_delete++;
            }*/
            if(index2[sort_index[i]] < countRootinos)
            {
                index_rootino_delete[nIndex_rootino_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_rootino_delete++;
            }
            else if(index2[sort_index[i]] < (countRootinos+countPhotons))
            {
                index_photon_delete[nIndex_photon_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_photon_delete++;
            }
/*            else
            {
                index_chargedPion_delete[nIndex_chargedPion_delete] = daughter_index[index2[sort_index[i]]];
                nIndex_chargedPion_delete++;
            }*/
        }
    }
    GetRootinos()->RemoveParticles(nIndex_rootino_delete, index_rootino_delete);
    GetPhotons()->RemoveParticles(nIndex_photon_delete, index_photon_delete);
//    GetChargedPions()->RemoveParticles(nIndex_chargedPion_delete, index_chargedPion_delete);


    return kTRUE;
}

void  GMesonReconstruction::ProcessEvent()
{
    if(!ProcessEventWithoutFilling())   return;

    GetEventParameters()->SetNReconstructed(GetNReconstructed());
    GetEventParameters()->Fill();

    GetNeutralPions()->Fill();
    GetEtas()->Fill();
    GetEtaPrimes()->Fill();
    FillReadList();
}
