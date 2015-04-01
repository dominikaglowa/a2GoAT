#ifndef __GTreeParticle_h__
#define __GTreeParticle_h__


#include "Rtypes.h"
#include "GTreeTrack.h"


class   GTreeMeson;

class  GTreeParticle    : public GTreeTrack
{
private:
    Int_t               nParticles;
    Double_t            mass[GTreeTrack_MAX];
    Int_t               trackIndex[GTreeTrack_MAX];  // index of the corresponding tack in the track list, -1 => No track


protected:
    virtual void    SetBranchAdresses();
    virtual void    SetBranches();
            void    PrintParticle(const Int_t i) const;

public:
    GTreeParticle(GTreeManager *Manager, const TString& _Name);
    virtual ~GTreeParticle();

            void            AddParticle(const Double_t _clusterEnergy = 0, const Double_t _theta = 0, const Double_t _phi = 0, const Double_t _mass = 0, const Double_t _time = 0, const Int_t _clusterSize = 0, const Int_t _centralCrystal = -1, const Int_t _centralVeto = -1, const Int_t _detectors = 0, const Double_t _vetoEnergy = 0, const Double_t _MWPC0Energy = 0, const Double_t _MWPC1Energy = 0, const Int_t _trackIndex = -1);
    virtual void            Clear() {nParticles = 0;}
            Int_t           GetNParticles()                     const	{return nParticles;}
            Double_t        GetKineticEnergy(const Int_t index) const   {return clusterEnergy[index];}
    const   Double_t*       GetKineticEnergy()                  const   {return clusterEnergy; }
            Double_t        GetTotalEnergy(const Int_t index)   const   {return (clusterEnergy[index]+mass[index]);}
            Double_t        GetMass(const Int_t index)          const   {return mass[index];}
    const   Double_t*       GetMass()                           const   {return mass; }
            Int_t           GetTrackIndex(const Int_t index)    const   {return trackIndex[index];}
    const   Int_t*          GetTrackIndex()                     const   {return trackIndex; }

    inline          TLorentzVector Particle(const Int_t particle);
    inline  const   TLorentzVector Particle(const Int_t particle) const;
    virtual void            Print() const;
            void            RemoveParticles(const Int_t nIndices, const Int_t* indices);
            void            RemoveAllParticles();
    virtual Bool_t          Write();

    friend  class GTreeMeson;
};


TLorentzVector         GTreeParticle::Particle(const Int_t particle)
{
    if(particle>=nParticles)
    {
        std::cout << "ERROR: Try to access particle " << particle << " in " << GetName() << ". Only " << nParticles << " in tree." << std::endl;
        static  TLorentzVector  helpVectorUnity(0,0,0,0);
        helpVectorUnity.SetPxPyPzE(0,0,0,0);
        return helpVectorUnity;
    }
    return GetVector(particle,GetMass(particle));
}

const   TLorentzVector GTreeParticle::Particle(const Int_t particle) const
{
    if(particle>=nParticles)
    {
        std::cout << "ERROR: Try to access particle " << particle << " in " << GetName() << ". Only " << nParticles << " in tree." << std::endl;
        static  TLorentzVector  helpVectorUnity(0,0,0,0);
        helpVectorUnity.SetPxPyPzE(0,0,0,0);
        return helpVectorUnity;
    }
    return GetVector(particle,GetMass(particle));
}


#endif
