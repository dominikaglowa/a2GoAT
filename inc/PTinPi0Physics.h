#ifndef __PTinPi0Physics_h__
#define __PTinPi0Physics_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysicsTin.h"

class	PTinPi0Physics  : public PPhysicsTin
{
private:
    GH1*	time;
    GH1*	time_cut;
    GH1*	time_2g;      
    GH1*	time_2g_cut;   
     
    GH1*	IM;
    GH1*	IM_2g;

    GH1*	MM;
    GH1*	MM_2g;

    //my added histograms
    GH1*	MMom;
    GH1*	MMom_2g;
    GH1*	MMomD;
    GH1*	MMomD_2g;
    GH1*	DeltaE_CM_D_2g;
    //GH2*       DeltaE_Missmom_BeamE[23];

    TH2F*       DeltaE_Missmom_BeamE_Prompt[23];
    TH2F*       DeltaE_Missmom_BeamE_Random[23];


    TH1F*       Test_histo[7000];


    //

    TH1*	TaggerAccScal;

protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PTinPi0Physics();
    virtual ~PTinPi0Physics();
    virtual Bool_t  Init();
    void FillDeltaE_Missmom(const GTreeMeson& tree);
    void FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index);
    void FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index);
	//


};
#endif
