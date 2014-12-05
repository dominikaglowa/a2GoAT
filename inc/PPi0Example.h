#ifndef __PPi0Example_h__
#define __PPi0Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

#include "GTreeManager.h"
#include "PPhysics.h"

class	PPi0Example  : public PPhysics
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
    //stuff from the old file 
    GH1*	MMom;
    GH1*	MMom_2g; 

    GH1*	MMomD;
    GH1*	MMomD_2g; 

    GH1*	DeltaE_CM_D_2g; 
    GH1*       DeltaE_Missmom_BeamE[7000];
    //
    TH1*	TaggerAccScal;

    string  	config;
protected:
    virtual Bool_t  Start();
    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
    virtual Bool_t    Write();
			
public:
    PPi0Example();
    virtual ~PPi0Example();
    virtual Bool_t  Init(const char* configfile);

};
#endif
