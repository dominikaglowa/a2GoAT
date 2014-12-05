#include "PPi0Example.h"

PPi0Example::PPi0Example()
{ 
    time 	= new GH1("time", 	"time", 	1400, -700, 700);
    time_cut 	= new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

    time_2g 	= new GH1("time_2g",	"time_2g", 	1400, -700, 700);
    time_2g_cut = new GH1("time_2g_cut","time_2g_cut", 	1400, -700, 700);

    IM 		= new GH1("IM", 	"IM", 		400,   0, 400);
    IM_2g 	= new GH1("IM_2g", 	"IM_2g", 	400,   0, 400);
  
    MM		= new GH1("MM", 	"MM", 	 	400,   800, 1200);     
    MM_2g	= new GH1("MM_2g", 	"MM_2g", 	400,   800, 1200);

    TaggerAccScal = new TH1D("TaggerAccScal","TaggerAccScal",352,0,352);
    //stuff from the old file
    MMom		= new GH1("MMom", 	"MMom;q(fm^{-1})", 	 	400, 0., 2.0);     
    MMom_2g	= new GH1("MMom_2g", 	"MMom_2g;q(fm^{-1})", 	400,   0, 2.0); 
    //
}

PPi0Example::~PPi0Example()
{
}

Bool_t	PPi0Example::Init(const char* configfile)
{
	cout << "Initialising physics analysis..." << endl;
	cout << "--------------------------------------------------" << endl << endl;
	if(configfile) SetConfigFile(configfile);

	if(!InitBackgroundCuts()) return kFALSE;
	if(!InitTargetMass()) return kFALSE;
	if(!InitTaggerChannelCuts()) return kFALSE;
	if(!InitTaggerScalers()) return kFALSE;
	cout << "--------------------------------------------------" << endl;
	return kTRUE;
}

Bool_t	PPi0Example::Start()
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

void	PPi0Example::ProcessEvent()
{
	// fill time diff (tagger - pi0), all pi0
	FillTime(*pi0,time);
	FillTimeCut(*pi0,time_cut);
	
	// fill missing mass, all pi0
	FillMissingMass(*pi0,MM);	
	
	// fill invariant mass, all pi0
	FillMass(*pi0,IM);

   //stuff from the old file
	// fill missing momentum, all pi0
	FillMissingMomentum(*pi0,MMom);
	//	
    // Some neutral decays
    for (Int_t i = 0; i < pi0->GetNParticles(); i++)
    {
        // Fill MM for 2 photon decay
        if ((pi0->GetNSubParticles(i) == 2) & (pi0->GetNSubPhotons(i) == 2))
        {
		// fill time diff (tagger - pi0), this pi0
		FillTime(*pi0,i,time_2g);
		FillTimeCut(*pi0,i,time_2g_cut);
			
		// fill missing mass, this pi0
            	FillMissingMass(*pi0,i,MM_2g);
            
		// fill invariant mass, this pi0
            FillMass(*pi0,i,IM_2g);
        }

    }

}

void	PPi0Example::ProcessScalerRead()
{
	// Fill Tagger Scalers
	FillScalers(GetTC_scaler_min(),GetTC_scaler_max(),TaggerAccScal);
}

Bool_t	PPi0Example::Write()
{
	// Write some TH1s
	GTreeManager::Write(TaggerAccScal);

	// Write all GH1's easily
	GTreeManager::Write();
}
