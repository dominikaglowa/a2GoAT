#ifndef __GMesonReconstruction_h__
#define __GMesonReconstruction_h__


#include "GTreeManager.h"

#define DEFAULT_WIDTH_NEUTRAL_PION 20.0
#define DEFAULT_WIDTH_ETA 44.0
#define DEFAULT_WIDTH_ETA_PRIME 60.0


class  GMesonReconstruction  : virtual public GTreeManager
{
private:
    Double_t	widthNeutralPion;
    Double_t	widthEta;
    Double_t	widthEtaPrime;

    Double_t	mesonThetaMin;
    Double_t	mesonThetaMax;

protected:

            Bool_t  ProcessEventWithoutFilling();
    virtual void    ProcessEvent();
    virtual Bool_t  Start();

public:
    GMesonReconstruction();
    virtual ~GMesonReconstruction();

            Bool_t  Init();
            void    SetNeutralPionWidth(const Double_t width)   {widthNeutralPion = width;}
            void    SetEtaWidth(const Double_t width)   {widthEta = width;}
            void    SetEtaPrimeWidth(const Double_t width)  {widthEtaPrime = width;}
};





#endif
