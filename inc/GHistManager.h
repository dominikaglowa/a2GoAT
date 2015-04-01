#ifndef __GHistManager_h__
#define __GHistManager_h__


#include <iostream>

#include <TBits.h>
#include <TObjArray.h>
#include <TDirectory.h>
#include <TH1.h>



class   GHistLinked;
class   GHistWriteList;
class   GHistWriteListEntry;

class   GHistManager
{
private:
    TObjArray   histList;

    void        AddHistogramToList(GHistLinked* hist);
    virtual     TDirectory* GetOutputDirectory() = 0;
    void        RemoveHistogramFromList(GHistLinked* hist);

protected:

public:
    GHistManager();
    virtual ~GHistManager();

    void    ClearLinkedHistograms();
    void    WriteLinkedHistograms(TDirectory* dir);

    friend  class   GHistLinked;
};





class   GHistLinked : public TObject
{
private:
    Bool_t  linked;

public:
    GHistLinked(Bool_t linkHistogram = kTRUE);
    virtual ~GHistLinked();

    virtual void        CalcResult() = 0;
    virtual Int_t       Fill(Double_t x) = 0;
    static  TDirectory* GetCreateDirectory(const char* name);
            void        Link();
    virtual void        PrepareWriteList(GHistWriteList* arr, const char* name = 0) = 0;
    virtual void        Reset(Option_t* option = "") = 0;
            void        Unlink();
    virtual Int_t       WriteWithoutCalcResult(const char* name = 0, Int_t option = 0, Int_t bufsize = 0) = 0;
    virtual Int_t       Write(const char* name = 0, Int_t option = 0, Int_t bufsize = 0)    {CalcResult(); WriteWithoutCalcResult(name, option, bufsize);}
};



class   GHistWriteList  : public TObjArray
{
private:
    GHistWriteListEntry* AddDirectory(const TString& _Name);

public:
    GHistWriteList()    {}
    ~GHistWriteList();

            void    AddHistogram(TH1* _Hist, const TString& _Name);
    GHistWriteList* GetDirectory(const TString& _Name);
    virtual void    Print();
    virtual Int_t	Write(const char* NotUsed = 0, Int_t option = 0, Int_t bufsize = 0);
};


class   GHistWriteListEntry  : public TObject
{
private:
    TObject*    obj;
    TString     name;
    Bool_t      isDirectory;

public:
    GHistWriteListEntry(TH1* _Hist, const TString& _Name)                              : obj(_Hist), name(_Name), isDirectory(kFALSE)   {}
    GHistWriteListEntry(const TString& _Name);
    ~GHistWriteListEntry();

    virtual Int_t       Write(const char* NotUsed = 0, Int_t option = 0, Int_t bufsize = 0);

    friend class    GHistWriteList;
};

#endif
