#include "GHistManager.h"



GHistManager*   gGHistManager = 0;


GHistManager::GHistManager()    :
    histList()
{
    gGHistManager   = this;
    histList.SetOwner(kFALSE);
}

GHistManager::~GHistManager()
{
    gGHistManager   = 0;
}

void    GHistManager::AddHistogramToList(GHistLinked* hist)
{
    histList.AddAtFree(hist);
}

void    GHistManager::RemoveHistogramFromList(GHistLinked* hist)
{
    histList.Remove(hist);
}

void GHistManager::ClearLinkedHistograms()
{
    TIter   iter(&histList);
    GHistLinked*    hist;
    while(hist=(GHistLinked*)iter.Next())
        hist->Reset();
}

void GHistManager::WriteLinkedHistograms(TDirectory* dir)
{
    std::cout << "Calc Result -->";
    {
        TIter   iter(&histList);
        GHistLinked*    hist;
        while(hist=(GHistLinked*)iter.Next())
            hist->CalcResult();
    }
    std::cout << " Done." << std::endl;

    std::cout << "Prepare Write List -->";
    GHistWriteList  writeList;
    {
        TIter   iter(&histList);
        GHistLinked*    hist;
        while(hist=(GHistLinked*)iter.Next())
            hist->PrepareWriteList(&writeList);
    }
    std::cout << " Done." << std::endl;

    std::cout << "Write -->";
    dir->cd();
    writeList.Write();
    std::cout << " Done." << std::endl;
}





GHistLinked::GHistLinked(Bool_t linkHistogram) :
    TObject(),
    linked(linkHistogram)
{
    if(linked == kTRUE)
        Link();
}

GHistLinked::~GHistLinked()
{
    Unlink();
}

TDirectory* GHistLinked::GetCreateDirectory(const char* name)
{
    TDirectory* ret = gDirectory->GetDirectory(name);
    if(!ret)
    {
        ret = gDirectory->mkdir(name);
        if(!ret)
            ret = gDirectory;
    }
    return ret;
}

void    GHistLinked::Link()
{
    if(linked)
        Unlink();
    if(gGHistManager)
    {
        gGHistManager->AddHistogramToList(this);
        linked = kTRUE;
    }
}

void    GHistLinked::Unlink()
{
    linked = kFALSE;
    if(gGHistManager)
        gGHistManager->RemoveHistogramFromList(this);
}





GHistWriteList::~GHistWriteList()
{
    TIter next(this);
    GHistWriteListEntry *obj;
    while ((obj = (GHistWriteListEntry*)next()))
        delete obj;
}

GHistWriteListEntry*    GHistWriteList::AddDirectory(const TString& _Name)
{
    GHistWriteListEntry*    entry = new GHistWriteListEntry(_Name);
    if(entry)
        AddAtFree(entry);
    return entry;
}

void    GHistWriteList::AddHistogram(TH1 *_Hist, const TString& _Name)
{
    if(!_Hist)
        return;
    GHistWriteListEntry*    entry = new GHistWriteListEntry(_Hist, _Name);
    if(entry)
        AddAtFree(entry);
}

GHistWriteList* GHistWriteList::GetDirectory(const TString& _Name)
{
    TIter next(this);
    GHistWriteListEntry *obj;
    while ((obj = (GHistWriteListEntry*)next()))
    {
        if(obj->isDirectory==kTRUE)
        {
            if(obj->name.EqualTo(_Name))
                return (GHistWriteList*)obj->obj;
        }
    }
    return (GHistWriteList*)AddDirectory(_Name)->obj;
}

void    GHistWriteList::Print()
{
    printf("Print GHistWriteList:\n");
    TIter next(this);
    GHistWriteListEntry *obj;
    while ((obj = (GHistWriteListEntry*)next()))
    {
        if(obj->isDirectory==kTRUE)
            printf("Directory: ");
        printf("%s\n", obj->name.Data());
    }
}

Int_t	GHistWriteList::Write(const char* NotUsed, Int_t option, Int_t bufsize)
{
    Int_t nbytes = 0;
    {
        TIter next(this);
        GHistWriteListEntry *obj;
        while ((obj = (GHistWriteListEntry*)next()))
        {
            if(obj->isDirectory==kFALSE)
                nbytes += obj->Write(NotUsed, option, bufsize);
        }
    }

    {
        TIter next(this);
        GHistWriteListEntry *obj;
        while ((obj = (GHistWriteListEntry*)next()))
        {
            if(obj->isDirectory==kTRUE)
            {
                TDirectory* parentDir   = gDirectory;
                GHistLinked::GetCreateDirectory(obj->name)->cd();
                nbytes += obj->Write(NotUsed, option, bufsize);
                parentDir->cd();
            }
        }
    }
    return nbytes;
}






GHistWriteListEntry::GHistWriteListEntry(const TString& _Name)  : obj(0), name(_Name), isDirectory(kTRUE)
{
    GHistWriteList* list = new GHistWriteList();
    obj = list;
}

GHistWriteListEntry::~GHistWriteListEntry()
{
    if(isDirectory==kTRUE)
        delete obj;
}

Int_t	GHistWriteListEntry::Write(const char* NotUsed, Int_t option, Int_t bufsize)
{
    if(!obj)
        return 0;
    Int_t   res = obj->Write(name, option, bufsize);
    if(isDirectory==kFALSE)
        ((TH1D*)obj)->SetDirectory(0);
    return  res;
}
