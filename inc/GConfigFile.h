#ifndef __GConfigFile_h__
#define __GConfigFile_h__


#include <string>
#include <iostream>
#include <TROOT.h>

class  GConfigFile
{
private:
    std::string	globalConfigFile;
    std::vector<std::string> inputFileList;
    std::vector<std::string> outputFileList;

protected:

public:
    GConfigFile();
    GConfigFile(const Char_t* configFile);
    virtual ~GConfigFile();

    virtual Bool_t	Init(const char* configFile);
            Bool_t	Init() {Init(globalConfigFile.c_str());}
            void    SetConfigFile(const Char_t* configFile)	{globalConfigFile = configFile;}
    const   Char_t* GetConfigFile() const {return globalConfigFile.c_str();}
            Bool_t  BaseConfig(const int argc, char* argv[], const std::string& defaultInputPrefix, const std::string& defaultOutputPrefix);
    const   Int_t   GetNFiles() {return inputFileList.size();}
    std::string GetInputFile(const Int_t i) {return inputFileList.at(i);}
    std::string GetOutputFile(const Int_t i) {return outputFileList.at(i);}

    std::string	ReadConfig(const std::string& inputKey, const Int_t instance, const Char_t* configName);
    std::string	ReadConfig(const std::string& inputKey, const Int_t instance)                             {return ReadConfig(inputKey, instance, globalConfigFile.c_str());}
    std::string	ReadConfig(const std::string& inputKey)                                                   {return ReadConfig(inputKey, 0, globalConfigFile.c_str());}
};

#endif
