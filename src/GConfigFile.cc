#include "GConfigFile.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include <algorithm>
#include <fstream>
#include <iostream>


GConfigFile::GConfigFile()  :
    globalConfigFile()
{
}

GConfigFile::GConfigFile(const Char_t* configFile)  :
    globalConfigFile(configFile)
{
}

GConfigFile::~GConfigFile()
{

}

Bool_t	GConfigFile::Init(const char* configFile)
{
    std::cout << "No Init function specified for this class." << std::endl;
    return kTRUE;
}

Bool_t	GConfigFile::BaseConfig(const int argc, char* argv[], const std::string& defaultInputPrefix, const std::string& defaultOutputPrefix)
{
    // Initialise strings
    std::string configFile = "";
    std::string serverFile = "";
    std::string inputDirectory = "";
    std::string outputDirectory = "";
    std::string inputFile = "";
    std::string outputFile = "";
    std::string inputPrefix = defaultInputPrefix;
    std::string outputPrefix = defaultOutputPrefix;

    Int_t length;
    std::string flag;

    Bool_t overwrite = kTRUE;
    if(argc == 1)
    {
        system("man ./documents/goat.man");
        return kFALSE;
    }
    else if(argc == 2)
    {
        flag = argv[1];
        flag.erase(0,1);
        if(strcmp(flag.c_str(), "help") == 0) return kFALSE;
        else configFile= argv[1];
    }
    else
    {
        for(Int_t i=1; i<argc; i++)
        {
            flag = argv[i];
            if(flag.find_first_of("-") == 0)
            {
                i++;
                flag.erase(0,1);
                if(strcmp(flag.c_str(), "s") == 0) serverFile = argv[i];
                else if(strcmp(flag.c_str(), "d") == 0) inputDirectory = argv[i];
                else if(strcmp(flag.c_str(), "D") == 0) outputDirectory = argv[i];
                else if(strcmp(flag.c_str(), "f") == 0) inputFile = argv[i];
                else if(strcmp(flag.c_str(), "F") == 0) outputFile = argv[i];
                else if(strcmp(flag.c_str(), "p") == 0) inputPrefix = argv[i];
                else if(strcmp(flag.c_str(), "P") == 0) outputPrefix = argv[i];
                else if(strcmp(flag.c_str(), "n") == 0)
                {
                        overwrite = kFALSE;
                    i--;
                }
                else
                {
                    std::cout << "Unknown flag " << flag << std::endl;
                    return kFALSE;
                }
            }
            else configFile = argv[i];
        }
    }

    // Check that config file exists:
    ifstream cfile(configFile.c_str());
    if(!cfile)
    {
        std::cout << "Config file '" << configFile << "' could not be found." << std::endl;
        return kFALSE;
    }
    SetConfigFile(configFile.c_str());

    // If server file is specified, check that it exists
    if(serverFile.length() > 0)
    {
        // Check that file exists:
        ifstream sfile(serverFile.c_str());
        if(!sfile)
        {
            std::cout << "Server file '" << serverFile << "' could not be found" << std::endl;
            return kFALSE;
        }
    }
    // If no server file is specified, allow for checking in the config file
    else serverFile = configFile;

    // If unset, scan server or config file for file settings
    if(inputDirectory.length() == 0)
    {
        flag = ReadConfig("Input-Directory",0,(Char_t*)serverFile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) inputDirectory = flag;
    }

    if(outputDirectory.length() == 0)
    {
        flag = ReadConfig("Output-Directory",0,(Char_t*)serverFile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) outputDirectory = flag;
    }

    if(inputFile.length() == 0)
    {
        flag = ReadConfig("Input-File",0,(Char_t*)serverFile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) inputFile = flag;
    }

    if(outputFile.length() == 0)
    {
        flag = ReadConfig("Output-File",0,(Char_t*)serverFile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) outputFile = flag;
    }

    if(inputPrefix.length() == 0)
    {
        flag = ReadConfig("Input-Prefix",0,(Char_t*)serverFile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) inputPrefix = flag;
    }

    if(outputPrefix.length() == 0)
    {
        flag = ReadConfig("Output-Prefix",0,(Char_t*)serverFile.c_str());
        flag.erase(0,flag.find_first_not_of(" "));
        if(strcmp(flag.c_str(),"nokey") != 0) outputPrefix = flag;
    }
    // Finished scanning for file settings

    // Fix directories to include final slash if not there
    if(inputDirectory.find_last_of("/") != (inputDirectory.length()-1)) inputDirectory += "/";
    if(outputDirectory.find_last_of("/") != (outputDirectory.length()-1)) outputDirectory += "/";

    // Output user settings (Set to defaults if still unspecified)
    std::cout << std::endl << "User inputs" << std::endl;
    std::cout << "Config file:      '" << configFile << "' chosen" << std::endl;
    if(inputDirectory.length() != 0)  std::cout << "Input directory:  '" << inputDirectory  << "' chosen" << std::endl;
    if(outputDirectory.length() != 0) std::cout << "Output directory: '" << outputDirectory << "' chosen" << std::endl;
    if(inputFile.length() != 0)  	  std::cout << "Input file:       '" << inputFile       << "' chosen" << std::endl;
    if(outputFile.length() != 0) 	  std::cout << "Output file:      '" << outputFile      << "' chosen" << std::endl;
    if(inputPrefix.length() != 0)  	  std::cout << "Input prefix:     '" << inputPrefix     << "' chosen" << std::endl;
    if(outputPrefix.length() != 0)    std::cout << "Output prefix:    '" << outputPrefix    << "' chosen" << std::endl;
    std::cout << std::endl;

    std::string file;
    std::string prefix;
    std::string suffix;

    Int_t filesFound = 0;
    // If input file is specified, use it
    if(inputFile.length() > 0)
    {
        std::cout << "Searching for input file(s)" << std::endl;
        file = inputFile;
        length = file.length();
        // File should at least have '.root' at the end
        if(length >= 5)
        {
            // Add input directory to it
            inputFile = inputDirectory+inputFile;
            std::cout << "Input file  '" << inputFile << "' chosen" << std::endl;

            // If output file is specified, use it
            if(outputFile.length() > 0) outputFile = outputDirectory+outputFile;
            // If output file is not specified, build it
            else
            {
                // If output directory is not specified, build it
                if(outputDirectory.length() == 0)
                {
                    prefix = file.substr(0,file.find_last_of("/")+1);
                    outputDirectory = inputDirectory+prefix;
                }
                // If input prefix doesn't match, simply prepend output prefix to the file name
                if(file.find(inputPrefix)>file.length()) suffix = ("_"+file.substr(file.find_last_of("/")+1,length-(file.find_last_of("/")+1)));
                // If input prefix does match, switch prefixes
                else suffix = file.substr(file.find_last_of("/")+1+inputPrefix.length(),length-(file.find_last_of("/")+1+inputPrefix.length()));
                // Build output file name
                outputFile = outputDirectory+outputPrefix+suffix;
            }

            std::cout << "Output file '" << outputFile << "' chosen" << std::endl << std::endl;
            inputFileList.push_back(inputFile);
            outputFileList.push_back(outputFile);
            filesFound++;
        }
    }
    // Otherwise scan input directory for matching files
    else
    {
        std::cout << "Searching input directory for files matching input prefix" << std::endl;

        // If output directory is not specified, use the input directory
        if(inputDirectory.length()  == 0) inputDirectory = "./";
        if(outputDirectory.length() == 0) outputDirectory = inputDirectory;

        // Create list of files in input directory
        TSystemFile *systemFile;
        TSystemDirectory *systemDirectory = new TSystemDirectory("files",inputDirectory.c_str());
        TList *fileList = systemDirectory->GetListOfFiles();
        fileList->Sort();
        TIter fileIterator(fileList);

        // Iterate over files
        while((systemFile=(TSystemFile*)fileIterator()))
        {
            file = systemFile->GetName();
            length = file.length();
            // File should at least have '.root' at the end
            if(length >= 5)
            {
                //Check that prefixes and suffixes match
                prefix = file.substr(0,inputPrefix.length());
                suffix = file.substr(length-5,5);
                if(((strcmp(prefix.c_str(),inputPrefix.c_str()) == 0)) && (strcmp(suffix.c_str(),".root") == 0))
                {
                    // Build input file name
                    inputFile = inputDirectory+file;
                    // Build output file name
                    suffix = file.substr(inputPrefix.length(),length-inputPrefix.length());
                    outputFile = outputDirectory+outputPrefix+suffix;

                    // Check for previously created output file
                    if((gSystem->IsFileInIncludePath(outputFile.c_str())) && !overwrite) continue;

                    inputFileList.push_back(inputFile);
                    outputFileList.push_back(outputFile);
                    filesFound++;
                }
            }
        }
    }
    if (filesFound == 0)
    {
        std::cout << "ERROR: No AcquRoot files found!" << std::endl;
        return kFALSE;
    }

    return kTRUE;
}

std::string GConfigFile::ReadConfig(const std::string& inputKey, const Int_t instance, const Char_t* configName)
{
    Int_t string_instance = 0;
    std::string key = inputKey;
    std::transform(key.begin(), key.end(),key.begin(), ::toupper);

    std::string str;
    std::string values;

    ifstream configFile;

    configFile.open(configName);
    //std::cout << "config: " << configName << std::endl;

    if (configFile.is_open())
    {
        while ( getline (configFile,str) )
        {
            std::string::size_type begin = str.find_first_not_of(" \f\t\v");
            if(begin == std::string::npos) continue;
            if(std::string("#").find(str[begin]) != std::string::npos)	continue;
            if(std::string("//").find(str[begin]) != std::string::npos)	continue;

            std::string firstWord;

            try {
                firstWord = str.substr(0,str.find(":"));
            }
            catch(std::exception& e) {
                firstWord = str.erase(str.find_first_of(":"),str.find_first_of(":"));
            }
            std::transform(firstWord.begin(),firstWord.end(),firstWord.begin(), ::toupper);

            values = str.substr(str.find(":")+1,str.length());

            if (strcmp(firstWord.c_str(),key.c_str()) == 0)
            {
                if (string_instance == instance)
                {
                    configFile.close();
                    return values;
                }
                else string_instance++;
            }
        }
        configFile.close();
    }

    return "nokey";
}
