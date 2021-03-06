// C++ HEADERS
#include <ctime>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <getopt.h>
#include <regex>

// ROOT INCLUDES
#include "TFile.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TNamed.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// ROOFIT INCLUDES
#include "RooAbsReal.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooDataSet.h"

// VEGAS includes
#include "VAException.h"
#include "VAModel3DData.h"
#include "VAArrayInfoFactoryLite.h"
#include "VARootIO.h"
#include "VARunHeader.h"
#include "VASkyMap.h"
#include "VAShowerData.h"
#include "VACoordinatePair.h"
#include "VATime.h"
#include "VASlalib.h"
#include "wcs.h"

using namespace VACoordinates;


class MLRooDataStore {

public:
    // Default constructor                                                                                                             
    MLRooDataStore();
    virtual ~MLRooDataStore();

    // Some constants to help us out                                                                                                   
    const double r2d = TMath::RadToDeg();
    const double d2r = TMath::DegToRad();

    
    // TODO : change these variables to RooRealProxy (Josh's comment)                                                                     
    // Create some variables needed by the datasets we're going to be creating.                                                          
    RooRealVar xCoord;  // Exagerate the range since we want to contain all events                                                        
    RooRealVar yCoord;  // Exagerate the range since we want to contain all events                                                        
    RooRealVar msw;
    //RooRealVar zenith("zenith","zenith",0.0,45.0);                                                                                      
    RooRealVar cosZ;
    RooRealVar azimuth;
    RooRealVar Erec;
    RooRealVar noise;
    RooRealVar Tels;
    RooRealVar Etrue;
    
    bool fillFromFile(std::string s5, std::string rooname);
    double getAvgAzimuth() ; // Due to the cyclic nature of azimuth, this has to be handled in a special way.
    double findAverageNoise(VAShowerData* showerdata, const VAQStatsData* pQStatsData, const VATime& time);
    void PrintInfo();
    void PrintInput(std::string input_name); 

protected:
    VACoordinatePair* _sourceCoordPair;  // Source position as read from run header                                                   
    VACoordinatePair* _trackCoordPair;   // tracking position calculated from source position and offset
    void saveObsInfo(VARootIO* file);
    int findEpoch(VATime time);
    int findATM_Data(VATime time);
    int findATM_Sims(VASimulationHeader* simHeader);
    double getAvgAzimuth(RooDataSet* ds) ; // Due to the cyclic nature of azimuth, this has to be handled in a special way. 

    RooDataSet* dataout;
    bool _isSim;                   // true if simulation header exists
    long long _runNum;
    double    _livetime;           // livetime in seconds
    double _avgZenith;
    double _avgAzimuth;
    double _avgNoise;
    double _trackRA;
    double _trackDec;
    double _srcRA;
    double _srcDec;
    double _offset{1.2f};
    double _xmean;
    double _ymean;
    double _epoch;
    double _ATM;
    std::string _SourceName;

};

