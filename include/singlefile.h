// C++ HEADERS
//#include <ctime>
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
#include "VAArrayInfoFactoryLite.h"
#include "VARootIO.h"
#include "VARunHeader.h"
#include "VASkyMap.h"
#include "VAShowerData.h"
#include "VATime.h"
#include "VASlalib.h"
//#include "/work/astropa/sw/veritas/wcstools-3.9.5/libwcs/wcs.h"
#include "wcs.h"

int singlefile(std::string s5);
