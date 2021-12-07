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
#include "VAArrayInfoFactoryLite.h"
#include "VARootIO.h"
#include "VARunHeader.h"
#include "VASkyMap.h"
#include "VAShowerData.h"
#include "VATime.h"
#include "VASlalib.h"

//#include "/work/astropa/sw/veritas/wcstools-3.9.5/libwcs/wcs.h"
#include "MLWCSHandler.h"
#include "wcs.h"
#include "MLRooDataStore.h"

int main(int argc, char * argv[])
{
  /* Arguments are as follows:
   *   1. string - list of files
   *   2. string - output file name that stores TTrees and RooDataSets
   */

  std::cout << "Name of file list: " << std::string(argv[1]) << std::endl;
  //Determine number of files in list
  std::string fileslist=std::string(argv[1]);
  std::string srcName=std::string(argv[2]);
  cout << "Finding list of root files in " << fileslist.c_str() << endl;
  int linenum=0,ii=0; std::string line1,line2;
  std::ifstream readfiles(fileslist);
  std::ifstream loopfiles(fileslist);

  while(std::getline(readfiles, line1)){
    linenum++;
  }
  cout << "Number of root files in " << fileslist.c_str() << ": " << linenum << std::endl;
  MLRooDataStore testfile[linenum];

  while(std::getline(loopfiles, line2)){
    std::string outputname(srcName+"_"+std::to_string(ii+1)+".root") ;
    //testfile[ii].fillFromFile(line2,ii);
    testfile[ii].fillFromFile(line2,outputname);
    ii++;
  }

  std::cout << "Done outputing RooDataSet" << std::endl;
  return 0;
}
