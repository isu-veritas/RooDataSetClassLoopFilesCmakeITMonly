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
#include "wcs.h"
#include "MLWCSHandler.h"

double findAverageNoise(VAShowerData* showerdata, const VAQStatsData* pQStatsData, const VATime& time) ;
//void Coordinate2Projected_Quick(double xcoord, double ycoord, WorldCoor* wcs, double *xproj, double *yproj) ;

//singlefile("../SourceAnalyses_MLM/Background/Segue1/stage5_MLM/60495.stage5.root")
int singlefile(std::string s5)
{

  // Load the run header and QStats for use in extracting the noise                                                                  
  //const VAQStatsData* qstats = file->loadTheQStatsData() ;

  // Some constants to help us out                                                                                                   
  const double r2d = TMath::RadToDeg();
  const double d2r = TMath::DegToRad();
  bool _isSim = false;
  // TODO : change these variables to RooRealProxy (Josh's comment)                                                                       
  // Create some variables needed by the datasets we're going to be creating.                                                   
  RooRealVar xCoord("xpos","xpos",-2.0, 2.0);  // Exagerate the range since we want to contain all events                        
  RooRealVar yCoord("ypos","ypos",-2.0, 2.0);  // Exagerate the range since we want to contain all events                        
  RooRealVar msw("msw","msw",0.0, 10.0);
  RooRealVar cosZ("cosZ","cos(zenith)",0.7, 1.0);
  RooRealVar zenith("zenith","zenith",0.0,45.0);
  RooRealVar azimuth("azimuth","azimuth", 0.0, 360.0);
  RooRealVar Erec("Erec","Erec", std::pow(10, -1.5), 100.0);
  RooRealVar noise("noise","noise", 0.0, 20.0);
  RooRealVar Tels("Tels","Tels", 1.5, 4.5);
  RooRealVar Etrue("Etrue","E_{true} [TeV]", 0.0, 100.0);

  // Define the x and y variables to be used for storing data values                                                                 
  double x=0.0, y=0.0, fMSW=0.0, zen=0.0, az=0.0, fErec=0.0, fNoise=0.0, fTels=0.0, fEtrue=0.0;

  // Setup the argument set containing a list of our RooFit observables                                                              
  // (you can specify up to 9 RooAbsArg's in the RooArgSet constructor)                                                              
  //RooArgSet* VariableArgSet = new RooArgSet(*xCoord,*yCoord,*msw,*zenith,*azimuth,*noise,*Erec,*Tels, "argset") ;
  RooArgSet* TestArgSet = new RooArgSet(xCoord,yCoord,msw,zenith,azimuth,noise,Tels,"argset");

  // Loop over all the entries and import necessary ones into the dataset                                                            
  double fMSL=0.0, eRA=0.0, eDec=0.0;
  int TelsParticipating = 0;                                                                                                       
//VACoordinatePair evntPos(eRA, eDec, VACoordinates::J2000, VACoordinates::Deg) ;
//VACoordinatePair instantTrackPos(eRA, eDec, VACoordinates::J2000, VACoordinates::Deg) ;

  TTree *dataTree = new TTree("datatree","datatree") ;
  dataTree->Branch(xCoord.GetName(), &x) ;                                                                                             
  dataTree->Branch(yCoord.GetName(), &y) ;                                                                                             
  dataTree->Branch(msw.GetName(), &fMSW) ;
  dataTree->Branch(zenith.GetName(), &zen) ;
  dataTree->Branch(azimuth.GetName(), &az) ;
  dataTree->Branch(noise.GetName(), &fNoise) ;
  //dataTree->Branch(Erec.GetName(), &fErec) ;                                                                                            
  dataTree->Branch(Tels.GetName(), &fTels) ;

  // Load the run header and QStats for use in extracting the noise 
  VARootIO* file = new VARootIO(s5, true);
  const VAQStatsData* qstats = file->loadTheQStatsData();
  // Check if there is a simulation header. If so, then treat this as a simulation file.                                                  
  VASimulationData* simdata = new VASimulationData() ;
  VASimulationHeader* simHeader = file->loadTheSimulationHeader(false) ;
  if (simHeader != nullptr) {
    _isSim = true ;
  }

  TFile* f = TFile::Open(s5.c_str(),"READ");
  TChain *ch = new TChain( "SelectedEvents/CombinedEventsTree" );
  ch->AddFile( s5.c_str() );
  VASimulationData *sim = 0;
  VAShowerData *sh = 0;

    
  // If we're dealing with a simulation file we need to do some special things with Etrue
  if (_isSim) {
    dataTree->Branch(Etrue.GetName(), &fEtrue) ;       // add 'Etrue' branch to TTree
    TestArgSet->add(Etrue) ;                                // add 'Etrue' to argset
    //ch->SetBranchAddress("Sim", &sim) ;     // set the address of the simulation branch
    // allowing us to read in true energy values
  }
  else {
    dataTree->Branch(Erec.GetName(), &fErec); //add 'Erec' branch to TTree if data                                                    
    TestArgSet->add(Erec) ;                  // add 'Erec' to argset if data 
    //if(itm) {ch->SetBranchAddress( "M3D", &sh );}
    //ch->SetBranchAddress( "S", &sh );
  }
  
  // Store information from the run header:
  // _ATM, _epoch, _livetime, _trackPos, _sourcePos, _runNum,
  //saveObsInfo(file) ;

  TTreeReader myReader("SelectedEvents/CombinedEventsTree",f);
  TTreeReaderValue<VAShowerData> showDat(myReader, "S");
  TTreeReaderValue<VASimulationData>* simDat = nullptr ; 
  // If we're dealing with simulations, set the simulation data value
  if (_isSim) {
    simDat = new TTreeReaderValue<VASimulationData>(myReader, "Sim") ;
    ch->SetBranchAddress( "Sim", &sim );
  }

  // Create a wcs objec that represents the tracking position                                                                             
  //WorldCoor* wcsobj = initWCS(getTrackingRA_Deg(), getTrackingDec_Deg(), "J2000") ;

  // Extract the actual data values                                                                                                  
  int count = 0, events=0, norm=0, convOK=0 ;
  //std::bitset<4> tmp_tel_multiplicity(0) ;
  while (myReader.Next()) {
    // Skip events which are not passing some of the basic requirements                                                            
    // This is just in case the user hasnt removed cut events                                                                      
    if (!showDat->fIsDirection || !showDat->fIsCorePosition || !showDat->fIsReconstructed) continue ;

    events++;
    fMSW = showDat->fMSW;
    // Find total number of telescopes participating in the event and make sure it passes                                          
    // the telescope requirement cut                                                                                               
    fTels = showDat->fTelUsedInReconstruction.at(0) +
    showDat->fTelUsedInReconstruction.at(1) +
    showDat->fTelUsedInReconstruction.at(2) +
    showDat->fTelUsedInReconstruction.at(3);
    count++;

    WorldCoor* wcsobj = initWCS((showDat->fArrayTrackingRA_J2000_Rad)*r2d, (showDat->fArrayTrackingDec_J2000_Rad)*r2d, "J2000");
    
    if (_isSim) {
      // Note this algorithm is necessary when doing simulations
      // otherwise the computed offset is always 0.5 degrees.
      Coordinate2Projected(showDat->fDirectionRA_J2000_Rad*r2d,
			   showDat->fDirectionDec_J2000_Rad*r2d,
			   showDat->fArrayTrackingRA_J2000_Rad*r2d,
			   showDat->fArrayTrackingDec_J2000_Rad*r2d,
			   &x, &y, "J2000") ;
    } else {
      // This is the correct algorithm for projecting the events based
      // on the central tracking position.
      Coordinate2Projected_Quick(showDat->fDirectionRA_J2000_Rad*r2d,
				 showDat->fDirectionDec_J2000_Rad*r2d,
				 wcsobj, &x, &y) ;
    }
    // Check that the event is inside the cut radius
    if ((x*x+y*y) >= 4.0) continue ;
    
    // Extract the event zenith angle                                                                                              
    zen = 90.0 - ( showDat->fDirectionElevation_Rad * r2d ) ;
    az = showDat->fDirectionAzimuth_Rad * r2d ;
    // Find the camera averaged noise for this event                                                                               
    fNoise = findAverageNoise(&(*showDat), qstats, showDat->fTime) ;
    // Check the reconstructed or true energy (in TeV) is within our limits                                                               
    // If we're reading in a simulation file get the true energy value (in TeV)                                                           
    if (_isSim && ((*simDat)->fEnergyGeV*0.001)<100) fEtrue = (*simDat)->fEnergyGeV*0.001 ;
    else { fErec = showDat->fEnergy_GeV*0.001; }
    // Fill the TTree with the appropriate values 
    dataTree->Fill();
    //tmp_tel_multiplicity.set( std::floor(fTels-0.5) ) ;
    norm++;

    wcsfree(wcsobj);
  }//finish looping over events
  std::cout << "Number of events looped over in first file: " << events << std::endl;
  //f->Close();
  //file->closeTheRootFile();
  std::cout << "closing root files" << std::endl;
  std::cout << "Number of events recorded in TTree: " << dataTree->GetEntries() << std::endl;
  delete f;
  delete file;
  delete simDat;

  // Print error message if no data is acceptable                                                                                         
  if (dataTree->GetEntries() == 0) {
    std::cout << "Number of events recorded in TTree: " << dataTree->GetEntries() << std::endl;
    delete dataTree;
    return 0;
  } else {
    TFile* rfoutRoo = new TFile("demoRoo.root", "recreate");
    // Fill a RooDataSet from the TTree  
    RooDataSet* dataout = new RooDataSet("dataout","dataout",dataTree,*TestArgSet);
    std::cout << "Number of events recorded in TTree: " << dataTree->GetEntries() << std::endl;
    dataTree->Write(); //write out data tree to root file

    //DIAGNOSTICS                                                                                                                         
    std::cout << "Read out RooDataSet before writing out to root file: " << std::endl;
    dataout->Print("v");
    std::cout << " " << std::endl;

    //TFile* rfoutRoo = new TFile("demoRoo.root", "recreate");
    dataout->Write(); //write out RooDataSet                                                                                            
    rfoutRoo->Close();
    /*
    TFile* rfinRoo = new TFile("demoRoo.root", "read");
    RooDataSet* d = (RooDataSet*)rfinRoo->Get("dataout");
    //DIAGNOSTICS                                                                                                                         
    std::cout << "Read in and print RooDataSet from root file: " << std::endl;
    d->Print("v");
    const RooArgSet* row = d->get();
    row->Print("v");
    row = d->get(90);
    row->Print("v");
    RooRealVar* mswrow = (RooRealVar*)row->find("msw");
    std::cout << mswrow->getVal() << std::endl;
    rfinRoo->Close();
    */
    delete dataTree;
    delete dataout;
    //delete d;
    delete rfoutRoo;
    //delete rfinRoo;
    return 1;
  }


}

double findAverageNoise(VAShowerData* showerdata, const VAQStatsData* pQStatsData, const VATime& time)
{
  // Function to read from the qstats what the noise is at a                                                                             
  // given time (based on VAStage6Analysis::getAverageNoise)                                                                             

  double meanNoise = 0.0 ;          // mean noise                                                                                        
  int nCameras = 0 ;                // number of cameras in the event                                                                    
  uint32_t fWindowSizeForNoise = 7 ;  // window size                                                                                     
  // (currently fixed at 7 samples)                                                                  

  // If the window size is 0 set it to 7 which is the default                                                                            
  if(fWindowSizeForNoise == 0) {
    cout << "MLM Analysis picked up a window size of 0, which likely means no EA is available.  Reverting to default size 7 for noise \
calc in EventStatsTree.\n";
    fWindowSizeForNoise = 7;
  }

  // Calculate the noise at each telescope which is participating in this event                                                          
  double noise=0.0 ;
  for(int telID = 0; telID < kMaxTels; telID++) {
    if((showerdata->fTelUsedInReconstruction.at(telID)) == 1) {
      //           std::cout << "   Tel " << telID << ": " << showerdata->fTelUsedInReconstruction.at(telID) << std::endl;                      
      // extract the noise for this telescope at this time                                                                           
      noise = pQStatsData->getCameraAverageTraceVarTimeIndpt(telID, fWindowSizeForNoise) ;

      // if the noise is greater than 0 add it and increment the number                                                              
      // of telescopes by 1                                                                                                          
      if(noise > 0.0) {
	meanNoise += noise;
	nCameras++;
      }
      //            noise = 0.0 ;                                                                                                                
    }
  }// end for loop over telID                                                                                                            

  // Divide the mean noise by the number of cameras                                                                                      
  // to get the average noise for this event                                                                                             
  if(nCameras > 0) {
    meanNoise /= nCameras;
  }

  return meanNoise ;

}

// Method for converting from world coordinates to gnomonic projected coordinates                                                                                               
// Make sure xcoord, ycoord, xTangentPoint, and yTangentPoint are in DEGREES coordinates                                                                                        
/*
void Coordinate2Projected_Quick(double xcoord, double ycoord,
				WorldCoor* wcs,
				double *xproj, double *yproj)
{
  // Fill the actual objects                                                                                                                                                   
  tnxpix(xcoord, ycoord, wcs, xproj, yproj) ;

  // Returned projected coordinates are in radians, so we must convert them to degrees                                                                                         
  *xproj *= (180.0 / TMath::Pi()) ;
  *yproj *= (180.0 / TMath::Pi()) ;

  return ;
}
*/
