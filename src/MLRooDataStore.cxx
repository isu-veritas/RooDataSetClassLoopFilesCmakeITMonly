// C++ HEADERS
#include <ctime>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
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
#include <TH2F.h>

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
#include "VATime.h"
#include "VASlalib.h"
#include "wcs.h"
#include "MLWCSHandler.h"
#include "MLRooDataStore.h"


// Default constructor, this doesnt do much of anything
MLRooDataStore::MLRooDataStore() : _sourceCoordPair(nullptr),_trackCoordPair(nullptr) {
  // initialize the source and tracking positions to some generic location                                                           
  _trackCoordPair  = new VACoordinatePair(0.0, 0.0, VACoordinates::J2000, VACoordinates::Deg) ;
  _sourceCoordPair = new VACoordinatePair(0.0, 0.0, VACoordinates::J2000, VACoordinates::Deg) ;
}
// Destructor
MLRooDataStore::~MLRooDataStore()
{
    // What do I need to delete so that memory management
    // doesnt become a problem?
    if (_sourceCoordPair != nullptr) delete _sourceCoordPair ;
    if (_trackCoordPair != nullptr) delete _trackCoordPair ;
}


//bool MLRooDataStore::fillFromFile(std::string s5, int nn) {
bool MLRooDataStore::fillFromFile(std::string s5, std::string rooname) {

  _isSim = false;
  std::cout << "Output root file: " << rooname << std::endl;

  // TODO : change these variables to RooRealProxy (Josh's comment)                                                                   
  // Create some variables needed by the datasets we're going to be creating.                                                   
  RooRealVar _xCoord("xpos","xpos",-2.0, 2.0);  xCoord = _xCoord;// Exagerate the range since we want to contain all events           
  RooRealVar _yCoord("ypos","ypos",-2.0, 2.0);  yCoord = _yCoord;// Exagerate the range since we want to contain all events           
  RooRealVar _msw("msw","msw",0.0, 10.0);  msw= _msw;
  //RooRealVar _zenith("zenith","zenith",0.0,45.0);
  RooRealVar _cosZ("cosZ","cos(zenith)", 0.7, 1.0);  cosZ = _cosZ;
  RooRealVar _azimuth("azimuth","azimuth", 0.0, 360.0);  azimuth = _azimuth;
  RooRealVar _Erec("Erec","Erec", std::pow(10, -1.5), 100.0);  Erec = _Erec;
  RooRealVar _noise("noise","noise", 0.0, 20.0);  noise = _noise;
  RooRealVar _Tels("Tels","Tels", 1.5, 4.5);  Tels = _Tels;  Tels = _Tels;
  //RooRealVar _Etrue("Etrue","E_{true} [TeV]", 0.0, 100.0);  Etrue = _Etrue;
  RooRealVar _Etrue("Etrue","E_{true} [TeV]", -1.5, 2.0);  Etrue = _Etrue; //Log(E_TeV) 30 GeV - 100 TeV
  TH2F *xvsy = new TH2F("xpos vs ypos", "xpos vs ypos",100,-1,1,100,-1,1); //all events
  // Define the x and y variables to be used for storing data values                                                                 
  double x=0.0, y=0.0, fMSW=0.0, zen=0.0, az=0.0, cosZen=0.0, fErec=0.0, fNoise=0.0, fTels=0.0, fEtrue=0.0;

  // Setup the argument set containing a list of our RooFit observables                                                              
  // (you can specify up to 9 RooAbsArg's in the RooArgSet constructor)                                                              
  //RooArgSet* VariableArgSet = new RooArgSet(*xCoord,*yCoord,*msw,*zenith,*azimuth,*noise,*Erec,*Tels, "argset") ;
  //RooArgSet* TestArgSet = new RooArgSet(xCoord,yCoord,msw,zenith,azimuth,noise,Tels,"argset");
  RooArgSet* TestArgSet = new RooArgSet(xCoord,yCoord,msw,cosZ,azimuth,noise,Tels,"argset"); 

  // Loop over all the entries and import necessary ones into the dataset                                                            
  double fMSL=0.0, eRA=0.0, eDec=0.0,  fTheta2_Deg2;
  double storeoff{1.2f}; storeoff=0.0;
  float OffsetBin[9] = {0.00,0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00};
  int TelsParticipating = 0;                                                                                                       
//VACoordinatePair evntPos(eRA, eDec, VACoordinates::J2000, VACoordinates::Deg) ;
//VACoordinatePair instantTrackPos(eRA, eDec, VACoordinates::J2000, VACoordinates::Deg) ;

  TTree *dataTree = new TTree("datatree","datatree") ;
  dataTree->Branch(xCoord.GetName(), &x) ;                                                                                            
  dataTree->Branch(yCoord.GetName(), &y) ;                                                                                            
  dataTree->Branch(msw.GetName(), &fMSW) ;
  //dataTree->Branch(zenith.GetName(), &zen);
  dataTree->Branch(cosZ.GetName(), &cosZen) ;
  dataTree->Branch(azimuth.GetName(), &az) ;
  dataTree->Branch(noise.GetName(), &fNoise) ;
  dataTree->Branch(Tels.GetName(), &fTels) ;

  std::vector<double> ObsInfoVect;
  std::vector<std::string> NameVect;
  std::vector<bool> Sim_Data;

  // Load the run header and QStats for use in extracting the noise 
  VARootIO* file = new VARootIO(s5, true);
  const VAQStatsData* qstats = file->loadTheQStatsData();

  // Check if there is a simulation header. If so, then treat this as a simulation file.                                                  
  VASimulationData* simdata = new VASimulationData() ;
  VASimulationHeader* simHeader = file->loadTheSimulationHeader(false) ;
  if (simHeader != nullptr) {
    _isSim = true ;
  }

  std::cout << "is simulation? " << _isSim << std::endl;

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
  }
  
  TTreeReader myReader("SelectedEvents/CombinedEventsTree",f);
  TTreeReaderValue<VAShowerData> showDat(myReader, "S");
  TTreeReaderValue<VAModel3DData> showDatM3D(myReader, "M3D");
  //ch->SetBranchAddress( "M3D", &sh );                                                                                        
  //ch->SetBranchAddress( "S", &sh ); 
  TTreeReaderValue<VASimulationData>* simDat = nullptr ; 
  // If we're dealing with simulations, set the simulation data value
  if (_isSim) {
    simDat = new TTreeReaderValue<VASimulationData>(myReader, "Sim") ;
    ch->SetBranchAddress( "Sim", &sim );
  }


  // Extract the actual data values                                                                                                  
  int count = 0, events=0, norm=0, convOK=0 ;
  //std::bitset<4> tmp_tel_multiplicity(0) ;
  while (myReader.Next()) {
    // Skip events which are not passing some of the basic requirements                                                            
    // This is just in case the user hasnt removed cut events                                                                      
    if (!showDatM3D->fIsDirection || !showDatM3D->fIsCorePosition || !showDatM3D->fIsReconstructed || !showDatM3D->fGoodness3D) continue;

    events++;
    storeoff=storeoff+((showDatM3D->fArrayTrackingDec_J2000_Rad)*r2d);
    fMSW = showDatM3D->fMSW;
    // Find total number of telescopes participating in the event and make sure it passes                                          
    // the telescope requirement cut                                                                                               
    fTels = showDat->fTelUsedInReconstruction.at(0) +
    showDat->fTelUsedInReconstruction.at(1) +
    showDat->fTelUsedInReconstruction.at(2) +
    showDat->fTelUsedInReconstruction.at(3);
    count++;

    WorldCoor* wcsobj = initWCS((showDatM3D->fArrayTrackingRA_J2000_Rad)*r2d, (showDatM3D->fArrayTrackingDec_J2000_Rad)*r2d, "J2000");
    
    if (_isSim) {
      // Note this algorithm is necessary when doing simulations
      // otherwise the computed offset is always 0.5 degrees.
      
      Coordinate2Projected(showDatM3D->fDirectionRA_J2000_Rad*r2d,
                           showDatM3D->fDirectionDec_J2000_Rad*r2d,
                           _srcRA,
                           _srcDec,
                           &x, &y, "J2000") ;
    } else {
      
      /*
      Coordinate2Projected(showDat->fDirectionRA_J2000_Rad*r2d,
			   showDat->fDirectionDec_J2000_Rad*r2d,
			   showDat->fArrayTrackingRA_J2000_Rad*r2d,
			   showDat->fArrayTrackingDec_J2000_Rad*r2d,
			   &x, &y, "J2000") ;
      
    } else {
      */
      // This is the correct algorithm for projecting the events based
      // on the central tracking position.
      Coordinate2Projected_Quick(showDatM3D->fDirectionRA_J2000_Rad*r2d,
				 showDatM3D->fDirectionDec_J2000_Rad*r2d,
				 wcsobj, &x, &y) ;
    }
    // Check that the event is inside the cut radius
    if ((x*x+y*y) >= 4.0) continue;
    //xvsy->Fill(x,y);
    
    // Extract the event zenith angle                                                                                              
    zen = 90.0 - ( showDatM3D->fDirectionElevation_Rad * r2d ) ;
    cosZen = std::cos(d2r*zen);
    az = showDatM3D->fDirectionAzimuth_Rad * r2d ;
    // Find the camera averaged noise for this event                                                                               
    fNoise = findAverageNoise(&(*showDat), qstats, showDat->fTime) ;
    // Check the reconstructed or true energy (in TeV) is within our limits                                                               
    // If we're reading in a simulation file get the true energy value (in TeV)                                                           
    //if (_isSim && ((*simDat)->fEnergyGeV*0.001)<100) fEtrue = ((*simDat)->fEnergyGeV*0.001) ;
    if (_isSim && ((*simDat)->fEnergyGeV*0.001)<100) fEtrue = TMath::Log10((*simDat)->fEnergyGeV*0.001) ;
    else { fErec = showDatM3D->fEnergy_GeV*0.001; }

    if((fEtrue<0.8)&&(fEtrue>0.64)) {xvsy->Fill(x,y);}

    // Fill the TTree with the appropriate values 
    dataTree->Fill();
    //tmp_tel_multiplicity.set( std::floor(fTels-0.5) ) ;
    norm++;
    wcsfree(wcsobj);
  }//finish looping over events  

  std::cout << "Number of events looped over in first file: " << events << std::endl;
  storeoff=storeoff/events;
  std::cout << "offset derived from showDat->fArrayTrackingDec_J2000_Rad)*r2d : " << storeoff << std::endl;
  f->Close();
  delete f;

  // Print error message if no data is acceptable                                                                                        
  if (dataTree->GetEntries() == 0) {
    std::cout << "Number of events recorded in TTree: " << dataTree->GetEntries() << std::endl;
    delete dataTree;
    return 0;
  } else {

    TFile* rfoutRoo = new TFile(rooname.c_str(), "recreate");
    // Fill a RooDataSet from the TTree                                                                                                     
    //RooDataSet* dataout = new RooDataSet("dataout","dataout",dataTree,*TestArgSet);                                                       
    dataout = new RooDataSet("dataout","dataout",dataTree,*TestArgSet);
    std::cout << "Number of events recorded in TTree: " << dataTree->GetEntries() << std::endl;
    dataTree->Write(); //write out data tree to root file                                                                                   
    dataout->Write(); //write out RooDataSet  
    xvsy->Write();

    _xmean = (dataout->mean(xCoord)) ; //unit degrees
    _ymean = (dataout->mean(yCoord)) ; //unite degrees
    std::cout << "xmean: " << _xmean << " ymean: " << _ymean << std::endl;

    // Now fill the offset value in the correct way if we're using simulations                                                              
    if (_isSim) {
      // The following uses the mean of the projected coordinates, which may not be 100%                                                    
      // accurate, but it is most likely within the errors of our pointing accuracy.                                                        
      _offset = std::sqrt((_xmean*_xmean)+(_ymean*_ymean)) ;
    } else {
      // Store the source offset in degrees                                                                                                 
      _offset = _sourceCoordPair->angularSeparation_Deg(*_trackCoordPair) ;
    }

    std::cout << "Filling ObsInfoVect" << std::endl;
    //START FILLING IN ObsInfoVect HERE
    //Set Observation Info and Source Name variable                                                                                        
    saveObsInfo(file);
    //push information to ObsInfoVector                                                                                                     
    ObsInfoVect.push_back(_epoch);
    ObsInfoVect.push_back(_livetime);
    ObsInfoVect.push_back(_ATM);
    ObsInfoVect.push_back(_trackRA);//tracking RA                                                                                           
    ObsInfoVect.push_back(_trackDec);//tracking Dec                                                                                         
    ObsInfoVect.push_back(_srcRA);//source RA                                                                                               
    ObsInfoVect.push_back(_srcDec);//source Dec                                                                                             
    NameVect.push_back(_SourceName);
    Sim_Data.push_back(_isSim);

    //_avgZenith=(dataout->mean(zenith));
    _avgZenith=( std::acos(dataout->mean(cosZ))*r2d );
    ObsInfoVect.push_back(_avgZenith);//avg zenith 
    //_avgAzimuth=(dataout->mean(azimuth));
    _avgAzimuth=getAvgAzimuth(dataout);
    _avgNoise=(dataout->mean(noise));

    ObsInfoVect.push_back(_avgAzimuth);//avg azimuth
    ObsInfoVect.push_back(_avgNoise);//avg noise
    //offset
    if (_isSim) { 
      //_offset=storeoff;
      for(int n=0; n<9; n++) {
	std::cout << "_offset: " << storeoff << " OffsetBin[n]: " << OffsetBin[n] << std::endl;
	//if(OffsetBin[n]==_offset) {rec_offset=Offsetstring[n]; continue;}
	if(abs(OffsetBin[n]-storeoff)<0.1) {_offset=storeoff; continue;}
      }
      ObsInfoVect.push_back(_offset); 
    }
    else {ObsInfoVect.push_back(_offset);}
    std::cout << "Size of ObsInfo Vector: " << ObsInfoVect.size() << std::endl;


    PrintInfo();
    //rfoutRoo->WriteObject(&_SourceName, "SourceName"); //write out source name
    rfoutRoo->WriteObject(&ObsInfoVect, "ObsInfoVect"); //write out ObsInfoVect
    rfoutRoo->WriteObject(&NameVect, "NameVect"); //write out strings
    rfoutRoo->WriteObject(&Sim_Data, "Sim_Data"); //write out boolean, defining sim or data
    rfoutRoo->Close();
    
    //PrintInput(rooname);

    //rfinRoo->GetListOfKeys()->Print();

    //std::cout << "Reading contents of ObsInfo Vector" << std::endl;
    //std::cout << "Reading Size of ObsInfo Vector: " << testvect->size() << std::endl;
    //std::cout << "Reading Size of Names Vector: " << testnamevect->size() << std::endl;
    //for(int aa=0; aa<testvect->size();aa++) {
    //std::cout << testvect->at(aa) << std::endl;
    //}
    //std::cout << "Name of source (after input): " << testnamevect->at(0) << std::endl;

    //const RooArgSet* row = d->get();
    //row->Print("v");
    //row = d->get(90);
    //row->Print("v");
    //RooRealVar* mswrow = (RooRealVar*)row->find("msw");
    //std::cout << mswrow->getVal() << std::endl;
    //rfinRoo->Close();
    file->closeTheRootFile();                                                                                                             
    std::cout << "closing root files" << std::endl; 
    delete file;
    delete simDat;
    delete dataTree;
    delete dataout;
    delete rfoutRoo;
    return 1;
  }


}


double MLRooDataStore::findAverageNoise(VAShowerData* showerdata, const VAQStatsData* pQStatsData, const VATime& time)
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

void MLRooDataStore::saveObsInfo(VARootIO* file)
{
  // Get the run header
  VARunHeader* RunHeader = file->loadTheRunHeader(true) ;
  
  // Get some information from the run header on source and offset positions
  // Data values are in radians
  double offsetRA  = RunHeader->fRunInfo.fOffsetRA ;
  double offsetDec = RunHeader->fRunInfo.fOffsetDec ;
  double srcRA     = RunHeader->fSourceInfo.fRA ;
  double srcDec    = RunHeader->fSourceInfo.fDec ;
  std::cout << "RunHeader->fRunInfo->fOffsetRA " << (RunHeader->fRunInfo.fOffsetRA) << std::endl; 
  std::cout << "RunHeader->fRunInfo->fOffsetDec " << (RunHeader->fRunInfo.fOffsetDec) << std::endl;
  std::cout << "RunHeader->fSourceInfo->fRA " << (RunHeader->fSourceInfo.fRA) << std::endl; 
  std::cout << "RunHeader->fSourceInfo->fDec " << (RunHeader->fSourceInfo.fDec) << std::endl;
  std::cout << "From run header:" << std::endl;
  std::cout << "offsetRA: " << offsetRA << std::endl;
  std::cout << "offsetDec: " << offsetDec << std::endl;
  std::cout << "srcRA: " << srcRA << std::endl;
  std::cout << "srcDec: " << srcDec << std::endl;
  
  std::vector<double> _sourcePos, _trackPos;
  //_sourcePos(std::vector<double>(2,0.0)); _trackPos(std::vector<double>(2,0.0));
  // Set information on the source and tracking position
  if (_isSim) {
    offsetRA  *= d2r ;
    offsetDec *= d2r ;
    srcRA     *= d2r ;
    srcDec    *= d2r ;
    _sourcePos = {srcRA, srcDec} ;
    _trackPos = {(_xmean+(_sourcePos[0]*r2d)), (_ymean+(_sourcePos[1]*r2d))} ;
    _sourceCoordPair->setCoordinates_Deg(_sourcePos[0]*r2d, _sourcePos[1]*r2d, VACoordinates::J2000 );
    _trackCoordPair->setCoordinates_Deg(_trackPos[0], _trackPos[1], VACoordinates::J2000) ;
    
    _srcRA=_sourcePos[0]*r2d;
    _srcDec=_sourcePos[1]*r2d;
    _trackRA=(_trackPos[0]);
    _trackDec=(_trackPos[1]);

  } 
  if (!_isSim) {
    _sourcePos = {srcRA, srcDec} ;
    _trackPos = {(srcRA+offsetRA), (srcDec+offsetDec)} ;
    _sourceCoordPair->setCoordinates_Deg(_sourcePos[0]*r2d, _sourcePos[1]*r2d, VACoordinates::J2000 );
    _trackCoordPair->setCoordinates_Deg((_trackPos[0]*r2d), (_trackPos[1]*r2d), VACoordinates::J2000) ;

    _srcRA=_sourcePos[0]*r2d;
    _srcDec=_sourcePos[1]*r2d;
    _trackRA=(_trackPos[0]*r2d);
    _trackDec=(_trackPos[1]*r2d);
  }

    // Store the run number
  _runNum = RunHeader->pfRunDetails->fRunNum;
    // Store the livetime
  _livetime = RunHeader->pfRunDetails->fRunCutLiveTimeSeconds;
  //std::cout << "livetime(sec): " << _livetime << std::endl;
    //_trackRA = (srcRA+offsetRA);
    //_trackDec = (srcDec+offsetDec);
    // Store the epoch
  _epoch = findEpoch( RunHeader->getStartTime() );
  //std::cout << "epoch: " << _epoch << std::endl;
  // Store the atmosphere
  if (!_isSim) {
    _ATM = findATM_Data( RunHeader->getStartTime() ) ;
    _SourceName = RunHeader->getSourceId();
  } 
  else {
    VASimulationHeader* simHeader = file->loadTheSimulationHeader(true);
    _ATM = findATM_Sims( simHeader ) ;
    _SourceName = "Simulation";
    delete simHeader;
  }
  //std::cout << "atmosphere: " << _ATM << std::endl;
  //std::cout << "Source Name: " << _SourceName << std::endl;
  delete RunHeader;
  return;
}

//_____________________________________________________________
int MLRooDataStore::findEpoch(VATime time)
{
    // Since this is based on the start time, it should be valid for both
    // actual data and simulations
    
    // Get year (y), month (m), and day (d)
    uint32_t y,m,d ;
    time.getCalendarDate(y, m, d) ;

    // This line constructs an integer of the form YYYYMMDD
    // I do this so I only have to make 1 integer comparison per epoch check
    int Date = 10000*y + 100*m + d ;
    
    // Get the Epoch based on the year and month
    // Pre-full-array epoch
    if (Date < 20070427) return 3 ;
    // V4 (old array)
    if (Date < 20090800) return 4 ;
    // V5 (new array)
    if (Date < 20120800) return 5 ;
    // V6 (upgrade array)
    if (Date < 20210800) return 6 ;
    
    // If we ever move to a V7 (8,9,etc...) array configuration, you should
    // modify the if-statement for V6 and add in an if-statement for V7 (8,9,etc...)
    
    // Should throw an error at this point
    return 0 ;
}

//_____________________________________________________________
int MLRooDataStore::findATM_Data(VATime time)
{
    // The following only works for real data.
    // For how to determine the atmosphere of a simulation
    // file, see VADataSet::findATM_Sims()
    
    // Get year (y), month (m), and day (d)
    uint32_t y,m,d ;
    time.getCalendarDate(y, m, d) ;
    
    // Get the atmosphere based on the year and month
    // return summer (22) if month is between May-October
    if ((m >= 5) && (m <= 10)) {
        return 22 ; // Summer
    } else {
        return 21 ; // Winter
    }
    // Should throw an error at this point
    return 0 ;
}


//_____________________________________________________________
int MLRooDataStore::findATM_Sims(VASimulationHeader* simHeader)
{
    // The following only works for simulation.
    // For how to determine the atmosphere of a real data
    // file, see VADataSet::findATM_Data()
    
    // In simulation files the "start date" is actually only
    // guaranteed to be from the same epoch, but not necessarily
    // the same season (i.e. atmosphere).
    
    // Setup a temporary variable for holding the atmosphere value we find.
    int tmpATM = 0 ;
    
    // We're going to have some fun with regular expressions here!
    // This is the string we want to parse
    std::string SimConfigStr = simHeader->fSimConfigFile ;
    
    // START: C++11 SOLUTION =================================================
    // Setup the regular expression that will find the line in the
    // string that we are looking for.
    std::regex rgxATM("(ATM)[[:blank:]]+(table)[[:blank:]]+([0-9]+)") ;
    std::regex rgxNMBR("[0-9]+") ;
    
    // If we find a match, we'll store it in this variable
    std::smatch matchstr ;
    
    // Now Check for a match to rgxATM
    if ( std::regex_search(SimConfigStr, matchstr, rgxATM) ) {
        // If we've found a match, parse the string for the date
        std::string match(matchstr.str()) ;
        if ( std::regex_search(match, matchstr, rgxNMBR) ){
            // Convert the number string into an integer
            tmpATM = std::stoi( matchstr.str() ) ;
        }
    }
    // END: C++11 SOLUTION ==================================================
 
    // Since there's a chance that we've used a different atmosphere,
    // (either updated or modified) I'll handle those special cases here
    if (tmpATM == 61) tmpATM = 21 ;  // special case of updated winter atmosphere
    if (tmpATM == 62) tmpATM = 22 ;  // special case of updated summer atmosphere
    
    // Return the found ATM value
    return tmpATM ;
}

void MLRooDataStore::PrintInfo()
{  
  // Print out some useful information about this run                                                                                 
  std::cout << "Print RooDataSet that has been output to root file: " << std::endl;
  dataout->Print("v");
  std::cout << "----------------------------------------------------" << std::endl;
  std::cout << " Printing out event params to be stored in root file " << std::endl;
  std::cout << "----------------------------------------------------" << std::endl ;
  std::cout << "INFORMATION FOR RUN    : " << _runNum << std::endl ;
  std::cout << "     Saved Name        : " << _SourceName << std::endl ;
  std::cout << "     Is Simulation     : " << _isSim  << std::endl ;
  std::cout << "     Livetime (Sec)    : " << _livetime << std::endl ;
  std::cout << "     Tracking RA       : " << _trackRA << std::endl ;
  std::cout << "     Tracking Dec      : " << _trackDec << std::endl ;
  std::cout << "     Source RA         : " << _srcRA << std::endl;
  std::cout << "     Source Dec        : " << _srcDec << std::endl;
  std::cout << "     Source Offset     : " << _offset << std::endl;
  //std::cout << "     Energy Range [TeV]: " << EnergyLowerBound << "-" << EnergyUpperBound << std::endl;
  std::cout << "     Events Saved      : " << dataout->numEntries() << std::endl ;
  std::cout << "     Epoch             : " << _epoch << std::endl ;
  std::cout << "     Season (ATM)      : " << _ATM << std::endl ;
  std::cout << "     Avg. Zenith       : " << _avgZenith << std::endl ;
  std::cout << "     Avg. Azimuth      : " << _avgAzimuth << std::endl ;
  std::cout << "     Avg. Noise        : " << _avgNoise << std::endl ;
  std::cout << "     Offset            : " << _offset << std::endl ;
  std::cout << "----------------------------------------------------" << std::endl ;

  return ;
}

//_____________________________________________________________
double MLRooDataStore::getAvgAzimuth(RooDataSet* ds)
{
  // Return value                                                                                                                                                   
  double avgAzimuth(0.0) ;

  // First get the range of variables in this dataset                                                                                                               
  double az_lowest(360.0), az_highest(0.0) ;
  ds->getRange(azimuth, az_lowest, az_highest) ;

  // If the highest and lowest azimuth values are more than 180 degrees apart,                                                                                      
  // we need to take special care when computing the average                                                                                                        
  if ((az_highest-az_lowest)>180.0) {
    // Sum all of the azimuth values in the dataset                                                                                                               
    double az_sum(0.0) ;
    for (int i=0; i<ds->numEntries(); i++) {
      // Get the current azimuth value                                                                                                                          
      double az_val = ds->get(i)->getRealValue(azimuth.GetName()) ;
      // Make sure the azimuth value is in the range -180 -> 180.                                                                                               
      if (az_val>180.0) az_val-=360.0 ;
      az_sum += az_val ;
    }
    // Compute the average value                                                                                                                                  
    avgAzimuth = az_sum / ds->numEntries() ;
    // Scale the value so that it's positive                                                                                                                      
    if (avgAzimuth < 0.0) avgAzimuth+=360.0 ;
  } else {
    // otherwise we can just compute the mean                                                                                                                     
    avgAzimuth = ds->mean(azimuth);
  }
  return avgAzimuth ;
}

void MLRooDataStore::PrintInput(std::string input_name) {

  TFile* rfinRoo = new TFile(input_name.c_str(), "read");
  RooDataSet* d = (RooDataSet*)rfinRoo->Get("dataout");
  std::vector<double> *testvect = (std::vector<double>*)rfinRoo->Get("ObsInfoVect");
  std::vector<std::string> *testnamevect = (std::vector<std::string>*)rfinRoo->Get("NameVect");

  //std::cout << "Reading contents of ObsInfo Vector" << std::endl;
  //std::cout << "Reading Size of ObsInfo Vector: " << testvect->size() << std::endl;
  //std::cout << "Reading Size of Names Vector: " << testnamevect->size() << std::endl;
  //for(int aa=0; aa<testvect->size();aa++) {
  //std::cout << testvect->at(aa) << std::endl;
  //}
  // Print out some useful information about this run                         
  std::cout << "Read in and print RooDataSet from root file: " << std::endl;
  d->Print("v");
  std::cout << "----------------------------------------------------" << std::endl;
  std::cout << " Printing out run params stored in root file " << std::endl;
  std::cout << "----------------------------------------------------" << std::endl ;
  std::cout << "     Saved Name        : " << testnamevect->at(0) << std::endl ;
  std::cout << "     Livetime (Sec)    : " << testvect->at(1) << std::endl ;
  std::cout << "     Tracking RA       : " << testvect->at(3) << std::endl ;
  std::cout << "     Tracking Dec      : " << testvect->at(4) << std::endl ;
  std::cout << "     Source RA         : " << testvect->at(5) << std::endl;
  std::cout << "     Source Dec        : " << testvect->at(6) << std::endl;
  std::cout << "     Events Saved      : " << d->numEntries() << std::endl ;
  std::cout << "     Epoch             : " << testvect->at(0) << std::endl ;
  std::cout << "     Season (ATM)      : " << testvect->at(2) << std::endl ;
  std::cout << "     Avg. Zenith       : " << testvect->at(7) << std::endl ;
  std::cout << "     Avg. Azimuth      : " << testvect->at(8) << std::endl ;
  std::cout << "     Avg. Noise        : " << testvect->at(9) << std::endl ;
  std::cout << "----------------------------------------------------" << std::endl ;

  rfinRoo->Close();

  delete testvect;
  delete testnamevect;
  delete d;
  delete rfinRoo;
  return;
}
