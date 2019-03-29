 
//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 27 15:49:34 2018 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/GluGluHToMuMu_M-125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/AA4BF847-3785-E811-B130-001E67DFFF5F.root
//////////////////////////////////////////////////////////

#ifndef HmmAnalyzer_h
#define HmmAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include "TRandom.h"
#include "MainEvent.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
#include "RoccoR.h"
#include "LeptonEfficiencyCorrector.h"
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"
#ifdef __CINT__

#pragma link C++ class vector<float>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<bool>+;
#endif
// Header file for the classes stored in the TTree if any.

class HmmAnalyzer : public MainEvent {
 public :
   HmmAnalyzer(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data",const char *isData="F");
   virtual ~HmmAnalyzer();
   void Analyze(bool isData, int option, string outputFileName, string label);
   TFile *oFile;
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *, const char *,double);
  //declare any specific function required
  TH1D *h_Zguess_pt = new TH1D("Zguess_pt","Z guess pt",100,0.,400.);
  TH1D *h_Zreal_pt = new TH1D("Zreal_pt","Z real pt",100,0.,400.);
  TH1D *h_Zguess_phi = new TH1D("Zguess_phi","Z guess #phi",50,-4.,4.);
  TH1D *h_Zreal_phi = new TH1D("Zreal_phi","Z real #phi",50,-4.,4.);
  TH1D *h_Zguess_eta = new TH1D("Zguess_eta","Z guess #eta",100,-3.,3.);
  TH1D *h_Zreal_eta = new TH1D("Zreal_eta","Z real #eta",100,-3.,3.);
};

#endif

#ifdef HmmAnalyzer_cxx
HmmAnalyzer::HmmAnalyzer(const TString &inputFileList, const char *outFileName, const char* dataset, const char *isData) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.


  TChain *tree = new TChain("Events");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
  
  MainEvent::Init(tree);
  oFile = new TFile(outFileName, "recreate");
  //TString histname(outFileName);
  //ohistFile = new TFile("hist_"+histname, "recreate");
  //BookTreeBranches();
}


bool HmmAnalyzer::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
}

 HmmAnalyzer::~HmmAnalyzer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   oFile->cd();
   h_Zguess_pt->Write();
   h_Zreal_pt->Write();
   h_Zguess_phi->Write();
   h_Zreal_phi->Write();
   h_Zguess_eta->Write();
   h_Zreal_eta->Write();
   oFile->Write();
   oFile->Close();
}

Long64_t HmmAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
#endif // #ifdef HmmAnalyzer_cxx
