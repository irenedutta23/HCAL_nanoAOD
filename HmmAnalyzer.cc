#define HmmAnalyzer_cxx
#include "HmmAnalyzer.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>

#ifdef MAKECINT
#pragma link C++ class vector<float>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<int>+;
#endif
#ifdef MAKECINT
#pragma link C++ class vector<bool>+;
#endif

int main(int argc, char* argv[])
{

  if(argc < 4) {
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
  
  map<TString,double> proc_scale;
  double lumi = 59.74*1000.;
  proc_scale["DYJetsToLL"]= 6225.42*lumi/(3258462550016+492179082112.0);
  proc_scale["ttTo2l2v"]=85.656*lumi/(623402174.0+4782395097.687500);
  proc_scale["ttTosemileptonic"]=6.871e+02*lumi/11784986264.000000;

  TString procname   = argv[3];
  HmmAnalyzer Hmm(inputFileList, outFileName, data, isData);
  cout << "dataset " << data << " " << endl;
  Hmm.EventLoop(data, isData,proc_scale[procname]);

  return 0;
}

void HmmAnalyzer::EventLoop(const char *data,const char *isData,double scale)
{ 
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  float muon_mass = 0.1056583745;
  Long64_t nbytes = 0, nb = 0;
  cout<<"Scale : "<<scale<<endl;
      
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%50000==0) cout <<"entry: "<<jentry<<endl;
      double evt_wt = genWeight*scale;
      std::vector<int> tagMuon_idx;
      tagMuon_idx.clear();
      for(int i=0;i<nMuon;i++){
	if(fabs(Muon_eta[i])<2.172 && fabs(Muon_eta[i])>1.740 && Muon_tightId[i] && Muon_pfRelIso04_all[i] < 0.15 && Muon_pt[i]>25.){
	  tagMuon_idx.push_back(i);
	}
      }
      if(!tagMuon_idx.size()) continue;
      int veto_flag=0;
      TLorentzVector dimu, mu1,mu2,mu_guess,mu_pr,dimu_guess;
      //reject events with two well reconstructed muons
      for(int i=0;i<nMuon;i++){
        if(tagMuon_idx[0]!=i && Muon_pt[i]>20. && fabs(Muon_eta[i])<2.5 && Muon_pfRelIso04_all[i] < 0.25  && Muon_charge[i]*Muon_charge[tagMuon_idx[0]]==-1){
          
	
	mu1.SetPtEtaPhiM(Muon_pt[tagMuon_idx[0]],(Muon_eta)[tagMuon_idx[0]],(Muon_phi)[tagMuon_idx[0]],muon_mass);
	mu2.SetPtEtaPhiM(Muon_pt[i],(Muon_eta)[i],(Muon_phi)[i],muon_mass);
	dimu=mu1+mu2;
	if(fabs(dimu.M()-91.18)<4. && fabs(Muon_eta[i])<2.172 && fabs(Muon_eta[i])>2.322){veto_flag =1;break;}
        
	}
      }
      if(veto_flag)continue;
      double phi_diff=9999.;
      int pr_mu_idx=-99;
      double phi_m2=DeltaPhi((Muon_phi)[tagMuon_idx[0]],M_PI);
      //select real probe muon closest in phi to candidate probe
      for(int i=0;i<nMuon;i++){
        if(fabs(Muon_eta[i])>2.172 && fabs(Muon_eta[i])<2.322 && (Muon_charge[i]*Muon_charge[tagMuon_idx[0]])==-1 && Muon_pfRelIso04_all[i] < 0.25 && Muon_pt[i]>10.){
	    if(abs(phi_diff)>abs(DeltaPhi((Muon_phi)[i],phi_m2))){
	      pr_mu_idx=i;
	      phi_diff = DeltaPhi(Muon_phi[i],phi_m2);
	    }
	    
        }
      }
      //setting real probe muon params
      if(pr_mu_idx!=-99){
	mu_pr.SetPtEtaPhiM(Muon_pt[pr_mu_idx],(Muon_eta)[pr_mu_idx],(Muon_phi)[pr_mu_idx],muon_mass);
	dimu=mu1+mu_pr;
      }
      //setting guess probe muon params
      if((Muon_eta)[tagMuon_idx[0]]>0)mu_guess.SetPtEtaPhiM(35.63,2.247,phi_m2,muon_mass);
      else mu_guess.SetPtEtaPhiM(35.63,-2.247,phi_m2,muon_mass);
      dimu_guess =mu1+mu_guess;
      
      veto_flag=0;
      //remove events with well reconstructed b jets
      for (int j =0;j<nJet;j++){
	if(Jet_pt[j]>30. && fabs(Jet_eta[j])<4.7 && Jet_jetId[j]>=2 && Jet_puId[j]>=1){
	    
	  if(Jet_btagDeepB[j]>0.4941){
	    veto_flag=1;
	    break;
	  }
 
	}
      }
      //remove events with well reconstructed electrons
      if(veto_flag)continue;
      for(int i=0;i<nElectron;i++){
	if(Electron_mvaFall17Iso_WP80){
	  veto_flag=1;
	  break;
	}  
      }
      if(veto_flag)continue;
      //remove events with well reconstructed photons
      for(int i=0;i<nPhoton;i++){
	if(Photon_mvaID_WP80){
	  veto_flag=1;
	  break;
	}  
      }
      if(veto_flag)continue;

      //remove events with high MET
      if(MET_pt>30.) continue;
      if(!(dimu_guess.M()>80))continue;
      h_Zguess_pt->Fill(dimu_guess.Pt(),evt_wt);
      h_Zguess_phi->Fill(dimu_guess.Phi(),evt_wt);
      h_Zguess_eta->Fill(dimu_guess.Eta(),evt_wt);
      if(pr_mu_idx!=-99){
	h_Zreal_pt->Fill(dimu.Pt(),evt_wt);
	h_Zreal_phi->Fill(dimu.Phi(),evt_wt);
	h_Zreal_eta->Fill(dimu.Eta(),evt_wt);
      }
      }
   }
