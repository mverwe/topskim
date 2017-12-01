#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "ZTree.h"
#include <string>
#include <vector>

const Float_t eleM = .000510998910;
const Float_t muM = .1056583715;
const float lepPtCut = 18.;

void anaZMuMuSkimTree(TString strIn = "emuskim_0.root", const std::string outFileName = "ZSkim_0.root") {

  const int iDebug = 0;
  
  TFile *f = TFile::Open(strIn.Data());
  TTree *tr = dynamic_cast<TTree*>(f->Get("skimTree"));

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  zTree_p = new TTree("zTree", "zTree");
  BookTree();

  //input tree variables
  // unsigned int run, lumi;
  // ULong64_t evt;
  int hiBin;
  float vz;
  
  int nLep;
  const int nLepMax = 4;
  int lepID[nLepMax];
  float lepPt[nLepMax];
  float lepPhi[nLepMax];
  float lepEta[nLepMax];
  int lepChg[nLepMax];
  float lepIso[nLepMax];
  std::vector<std::vector<float> > *lepIsoCones = 0;
  
  const int nMaxJets = 500;
  int nJt;
  float jtPt[nMaxJets];
  float jtPhi[nMaxJets];
  float jtEta[nMaxJets];
  float jtM[nMaxJets];
  float discr_csvV1[nMaxJets];

  tr->SetBranchAddress("run", &run_);
  tr->SetBranchAddress("evt", &evt_);
  tr->SetBranchAddress("lumi", &lumi_);
  tr->SetBranchAddress("hiBin", &hiBin);
  tr->SetBranchAddress("vz", &vz);
  tr->SetBranchAddress("nLep", &nLep);
  tr->SetBranchAddress("lepID", lepID);
  tr->SetBranchAddress("lepPt", lepPt);
  tr->SetBranchAddress("lepPhi", lepPhi);
  tr->SetBranchAddress("lepEta", lepEta);
  tr->SetBranchAddress("lepChg", lepChg);
  tr->SetBranchAddress("lepIso", lepIso);
  tr->SetBranchAddress("lepIsoCones", &lepIsoCones);
  tr->SetBranchAddress("nJt", &nJt);
  tr->SetBranchAddress("jtPt", jtPt);
  tr->SetBranchAddress("jtPhi", jtPhi);
  tr->SetBranchAddress("jtEta", jtEta);
  tr->SetBranchAddress("jtM", jtM);
  tr->SetBranchAddress("discr_csvV1", discr_csvV1);

  int dilepPairs = 0;
    
  int nEntries = (int)tr->GetEntries();
  for(int entry = 0; entry < nEntries; entry++){
    tr->GetEntry(entry);

    hiBin_ = hiBin;
    vz_    = vz;

    for(int ilep = 0; ilep<nDilepMax; ++ilep) {
      dilepPt_[ilep] = -999.;
      dilepEta_[ilep] = -999.;
      dilepPhi_[ilep] = -999.;
      dilepM_[ilep] = -999.;

      dilepDeltaPhi_[ilep] = -999.;
        
      lep1Pt_[ilep] = -999.; 
      lep2Pt_[ilep] = -999.; 
      lep1Eta_[ilep] = -999.;
      lep2Eta_[ilep] = -999.;
      lep1Phi_[ilep] = -999.;
      lep2Phi_[ilep] = -999.;
      lep1Iso_[ilep] = -999.;
      lep2Iso_[ilep] = -999.;

      lep1IsoCones_.clear();
      lep2IsoCones_.clear();
    }
    
    if(nLep<2) continue;

    dilepPairs = 0;
    
    std::vector<TLorentzVector> dileptons;
    for(int ilep = 0; ilep<nLep;  ilep++) {
      for(int jlep = ilep+1; jlep<nLep;  jlep++) {
        if(lepID[ilep]!=lepID[jlep]) continue;   //reject different flavor
        if(lepChg[ilep]*lepChg[jlep]>0) continue; //reject same sign
        if(lepPt[ilep]<lepPtCut) continue;

        //make dilepton pair
        double lepM = muM;
        dilepChannel_[dilepPairs] = 13;
        if(lepID[ilep]==11) {
          lepM = eleM;
          dilepChannel_[dilepPairs] = 11;
        }
        TLorentzVector lep1;
        lep1.SetPtEtaPhiM(lepPt[ilep],lepEta[ilep],lepPhi[ilep],lepM);
        TLorentzVector lep2;
        lep2.SetPtEtaPhiM(lepPt[jlep],lepEta[jlep],lepPhi[jlep],lepM);

        TLorentzVector dilepton = lep1 + lep2;
        
        if(dilepton.M()<20.) continue;

        dileptons.push_back(dilepton);

        dilepPt_[dilepPairs] = dilepton.Pt();
        dilepEta_[dilepPairs] = dilepton.Eta();
        dilepPhi_[dilepPairs] = dilepton.Phi();
        dilepM_[dilepPairs] = dilepton.M();

        dilepDeltaPhi_[ilep] = fabs(TVector2::Phi_mpi_pi(lepPhi[ilep]-lepPhi[jlep]))/TMath::Pi();
        
        lep1Pt_[dilepPairs] = lepPt[ilep];
        lep2Pt_[dilepPairs] = lepPt[jlep];
        lep1Eta_[dilepPairs] = lepEta[ilep];
        lep2Eta_[dilepPairs] = lepEta[jlep];
        lep1Phi_[dilepPairs] = lepPhi[ilep];
        lep2Phi_[dilepPairs] = lepPhi[jlep];
        lep1Iso_[dilepPairs] = lepIso[ilep];
        lep2Iso_[dilepPairs] = lepIso[jlep];

        lep1IsoCones_.push_back(lepIsoCones->at(ilep));
        lep2IsoCones_.push_back(lepIsoCones->at(jlep));
          
        if(iDebug) Printf("dilepton pair %d mass: %f pt: %f",ilep,dilepton.M(),dilepton.Pt());
        ++dilepPairs;
        
      }
    }
    
    nDilep_ = dilepPairs;

    zTree_p->Fill();
  }
  

  Printf("Total dilep pairs: %d ",dilepPairs);

  outFile_p->cd();
  zTree_p->Write();
  outFile_p->Close();
  delete outFile_p;
  
}
