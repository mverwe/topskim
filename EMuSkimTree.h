#ifndef EMuSkimTree_h
#define EMuSkimTree_h

#include "TTree.h"

#include <iostream>

TTree* skimTree_p = 0;

const Int_t nLep = 4;
const Int_t eleID = 11;
const Int_t muID = 13;
const Float_t eleM = .000510998910;
const Float_t muM = .1056583715;

UInt_t run_, lumi_;
ULong64_t evt_;
Int_t hiBin_;
Float_t hiHF_;
Float_t vz_;
Float_t rho_;

Int_t nLep_;
Int_t lepID_[nLep];
Float_t lepPt_[nLep];
Float_t lepPhi_[nLep];
Float_t lepEta_[nLep];
Int_t lepChg_[nLep];
Float_t lepIso_[nLep];
std::vector<std::vector<float>> lepIsoCones_;

const Int_t nMaxJets = 500;
Int_t nJt_;
Float_t jtPt_[nMaxJets];
Float_t jtPhi_[nMaxJets];
Float_t jtEta_[nMaxJets];
Float_t jtM_[nMaxJets];
Float_t discr_csvV1_[nMaxJets];
int     refparton_flavorForB_[nMaxJets];

// List of branches
TBranch        *b_run;   //!
TBranch        *b_evt;   //!
TBranch        *b_lumi;   //!
TBranch        *b_hiBin;   //!
TBranch        *b_hiHF;   //!
TBranch        *b_vz;   //!
TBranch        *b_rho;   //!
TBranch        *b_nLep;   //!
TBranch        *b_lepID;   //!
TBranch        *b_lepPt;   //!
TBranch        *b_lepPhi;   //!
TBranch        *b_lepEta;   //!
TBranch        *b_lepChg;   //!
TBranch        *b_lepIso; //!
TBranch        *b_lepIsoCones; //!
TBranch        *b_nJt;   //!
TBranch        *b_jtPt;   //!
TBranch        *b_jtPhi;   //!
TBranch        *b_jtEta;   //!
TBranch        *b_jtM;   //!
TBranch        *b_discr_csvV1;//!
TBranch        *b_refparton_flavorForB; //!

void BookTree()
{
  if(skimTree_p == NULL){
    std::cout << "BOOKTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }

  skimTree_p->Branch("run", &run_, "run/i");
  skimTree_p->Branch("evt", &evt_, "evt/l");
  skimTree_p->Branch("lumi", &lumi_, "lumi/i");
  skimTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  skimTree_p->Branch("hiHF", &hiHF_, "hiHF/F");
  skimTree_p->Branch("vz", &vz_, "vz/F");
  skimTree_p->Branch("rho", &rho_, "rho/F");

  skimTree_p->Branch("nLep", &nLep_, "nLep/I");
  skimTree_p->Branch("lepID", lepID_, Form("lepID[%d]/I",nLep));
  skimTree_p->Branch("lepPt", lepPt_, Form("lepPt[%d]/F", nLep));
  skimTree_p->Branch("lepPhi", lepPhi_, Form("lepPhi[%d]/F", nLep));
  skimTree_p->Branch("lepEta", lepEta_, Form("lepEta[%d]/F", nLep));
  skimTree_p->Branch("lepChg", lepChg_, Form("lepChg[%d]/I", nLep));
  skimTree_p->Branch("lepIso", lepIso_, Form("lepIso[%d]/F", nLep));
  skimTree_p->Branch("lepIsoCones", &lepIsoCones_); 
 
  skimTree_p->Branch("nJt", &nJt_, "nJt/I");
  skimTree_p->Branch("jtPt", jtPt_, "jtPt[nJt]/F");
  skimTree_p->Branch("jtPhi", jtPhi_, "jtPhi[nJt]/F");
  skimTree_p->Branch("jtEta", jtEta_, "jtEta[nJt]/F");
  skimTree_p->Branch("jtM", jtM_, "jtM[nJt]/F");
  skimTree_p->Branch("discr_csvV1", discr_csvV1_, "discr_csvV1[nJt]/F");
  skimTree_p->Branch("refparton_flavorForB", refparton_flavorForB_, "refparton_flavorForB[nJt]/I");
  
  return;
}


void ReadTree()
{
  if(skimTree_p == NULL){
    std::cout << "READTREE error; skimTree_p is NULL. return" << std::endl;
    return;
  }

  /* skimTree_p->SetBranchAddress("run", &run, &b_run); */
  /* skimTree_p->SetBranchAddress("evt", &evt, &b_evt); */
  /* skimTree_p->SetBranchAddress("lumi", &lumi, &b_lumi); */
  /* skimTree_p->SetBranchAddress("hiBin", &hiBin, &b_hiBin); */
  /* skimTree_p->SetBranchAddress("vz", &vz, &b_vz); */
  /* skimTree_p->SetBranchAddress("nLep", &nLep, &b_nLep); */
  /* skimTree_p->SetBranchAddress("lepID", lepID, &b_lepID); */
  /* skimTree_p->SetBranchAddress("lepPt", lepPt, &b_lepPt); */
  /* skimTree_p->SetBranchAddress("lepPhi", lepPhi, &b_lepPhi); */
  /* skimTree_p->SetBranchAddress("lepEta", lepEta, &b_lepEta); */
  /* skimTree_p->SetBranchAddress("lepChg", lepChg, &b_lepChg); */
  /* skimTree_p->SetBranchAddress("nJt", &nJt, &b_nJt); */
  /* skimTree_p->SetBranchAddress("jtPt", jtPt, &b_jtPt); */
  /* skimTree_p->SetBranchAddress("jtPhi", jtPhi, &b_jtPhi); */
  /* skimTree_p->SetBranchAddress("jtEta", jtEta, &b_jtEta); */
  /* skimTree_p->SetBranchAddress("jtM", jtM, &b_jtM); */


  return;
}

#endif
