#ifndef IsoEMuTree_h
#define IsoEMuTree_h

#include "TTree.h"

#include <iostream>

TTree* isoEMuTree_p = 0;

UInt_t run_, lumi_;
ULong64_t evt_;
Int_t hiBin_;
Float_t vz_;

const int nLepMax = 500;
Int_t nLep_;
Int_t   lepID_[nLepMax];
Float_t lepPt_[nLepMax];
Float_t lepPhi_[nLepMax];
Float_t lepEta_[nLepMax];
Int_t   lepChg_[nLepMax];
Float_t lepIso_[nLepMax];
std::vector<std::vector<float>> lepIsoCones_;
Int_t   lepFromB_[nLepMax];

// List of branches
TBranch        *b_run;   //!
TBranch        *b_evt;   //!
TBranch        *b_lumi;   //!
TBranch        *b_hiBin;   //!
TBranch        *b_vz;   //!

TBranch        *b_nLep;   //!
TBranch        *b_lepID;   //!
TBranch        *b_lepPt;   //!
TBranch        *b_lepPhi;   //!
TBranch        *b_lepEta;   //!
TBranch        *b_lepChg;   //!
TBranch        *b_lepIso; //!
TBranch        *b_lepIsoCones; //!
TBranch        *b_lepFromB; //!


void BookTree()
{
  if(isoEMuTree_p == NULL){
    std::cout << "BOOKTREE error; isoEMuTree_p is NULL. return" << std::endl;
    return;
  }

  isoEMuTree_p->Branch("run", &run_, "run/i");
  isoEMuTree_p->Branch("evt", &evt_, "evt/l");
  isoEMuTree_p->Branch("lumi", &lumi_, "lumi/i");
  isoEMuTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  isoEMuTree_p->Branch("vz", &vz_, "vz/F");

  isoEMuTree_p->Branch("nLep", &nLep_, "nLep/I");
  isoEMuTree_p->Branch("lepID", lepID_, Form("lepID[%d]/I",nLepMax));
  isoEMuTree_p->Branch("lepPt", lepPt_, Form("lepPt[%d]/F", nLepMax));
  isoEMuTree_p->Branch("lepPhi", lepPhi_, Form("lepPhi[%d]/F", nLepMax));
  isoEMuTree_p->Branch("lepEta", lepEta_, Form("lepEta[%d]/F", nLepMax));
  isoEMuTree_p->Branch("lepChg", lepChg_, Form("lepChg[%d]/I", nLepMax));
  isoEMuTree_p->Branch("lepIso", lepIso_, Form("lepIso[%d]/F", nLepMax));
  isoEMuTree_p->Branch("lepIsoCones", &lepIsoCones_); 
  isoEMuTree_p->Branch("lepFromB", lepFromB_, Form("lepFromB[%d]/I",nLepMax));
  
  return;
}


void ReadTree()
{
  if(isoEMuTree_p == NULL){
    std::cout << "READTREE error; isoEMuTree_p is NULL. return" << std::endl;
    return;
  }

  return;
}

#endif
