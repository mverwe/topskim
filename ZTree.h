#ifndef ZTree_h
#define ZTree_h

#include "TTree.h"

#include <iostream>

TTree* zTree_p = 0;

UInt_t run_, lumi_;
ULong64_t evt_;
Int_t hiBin_;
Float_t vz_;

const int nDilepMax = 500;
Int_t nDilep_;
Int_t   dilepChannel_[nDilepMax];
Float_t dilepPt_[nDilepMax];
Float_t dilepPhi_[nDilepMax];
Float_t dilepEta_[nDilepMax];
Float_t dilepM_[nDilepMax];
Float_t dilepDeltaPhi_[nDilepMax];

Float_t lep1Pt_[nDilepMax];
Float_t lep2Pt_[nDilepMax];
Float_t lep1Eta_[nDilepMax];
Float_t lep2Eta_[nDilepMax];
Float_t lep1Phi_[nDilepMax];
Float_t lep2Phi_[nDilepMax];
Float_t lep1Iso_[nDilepMax];
Float_t lep2Iso_[nDilepMax];
std::vector<std::vector<float>> lep1IsoCones_;
std::vector<std::vector<float>> lep2IsoCones_;

// List of branches
TBranch        *b_run;   //!
TBranch        *b_evt;   //!
TBranch        *b_lumi;   //!
TBranch        *b_hiBin;   //!
TBranch        *b_vz;   //!

TBranch        *b_nDilep;   //!
TBranch        *b_dilepChannel;   //!
TBranch        *b_dilepPt;   //!
TBranch        *b_dilepPhi;   //!
TBranch        *b_dilepEta;   //!
TBranch        *b_dilepM;   //!
TBranch        *b_dilepDeltaPhi;   //!
TBranch        *b_lep1Pt; //!
TBranch        *b_lep2Pt; //!
TBranch        *b_lep1Eta; //!
TBranch        *b_lep2Eta; //!
TBranch        *b_lep1Phi; //!
TBranch        *b_lep2Phi; //!
TBranch        *b_lep1Iso; //!
TBranch        *b_lep2Iso; //!
TBranch        *b_lep1IsoCones; //!
TBranch        *b_lep2IsoCones; //!

void BookTree()
{
  if(zTree_p == NULL){
    std::cout << "BOOKTREE error; zTree_p is NULL. return" << std::endl;
    return;
  }

  zTree_p->Branch("run", &run_, "run/i");
  zTree_p->Branch("evt", &evt_, "evt/l");
  zTree_p->Branch("lumi", &lumi_, "lumi/i");
  zTree_p->Branch("hiBin", &hiBin_, "hiBin/I");
  zTree_p->Branch("vz", &vz_, "vz/F");

  zTree_p->Branch("nDilep", &nDilep_, "nDilep/I");
  zTree_p->Branch("dilepChannel", &dilepChannel_, "dilepChannel[nDilep]/I");
  zTree_p->Branch("dilepPt", dilepPt_, "dilepPt[nDilep]/F");
  zTree_p->Branch("dilepPhi", dilepPhi_, "dilepPhi[nDilep]/F");
  zTree_p->Branch("dilepEta", dilepEta_, "dilepEta[nDilep]/F");
  zTree_p->Branch("dilepM", dilepM_, "dilepM[nDilep]/F");
  zTree_p->Branch("dilepDeltaPhi", dilepDeltaPhi_, "dilepDeltaPhi[nDilep]/F");
  zTree_p->Branch("lep1Pt",lep1Pt_,"lep1Pt[nDilep]/F");
  zTree_p->Branch("lep2Pt",lep2Pt_,"lep2Pt[nDilep]/F");
  zTree_p->Branch("lep1Eta",lep1Eta_,"lep1Eta[nDilep]/F");
  zTree_p->Branch("lep2Eta",lep2Eta_,"lep2Eta[nDilep]/F");
  zTree_p->Branch("lep1Phi",lep1Phi_,"lep1Phi[nDilep]/F");
  zTree_p->Branch("lep2Phi",lep2Phi_,"lep2Phi[nDilep]/F");
  zTree_p->Branch("lep1Iso",lep1Iso_,"lep1Iso[nDilep]/F");
  zTree_p->Branch("lep2Iso",lep2Iso_,"lep2Iso[nDilep]/F");
  zTree_p->Branch("lep1IsoCones", &lep1IsoCones_);
  zTree_p->Branch("lep2IsoCones", &lep2IsoCones_); 
 
  return;
}


void ReadTree()
{
  if(zTree_p == NULL){
    std::cout << "READTREE error; zTree_p is NULL. return" << std::endl;
    return;
  }

  return;
}

#endif
