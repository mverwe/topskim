#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "IsoEMuTree.h"
#include <string>
#include <vector>

const Float_t eleM = .000510998910;
const Float_t muM = .1056583715;
const float lepPtCut = 18.;

void anaIsoEMuSkimTree(TString strIn = "emuskim_0.root", const std::string outFileName = "topEmuSkim_0.root") {

  const int iDebug = 0;
  
  TFile *f = TFile::Open(strIn.Data());
  TTree *tr = dynamic_cast<TTree*>(f->Get("skimTree"));

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  isoEMuTree_p = new TTree("isoEMuTree", "isoEMuTree");
  BookTree();

  //input tree variables
  // unsigned int run, lumi;
  // ULong64_t evt;
  int hiBin;
  float vz;
  
  int nLep;
  const int nLepMax = 4;
  int   lepID[nLepMax];
  float lepPt[nLepMax];
  float lepPhi[nLepMax];
  float lepEta[nLepMax];
  int   lepChg[nLepMax];
  float lepIso[nLepMax];
  std::vector<std::vector<float> > *lepIsoCones = 0;
  
  const int nMaxJets = 500;
  int nJt;
  float jtPt[nMaxJets];
  float jtPhi[nMaxJets];
  float jtEta[nMaxJets];
  float jtM[nMaxJets];
  float discr_csvV1[nMaxJets];
  int   refparton_flavorForB[nMaxJets];
  
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
  tr->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);
  
  int emuPairs = 0;
  int emuPairsMassCut = 0;
  int emuPairsMassCutNJet2 = 0;
  
  int nEntries = (int)tr->GetEntries();
  for(int entry = 0; entry < nEntries; entry++){
    tr->GetEntry(entry);

    hiBin_ = hiBin;
    vz_    = vz;

    nLep_ = 0;
    lepIsoCones_.clear();
        
    for(int ilep = 0; ilep<nLep;  ilep++) {
      if(lepPt[ilep]<lepPtCut) continue;

      //find closest jet
      double minDR = 999.;
      int    ijMatch = -1;
      for(int ij = 0; ij<nJt; ++ij) {
        double drToLep = sqrt(pow(acos(cos(jtPhi[ij]-lepPhi[ilep])),2)+pow(jtEta[ij]-lepEta[ilep],2));
        if(drToLep<minDR) {
          ijMatch = ij;
          minDR = drToLep;
        }
      }

      //assign values to variables for output tree
      lepID_[nLep_] = lepID[ilep];
      lepPt_[nLep_] = lepPt[ilep];
      lepEta_[nLep_] = lepEta[ilep];
      lepPhi_[nLep_] = lepPhi[ilep];
      lepChg_[nLep_] = lepChg[ilep];
      lepIso_[nLep_] = lepIso[ilep];
      lepIsoCones_.push_back(lepIsoCones->at(ilep));

      if(minDR<0.3 && abs(refparton_flavorForB[ijMatch])==5)
        lepFromB_[nLep_] = 1;
      else
        lepFromB_[nLep_] = 0;
      
      ++nLep_;
    }

    isoEMuTree_p->Fill();
  }

  outFile_p->cd();
  isoEMuTree_p->Write();
  outFile_p->Close();
  delete outFile_p;
  
}
