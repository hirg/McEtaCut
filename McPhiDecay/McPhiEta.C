#include <iostream>
#include <string> 
#include <map>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TVector3.h"
#include "../Utility/functions.h"
#include "../Utility/StSpinAlignmentCons.h"

using namespace std;

TH1F* readpt(int energy, int pid, int centrality);
TH1F* readeta(int energy, int pid, int centrality);
TH1F* readphi(int energy, int pid, int centrality);
void getKinematics(TLorentzVector& lPhi, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters);
void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus);
void write(int energy);
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
bool passEtaCut(float eta, int BinEta);

// histograms
TH3F *h_Tracks;
TH3F *h_Eta;
TH2F *h_phiRP, *h_cosRP;
TH2F *h_CosEtaKaon[20], *h_CosEtaPhi[20], *h_CosEtaKOnly[20];

TFile *File_InPut;
// sampling functions
TH1F *h_pt, *h_eta, *h_phi;

TPythia6Decayer* pydecay;

void McPhiEta(int energy = 6, int pid = 0, int cent = 0, int const NMax = 10000000)
{
  int   const BinPt    = vmsa::BinPt;
  int   const BinY     = vmsa::BinY;
  int   const BinPhi   = vmsa::BinPhi;

  h_Tracks = new TH3F("h_Tracks","h_Tracks",BinPt,vmsa::ptMin,vmsa::ptMax,10.0*BinY,-10.0,10.0,BinPhi,-TMath::Pi(),TMath::Pi());
  h_Eta = new TH3F("h_Eta","h_Eta",10*BinY,-10.0,10.0,10*BinY,-10.0,10.0,10*BinY,-10.0,10.0); // eta for phi, K+ and K-

  h_phiRP = new TH2F("h_phiRP","h_phiRP",BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());
  h_cosRP = new TH2F("h_cosRP","h_cosRP",BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);

  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    string HistName;
    HistName = Form("h_CosEtaKaon_%d",i_eta);
    h_CosEtaKaon[i_eta] = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);
    HistName = Form("h_CosEtaPhi_%d",i_eta);
    h_CosEtaPhi[i_eta] = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);
    HistName = Form("h_CosEtaKOnly_%d",i_eta);
    h_CosEtaKOnly[i_eta] = new TH2F(HistName.c_str(),HistName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax,BinY,-1.0,1.0);
  }

  string InPutKinematics = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Kinematics.root",vmsa::mBeamEnergy[energy].c_str());
  File_InPut = TFile::Open(InPutKinematics.c_str());

  h_pt  = readpt(energy,pid,cent);
  h_eta = readeta(energy,pid,cent);
  h_phi = readphi(energy,pid,cent);

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

  pydecay = TPythia6Decayer::Instance();
  pydecay->Init();
  setDecayChannels(pid); // phi--> K+K-

  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lPhi = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;

    getKinematics(*lPhi,vmsa::InvMass[pid]);
    decayAndFill(vmsa::decayMother[pid],lPhi,ptl);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write(energy);

  stopWatch->Stop();   
  stopWatch->Print();
}

TH1F* readpt(int energy, int pid, int centrality)
{
  TH3F *h_Kinematics = (TH3F*)File_InPut->Get("h_phi");
  TH1F *h_pt = (TH1F*)h_Kinematics->Project3D("x")->Clone();

  return h_pt;
}

TH1F* readeta(int energy, int pid, int centrality)
{
  TH3F *h_Kinematics = (TH3F*)File_InPut->Get("h_phi");
  TH1F *h_eta = (TH1F*)h_Kinematics->Project3D("y")->Clone();

  return h_eta;
}

TH1F* readphi(int energy, int pid, int centrality)
{
  TH3F *h_Kinematics = (TH3F*)File_InPut->Get("h_phi");
  TH1F *h_phi = (TH1F*)h_Kinematics->Project3D("z")->Clone();

  return h_phi;
}

void getKinematics(TLorentzVector& lPhi, double const mass)
{
  double const pt = h_pt->GetRandom();
  double const eta = h_eta->GetRandom();
  double const phi = h_phi->GetRandom();

  lPhi.SetPtEtaPhiM(pt,eta,phi,mass);
}

void setDecayChannels(int const pid)
{
  int const mdme = vmsa::decayChannels[pid];
  cout << "mdme = " << mdme << endl;
  for (int idc = vmsa::decayChannelsFirst[pid]; idc < vmsa::decayChannelsSecond[pid] + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); // close all decay channel
  TPythia6::Instance()->SetMDME(mdme, 1, 1); // open the one we need
  int *PYSeed = new int;
  TPythia6::Instance()->SetMRPY(1,(int)PYSeed); // Random seed
}

void decayAndFill(int const kf, TLorentzVector* lPhi, TClonesArray& daughters)
{
  pydecay->Decay(kf, lPhi);
  pydecay->ImportParticles(&daughters);

  TLorentzVector lKplus;
  TLorentzVector lKminus;

  int nTrk = daughters.GetEntriesFast();
  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

    switch (ptl0->GetPdgCode())
    {
      case 321:
	ptl0->Momentum(lKplus);
	break;
      case -321:
	ptl0->Momentum(lKminus);
	break;
      default:
	break;
    }
  }
  daughters.Clear("C");

  fill(lPhi,lKplus,lKminus);
}

void fill(TLorentzVector* lPhi, TLorentzVector const& lKplus, TLorentzVector const& lKminus)
{
  TVector3 vMcKpBoosted = CalBoostedVector(lKplus,lPhi); // boost Kplus back to phi-meson rest frame

  float Pt_lPhi = lPhi->Pt();
  float Eta_lPhi = lPhi->Eta();
  float Eta_lKplus = lKplus.Eta();
  float Eta_lKminus = lKminus.Eta();

  TVector3 nQ(0.0,-1.0,0.0); // direction of angular momentum with un-smeared EP
  float CosThetaStarRP = vMcKpBoosted.Dot(nQ);

  h_phiRP->Fill(Pt_lPhi,lPhi->Phi());
  h_cosRP->Fill(Pt_lPhi,CosThetaStarRP);
  h_Tracks->Fill(Pt_lPhi,Eta_lPhi,lPhi->Phi());
  h_Eta->Fill(Eta_lPhi,Eta_lKplus,Eta_lKminus);

  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    if( passEtaCut(Eta_lPhi,i_eta) ) h_CosEtaPhi[i_eta]->Fill(Pt_lPhi,CosThetaStarRP);

    if( passEtaCut(Eta_lKplus,i_eta) && passEtaCut(Eta_lKminus,i_eta) && passEtaCut(Eta_lPhi,i_eta) )
    {
      h_CosEtaKaon[i_eta]->Fill(Pt_lPhi,CosThetaStarRP);
    }

    if( passEtaCut(Eta_lKplus,i_eta) && passEtaCut(Eta_lKminus,i_eta) )
    {
      h_CosEtaKOnly[i_eta]->Fill(Pt_lPhi,CosThetaStarRP);
    }
  }
}

TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec)
{
  TVector3 vMcBeta = -1.0*lMcVec->BoostVector(); // boost vector

  TLorentzVector lKaon = lMcDau;
  lKaon.Boost(vMcBeta); // boost Kplus back to phi-meson rest frame
  TVector3 vMcDauStar = lKaon.Vect().Unit(); // momentum direction of Kplus in phi-meson rest frame

  return vMcDauStar;
}

bool passEtaCut(float eta, int BinEta)
{
  if(TMath::Abs(eta) >= vmsa::McEtaBin[BinEta]) return kFALSE;

  return kTRUE;
}

void write(int energy)
{
  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPhiEta.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();

  h_Tracks->Write();
  h_phiRP->Write();
  h_cosRP->Write();
  h_Eta->Write();

  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    h_CosEtaKaon[i_eta]->Write();
    h_CosEtaPhi[i_eta]->Write();
    h_CosEtaKOnly[i_eta]->Write();
  }

  File_OutPut->Close();
}

