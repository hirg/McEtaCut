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
#include "TProfile.h"
#include "../Utility/functions.h"
#include "../Utility/StSpinAlignmentCons.h"

using namespace std;

TH1F* readpt(int energy, int pid, int centrality);
TH1F* readeta(int energy, int pid, int centrality);
TH1F* readphi(int energy, int pid, int centrality);
void getKinematics(TLorentzVector& lLambda, double const mass);
void setDecayChannels(int const pid);
void decayAndFill(int const pid, TLorentzVector* lLambda, TClonesArray& daughters);
void fill(int const pid, TLorentzVector* lLambda, TLorentzVector const& lProton, TLorentzVector const& lPion);
void write(int energy,int pid,int counter);
TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec);
bool passEtaCut(float eta, int BinEta);
bool passPtCut(float Pt_Proton, float Pt_Pion, float Pt_Lambda);
bool passSTARCut(TLorentzVector lProton, TLorentzVector lPion, TLorentzVector* lLambda);
bool Sampling(int const pid, TF1 *f_pHPhy,float CosThetaStar);

double polarization(double *x_val, double *par)
{
  double x = x_val[0];
  double pH = par[0];
  double Norm = par[1];
  double alpha = par[2];

  double dNdCosThetaStar = Norm*(1+alpha*pH*x);

  return dNdCosThetaStar;
}

int const decayMother[2] = {3122,-3122};
int const decayChannelsFirst = 1058;
int const decayChannelsSecond = 1061;
int const decayChannels = 1058; // 0: Lambda->p+pi-, 1: Lambdabar->pbar+pi+
float const spinDirection[2] = {1.0,-1.0}; // pbar's momentum is opposite to anti-Lambda spin
float const alphaH = 0.642;
float const invMass = 1.11568;

// no STAR cuts
TH3F *h_Tracks, *h_TracksProton, *h_TracksPion;
TH3F *h_Eta;
TH2F *h_phiRP, *h_cosRP;
TProfile *p_cosRP, *p_sinRP;
TProfile *p_cosDau[20],     *p_cosLambda[20],     *p_cosDauOnly[20];
TProfile *p_cosInteDau[20], *p_cosInteLambda[20], *p_cosInteDauOnly[20];
TProfile *p_sinDau[20],     *p_sinLambda[20],     *p_sinDauOnly[20];
TProfile *p_sinInteDau[20], *p_sinInteLambda[20], *p_sinInteDauOnly[20];

// STAR pT cuts
TProfile *p_cosPt, *p_sinPt;
TProfile *p_cosInteDauPt[20];
TProfile *p_sinInteDauPt[20];

// STAR cuts
TProfile *p_cosSTAR, *p_sinSTAR;

TFile *File_InPut;
// sampling functions
TF1 *f_pHPhy;
TH1F *h_pt, *h_eta, *h_phi;


TPythia6Decayer* pydecay;

void McLambdaEta(int energy = 6, int pid = 0, int cent = 0, int const NMax = 10000000, int counter = 0) // pid = 0 for Lambda, 1 for anti-Lambda
{
  int const BinPt    = vmsa::BinPt;
  int const BinY     = vmsa::BinY;
  int const BinPhi   = vmsa::BinPhi;

  h_Tracks = new TH3F("h_Tracks","h_Tracks",BinPt,vmsa::ptMin,vmsa::ptMax,10.0*BinY,-10.0,10.0,BinPhi,-TMath::Pi(),TMath::Pi());
  h_TracksProton = new TH3F("h_TracksProton","h_TracksProton",BinPt,vmsa::ptMin,vmsa::ptMax,10.0*BinY,-10.0,10.0,BinPhi,-TMath::Pi(),TMath::Pi());
  h_TracksPion = new TH3F("h_TracksPion","h_TracksPion",BinPt,vmsa::ptMin,vmsa::ptMax,10.0*BinY,-10.0,10.0,BinPhi,-TMath::Pi(),TMath::Pi());
  h_Eta = new TH3F("h_Eta","h_Eta",10*BinY,-10.0,10.0,10*BinY,-10.0,10.0,10*BinY,-10.0,10.0); // eta for phi, K+ and K-

  h_phiRP = new TH2F("h_phiRP","h_phiRP",BinPt,vmsa::ptMin,vmsa::ptMax,BinPhi,-TMath::Pi(),TMath::Pi());
  h_cosRP = new TH2F("h_cosRP","h_cosRP",BinPt,vmsa::ptMin,vmsa::ptMax,2.0*BinY,-1.0,1.0);

  p_cosRP = new TProfile("p_cosRP","p_cosRP",BinPt,vmsa::ptMin,vmsa::ptMax);
  p_sinRP = new TProfile("p_sinRP","p_sinRP",BinPt,vmsa::ptMin,vmsa::ptMax);


  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    string ProName;
    ProName = Form("p_cosDau_%d",i_eta);
    p_cosDau[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax);
    ProName = Form("p_cosInteDau_%d",i_eta);
    p_cosInteDau[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);

    ProName = Form("p_cosLambda_%d",i_eta);
    p_cosLambda[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax);
    ProName = Form("p_cosInteLambda_%d",i_eta);
    p_cosInteLambda[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);

    ProName = Form("p_sinDau_%d",i_eta);
    p_sinDau[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax);
    ProName = Form("p_sinInteDau_%d",i_eta);
    p_sinInteDau[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);

    ProName = Form("p_sinLambda_%d",i_eta);
    p_sinLambda[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax);
    ProName = Form("p_sinInteLambda_%d",i_eta);
    p_sinInteLambda[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);

    ProName = Form("p_cosDauOnly_%d",i_eta);
    p_cosDauOnly[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax);
    ProName = Form("p_cosInteDauOnly_%d",i_eta);
    p_cosInteDauOnly[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);

    ProName = Form("p_sinDauOnly_%d",i_eta);
    p_sinDauOnly[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),BinPt,vmsa::ptMin,vmsa::ptMax);
    ProName = Form("p_sinInteDauOnly_%d",i_eta);
    p_sinInteDauOnly[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);
  }

  p_cosPt = new TProfile("p_cosPt","p_cosPt",BinPt,vmsa::ptMin,vmsa::ptMax);
  p_sinPt = new TProfile("p_sinPt","p_sinPt",BinPt,vmsa::ptMin,vmsa::ptMax);

  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    string ProName;
    ProName = Form("p_cosInteDauPt_%d",i_eta);
    p_cosInteDauPt[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);

    ProName = Form("p_sinInteDauPt_%d",i_eta);
    p_sinInteDauPt[i_eta] = new TProfile(ProName.c_str(),ProName.c_str(),1,vmsa::McEtaBin[i_eta]-0.1,vmsa::McEtaBin[i_eta]+0.1);
  }

  p_cosSTAR = new TProfile("p_cosSTAR","p_cosSTAR",1,-1,1);
  p_sinSTAR = new TProfile("p_sinSTAR","p_sinSTAR",1,-1,1);

  string InPutKinematics = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/Data/Kinematics.root",vmsa::mBeamEnergy[energy].c_str());
  File_InPut = TFile::Open(InPutKinematics.c_str());

  h_pt  = readpt(energy,pid,cent);
  h_eta = readeta(energy,pid,cent);
  h_phi = readphi(energy,pid,cent);

  float pHPhy = 0.02; // input pH value
  f_pHPhy = new TF1("f_pHPhy",polarization,-1.0,1.0,3);
  f_pHPhy->FixParameter(0,pHPhy);
  f_pHPhy->FixParameter(1,1.0);
  f_pHPhy->FixParameter(2,alphaH*spinDirection[pid]);

  TStopwatch* stopWatch = new TStopwatch();
  stopWatch->Start();
  if(gRandom) delete gRandom;
  gRandom = new TRandom3();
  gRandom->SetSeed();

  pydecay = TPythia6Decayer::Instance();
  pydecay->Init();
  setDecayChannels(pid); // Lambda->p+pi- | Lambdabar->pbar+pi+

  TClonesArray ptl("TParticle", 10);
  TLorentzVector *lLambda = new TLorentzVector();
  for(int i_ran = 0; i_ran < NMax; ++i_ran)
  {
    if (floor(10.0*i_ran/ static_cast<float>(NMax)) > floor(10.0*(i_ran-1)/ static_cast<float>(NMax)))
    cout << "=> processing data: " << 100.0*i_ran/ static_cast<float>(NMax) << "%" << endl;

    getKinematics(*lLambda,invMass);
    decayAndFill(pid,lLambda,ptl);
  }
  cout << "=> processing data: 100%" << endl;
  cout << "work done!" << endl;

  write(energy,pid,counter);

  stopWatch->Stop();   
  stopWatch->Print();
}

TH1F* readpt(int energy, int pid, int centrality)
{
  TH3F *h_Kinematics = (TH3F*)File_InPut->Get("h_Lambda_pt_eta_phi");
  TH1F *h_pt = (TH1F*)h_Kinematics->Project3D("x")->Clone();

  return h_pt;
}

TH1F* readeta(int energy, int pid, int centrality)
{
  TH3F *h_Kinematics = (TH3F*)File_InPut->Get("h_Lambda_pt_eta_phi");
  TH1F *h_eta = (TH1F*)h_Kinematics->Project3D("y")->Clone();

  return h_eta;
}

TH1F* readphi(int energy, int pid, int centrality)
{
  TH3F *h_Kinematics = (TH3F*)File_InPut->Get("h_Lambda_pt_eta_phi");
  TH1F *h_phi = (TH1F*)h_Kinematics->Project3D("z")->Clone();

  return h_phi;
}

void getKinematics(TLorentzVector& lLambda, double const mass)
{
  double const pt = h_pt->GetRandom();
  double const eta = h_eta->GetRandom();
  double const phi = h_phi->GetRandom();

  lLambda.SetPtEtaPhiM(pt,eta,phi,mass);
}

void setDecayChannels(int const pid)
{
  int const mdme = decayChannels;
  cout << "mdme = " << mdme << endl;
  for (int idc = decayChannelsFirst; idc < decayChannelsSecond + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0); // close all decay channel
  TPythia6::Instance()->SetMDME(mdme, 1, 1); // open the one we need
  int *PYSeed = new int;
  TPythia6::Instance()->SetMRPY(1,(int)PYSeed); // Random seed
}

void decayAndFill(int const pid, TLorentzVector* lLambda, TClonesArray& daughters)
{
  pydecay->Decay(decayMother[pid], lLambda);
  pydecay->ImportParticles(&daughters);

  TLorentzVector lProton;
  TLorentzVector lPion;

  int nTrk = daughters.GetEntriesFast();
  // cout << "nTrk = " << nTrk << endl;
  for (int iTrk = 0; iTrk < nTrk; ++iTrk)
  {
    TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
    // cout << "PdgCode = " << ptl0->GetPdgCode() << endl;

    switch (TMath::Abs(ptl0->GetPdgCode()))
    {
      case 2212:
	ptl0->Momentum(lProton);
	break;
      case 211:
	ptl0->Momentum(lPion);
	break;
      default:
	break;
    }
  }
  daughters.Clear("C");
  // cout << "lLambda.M() = " << lLambda->M() << endl;
  // cout << "lProton.M() = " << lProton.M() << endl;
  // cout << "lPion.M() = " << lPion.M() << endl;

  fill(pid,lLambda,lProton,lPion);
}

void fill(int const pid, TLorentzVector* lLambda, TLorentzVector const& lProton, TLorentzVector const& lPion)
{
  TVector3 vMcKpBoosted = spinDirection[pid]*CalBoostedVector(lProton,lLambda); // boost Lambda back to Lambda rest frame

  float Pt_Lambda = lLambda->Pt();
  float Eta_Lambda = lLambda->Eta();
  float Phi_Lambda = lLambda->Phi();

  float Pt_Proton = lProton.Pt();
  float Eta_Proton = lProton.Eta();
  float Phi_Proton = lProton.Phi();

  float Pt_Pion = lPion.Pt();
  float Eta_Pion = lPion.Eta();
  float Phi_Pion = lPion.Phi();

  float Psi = 0.0;
  TVector3 nQ(TMath::Sin(Psi),-TMath::Cos(Psi),0.0); // direction of angular momentum with un-smeared EP
  float CosThetaStarSimple = vMcKpBoosted.Dot(nQ);
  if(!Sampling(pid,f_pHPhy,CosThetaStarSimple)) return;

  // float SinPhiStarRP = TMath::Sin(vMcKpBoosted.Theta())*TMath::Sin(Psi-vMcKpBoosted.Phi());

  float CosThetaStarRP = (3.0*spinDirection[pid]/alphaH)*CosThetaStarSimple;
  float SinPhiStarRP = (8.0*spinDirection[pid]/(alphaH*TMath::Pi()))*TMath::Sin(Psi-vMcKpBoosted.Phi());

  h_phiRP->Fill(Pt_Lambda,Phi_Lambda);
  h_cosRP->Fill(Pt_Lambda,CosThetaStarRP);
  h_Tracks->Fill(Pt_Lambda,Eta_Lambda,Phi_Lambda);
  h_TracksProton->Fill(Pt_Proton,Eta_Proton,Phi_Proton);
  h_TracksPion->Fill(Pt_Pion,Eta_Pion,Phi_Pion);
  h_Eta->Fill(Eta_Lambda,Eta_Proton,Eta_Pion);
  p_cosRP->Fill(Pt_Lambda,CosThetaStarRP);
  p_sinRP->Fill(Pt_Lambda,SinPhiStarRP);

  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    if( passEtaCut(Eta_Lambda,i_eta) ) 
    {
      p_cosLambda[i_eta]->Fill(Pt_Lambda,CosThetaStarRP);
      p_cosInteLambda[i_eta]->Fill(vmsa::McEtaBin[i_eta],CosThetaStarRP);
      p_sinLambda[i_eta]->Fill(Pt_Lambda,SinPhiStarRP);
      p_sinInteLambda[i_eta]->Fill(vmsa::McEtaBin[i_eta],SinPhiStarRP);
    }

    if( passEtaCut(Eta_Proton,i_eta) && passEtaCut(Eta_Pion,i_eta) && passEtaCut(Eta_Lambda,i_eta) )
    {
      p_cosDau[i_eta]->Fill(Pt_Lambda,CosThetaStarRP);
      p_cosInteDau[i_eta]->Fill(vmsa::McEtaBin[i_eta],CosThetaStarRP);
      p_sinDau[i_eta]->Fill(Pt_Lambda,SinPhiStarRP);
      p_sinInteDau[i_eta]->Fill(vmsa::McEtaBin[i_eta],SinPhiStarRP);
    }

    if( passEtaCut(Eta_Proton,i_eta) && passEtaCut(Eta_Pion,i_eta) )
    {
      p_cosDauOnly[i_eta]->Fill(Pt_Lambda,CosThetaStarRP);
      p_cosInteDauOnly[i_eta]->Fill(vmsa::McEtaBin[i_eta],CosThetaStarRP);
      p_sinDauOnly[i_eta]->Fill(Pt_Lambda,SinPhiStarRP);
      p_sinInteDauOnly[i_eta]->Fill(vmsa::McEtaBin[i_eta],SinPhiStarRP);
    }
  }

  if( passPtCut(Pt_Proton,Pt_Pion,Pt_Lambda) )
  { // apply STAR pT cuts
    p_cosRP->Fill(Pt_Lambda,CosThetaStarRP);
    p_sinRP->Fill(Pt_Lambda,SinPhiStarRP);

    if( passEtaCut(Eta_Proton,i_eta) && passEtaCut(Eta_Pion,i_eta) && passEtaCut(Eta_Lambda,i_eta) )
    {
      p_cosInteDauPt[i_eta]->Fill(vmsa::McEtaBin[i_eta],CosThetaStarRP);
      p_sinInteDauPt[i_eta]->Fill(vmsa::McEtaBin[i_eta],SinPhiStarRP);
    }
  }

  if( passSTARCut(lProton,lPion,lLambda) )
  { // apply STAR cuts
    p_cosSTAR->Fill(0.0,CosThetaStarRP);
    p_sinSTAR->Fill(0.0,SinPhiStarRP);
  }
}

TVector3 CalBoostedVector(TLorentzVector const lMcDau, TLorentzVector *lMcVec)
{
  TVector3 vMcBeta = -1.0*lMcVec->BoostVector(); // boost vector

  TLorentzVector lProton = lMcDau;
  lProton.Boost(vMcBeta); // boost proton back to Lambda rest frame
  TVector3 vMcDauStar = lProton.Vect().Unit(); // momentum direction of proton in Lambda rest frame

  return vMcDauStar;
}

bool passEtaCut(float eta, int BinEta)
{
  if(TMath::Abs(eta) >= vmsa::McEtaBin[BinEta]) return kFALSE;

  return kTRUE;
}

bool passPtCut(float Pt_Proton, float Pt_Pion, float Pt_Lambda)
{
  if( !(Pt_Proton > 0.15 && Pt_Pion > 0.15 && Pt_Lambda > 0.4 && Pt_Lambda < 3.0) ) return kFALSE; // pT cut on daughter particles to be consistent with STAR

  return kTRUE;
}

bool passSTARCut(TLorentzVector lProton, TLorentzVector lPion, TLorentzVector* lLambda)
{
  if( !(lProton.Pt() > 0.15 && lProton.Eta() > -1.0 && lProton.Eta() < 1.0) ) 
    return kFALSE; // STAR cut for proton

  if( !(lPion.Pt() > 0.15 && lPion.Eta() > -1.0 && lPion.Eta() < 1.0) ) 
    return kFALSE; // STAR cut for pion

  if( !(lLambda->Pt() > 0.4 && lLambda->Pt() < 3.0 && lLambda->Rapidity() > -1.0 && lLambda->Rapidity() < 1.0) ) 
    return kFALSE; // STAR cut for Lambda

  return kTRUE;
}

bool Sampling(int const pid, TF1 *f_pHPhy,float CosThetaStar)
{
  float wMax;
  if(pid == 0) wMax = f_pHPhy->Eval(1.0);
  else wMax = f_pHPhy->Eval(-1.0);
  return !(gRandom->Rndm() > f_pHPhy->Eval(CosThetaStar)/wMax);
}

void write(int energy, int pid, int counter)
{
  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McLambdaEta_%d_%d.root",vmsa::mBeamEnergy[energy].c_str(),pid,counter);
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();

  h_Tracks->Write();
  h_TracksProton->Write();
  h_TracksPion->Write();
  h_phiRP->Write();
  h_cosRP->Write();
  h_Eta->Write();

  p_cosRP->Write();
  p_sinRP->Write();

  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    p_cosDau[i_eta]->Write();
    p_cosInteDau[i_eta]->Write();
    p_cosLambda[i_eta]->Write();
    p_cosInteLambda[i_eta]->Write();
    p_cosDauOnly[i_eta]->Write();
    p_cosInteDauOnly[i_eta]->Write();

    p_sinDau[i_eta]->Write();
    p_sinInteDau[i_eta]->Write();
    p_sinLambda[i_eta]->Write();
    p_sinInteLambda[i_eta]->Write();
    p_sinDauOnly[i_eta]->Write();
    p_sinInteDauOnly[i_eta]->Write();
  }

    p_cosRP->Write()
    p_sinRP->Write()

    if( passEtaCut(Eta_Proton,i_eta) && passEtaCut(Eta_Pion,i_eta) && passEtaCut(Eta_Lambda,i_eta) )
    {
      p_cosInteDauPt[i_eta]->Write();
      p_sinInteDauPt[i_eta]->Write();
    }

  p_cosSTAR->Write();
  p_sinSTAR->Write();

  File_OutPut->Close();
}

