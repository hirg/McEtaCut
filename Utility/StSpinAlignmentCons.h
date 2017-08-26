#ifndef StSpinAlignmentCons_h
#define StSpinAlignmentCons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace vmsa
{
  //--------------------------------------------------
  // used in TreeProduction and FillSpinAlginment
  int const NumBeamEnergy = 7;
  // event cut
  float const mVzMaxMap[NumBeamEnergy] = {70.0,50.0,70.0,70.0,40.0,40.0,30.0}; // 7.7 -- 200 GeV
  float const mVrMax = 2.0;
  float const mVzVpdDiffMax = 3.0;
  int const mMatchedToFMin = 2;

  // track cut
  float const mSigScaleMap[NumBeamEnergy] = {1.0,1.0,1.0,1.9,1.0,1.0,1.0}; // 7.7 -- 200 GeV
  float const mDcaEPMax[NumBeamEnergy] = {1.0,1.0,1.0,1.0,1.0,1.0,3.0}; // for event plane reconstruction: 1.0 for BES, 3.0 for 200GeV
  float const mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow
  float const mDcaTrMax_phi = 2.0; // for phi meson mDcaTrMax = 2.0 to fill a tree and apply an additional cut
  int const mHitsDedxMin = 5;
  int const mHitsFitTPCMin = 15;
  int const mHitsMaxTPCMin = 0;
  float const mHitsRatioTPCMin = 0.51;
  float const mEtaMax = 1.0;
  float const mPrimPtMin[NumBeamEnergy] = {0.2,0.2,0.2,0.2,0.2,0.2,0.15}; // for event plane reconstruction and for pion, kaon, proton: 0.15 for 200 GeV, 0.2 for BES
  float const mGlobPtMin = 0.1; // for phi, Lambda, K0s
  float const mPrimPtMax = 2.0;
  float const mPrimPtWeight = 2.0;
  float const mPrimMomMax = 10.0; // also use for gMom
  float const mMass2Min = -10.0;
  // double const MAGFIELDFACTOR = kilogauss;
  int const mTrackMin = 2;
  int const mTrackMin_Full = 4;
  float const mToFYLocalMax = 1.8;
  float const mToFZLocalMax = 1.8;
  float const mNSigmaElectronMax = 2.5;
  float const mNSigmaPionMax = 2.5;
  float const mNSigmaKaonMax = 2.5;
  float const mNSigmaProtonMax = 2.5;
  float const mMassPion = 0.13957;
  float const mMassKaon = 0.49368;
  float const mMassProton = 0.93827;
  float const mSigKaon = 2.5;

  // used constant
  // float const mEta_Gap[4] = {0.05,0.10,0.20,0.50};
  float const mEta_Gap = 0.05;
  float const mShiftOrder[5] = {2.0, 4.0, 6.0, 8.0, 10.0};

  int const pt_total = 25; // pt bin
  int const pt_start = 0;
  int const pt_stop  = 25;
  float const ptRawStart[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
  float const ptRawStop[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};
  double const pt_bin[pt_total+1] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

  // mix event
  int const Bin_Centrality = 9;
  int const Bin_VertexZ = 10;
  int const Bin_Phi_Psi = 5;
  int const Buffer_depth = 5;
  TString const MixEvent[2] = {"SE","ME"};

  TString const vm_tree[2]  = {"PhiMesonEvent","KStarEvent"};
  TString const vm_branch[2] = {"phi_SpinAlignment_branch","KStar_SpinAlignment_branch"};
  int const mList_Delta = 50;
  //--------------------------------------------------

  // used in CalSpinAlginment
  int const pt_rebin = 9; // maximum pt binning
  float const pt_low[NumBeamEnergy][pt_rebin] = {
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2},
    {0.4,0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2}
  };
  float const pt_up[NumBeamEnergy][pt_rebin]  = { 
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0},
    {0.8,1.2,1.8,2.4,3.0,4.2,5.4,7.2,8.0}
  };
  int const pt_rebin_start[NumBeamEnergy][pt_rebin] = {
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24},
    {1,3,5, 8,11,14,17,20,24}
  };
  int const pt_rebin_stop[NumBeamEnergy][pt_rebin]  = {
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24},
    {2,4,7,10,13,16,19,23,24}
  };
  int const pt_rebin_first[NumBeamEnergy] = {0,0,0,0,0,0,0};
  int const pt_rebin_last[NumBeamEnergy]  = {8,8,6,6,6,6,8};
  int const pt_QA[NumBeamEnergy]    = {1,1,2,2,3,3,0};
  int const pt_RawQA[NumBeamEnergy]    = {2,4,6,3,10,12,14};

  std::string const Centrality[9] = {"70%-80%","60%-70%","50%-60%","40%-50%","30%-40%","20%-30%","10%-20%","5%-10%","0%-5%"}; // Centrality bin
  int const Cent_Total = 5;
  int const Cent_start = 0;
  int const Cent_stop  = 5;
  int const cent_low[5] = {2,0,7,4,0}; // 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
  int const cent_up[5]  = {5,8,8,6,3}; // 0 = 20%-60%, 1 = 0-80%, 2 = 0-10%, 3 = 10-40%, 4 = 40-80%
  int const Cent_QA    = 0;

  int const CTS_total = 7; // cos(theta*) bin
  int const CTS_start = 0;
  int const CTS_stop  = 7;
  float const CTS_low[7] = {0.0/7.0,1.0/7.0,2.0/7.0,3.0/7.0,4.0/7.0,5.0/7.0,6.0/7.0};
  float const CTS_up[7]  = {1.0/7.0,2.0/7.0,3.0/7.0,4.0/7.0,5.0/7.0,6.0/7.0,7.0/7.0};

  int const EtaGap_total = 4; // EtaGap bin
  int const Eta_start = 0; // EtaGap bin
  int const Eta_stop  = 1;
  int const Eta_QA    = 0;

  int const Charge_total = 2;
  int const Charge_start = 0;
  int const Charge_stop  = 2;

  int const Sys_start = 0; // SysError bin
  int const Sys_stop  = 1;
  int const Sys_QA    = 0;

  int const Norm_start = 0;
  int const Norm_stop  = 3;
  int const Norm_QA    = 2;

  std::string const mInteMethod[2] = {"Count","Inte"};
  int const Method_start = 0;
  int const Method_stop  = 2;
  int const Method_QA    = 0;

  float const nSigVecSys[3] = {2.0,2.5,3.0};
  int const Sig_start = 0;
  int const Sig_stop  = 3;
  int const Sig_QA    = 2;

  // used for systematic errors
  int const FuncParNum[4] = {5,6,6,6};
  int const Func_start = 0;
  int const Func_stop  = 1;
  int const Func_QA    = 0;

  // shared constant
  std::string const mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  float const mEnergyValue[NumBeamEnergy] = {7.7,11.5,19.6,27.0,39.0,62.4,200.0};
  int const mBeamYear[NumBeamEnergy] = {2010,2010,2011,2011,2010,2010,2011};

  std::string const mPID[3]   = {"Phi","KStar","K0S"};
  float const Norm_Start[3][2]  = {{1.04,0.99},{0.41,0.30},{0.41,0.30}}; // normalise to right and left
  float const Norm_Stop[3][2]   = {{1.05,1.00},{0.46,0.31},{0.46,0.31}};
  float const BW_Start[3]     = {0.994,1.0,1.0};
  // float const BW_Start[3]     = {0.99,1.0,1.0}; // for RooFit
  float const BW_Stop[3]      = {1.050,1.0,1.0};
  float const Width[3]        = {0.00426,0.0487,0.0487};
  float const InvMass_low[3]  = {0.98,0.4,0.4};
  float const InvMass_high[3] = {1.05,0.6,0.6};
  float const nSigVec = 2.0;

  float const ptEffMax = 8.0;
  float const ptMin = 0.2;
  float const ptMax = 5.0;
  // int const BinPt  = 80; // for efficiency
  int const BinPt  = 20;
  int const BinEta = 10;
  int const BinY = 20;
  int const BinPhi = 24;
  int const BinCos = 7;

  // used in McPhiResCorr
  double const acceptanceRapidity = 1.0;
  float const rhoDelta = 0.01;
  float const InvMass[3] = {1.01946,0.89594,0.49761}; // 0: phi, 1: K*, 2 K0S
  int const decayChannelsFirst[3] = {656,617,613};
  int const decayChannelsSecond[3] = {666,619,614};
  int const decayMother[3] = {333,313,310};
  int const decayChannels[3] = {656,617,613}; // 0: phi->K+K-, 1: K*->Kpi, 2 K0S->pi+pi-
  float const McEtaBin[20] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,3.5,4.0};

  // used in RcPhiEffCorr
  std::string const mParType[2] = {"Kplus","Kminus"};
  std::string const mYear[2] = {"run11","run10"};
  std::string const mCuts[2] = {"pr","gl"};
  int const NCentMax = 9; 
  float const weight[NCentMax] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.5};

  // plotting
  int const Color[pt_rebin] = {kGray+2,kBlack,kRed,kCyan,kMagenta,kAzure,kViolet,kBlue,kRed};
  int const Style[pt_rebin] = {20,21,22,23,24,25,26,28,29};

  // ZDC-SMD constant
  std::string const mEastWest[2] = {"East","West"};
  std::string const mVertHori[2] = {"Vertical","Horizontal"};
}

#endif
