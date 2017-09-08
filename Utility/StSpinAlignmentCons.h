#ifndef StSpinAlignmentCons_h
#define StSpinAlignmentCons_h

#include <string>
#include "TString.h"
// #include "StarClassLibrary/SystemOfUnits.h"

namespace vmsa
{
  // shared constant
  int const NumBeamEnergy = 7;
  std::string const mBeamEnergy[NumBeamEnergy] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  float const mEnergyValue[NumBeamEnergy] = {7.7,11.5,19.6,27.0,39.0,62.4,200.0};
  int const mBeamYear[NumBeamEnergy] = {2010,2010,2011,2011,2010,2010,2011};

  std::string const mPID[3]   = {"Phi","KStar","K0S"};

  float const ptEffMax = 8.0;
  float const ptMin = 0.0;
  float const ptMax = 5.0;
  int const BinPt  = 50;
  int const BinEta = 10;
  int const BinY = 20;
  int const BinPhi = 24;
  int const BinCos = 7;

  // used in McPhiResCorr
  double const acceptanceRapidity = 1.0;
  float const InvMass[3] = {1.01946,0.89594,0.49761}; // 0: phi, 1: K*, 2 K0S
  int const decayChannelsFirst[3] = {656,617,613};
  int const decayChannelsSecond[3] = {666,619,614};
  int const decayMother[3] = {333,313,310};
  int const decayChannels[3] = {656,617,613}; // 0: phi->K+K-, 1: K*->Kpi, 2 K0S->pi+pi-
  float const McEtaBin[20] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,2.0,2.5,3.0,3.5,4.0};
}

#endif
