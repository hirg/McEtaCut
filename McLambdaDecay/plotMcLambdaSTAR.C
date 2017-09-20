#include <string>
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/draw.h"

void plotMcLambdaSTAR(int pid = 0)
{
  TGraphAsymmErrors *g_pH = new TGraphAsymmErrors();
  int i_point = 0;
  for(int i_energy = 0; i_energy < 7; ++i_energy)
  {
    if(i_energy == 0 || i_energy == 5 || i_energy == 6)
    {
      // string InPutHist = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McLambdaEta_%d.root",vmsa::mBeamEnergy[energy].c_str(),pid);
      string InPutHist = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/MonteCarlo/McLambdaEta_%d.root",vmsa::mBeamEnergy[i_energy].c_str(),pid);
      TFile *File_InPut = TFile::Open(InPutHist.c_str());
      TProfile *p_cosSTAR = (TProfile*)File_InPut->Get("p_cosSTAR");
      float pH = p_cosSTAR->GetBinContent(1);
      float err = p_cosSTAR->GetBinError(1);
      g_pH->SetPoint(i_point,vmsa::mEnergyValue[i_energy],pH);
      g_pH->SetPointError(i_point,0.0,0.0,err,err);
      i_point++;
    }
  }

  TCanvas *c_pH = new TCanvas("c_pH","c_pH",10,10,800,800);
  c_pH->cd()->SetLeftMargin(0.15);
  c_pH->cd()->SetBottomMargin(0.15);
  c_pH->cd()->SetTicks(1,1);
  c_pH->cd()->SetGrid(0,0);

  TH1F *h_play = new TH1F("h_play","h_play",1000,0.0,1000.0);
  for(int i_bin = 0; i_bin < 1000; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetXaxis()->SetNdivisions(505);
  h_play->GetXaxis()->SetRangeUser(0.1,210);

  h_play->GetYaxis()->SetTitle("P_{H}");
  h_play->GetYaxis()->SetTitleSize(0.04);
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(-0.005,0.005);
  h_play->Draw("pE");
  PlotLine(0.1,210.0,0.0,0.0,1,2,2);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_pH,24,2,1.4);
  c_pH->SaveAs("../figures/c_pH.eps");
}
