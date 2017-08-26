#include <string>
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLatex.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "../Utility/functions.h"
#include "../Utility/draw.h"
#include "../Utility/StSpinAlignmentCons.h"

float const McEtaBinFake[25] = {0.001,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0};
float const McRhoFake[25] = {0.75,0.72501,0.66102,0.58191,0.50832,0.45017,0.40863,0.38076,0.36277,0.35144,0.34441,0.34009,0.33744,0.33583,0.33485,0.33425,0.33389,0.33367,0.33354,0.33346,0.33341,0.33338,0.33336,0.33335,0.33334}; 

void plotMcPhiEtaBoost()
{
  string InPutHist = "/Users/xusun/Data/SpinAlignment/AuAu200GeV/MonteCarlo/McPhiEtaBoost.root";
  TFile *File_InPut = TFile::Open(InPutHist.c_str());

  TProfile *p_Cos2EtaKaon[25];
  TGraphAsymmErrors *g_phiMc = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_phiAna = new TGraphAsymmErrors();
  for(int i_eta = 0; i_eta < 25; ++i_eta)
  {
    string HistName = Form("p_Cos2EtaKaon_%d",i_eta);
    p_Cos2EtaKaon[i_eta] = (TProfile*)File_InPut->Get(HistName.c_str());

    g_phiMc->SetPoint(i_eta,McEtaBinFake[i_eta],p_Cos2EtaKaon[i_eta]->GetBinContent(1));
    g_phiMc->SetPointError(i_eta,0.0,0.0,p_Cos2EtaKaon[i_eta]->GetBinError(1),p_Cos2EtaKaon[i_eta]->GetBinError(1));
    g_phiAna->SetPoint(i_eta,McEtaBinFake[i_eta],McRhoFake[i_eta]);
  }


  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd()->SetLeftMargin(0.15);
  c_play->cd()->SetBottomMargin(0.15);
  c_play->cd()->SetTicks(1,1);
  c_play->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,-0.5,9.5);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetTitle("|#eta| <");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetXaxis()->SetLabelSize(0.04);
  h_play->GetXaxis()->SetNdivisions(505);

  h_play->GetYaxis()->SetTitle("#rho_{00}");
  h_play->GetYaxis()->SetTitleSize(0.04);
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(0.25,0.85);
  h_play->Draw("pE");
  PlotLine(-0.5,9.5,1.0/3.0,1.0/3.0,1,2,2);

  g_phiMc->SetMarkerStyle(24);
  g_phiMc->SetMarkerSize(1.4);
  g_phiMc->SetMarkerColor(4);
  g_phiMc->Draw("pE same");


  g_phiAna->SetMarkerStyle(24);
  g_phiAna->SetMarkerSize(1.4);
  g_phiAna->SetMarkerColor(2);
  g_phiAna->SetLineColor(2);
  g_phiAna->SetLineStyle(2);
  g_phiAna->SetLineWidth(2);
  g_phiAna->Draw("l same");

  TLegend *leg = new TLegend(0.4,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(g_phiMc,"MonteCarlo Simulation","p");
  leg->AddEntry(g_phiAna,"Analytic Calculation","l");
  leg->Draw("same");

  c_play->SaveAs("../figures/c_phiEtaBoost.eps");

  TFile *File_OutPut = new TFile("/Users/xusun/Data/SpinAlignment/AuAu200GeV/MonteCarlo/McRhoBoost.root","RECREATE");
  File_OutPut->cd();
  g_phiMc->SetName("g_phiMc");
  g_phiMc->Write();
  g_phiAna->SetName("g_phiAna");
  g_phiAna->Write();
  File_OutPut->Close();
}
