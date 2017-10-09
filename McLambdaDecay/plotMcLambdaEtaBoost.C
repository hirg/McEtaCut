#include <string>
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLegend.h"
#include "TF1.h"
#include "../Utility/StSpinAlignmentCons.h"
#include "../Utility/draw.h"

float const McEtaBinFake[17] = {0.001,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0};	
float const McSinFake[17] = {1.2732,1.2604,1.2263,1.1815,1.1363,1.0972,1.0666,1.0443,1.0289,1.0185,1.0117,1.0073,1.0046,1.0028,1.0017,1.0011,1.0007};
float const McCosFake[17] = {1.5000,1.4700,1.3932,1.2983,1.21  ,1.1402,1.0904,1.0569,1.0353,1.0217,1.0133,1.0081,1.0049,1.0030,1.0018,1.0011,1.0007};

void plotMcLambdaEtaBoost(int energy = 6, int pid = 0)
{

  string PID[2] = {"L","Lbar"};
  // string InPutHist = Form("/Users/xusun/Data/SpinAlignment/AuAu200GeV/MonteCarlo/McLambdaEtaBoost_%d.root",pid);
  string InPutHist = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McLambdaEtaBoost_%d.root",vmsa::mBeamEnergy[energy].c_str(),pid);
  TFile *File_InPut = TFile::Open(InPutHist.c_str());

  TProfile *p_cosRP, *p_sinRP;
  p_cosRP = (TProfile*)File_InPut->Get("p_cosRP");
  p_sinRP = (TProfile*)File_InPut->Get("p_sinRP");

  float pHPhy = 0.01;
  TGraphAsymmErrors *g_LambdaMcSin = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_LambdaMcCos = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_LambdaAnaSin = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_LambdaAnaCos = new TGraphAsymmErrors();
  TProfile *p_cosInteDau[17];
  TProfile *p_sinInteDau[17];
  for(int i_eta = 0; i_eta < 17; ++i_eta)
  {
    string ProName;
    ProName = Form("p_cosInteDau_%d",i_eta);
    p_cosInteDau[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());
    g_LambdaMcCos->SetPoint(i_eta,McEtaBinFake[i_eta],p_cosInteDau[i_eta]->GetBinContent(1));
    g_LambdaMcCos->SetPointError(i_eta,0.0,0.0,p_cosInteDau[i_eta]->GetBinError(1),p_cosInteDau[i_eta]->GetBinError(1));
    g_LambdaAnaCos->SetPoint(i_eta,McEtaBinFake[i_eta],McCosFake[i_eta]*pHPhy);

    ProName = Form("p_sinInteDau_%d",i_eta);
    p_sinInteDau[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());
    g_LambdaMcSin->SetPoint(i_eta,McEtaBinFake[i_eta],p_sinInteDau[i_eta]->GetBinContent(1));
    g_LambdaMcSin->SetPointError(i_eta,0.0,0.0,p_sinInteDau[i_eta]->GetBinError(1),p_sinInteDau[i_eta]->GetBinError(1));

    g_LambdaAnaSin->SetPoint(i_eta,McEtaBinFake[i_eta],McSinFake[i_eta]*pHPhy);
  }

  TCanvas *c_PolaPt = new TCanvas("c_PolaPt","c_PolaPt",10,10,800,800);
  c_PolaPt->cd()->SetLeftMargin(0.15);
  c_PolaPt->cd()->SetBottomMargin(0.15);
  c_PolaPt->cd()->SetTicks(1,1);
  c_PolaPt->cd()->SetGrid(0,0);

  p_cosRP->SetTitle("");
  p_cosRP->SetStats(0);
  p_cosRP->GetXaxis()->SetTitle("no #eta cut");
  p_cosRP->GetXaxis()->CenterTitle();
  p_cosRP->GetXaxis()->SetLabelSize(0.00);
  p_cosRP->GetXaxis()->SetNdivisions(500);

  p_cosRP->GetYaxis()->SetTitle("P_{H}");
  p_cosRP->GetYaxis()->SetTitleSize(0.04);
  p_cosRP->GetYaxis()->CenterTitle();
  p_cosRP->GetYaxis()->SetLabelSize(0.04);
  p_cosRP->GetYaxis()->SetNdivisions(505);
  p_cosRP->GetYaxis()->SetRangeUser(pHPhy-0.01,pHPhy+0.03);
  p_cosRP->SetMarkerStyle(20);
  p_cosRP->SetMarkerSize(1.4);
  p_cosRP->SetMarkerColor(kGray+2);
  p_cosRP->Draw("pE");
  p_sinRP->SetMarkerStyle(24);
  p_sinRP->SetMarkerSize(1.4);
  p_sinRP->SetMarkerColor(2);
  p_sinRP->Draw("pE same");
  TLegend *leg = new TLegend(0.4,0.5,0.8,0.7);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(p_cosRP,"#frac{3}{#alpha_{H}}<cos(#theta*)>","p");
  leg->AddEntry(p_sinRP,"#frac{#pi}{8#alpha_{H}}<sin(#Psi-#phi_{p}*)>","p");
  leg->Draw("same");

  string outputPolaPt = Form("../figures/c_PolaPt_%s.eps",PID[pid].c_str());
  c_PolaPt->SaveAs(outputPolaPt.c_str());

  TCanvas *c_PolaCom = new TCanvas("c_PolaCom","c_PolaCom",10,10,800,800);
  c_PolaCom->cd()->SetLeftMargin(0.15);
  c_PolaCom->cd()->SetBottomMargin(0.15);
  c_PolaCom->cd()->SetTicks(1,1);
  c_PolaCom->cd()->SetGrid(0,0);
  TH1F *h_playCom = new TH1F("h_playCom","h_playCom",100,-0.5,9.5);
  for(int i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_playCom->SetBinContent(i_bin+1,-10.0);
    h_playCom->SetBinError(i_bin+1,1.0);
  }
  h_playCom->SetTitle("");
  h_playCom->SetStats(0);
  h_playCom->GetXaxis()->SetTitle("|#eta| <");
  h_playCom->GetXaxis()->CenterTitle();
  h_playCom->GetXaxis()->SetLabelSize(0.04);
  h_playCom->GetXaxis()->SetNdivisions(505);
  h_playCom->GetXaxis()->SetRangeUser(-0.05,4.1);

  h_playCom->GetYaxis()->SetTitle("P_{H}");
  h_playCom->GetYaxis()->SetTitleSize(0.04);
  h_playCom->GetYaxis()->CenterTitle();
  h_playCom->GetYaxis()->SetLabelSize(0.04);
  h_playCom->GetYaxis()->SetNdivisions(505);
  h_playCom->GetYaxis()->SetRangeUser(pHPhy-0.01,pHPhy+0.03);
  h_playCom->Draw("pE");

  for(int i_eta = 0; i_eta < 17; ++i_eta)
  {
    p_cosInteDau[i_eta]->SetMarkerStyle(20);
    p_cosInteDau[i_eta]->SetMarkerSize(1.4);
    p_cosInteDau[i_eta]->SetMarkerColor(kGray+2);
    p_cosInteDau[i_eta]->Draw("pE same");

    p_sinInteDau[i_eta]->SetMarkerStyle(24);
    p_sinInteDau[i_eta]->SetMarkerSize(1.4);
    p_sinInteDau[i_eta]->SetMarkerColor(2);
    p_sinInteDau[i_eta]->Draw("pE same");
  }

  TLegend *legEta = new TLegend(0.4,0.6,0.8,0.8);
  legEta->SetBorderSize(0);
  legEta->SetFillColor(0);
  legEta->AddEntry(p_cosInteDau[0],"#frac{3}{#alpha_{H}}<cos(#theta*)>","p");
  legEta->AddEntry(p_sinInteDau[0],"#frac{#pi}{8#alpha_{H}}<sin(#Psi-#phi_{p}*)>","p");
  legEta->Draw("same");
  PlotLine(-0.05,4.1,pHPhy,pHPhy,1,2,2);
  string outputPolaCom = Form("../figures/c_PolaCom_%s.eps",PID[pid].c_str());
  c_PolaCom->SaveAs(outputPolaCom.c_str());

  TCanvas *c_PolaEta = new TCanvas("c_PolaEta","c_PolaEta",10,10,800,800);
  c_PolaEta->cd()->SetLeftMargin(0.15);
  c_PolaEta->cd()->SetBottomMargin(0.15);
  c_PolaEta->cd()->SetTicks(1,1);
  c_PolaEta->cd()->SetGrid(0,0);
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
  h_play->GetXaxis()->SetRangeUser(-0.05,4.1);

  h_play->GetYaxis()->SetTitle("P_{H}");
  h_play->GetYaxis()->SetTitleSize(0.04);
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(pHPhy-0.01,pHPhy+0.03);
  h_play->Draw("pE");
  PlotLine(-0.05,4.1,pHPhy,pHPhy,1,2,2);

  g_LambdaMcCos->SetMarkerStyle(20);
  g_LambdaMcCos->SetMarkerSize(1.4);
  g_LambdaMcCos->SetMarkerColor(kGray+2);
  g_LambdaMcCos->Draw("pE same");

  g_LambdaAnaCos->SetMarkerStyle(24);
  g_LambdaAnaCos->SetMarkerSize(1.4);
  g_LambdaAnaCos->SetMarkerColor(kGray+2);
  g_LambdaAnaCos->SetLineColor(kGray+2);
  g_LambdaAnaCos->SetLineWidth(2);
  g_LambdaAnaCos->SetLineStyle(2);
  g_LambdaAnaCos->Draw("l same");

  g_LambdaMcSin->SetMarkerStyle(20);
  g_LambdaMcSin->SetMarkerSize(1.4);
  g_LambdaMcSin->SetMarkerColor(2);
  g_LambdaMcSin->Draw("pE same");

  g_LambdaAnaSin->SetMarkerStyle(24);
  g_LambdaAnaSin->SetMarkerSize(1.4);
  g_LambdaAnaSin->SetMarkerColor(2);
  g_LambdaAnaSin->SetLineColor(2);
  g_LambdaAnaSin->SetLineWidth(2);
  g_LambdaAnaSin->SetLineStyle(2);
  g_LambdaAnaSin->Draw("l same");

  TLegend *legEta = new TLegend(0.4,0.6,0.8,0.8);
  legEta->SetBorderSize(0);
  legEta->SetFillColor(0);
  legEta->AddEntry(g_LambdaMcCos,"MC cos","p");
  legEta->AddEntry(g_LambdaAnaCos,"Analytic cos","l");
  legEta->AddEntry(g_LambdaMcSin,"MC sin","p");
  legEta->AddEntry(g_LambdaAnaSin,"Analytic sin","l");
  legEta->Draw("same");
  string outputPolaEta = Form("../figures/c_PolaEta_%s.eps",PID[pid].c_str());
  c_PolaEta->SaveAs(outputPolaEta.c_str());

  // string OutPutFile = Form("/Users/xusun/Data/SpinAlignment/AuAu200GeV/MonteCarlo/McPHBoost_%d.root,pid");
  string OutPutFile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPHBoost_%d.root",vmsa::mBeamEnergy[energy].c_str(),pid);
  TFile *File_OutPut = new TFile(OutPutFile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_play->Write();
  g_LambdaMcCos->SetName("g_LambdaMcCos");
  g_LambdaMcCos->Write();
  g_LambdaAnaCos->SetName("g_LambdaAnaCos");
  g_LambdaAnaCos->Write();
  g_LambdaMcSin->SetName("g_LambdaMcSin");
  g_LambdaMcSin->Write();
  g_LambdaAnaSin->SetName("g_LambdaAnaSin");
  g_LambdaAnaSin->Write();
  File_OutPut->Close();
}
