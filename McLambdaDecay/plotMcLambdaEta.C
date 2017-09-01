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

void plotMcLambdaEta(int energy = 6, int pid = 0)
{
  string PID[2] = {"L","Lbar"};
  // string InPutHist = Form("/Users/xusun/Data/SpinAlignment/AuAu200GeV/MonteCarlo/McLambdaEta_%d.root",pid);
  string InPutHist = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McLambdaEta_%d.root",vmsa::mBeamEnergy[energy].c_str(),pid);
  TFile *File_InPut = TFile::Open(InPutHist.c_str());

  TProfile *p_cosRP, *p_sinRP;
  p_cosRP = (TProfile*)File_InPut->Get("p_cosRP");
  p_sinRP = (TProfile*)File_InPut->Get("p_sinRP");

  TGraphAsymmErrors *g_Dau = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_Lambda = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_DauOnly = new TGraphAsymmErrors();
  TProfile *p_cosInteDau[20], *p_cosInteLambda[20], *p_cosInteDauOnly[20];
  TProfile *p_sinInteDau[20], *p_sinInteLambda[20], *p_sinInteDauOnly[20];
  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    string ProName;
    ProName = Form("p_cosInteDau_%d",i_eta);
    p_cosInteDau[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());
    g_Dau->SetPoint(i_eta,vmsa::McEtaBin[i_eta],p_cosInteDau[i_eta]->GetBinContent(1));
    g_Dau->SetPointError(i_eta,0.0,0.0,p_cosInteDau[i_eta]->GetBinError(1),p_cosInteDau[i_eta]->GetBinError(1));

    ProName = Form("p_cosInteLambda_%d",i_eta);
    p_cosInteLambda[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());
    g_Lambda->SetPoint(i_eta,vmsa::McEtaBin[i_eta],p_cosInteLambda[i_eta]->GetBinContent(1));
    g_Lambda->SetPointError(i_eta,0.0,0.0,p_cosInteLambda[i_eta]->GetBinError(1),p_cosInteLambda[i_eta]->GetBinError(1));

    ProName = Form("p_cosInteDauOnly_%d",i_eta);
    p_cosInteDauOnly[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());
    g_DauOnly->SetPoint(i_eta,vmsa::McEtaBin[i_eta],p_cosInteDauOnly[i_eta]->GetBinContent(1));
    g_DauOnly->SetPointError(i_eta,0.0,0.0,p_cosInteDauOnly[i_eta]->GetBinError(1),p_cosInteDauOnly[i_eta]->GetBinError(1));

    ProName = Form("p_sinInteDau_%d",i_eta);
    p_sinInteDau[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());

    ProName = Form("p_sinInteLambda_%d",i_eta);
    p_sinInteLambda[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());

    ProName = Form("p_sinInteDauOnly_%d",i_eta);
    p_sinInteDauOnly[i_eta] = (TProfile*)File_InPut->Get(ProName.c_str());
  }

  TCanvas *c_PolaPt = new TCanvas("c_PolaPt","c_PolaPt",10,10,800,800);
  c_PolaPt->cd()->SetLeftMargin(0.15);
  c_PolaPt->cd()->SetBottomMargin(0.15);
  c_PolaPt->cd()->SetTicks(1,1);
  c_PolaPt->cd()->SetGrid(0,0);

  p_cosRP->SetTitle("");
  p_cosRP->SetStats(0);
  p_cosRP->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_cosRP->GetXaxis()->CenterTitle();
  p_cosRP->GetXaxis()->SetLabelSize(0.04);
  p_cosRP->GetXaxis()->SetNdivisions(505);

  p_cosRP->GetYaxis()->SetTitle("P_{H}");
  p_cosRP->GetYaxis()->SetTitleSize(0.04);
  p_cosRP->GetYaxis()->CenterTitle();
  p_cosRP->GetYaxis()->SetLabelSize(0.04);
  p_cosRP->GetYaxis()->SetNdivisions(505);
  p_cosRP->GetYaxis()->SetRangeUser(0.00,0.05);
  p_cosRP->SetMarkerStyle(20);
  p_cosRP->SetMarkerSize(1.4);
  p_cosRP->SetMarkerColor(kGray+2);
  p_cosRP->Draw("pE");
  PlotLine(vmsa::ptMin,vmsa::ptMax,0.02,0.02,1,2,2);
  TF1 *f_pol0 = new TF1("f_pol0","pol0",vmsa::ptMin,vmsa::ptMax);
  f_pol0->SetParameter(0,0.02);
  p_cosRP->Fit(f_pol0,"NQ");
  f_pol0->SetLineColor(2);
  f_pol0->SetLineWidth(2);
  f_pol0->SetLineStyle(2);
  f_pol0->Draw("l same");
  plotTopLegend("input P_{H} = 0.02",3.0,0.03,0.03,1,0.0,42,0,1);
  string PHfit = Form("extract P_{H} = %0.4f #pm %0.4f",f_pol0->GetParameter(0),f_pol0->GetParError(0));
  plotTopLegend((char*)PHfit.c_str(),3.0,0.04,0.03,1,0.0,42,0,1);

  string outputPolaPt = Form("../figures/c_PolaPt_%s.eps",PID[pid].c_str());
  c_PolaPt->SaveAs(outputPolaPt.c_str());


  TCanvas *c_PolaEta = new TCanvas("c_PolaEta","c_PolaEta",10,10,800,800);
  c_PolaEta->cd()->SetLeftMargin(0.15);
  c_PolaEta->cd()->SetBottomMargin(0.15);
  c_PolaEta->cd()->SetTicks(1,1);
  c_PolaEta->cd()->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
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

  h_play->GetYaxis()->SetTitle("P_{H}");
  h_play->GetYaxis()->SetTitleSize(0.04);
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetLabelSize(0.04);
  h_play->GetYaxis()->SetNdivisions(505);
  h_play->GetYaxis()->SetRangeUser(0.01,0.03);
  h_play->Draw("pE");
  PlotLine(0.0,10.0,0.0,0.0,1,2,2);

  g_Dau->SetMarkerStyle(24);
  g_Dau->SetMarkerSize(1.4);
  g_Dau->SetMarkerColor(4);
  g_Dau->Draw("pE same");

  g_DauOnly->SetMarkerStyle(26);
  g_DauOnly->SetMarkerSize(1.4);
  g_DauOnly->SetMarkerColor(2);
  g_DauOnly->Draw("pE same");

  g_Lambda->SetMarkerStyle(30);
  g_Lambda->SetMarkerSize(1.4);
  g_Lambda->SetMarkerColor(1);
  g_Lambda->Draw("pE same");

  TLegend *legEta = new TLegend(0.4,0.6,0.8,0.8);
  legEta->SetBorderSize(0);
  legEta->SetFillColor(0);
  legEta->AddEntry(g_Dau,"cut on daughters and #Lambda","p");
  legEta->AddEntry(g_DauOnly,"cut on daughters","p");
  legEta->AddEntry(g_Lambda,"cut on #Lambda only","p");
  legEta->Draw("same");
  string outputPolaEta = Form("../figures/c_PolaEta_%s.eps",PID[pid].c_str());
  c_PolaEta->SaveAs(outputPolaEta.c_str());

  TH3F *h_Eta = (TH3F*)File_InPut->Get("h_Eta");
  TH2F *h_EtaProton = (TH2F*)h_Eta->Project3D("yx")->Clone("h_EtaProton");
  TH2F *h_EtaPion   = (TH2F*)h_Eta->Project3D("zx")->Clone("h_EtaPion");
  TH2F *h_EtaPPion  = (TH2F*)h_Eta->Project3D("zy")->Clone("h_EtaPPion");
  TCanvas *c_EtaPL = new TCanvas("c_EtaPL","c_EtaPL",10,10,800,800);
  c_EtaPL->cd()->SetLeftMargin(0.15);
  c_EtaPL->cd()->SetBottomMargin(0.15);
  c_EtaPL->cd()->SetTicks(1,1);
  c_EtaPL->cd()->SetGrid(0,0);
  h_EtaProton->SetTitle("");
  h_EtaProton->SetStats(0);
  h_EtaProton->GetXaxis()->SetTitle("#eta_{p}");
  h_EtaProton->GetXaxis()->CenterTitle();
  h_EtaProton->GetXaxis()->SetLabelSize(0.04);
  h_EtaProton->GetXaxis()->SetNdivisions(505);

  h_EtaProton->GetYaxis()->SetTitle("#eta_{#Lambda}");
  h_EtaProton->GetYaxis()->SetTitleSize(0.04);
  h_EtaProton->GetYaxis()->CenterTitle();
  h_EtaProton->GetYaxis()->SetLabelSize(0.04);
  h_EtaProton->GetYaxis()->SetNdivisions(505);
  h_EtaProton->Draw("colz");
  string outputEtaPL = Form("../figures/c_EtaPL_%s.eps",PID[pid].c_str());
  c_EtaPL->SaveAs(outputEtaPL.c_str());

  TCanvas *c_EtaPiL = new TCanvas("c_EtaPiL","c_EtaPiL",10,10,800,800);
  c_EtaPiL->cd()->SetLeftMargin(0.15);
  c_EtaPiL->cd()->SetBottomMargin(0.15);
  c_EtaPiL->cd()->SetTicks(1,1);
  c_EtaPiL->cd()->SetGrid(0,0);
  h_EtaPion->SetTitle("");
  h_EtaPion->SetStats(0);
  h_EtaPion->GetXaxis()->SetTitle("#eta_{#pi}");
  h_EtaPion->GetXaxis()->CenterTitle();
  h_EtaPion->GetXaxis()->SetLabelSize(0.04);
  h_EtaPion->GetXaxis()->SetNdivisions(505);

  h_EtaPion->GetYaxis()->SetTitle("#eta_{#Lambda}");
  h_EtaPion->GetYaxis()->SetTitleSize(0.04);
  h_EtaPion->GetYaxis()->CenterTitle();
  h_EtaPion->GetYaxis()->SetLabelSize(0.04);
  h_EtaPion->GetYaxis()->SetNdivisions(505);
  h_EtaPion->Draw("colz");
  string outputEtaPiL = Form("../figures/c_EtaPiL_%s.eps",PID[pid].c_str());
  c_EtaPiL->SaveAs(outputEtaPiL.c_str());

  TCanvas *c_EtaPPi = new TCanvas("c_EtaPPi","c_EtaPPi",10,10,800,800);
  c_EtaPPi->cd()->SetLeftMargin(0.15);
  c_EtaPPi->cd()->SetBottomMargin(0.15);
  c_EtaPPi->cd()->SetTicks(1,1);
  c_EtaPPi->cd()->SetGrid(0,0);
  h_EtaPPion->SetTitle("");
  h_EtaPPion->SetStats(0);
  h_EtaPPion->GetXaxis()->SetTitle("#eta_{p}");
  h_EtaPPion->GetXaxis()->CenterTitle();
  h_EtaPPion->GetXaxis()->SetLabelSize(0.04);
  h_EtaPPion->GetXaxis()->SetNdivisions(505);

  h_EtaPPion->GetYaxis()->SetTitle("#eta_{#pi}");
  h_EtaPPion->GetYaxis()->SetTitleSize(0.04);
  h_EtaPPion->GetYaxis()->CenterTitle();
  h_EtaPPion->GetYaxis()->SetLabelSize(0.04);
  h_EtaPPion->GetYaxis()->SetNdivisions(505);
  h_EtaPPion->Draw("colz");
  string outputEtaPPi = Form("../figures/c_EtaPPi_%s.eps",PID[pid].c_str());
  c_EtaPPi->SaveAs(outputEtaPPi.c_str());
  
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/MonteCarlo/McPH_%d.root",vmsa::mBeamEnergy[energy].c_str(),pid);
  string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPH_%d.root",vmsa::mBeamEnergy[energy].c_str(),pid);
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_play->Write();
  g_Dau->SetName("g_Dau");
  g_Dau->Write();
  g_Lambda->SetName("g_Lambda");
  g_Lambda->Write();
  g_DauOnly->SetName("g_DauOnly");
  g_DauOnly->Write();
  File_OutPut->Close();
}
