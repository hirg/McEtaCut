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

void plotMcPhiEta(int energy = 6)
{
  // string InPutHist = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/MonteCarlo/McPhiEta.root",vmsa::mBeamEnergy[energy].c_str());
  string InPutHist = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McPhiEta.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_InPut = TFile::Open(InPutHist.c_str());

  TH3F *h_Eta = (TH3F*)File_InPut->Get("h_Eta");
  TH2F *h_eta = (TH2F*)h_Eta->Project3D("yx");
  h_eta->SetTitle("");
  h_eta->SetStats(0);
  h_eta->GetXaxis()->SetTitle("#eta(#phi-meson)");
  h_eta->GetXaxis()->CenterTitle();
  h_eta->GetXaxis()->SetLabelSize(0.04);
  h_eta->GetXaxis()->SetNdivisions(505);

  h_eta->GetYaxis()->SetTitle("#eta(K^{+})");
  h_eta->GetYaxis()->SetTitleSize(0.04);
  h_eta->GetYaxis()->CenterTitle();
  h_eta->GetYaxis()->SetLabelSize(0.04);
  h_eta->GetYaxis()->SetNdivisions(505);
  TCanvas *c_eta = new TCanvas("c_eta","c_eta",10,10,800,800);
  c_eta->cd()->SetLeftMargin(0.15);
  c_eta->cd()->SetBottomMargin(0.15);
  c_eta->cd()->SetTicks(1,1);
  c_eta->cd()->SetGrid(0,0);
  h_eta->Draw("colz");
  c_eta->SaveAs("../figures/c_eta.eps");

  TH2F *h_cosRP = (TH2F*)File_InPut->Get("h_cosRP");
  TH1F *h_Cos = (TH1F*)h_cosRP->ProjectionY();
  TF1 *f_rho = new TF1("f_rho",SpinDensity,-1.0,1.0,2);
  f_rho->SetParameter(0,0.33);
  f_rho->SetParameter(1,100);
  h_Cos->Fit(f_rho,"NQ");
  float rho = f_rho->GetParameter(0);
  float err = f_rho->GetParError(0);
  float Norm = f_rho->GetParameter(1);

  TH2F *h_CosEtaKaon[20], *h_CosEtaPhi[20], *h_CosEtaKOnly[20];
  TH1F *h_CosKaon[20], *h_CosPhi[20], *h_CosKOnly[20];
  float Norm_Kaon[20], Norm_Phi[20], Norm_KOnly[20];
  float rho_Kaon[20], rho_Phi[20], rho_KOnly[20];
  float err_Kaon[20], err_Phi[20], err_KOnly[20];
  TGraphAsymmErrors *g_Kaon = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_Phi = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_KOnly = new TGraphAsymmErrors();
  for(int i_eta = 0; i_eta < 20; ++i_eta)
  {
    string HistName;
    HistName = Form("h_CosEtaKaon_%d",i_eta);
    h_CosEtaKaon[i_eta] = (TH2F*)File_InPut->Get(HistName.c_str());
    HistName = Form("h_CosKaon_%d",i_eta);
    h_CosKaon[i_eta] = (TH1F*)h_CosEtaKaon[i_eta]->ProjectionY(HistName.c_str());
    h_CosKaon[i_eta]->SetTitle("");
    h_CosKaon[i_eta]->SetStats(0);
    h_CosKaon[i_eta]->GetXaxis()->SetTitle("cos(#theta*)");
    h_CosKaon[i_eta]->GetXaxis()->CenterTitle();
    h_CosKaon[i_eta]->GetXaxis()->SetLabelSize(0.04);
    h_CosKaon[i_eta]->GetXaxis()->SetNdivisions(505);

    h_CosKaon[i_eta]->GetYaxis()->SetTitle("Counts");
    h_CosKaon[i_eta]->GetYaxis()->SetTitleSize(0.04);
    h_CosKaon[i_eta]->GetYaxis()->CenterTitle();
    h_CosKaon[i_eta]->GetYaxis()->SetLabelSize(0.04);
    h_CosKaon[i_eta]->GetYaxis()->SetNdivisions(505);
    h_CosKaon[i_eta]->GetYaxis()->SetRangeUser(0.9*h_CosKaon[i_eta]->GetMinimum(),1.1*h_CosKaon[i_eta]->GetMaximum());

    h_CosKaon[i_eta]->SetMarkerStyle(24);
    h_CosKaon[i_eta]->SetMarkerSize(1.1);
    h_CosKaon[i_eta]->SetMarkerColor(2);
    h_CosKaon[i_eta]->SetLineColor(2);
    TF1 *f_rhoKaon = new TF1("f_rhoKaon",SpinDensity,-1.0,1.0,2);
    f_rhoKaon->SetParameter(0,0.33);
    f_rhoKaon->SetParameter(1,100);
    h_CosKaon[i_eta]->Fit(f_rhoKaon,"NQ");
    rho_Kaon[i_eta] = f_rhoKaon->GetParameter(0);
    err_Kaon[i_eta] = f_rhoKaon->GetParError(0);
    Norm_Kaon[i_eta] = f_rhoKaon->GetParameter(1);
    g_Kaon->SetPoint(i_eta,vmsa::McEtaBin[i_eta],rho_Kaon[i_eta]);
    g_Kaon->SetPointError(i_eta,0.0,0.0,err_Kaon[i_eta],err_Kaon[i_eta]);

    HistName = Form("h_CosEtaPhi_%d",i_eta);
    h_CosEtaPhi[i_eta] = (TH2F*)File_InPut->Get(HistName.c_str());
    HistName = Form("h_CosPhi_%d",i_eta);
    h_CosPhi[i_eta] = (TH1F*)h_CosEtaPhi[i_eta]->ProjectionY(HistName.c_str());
    h_CosPhi[i_eta]->SetTitle("");
    h_CosPhi[i_eta]->SetStats(0);
    h_CosPhi[i_eta]->GetXaxis()->SetTitle("cos(#theta*)");
    h_CosPhi[i_eta]->GetXaxis()->CenterTitle();
    h_CosPhi[i_eta]->GetXaxis()->SetLabelSize(0.04);
    h_CosPhi[i_eta]->GetXaxis()->SetNdivisions(505);

    h_CosPhi[i_eta]->GetYaxis()->SetTitle("Counts");
    h_CosPhi[i_eta]->GetYaxis()->SetTitleSize(0.04);
    h_CosPhi[i_eta]->GetYaxis()->CenterTitle();
    h_CosPhi[i_eta]->GetYaxis()->SetLabelSize(0.04);
    h_CosPhi[i_eta]->GetYaxis()->SetNdivisions(505);
    h_CosPhi[i_eta]->GetYaxis()->SetRangeUser(0.9*h_CosPhi[i_eta]->GetMinimum(),1.1*h_CosPhi[i_eta]->GetMaximum());

    h_CosPhi[i_eta]->SetMarkerStyle(24);
    h_CosPhi[i_eta]->SetMarkerSize(1.1);
    h_CosPhi[i_eta]->SetMarkerColor(4);
    h_CosPhi[i_eta]->SetLineColor(4);
    TF1 *f_rhoPhi = new TF1("f_rhoPhi",SpinDensity,-1.0,1.0,2);
    f_rhoPhi->SetParameter(0,0.33);
    f_rhoPhi->SetParameter(1,100);
    h_CosPhi[i_eta]->Fit(f_rhoPhi,"NQ");
    rho_Phi[i_eta] = f_rhoPhi->GetParameter(0);
    err_Phi[i_eta] = f_rhoPhi->GetParError(0);
    Norm_Phi[i_eta] = f_rhoPhi->GetParameter(1);
    g_Phi->SetPoint(i_eta,vmsa::McEtaBin[i_eta],rho_Phi[i_eta]);
    g_Phi->SetPointError(i_eta,0.0,0.0,err_Phi[i_eta],err_Phi[i_eta]);

    HistName = Form("h_CosEtaKOnly_%d",i_eta);
    h_CosEtaKOnly[i_eta] = (TH2F*)File_InPut->Get(HistName.c_str());
    HistName = Form("h_CosKOnly_%d",i_eta);
    h_CosKOnly[i_eta] = (TH1F*)h_CosEtaKOnly[i_eta]->ProjectionY(HistName.c_str());
    h_CosKOnly[i_eta]->SetTitle("");
    h_CosKOnly[i_eta]->SetStats(0);
    h_CosKOnly[i_eta]->GetXaxis()->SetTitle("cos(#theta*)");
    h_CosKOnly[i_eta]->GetXaxis()->CenterTitle();
    h_CosKOnly[i_eta]->GetXaxis()->SetLabelSize(0.04);
    h_CosKOnly[i_eta]->GetXaxis()->SetNdivisions(505);

    h_CosKOnly[i_eta]->GetYaxis()->SetTitle("Counts");
    h_CosKOnly[i_eta]->GetYaxis()->SetTitleSize(0.04);
    h_CosKOnly[i_eta]->GetYaxis()->CenterTitle();
    h_CosKOnly[i_eta]->GetYaxis()->SetLabelSize(0.04);
    h_CosKOnly[i_eta]->GetYaxis()->SetNdivisions(505);
    h_CosKOnly[i_eta]->GetYaxis()->SetRangeUser(0.9*h_CosKOnly[i_eta]->GetMinimum(),1.1*h_CosKOnly[i_eta]->GetMaximum());

    h_CosKOnly[i_eta]->SetMarkerStyle(24);
    h_CosKOnly[i_eta]->SetMarkerSize(1.1);
    h_CosKOnly[i_eta]->SetMarkerColor(2);
    h_CosKOnly[i_eta]->SetLineColor(2);
    TF1 *f_rhoKOnly = new TF1("f_rhoKOnly",SpinDensity,-1.0,1.0,2);
    f_rhoKOnly->SetParameter(0,0.33);
    f_rhoKOnly->SetParameter(1,100);
    h_CosKOnly[i_eta]->Fit(f_rhoKOnly,"NQ");
    rho_KOnly[i_eta] = f_rhoKOnly->GetParameter(0);
    err_KOnly[i_eta] = f_rhoKOnly->GetParError(0);
    Norm_KOnly[i_eta] = f_rhoKOnly->GetParameter(1);
    g_KOnly->SetPoint(i_eta,vmsa::McEtaBin[i_eta],rho_KOnly[i_eta]);
    g_KOnly->SetPointError(i_eta,0.0,0.0,err_KOnly[i_eta],err_KOnly[i_eta]);
  }
  g_Kaon->SetPoint(20,8.0,rho);
  g_Kaon->SetPointError(20,0.0,0.0,err,err);
  g_Phi->SetPoint(20,8.0,rho);
  g_Phi->SetPointError(20,0.0,0.0,err,err);
  g_KOnly->SetPoint(20,8.0,rho);
  g_KOnly->SetPointError(20,0.0,0.0,err,err);

  TCanvas *c_EtaCutEff = new TCanvas("c_EtaCutEff","c_EtaCutEff",10,10,1600,800);
  c_EtaCutEff->Divide(2,1);
  for(int i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_EtaCutEff->cd(i_pad+1)->SetLeftMargin(0.15);
    c_EtaCutEff->cd(i_pad+1)->SetBottomMargin(0.15);
    c_EtaCutEff->cd(i_pad+1)->SetTicks(1,1);
    c_EtaCutEff->cd(i_pad+1)->SetGrid(0,0);
  }
  c_EtaCutEff->cd(1);
  h_CosKaon[10]->Draw("pE");
  TF1 *f_rhoKaonQA = new TF1("f_rhoKaonQA",SpinDensity,-1.0,1.0,2);
  f_rhoKaonQA->SetParameter(0,rho_Kaon[10]);
  f_rhoKaonQA->SetParameter(1,Norm_Kaon[10]);
  f_rhoKaonQA->SetLineColor(2);
  f_rhoKaonQA->SetLineWidth(2);
  f_rhoKaonQA->SetLineStyle(1);
  f_rhoKaonQA->Draw("l same");
  plotTopLegend("cut on K and #phi-meson",0.2,0.80,0.04,1,0.0,42,1);
  plotTopLegend("|#eta| < 1.0",0.2,0.75,0.04,1,0.0,42,1);
  string printRhoKaon = Form("#rho_{00} = %1.3f #pm %1.3f",rho_Kaon[10],err_Kaon[10]);
  plotTopLegend((char*)printRhoKaon.c_str(),0.2,0.70,0.04,1,0.0,42,1);

  c_EtaCutEff->cd(2);
  h_CosPhi[10]->Draw("pE");
  TF1 *f_rhoPhiQA = new TF1("f_rhoPhiQA",SpinDensity,-1.0,1.0,2);
  f_rhoPhiQA->SetParameter(0,rho_Phi[10]);
  f_rhoPhiQA->SetParameter(1,Norm_Phi[10]);
  f_rhoPhiQA->SetLineColor(4);
  f_rhoPhiQA->SetLineWidth(2);
  f_rhoPhiQA->SetLineStyle(1);
  f_rhoPhiQA->Draw("l same");
  plotTopLegend("cut on #phi-meson only",0.2,0.80,0.04,1,0.0,42,1);
  plotTopLegend("|#eta| < 1.0",0.2,0.75,0.04,1,0.0,42,1);
  string printRhoPhi = Form("#rho_{00} = %1.3f #pm %1.3f",rho_Phi[10],err_Phi[10]);
  plotTopLegend((char*)printRhoPhi.c_str(),0.2,0.70,0.04,1,0.0,42,1);
  c_EtaCutEff->SaveAs("../figures/c_phiEtaCutEff.eps");

  TCanvas *c_rhoEta = new TCanvas("c_rhoEta","c_rhoEta",10,10,800,800);
  c_rhoEta->cd()->SetLeftMargin(0.15);
  c_rhoEta->cd()->SetBottomMargin(0.15);
  c_rhoEta->cd()->SetTicks(1,1);
  c_rhoEta->cd()->SetGrid(0,0);

  TH1F *h_rho = new TH1F("h_rho","h_rho",100,0.0,10.0);
  for(int i_bin = 1; i_bin < 101; ++i_bin)
  {
    h_rho->SetBinContent(i_bin,-10.0);
    h_rho->SetBinError(i_bin,1.0);
  }
  h_rho->SetTitle("");
  h_rho->SetStats(0);
  h_rho->GetXaxis()->SetTitle("|#eta| <");
  h_rho->GetXaxis()->CenterTitle();
  h_rho->GetXaxis()->SetRangeUser(0.0,10.0);
  h_rho->GetYaxis()->SetTitle("#rho_{00}");
  h_rho->GetYaxis()->SetTitleSize(0.04);
  h_rho->GetXaxis()->SetNdivisions(505);

  h_rho->GetYaxis()->CenterTitle();
  h_rho->GetYaxis()->SetRangeUser(0.25,0.7);
  h_rho->GetYaxis()->SetNdivisions(505);
  h_rho->Draw("pE");
  PlotLine(0.0,10.0,1.0/3.0,1.0/3.0,1,2,2);

  g_Kaon->SetMarkerStyle(24);
  g_Kaon->SetMarkerColor(4);
  g_Kaon->SetMarkerSize(1.4);
  g_Kaon->Draw("pE same");

  g_Phi->SetMarkerStyle(30);
  g_Phi->SetMarkerColor(1);
  g_Phi->SetMarkerSize(1.4);
  g_Phi->Draw("pE same");

  g_KOnly->SetMarkerStyle(26);
  g_KOnly->SetMarkerColor(2);
  g_KOnly->SetMarkerSize(1.4);
  g_KOnly->Draw("pE same");

  TLegend *leg = new TLegend(0.4,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(g_Kaon,"cut on K and #phi-meson","p");
  leg->AddEntry(g_KOnly,"cut on K only","p");
  leg->AddEntry(g_Phi,"cut on #phi-meson only","p");
  leg->Draw("same");
  c_rhoEta->SaveAs("../figures/c_phiRhoEta.eps");
  
  // string outputfile = Form("/Users/xusun/Data/SpinAlignment/AuAu%s/MonteCarlo/McRho.root",vmsa::mBeamEnergy[energy].c_str());
  string outputfile = Form("/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/SpinAlignment/Phi/MonteCarlo/McRho.root",vmsa::mBeamEnergy[energy].c_str());
  TFile *File_OutPut = new TFile(outputfile.c_str(),"RECREATE");
  File_OutPut->cd();
  h_rho->Write();
  g_Kaon->SetName("g_Kaon");
  g_Kaon->Write();
  g_Phi->SetName("g_Phi");
  g_Phi->Write();
  g_KOnly->SetName("g_KOnly");
  g_KOnly->Write();
  File_OutPut->Close();
}
