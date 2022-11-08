#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGaxis.h>

using namespace std;

void hist_draw()
{
  TFile *file_in;
  TFile *file_in1;

  file_in = new TFile("data_e2a_ep_C12_2261_neutrino6_united4_radphot_test.root");
  file_in1 = new TFile("C12_2261.root");


  TH1F *old_1p2pi = (TH1F*)file_in1->Get("2pi_1p_rot_cal_pimi");
  TH1F *old_2p1pi = (TH1F*)file_in1->Get("1pi_2p_rot_cal_pimi");
  TH1F *old_2p2pi = (TH1F*)file_in1->Get("2pi_2p_rot_cal_pimi");
  TH1F *old_1p3pi = (TH1F*)file_in1->Get("3pi_1p_rot_cal_pimi");
  TH1F *old_3p1pi = (TH1F*)file_in1->Get("1pi_3p_rot_cal_pimi");
  TH1F *old_E_cal_sub = (TH1F*)file_in1->Get("sub_cal_all_pimi");
  TH1F *new_1p2pi = (TH1F*)file_in->Get("h1_E_tot_1p2pi_pimi");
  TH1F *new_2p1pi = (TH1F*)file_in->Get("h1_E_tot_2p1pi_1p1pi_pimi");
  TH1F *new_2p2pi = (TH1F*)file_in->Get("h1_E_tot_2p2pi_pimi");
  TH1F *new_1p3pi = (TH1F*)file_in->Get("h1_E_tot_1p3pi_pimi");
  TH1F *new_3p1pi = (TH1F*)file_in->Get("h1_E_tot_3p1pi_pimi");
  TH1F *new_E_cal_sub = (TH1F*)file_in->Get("h1_E_cal_pimi_sub");

  double N_E_bins = old_1p2pi->GetNbinsX();

  for(int i = 1; i<=N_E_bins;i++){
    old_1p2pi->SetBinContent(i,old_1p2pi->GetBinContent(i)/old_1p2pi->GetBinWidth(i));
  }

  N_E_bins = old_2p1pi->GetNbinsX();

  for(int i = 1; i<=N_E_bins;i++){
    old_2p1pi->SetBinContent(i,old_2p1pi->GetBinContent(i)/old_2p1pi->GetBinWidth(i));
  }

  N_E_bins = old_2p2pi->GetNbinsX();

  for(int i = 1; i<=N_E_bins;i++){
    old_2p2pi->SetBinContent(i,old_2p2pi->GetBinContent(i)/old_2p2pi->GetBinWidth(i));
  }

  N_E_bins = old_1p3pi->GetNbinsX();

  for(int i = 1; i<=N_E_bins;i++){
    old_1p3pi->SetBinContent(i,old_1p3pi->GetBinContent(i)/old_1p3pi->GetBinWidth(i));
  }

  N_E_bins = old_3p1pi->GetNbinsX();

  for(int i = 1; i<=N_E_bins;i++){
    old_3p1pi->SetBinContent(i,old_3p1pi->GetBinContent(i)/old_3p1pi->GetBinWidth(i));
  }

  N_E_bins = old_E_cal_sub->GetNbinsX();

  for(int i = 1; i<=N_E_bins;i++){
    old_E_cal_sub->SetBinContent(i,old_E_cal_sub->GetBinContent(i)/old_E_cal_sub->GetBinWidth(i));
  }

  double N_E_bins1 = new_1p2pi->GetNbinsX();

  for(int j = 1; j<=N_E_bins1;j++){
    new_1p2pi->SetBinContent(j,new_1p2pi->GetBinContent(j)/new_1p2pi->GetBinWidth(j));
  }

  N_E_bins1 = new_2p1pi->GetNbinsX();

  for(int j = 1; j<=N_E_bins1;j++){
    new_2p1pi->SetBinContent(j,new_2p1pi->GetBinContent(j)/new_2p1pi->GetBinWidth(j));
  }

  N_E_bins1 = new_2p2pi->GetNbinsX();

  for(int j = 1; j<=N_E_bins1;j++){
    new_2p2pi->SetBinContent(j,new_2p2pi->GetBinContent(j)/new_2p2pi->GetBinWidth(j));
  }

  N_E_bins1 = new_1p3pi->GetNbinsX();

  for(int j = 1; j<=N_E_bins1;j++){
    new_1p3pi->SetBinContent(j,new_1p3pi->GetBinContent(j)/new_1p3pi->GetBinWidth(j));
  }

  N_E_bins1 = new_3p1pi->GetNbinsX();

  for(int j = 1; j<=N_E_bins1;j++){
    new_3p1pi->SetBinContent(j,new_3p1pi->GetBinContent(j)/new_3p1pi->GetBinWidth(j));
  }

  double N_E_bins2 = new_E_cal_sub->GetNbinsX();

  for(int k = 1; k<=N_E_bins2;k++){
    new_E_cal_sub->SetBinContent(k,new_E_cal_sub->GetBinContent(k)/new_E_cal_sub->GetBinWidth(k));
  }

  TCanvas *c1 = new TCanvas("c1","",567,370);

  old_1p2pi->GetXaxis()->SetTitle("1p2pi Old");
  old_1p2pi->GetXaxis()->SetTitleSize(0.05);
  old_1p2pi->GetXaxis()->SetLabelSize(0.05);
  old_1p2pi->GetYaxis()->SetLabelSize(0.05);
  old_1p2pi->GetXaxis()->SetTitleOffset(0.7);
  old_1p2pi->Draw("HIST");

  c1->Print("/u/home/ltracy/Root_work/e4nu/histograms/old_1p2pi.svg");

  TCanvas *c2 = new TCanvas("c2","",567,370);

  old_2p1pi->GetXaxis()->SetTitle("2p1pi Old");
  old_2p1pi->GetXaxis()->SetTitleSize(0.05);
  old_2p1pi->GetXaxis()->SetLabelSize(0.05);
  old_2p1pi->GetYaxis()->SetLabelSize(0.05);
  old_2p1pi->GetXaxis()->SetTitleOffset(0.7);
  old_2p1pi->Draw("HIST");

  c2->Print("/u/home/ltracy/Root_work/e4nu/histograms/old_2p1pi.svg");

  TCanvas *c3 = new TCanvas("c3","",567,370);

  old_2p2pi->GetXaxis()->SetTitle("2p2pi Old");
  old_2p2pi->GetXaxis()->SetTitleSize(0.05);
  old_2p2pi->GetXaxis()->SetLabelSize(0.05);
  old_2p2pi->GetYaxis()->SetLabelSize(0.05);
  old_2p2pi->GetXaxis()->SetTitleOffset(0.7);
  old_2p2pi->Draw("HIST");

  c3->Print("/u/home/ltracy/Root_work/e4nu/histograms/old_2p2pi.svg");

  TCanvas *c4 = new TCanvas("c4","",567,370);

  old_1p3pi->GetXaxis()->SetTitle("1p3pi Old");
  old_1p3pi->GetXaxis()->SetTitleSize(0.05);
  old_1p3pi->GetXaxis()->SetLabelSize(0.05);
  old_1p3pi->GetYaxis()->SetLabelSize(0.05);
  old_1p3pi->GetXaxis()->SetTitleOffset(0.7);
  old_1p3pi->Draw("HIST");

  c4->Print("/u/home/ltracy/Root_work/e4nu/histograms/old_1p3pi.svg");

  TCanvas *c5 = new TCanvas("c5","",567,370);

  old_3p1pi->GetXaxis()->SetTitle("3p1pi Old");
  old_3p1pi->GetXaxis()->SetTitleSize(0.05);
  old_3p1pi->GetXaxis()->SetLabelSize(0.05);
  old_3p1pi->GetYaxis()->SetLabelSize(0.05);
  old_3p1pi->GetXaxis()->SetTitleOffset(0.7);
  old_3p1pi->Draw("HIST");

  c5->Print("/u/home/ltracy/Root_work/e4nu/histograms/old_3p1pi.svg");

  TCanvas *c6 = new TCanvas("c6","",567,370);

  new_1p2pi->GetXaxis()->SetTitle("1p2pi New");
  new_1p2pi->GetXaxis()->SetTitleSize(0.05);
  new_1p2pi->GetXaxis()->SetLabelSize(0.05);
  new_1p2pi->GetYaxis()->SetLabelSize(0.05);
  new_1p2pi->GetXaxis()->SetTitleOffset(0.7);
  new_1p2pi->Draw("HIST");

  c6->Print("/u/home/ltracy/Root_work/e4nu/histograms/new_1p2pi.svg");

  TCanvas *c7 = new TCanvas("c7","",567,370);

  new_2p1pi->GetXaxis()->SetTitle("2p1pi New");
  new_2p1pi->GetXaxis()->SetTitleSize(0.05);
  new_2p1pi->GetXaxis()->SetLabelSize(0.05);
  new_2p1pi->GetYaxis()->SetLabelSize(0.05);
  new_2p1pi->GetXaxis()->SetTitleOffset(0.7);
  new_2p1pi->Draw("HIST");

  c7->Print("/u/home/ltracy/Root_work/e4nu/histograms/new_2p1pi.svg");

  TCanvas *c8 = new TCanvas("c8","",567,370);

  new_2p2pi->GetXaxis()->SetTitle("2p2pi New");
  new_2p2pi->GetXaxis()->SetTitleSize(0.05);
  new_2p2pi->GetXaxis()->SetLabelSize(0.05);
  new_2p2pi->GetYaxis()->SetLabelSize(0.05);
  new_2p2pi->GetXaxis()->SetTitleOffset(0.7);
  new_2p2pi->Draw("HIST");

  c8->Print("/u/home/ltracy/Root_work/e4nu/histograms/new_2p2pi.svg");

  TCanvas *c9 = new TCanvas("c9","",567,370);

  new_1p3pi->GetXaxis()->SetTitle("1p3pi New");
  new_1p3pi->GetXaxis()->SetTitleSize(0.05);
  new_1p3pi->GetXaxis()->SetLabelSize(0.05);
  new_1p3pi->GetYaxis()->SetLabelSize(0.05);
  new_1p3pi->GetXaxis()->SetTitleOffset(0.7);
  new_1p3pi->Draw("HIST");

  c9->Print("/u/home/ltracy/Root_work/e4nu/histograms/new_1p3pi.svg");

  TCanvas *c10 = new TCanvas("c10","",567,370);

  new_3p1pi->GetXaxis()->SetTitle("3p1pi New");
  new_3p1pi->GetXaxis()->SetTitleSize(0.05);
  new_3p1pi->GetXaxis()->SetLabelSize(0.05);
  new_3p1pi->GetYaxis()->SetLabelSize(0.05);
  new_3p1pi->GetXaxis()->SetTitleOffset(0.7);
  new_3p1pi->Draw("HIST");

  c10->Print("/u/home/ltracy/Root_work/e4nu/histograms/new_3p1pi.svg");

  TCanvas*c11 = new TCanvas("c11","",567,370);

  old_E_cal_sub->GetXaxis()->SetTitle("E_{cal} Old");
  old_E_cal_sub->GetXaxis()->SetTitleSize(0.05);
  old_E_cal_sub->GetXaxis()->SetLabelSize(0.05);
  old_E_cal_sub->GetYaxis()->SetLabelSize(0.05);
  old_E_cal_sub->GetXaxis()->SetTitleOffset(0.7);
  old_E_cal_sub->Draw("HIST");

  c11->Print("/u/home/ltracy/Root_work/e4nu/histograms/CalOld.svg");

  TCanvas *c12 = new TCanvas("c12","",567,370);

  new_E_cal_sub->GetXaxis()->SetTitle("E_{cal} New");
  new_E_cal_sub->GetXaxis()->SetTitleSize(0.05);
  new_E_cal_sub->GetXaxis()->SetLabelSize(0.05);
  new_E_cal_sub->GetYaxis()->SetLabelSize(0.05);
  new_E_cal_sub->GetXaxis()->SetTitleOffset(0.7);
  new_E_cal_sub->Draw("HIST");

  c12->Print("/u/home/ltracy/Root_work/e4nu/histograms/CalNew.svg");
}
