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

double N_E_bins;

void divideByBinWidth(double NumBins, TH1F* h1);
void reflectOverX(double NumBins, TH1F* h1);
void normalizeHisto(TH1F* h1);
int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    cout<< "Incorrect Use of file. Please use ./ProgramName target TargetEnergy first number" << endl;
  }

  TFile *data_file;
  TFile *genie_file;

  char* e2;
  *e2 = '2';
  char* e4;
  *e4 = '4';

  //char* dfileName = Form("/u/home/amand/RootWork/GenieAnalysis/data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  char* dfileName = Form("/mnt/c/Users/alici/Documents/Git/WorkingCode/PresentationStuff/drawFunctions/data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  data_file = new TFile(dfileName);

  //char* gfileName = Form("/u/home/amand/RootWork/GenieAnalysis/genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  char* gfileName = Form("/mnt/c/Users/alici/Documents/Git/WorkingCode/PresentationStuff/drawFunctions/genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  genie_file = new TFile(gfileName);

  //Pull in all of our histograms :)
  TH1F *dh1_E_cal_pimi_sub = (TH1F*)data_file->Get("h1_E_cal_pimi_sub");
  TH1F *gh1_E_cal_pimi_sub = (TH1F*)genie_file->Get("h1_E_cal_pimi_sub");
  TH1F *dh1_E_cal_pipl_sub = (TH1F*)data_file->Get("h1_E_cal_pipl_sub");
  TH1F *gh1_E_cal_pipl_sub = (TH1F*)genie_file->Get("h1_E_cal_pipl_sub");
  //1p1pi Histogram stuff (data)
  TH1F *dh1_E_cal_1p1pi_pimi_tot   = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pimi_tot");
  TH1F *dh1_E_cal_1p1pi_pimi_ONLY  = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pimi_ONLY");
  TH1F *dh1_E_cal_1p1pi_pipl_tot   = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pipl_tot");
  TH1F *dh1_E_cal_1p1pi_pipl_ONLY  = (TH1F*)data_file->Get("h1_E_cal_1p1pi_pipl_ONLY");
  //genie
  TH1F *gh1_E_cal_1p1pi_pimi_tot   = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pimi_tot");
  TH1F *gh1_E_cal_1p1pi_pimi_ONLY  = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pimi_ONLY");
  TH1F *gh1_E_cal_1p1pi_pipl_tot   = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pipl_tot");
  TH1F *gh1_E_cal_1p1pi_pipl_ONLY  = (TH1F*)genie_file->Get("h1_E_cal_1p1pi_pipl_ONLY");

  //Divide by bin width for all histograms
  divideByBinWidth(N_E_bins, dh1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pimi_tot);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pimi_ONLY);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pipl_tot);
  divideByBinWidth(N_E_bins, dh1_E_cal_1p1pi_pipl_ONLY);
  //genie stuff
  divideByBinWidth(N_E_bins, gh1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pimi_tot);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pimi_ONLY);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pipl_tot);
  divideByBinWidth(N_E_bins, gh1_E_cal_1p1pi_pipl_ONLY);

  //Normalize all of our histograms
  normalizeHisto(dh1_E_cal_pimi_sub);
  normalizeHisto(dh1_E_cal_1p1pi_pimi_tot);
  normalizeHisto(dh1_E_cal_1p1pi_pimi_ONLY);
  normalizeHisto(dh1_E_cal_1p1pi_pipl_tot);
  normalizeHisto(dh1_E_cal_1p1pi_pipl_ONLY);
  //genie stuff
  normalizeHisto(gh1_E_cal_pimi_sub);
  normalizeHisto(gh1_E_cal_1p1pi_pimi_tot);
  normalizeHisto(gh1_E_cal_1p1pi_pimi_ONLY);
  normalizeHisto(gh1_E_cal_1p1pi_pipl_tot);
  normalizeHisto(gh1_E_cal_1p1pi_pipl_ONLY);

  gStyle->SetOptStat(0);

  if(strcmp(argv[3], e2)==0)
  {
    TCanvas *c1 = new TCanvas("c1","",567,370);
    gh1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,3.5,"X");
    gh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
    gh1_E_cal_pimi_sub->SetLineColor(3); //green
    //data stuff
    dh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
    dh1_E_cal_pimi_sub->SetLineColor(3); //green
    dh1_E_cal_1p1pi_pimi_tot->SetLineStyle(2); //black
    dh1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(2); //red
    dh1_E_cal_pimi_sub->SetLineStyle(2); //green
    //
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pimi_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    gh1_E_cal_pimi_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    dh1_E_cal_pimi_sub->Draw("HIST SAME");
    c1->Print("genie_1p1pi_pimiStuff_2GeV.png");

    TCanvas *c2 = new TCanvas("c2","",567,370);
    gh1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,3.5,"X");
    gh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //black
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //red
    gh1_E_cal_pipl_sub->SetLineColor(3); //green
    //data stuff
    dh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //red
    dh1_E_cal_pipl_sub->SetLineColor(3); //green
    dh1_E_cal_1p1pi_pipl_tot->SetLineStyle(2); //black
    dh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2); //red
    dh1_E_cal_pipl_sub->SetLineStyle(2); //green
    //
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pipl_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    gh1_E_cal_pipl_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    dh1_E_cal_pipl_sub->Draw("HIST SAME");
    c2->Print("genie_1p1pi_piplStuff_2GeV.png");
  }
  else if(strcmp(argv[3], e4)==0)
  {
    TCanvas *c1 = new TCanvas("c1","",567,370);
    gh1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,5.5,"X");
    gh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
    gh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
    gh1_E_cal_pimi_sub->SetLineColor(3); //green
    //data stuff
    dh1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
    dh1_E_cal_pimi_sub->SetLineColor(3); //green
    dh1_E_cal_1p1pi_pimi_tot->SetLineStyle(2); //black
    dh1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(2); //red
    dh1_E_cal_pimi_sub->SetLineStyle(2); //green
    //
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pimi_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    gh1_E_cal_pimi_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
    dh1_E_cal_pimi_sub->Draw("HIST SAME");
    c1->Print("genie_1p1pi_pimiStuff_4GeV.png");

    TCanvas *c2 = new TCanvas("c2","",567,370);
    gh1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,5.5,"X");
    gh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //black
    gh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //red
    gh1_E_cal_pipl_sub->SetLineColor(3); //green
    //data stuff
    dh1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //black
    dh1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //red
    dh1_E_cal_pipl_sub->SetLineColor(3); //green
    dh1_E_cal_1p1pi_pipl_tot->SetLineStyle(2); //black
    dh1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2); //red
    dh1_E_cal_pipl_sub->SetLineStyle(2); //green
    //
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
    gh1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
    gh1_E_cal_1p1pi_pipl_tot->Draw("HIST");
    gh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    gh1_E_cal_pipl_sub->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_tot->Draw("HIST SAME");
    dh1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
    dh1_E_cal_pipl_sub->Draw("HIST SAME");
    c2->Print("genie_1p1pi_piplStuff_4GeV.png");
  }
  else{
    return 0;
  }
  return 0;
}

//Helper functions
//This divides by the bin width
void divideByBinWidth(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  for(int i = 1; i<=NumBins; i++){
 	 h1->SetBinContent(i,h1->GetBinContent(i)/h1->GetBinWidth(i));
  }
}
//This function reflects over the x axis so it looks better
void reflectOverX(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  double binCont;
  for(int i = 1; i<NumBins; i++)
  {
    binCont = h1->GetBinContent(i);
    h1->SetBinContent(i,-1*binCont);
  }
}
//Normalize Histogram function
void normalizeHisto(TH1F* h1)
{
  Double_t factor = 1.;
  h1->Scale(factor/h1->Integral());
}
