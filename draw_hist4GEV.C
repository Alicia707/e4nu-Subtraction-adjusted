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

int main(int argc, char* argv[]){
  if(argc != 2)
  {
    cout << "Incorrect Use of file. Please use ./ProgramName filePath" << endl;
  }
  TFile *file_in;

  file_in = new TFile(argv[1]);

  //Subtraction histograms and all that jazz
  TH1F *h1_E_cal_pimi_sub         = (TH1F*)file_in->Get("h1_E_cal_pimi_sub");
  TH1F *h1_E_cal_pipl_sub         = (TH1F*)file_in->Get("h1_E_cal_pipl_sub");
  TH1F *h1_E_tot_2p1pi_1p1pi_pimi = (TH1F*)file_in->Get("h1_E_tot_2p1pi_1p1pi_pimi");
  TH1F *h1_E_tot_1p2pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_1p2pi_pimi");
  TH1F *h1_E_tot_2p2pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_2p2pi_pimi");
  TH1F *h1_E_tot_1p3pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_1p3pi_pimi");
  TH1F *h1_E_tot_3p1pi_pimi       = (TH1F*)file_in->Get("h1_E_tot_3p1pi_pimi");
  TH1F *h1_E_tot_2p1pi_1p1pi_pipl = (TH1F*)file_in->Get("h1_E_tot_2p1pi_1p1pi_pipl");
  TH1F *h1_E_tot_1p2pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_1p2pi_pipl");
  TH1F *h1_E_tot_2p2pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_2p2pi_pipl");
  TH1F *h1_E_tot_1p3pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_1p3pi_pipl");
  TH1F *h1_E_tot_3p1pi_pipl       = (TH1F*)file_in->Get("h1_E_tot_3p1pi_pipl");

  //1p1pi Histogram stuff
  TH1F *h1_E_cal_1p1pi_pimi_tot   = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pimi_tot");
  TH1F *h1_E_cal_1p1pi_pimi_ONLY  = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pimi_ONLY");
  TH1F *h1_E_cal_1p1pi_pipl_tot   = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pipl_tot");
  TH1F *h1_E_cal_1p1pi_pipl_ONLY  = (TH1F*)file_in->Get("h1_E_cal_1p1pi_pipl_ONLY");

  //variable to be used later ;)
  int min = 0;

  gStyle->SetOptStat(0);

  divideByBinWidth(N_E_bins, h1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, h1_E_cal_pipl_sub);
  divideByBinWidth(N_E_bins, h1_E_tot_2p1pi_1p1pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_1p2pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_1p3pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_3p1pi_pimi);
  divideByBinWidth(N_E_bins, h1_E_tot_2p1pi_1p1pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_1p2pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_2p2pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_1p3pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_tot_3p1pi_pipl);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_tot);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pimi_ONLY);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_tot);
  divideByBinWidth(N_E_bins, h1_E_cal_1p1pi_pipl_ONLY);

  reflectOverX(N_E_bins, h1_E_tot_2p1pi_1p1pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_1p2pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_2p2pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_1p3pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_3p1pi_pimi);
  reflectOverX(N_E_bins, h1_E_tot_2p1pi_1p1pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_1p2pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_2p2pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_1p3pi_pipl);
  reflectOverX(N_E_bins, h1_E_tot_3p1pi_pipl);

  TCanvas *c1 = new TCanvas("c1","",567,370);
  h1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,5.5,"X");
  h1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
  h1_E_cal_pimi_sub->SetLineColor(3); //green
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pimi_tot->Draw("HIST");
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  c1->Print("1p1pi_pimiStuff_4GeV.png");

  TCanvas *c2 = new TCanvas("c2","",567,370);

  h1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,5.5,"X");
  h1_E_cal_1p1pi_pipl_tot->SetLineColor(1);
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2);
  h1_E_cal_pipl_sub->SetLineColor(3);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
  h1_E_cal_1p1pi_pipl_tot->Draw("HIST");
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  c2->Print("1p1pi_piplStuff_4GeV.png");


  TCanvas *c3 = new TCanvas("c3","",567,370);

  h1_E_cal_1p1pi_pimi_tot->SetAxisRange(0,5.5,"X");
  //Make the minimum prettier :-)
  min = h1_E_cal_1p1pi_pimi_tot->GetMaximum()*-.35;
  h1_E_cal_1p1pi_pimi_tot->SetMinimum(min);
  //Set colors and line styles for ALL 1p1pi, ONLY 1p1pi, and subtraction
  h1_E_cal_1p1pi_pimi_tot->SetLineColor(1); //black
  h1_E_cal_1p1pi_pimi_ONLY->SetLineStyle(2);
  h1_E_cal_1p1pi_pimi_ONLY->SetLineColor(2); //red
  h1_E_cal_pimi_sub->SetLineStyle(2);
  h1_E_cal_pimi_sub->SetLineColor(3); //green
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pimi->SetLineColor(4); //blue
  h1_E_tot_1p2pi_pimi->SetLineColor(40); //purple
  h1_E_tot_2p2pi_pimi->SetLineColor(6); //magenta
  h1_E_tot_1p3pi_pimi->SetLineColor(7); //cyan
  h1_E_tot_3p1pi_pimi->SetLineColor(8); //dirt green
  //Make the axes prettier
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pimi_tot->GetXaxis()->SetTitleOffset(0.7);
  //Plot 1p1pi, 1p1pi AT LEAST, and sub
  h1_E_cal_1p1pi_pimi_tot->Draw("HIST");
  h1_E_cal_pimi_sub->Draw("HIST SAME");
  h1_E_cal_1p1pi_pimi_ONLY->Draw("HIST SAME");
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pimi->Draw("HIST SAME");
  h1_E_tot_1p2pi_pimi->Draw("HIST SAME");
  h1_E_tot_2p2pi_pimi->Draw("HIST SAME");
  h1_E_tot_1p3pi_pimi->Draw("HIST SAME");
  h1_E_tot_3p1pi_pimi->Draw("HIST SAME");
  //Create png
  c3->Print("pimiSubtraction_4GeV.png");

  TCanvas *c4 = new TCanvas("c4","",567,370);

  h1_E_cal_1p1pi_pipl_tot->SetAxisRange(0,5.5,"X");
  //Make the minimum prettier
  min = h1_E_cal_1p1pi_pipl_tot->GetMaximum()*-.35;
  h1_E_cal_1p1pi_pipl_tot->SetMinimum(min);
  //Set colors and line styles for ALL 1p1pi, ONLY 1p1pi, and subtraction
  h1_E_cal_1p1pi_pipl_tot->SetLineColor(1); //Black
  h1_E_cal_1p1pi_pipl_ONLY->SetLineStyle(2);  //Dotted Line
  h1_E_cal_1p1pi_pipl_ONLY->SetLineColor(2); //Red
  h1_E_cal_pipl_sub->SetLineStyle(2);
  h1_E_cal_pipl_sub->SetLineColor(3); //green
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pipl->SetLineColor(4); //blue
  h1_E_tot_1p2pi_pipl->SetLineColor(40); //purple
  h1_E_tot_2p2pi_pipl->SetLineColor(6); //magenta
  h1_E_tot_1p3pi_pipl->SetLineColor(7); //cyan
  h1_E_tot_3p1pi_pipl->SetLineColor(8); //dirt green
  //Axes stuff
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitle("E_{cal}");
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetYaxis()->SetLabelSize(0.05);
  h1_E_cal_1p1pi_pipl_tot->GetXaxis()->SetTitleOffset(0.7);
  //Plot 1p1pi, 1p1pi AT LEAST, and sub
  h1_E_cal_1p1pi_pipl_tot->Draw("HIST");
  h1_E_cal_1p1pi_pipl_ONLY->Draw("HIST SAME");
  h1_E_cal_pipl_sub->Draw("HIST SAME");
  //Plot members
  h1_E_tot_2p1pi_1p1pi_pipl->Draw("HIST SAME");
  h1_E_tot_1p2pi_pipl->Draw("HIST SAME");
  h1_E_tot_2p2pi_pipl->Draw("HIST SAME");
  h1_E_tot_1p3pi_pipl->Draw("HIST SAME");
  h1_E_tot_3p1pi_pipl->Draw("HIST SAME");
  //Create png
  c4->Print("piplSubtraction_4GeV.png");

  return 0;

}

void divideByBinWidth(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  for(int i = 1; i<=NumBins; i++){
 	 h1->SetBinContent(i,h1->GetBinContent(i)/h1->GetBinWidth(i));
  }
}

void reflectOverX(double NumBins, TH1F* h1)
{
  NumBins = h1->GetNbinsX();
  for(int i = 1; i<NumBins; i++)
  {
    h1->SetBinContent(i,-(h1->GetBinContent(i)));
  }
}

/*
 - Use subprograms to make plots with more lines on each plot
  - Another program that reads in the histogram and then makes pretty pictures
  - Make sure to label the signs the contributions come in as (added/subtracted)
    - If subtracted -> Plot it as negative for ease of understanding
      - Will be able to see graphically that blank + blank2 sums to subtraction histogram
      -
*/
