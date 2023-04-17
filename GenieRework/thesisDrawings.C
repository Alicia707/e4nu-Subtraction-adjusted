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
double findMaximum(TH1F* h1, TH1F *h2);
int main(int argc, char* argv[])
{
  if(argc != 4)
  {
    cout<< "Incorrect Use of file. Please use ./thesisDrawings target [TargetEnergy] [NRG first number]" << endl;
  }
  TFile *data_file;
  TFile *genie_file;
  char* e2 = new char;
  *e2 = '2';
  char* e4 = new char;
  *e4 = '4';
  char* dfileName = Form("/u/home/amand/RootWork/adjustedSubtraction/data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  //char* dfileName = Form("/mnt/c/Users/alici/Documents/Git/WorkingCode/PresentationStuff/drawFunctions/data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  data_file = new TFile(dfileName);

  char* gfileName = Form("/u/home/amand/RootWork/adjustedSubtraction/genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  //char* gfileName = Form("/mnt/c/Users/alici/Documents/Git/WorkingCode/PresentationStuff/drawFunctions/genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root", argv[1], argv[2]);
  genie_file = new TFile(gfileName);

  //Pull in all of our histograms :)
  TH1F *dh1_E_cal_pimi_sub = (TH1F*)data_file->Get("h1_E_cal_pimi_sub_thesis");
  TH1F *gh1_E_cal_pimi_sub = (TH1F*)genie_file->Get("h1_E_cal_pimi_sub_thesis");
  TH1F *dh1_E_cal_pipl_sub = (TH1F*)data_file->Get("h1_E_cal_pipl_sub_thesis");
  TH1F *gh1_E_cal_pipl_sub = (TH1F*)genie_file->Get("h1_E_cal_pipl_sub_thesis");


  //Divide by bin width for all histograms
  divideByBinWidth(N_E_bins, dh1_E_cal_pimi_sub); 
  divideByBinWidth(N_E_bins, dh1_E_cal_pipl_sub);
  //genie stuff
  divideByBinWidth(N_E_bins, gh1_E_cal_pimi_sub);
  divideByBinWidth(N_E_bins, gh1_E_cal_pipl_sub);
  //Normalize all of our histograms
  normalizeHisto(dh1_E_cal_pimi_sub); 
  normalizeHisto(dh1_E_cal_pipl_sub);
  //genie stuff
  normalizeHisto(gh1_E_cal_pimi_sub);
  normalizeHisto(gh1_E_cal_pipl_sub);

  gStyle->SetOptStat(0);

  //Set Labels
  gh1_E_cal_pimi_sub->GetXaxis()->SetTitle("E [GeV]");
  gh1_E_cal_pimi_sub->GetYaxis()->SetTitle("Fraction of Total Events");
  gh1_E_cal_pimi_sub->GetXaxis()->SetLabelSize(0.05);
  gh1_E_cal_pimi_sub->GetYaxis()->SetLabelSize(0.05);
  gh1_E_cal_pimi_sub->GetXaxis()->SetTitleSize(0.05);
  gh1_E_cal_pimi_sub->GetYaxis()->SetTitleSize(0.05);
  gh1_E_cal_pimi_sub->GetXaxis()->SetTitleOffset(0.95);
  gh1_E_cal_pimi_sub->GetYaxis()->SetTitleOffset(1.);
  gh1_E_cal_pimi_sub->GetYaxis()->CenterTitle(true);

  gh1_E_cal_pipl_sub->GetXaxis()->SetTitle("E [GeV]");
  gh1_E_cal_pipl_sub->GetYaxis()->SetTitle("Fraction of Total Events");
  gh1_E_cal_pipl_sub->GetXaxis()->SetLabelSize(0.05);
  gh1_E_cal_pipl_sub->GetYaxis()->SetLabelSize(0.05);
  gh1_E_cal_pipl_sub->GetXaxis()->SetTitleSize(0.05);
  gh1_E_cal_pipl_sub->GetYaxis()->SetTitleSize(0.05);
  gh1_E_cal_pipl_sub->GetXaxis()->SetTitleOffset(0.95);
  gh1_E_cal_pipl_sub->GetYaxis()->SetTitleOffset(1.);
  gh1_E_cal_pipl_sub->GetYaxis()->CenterTitle(true);

  dh1_E_cal_pimi_sub->GetXaxis()->SetTitle("E [GeV]");
  dh1_E_cal_pimi_sub->GetYaxis()->SetTitle("Fraction of Total Events");
  dh1_E_cal_pimi_sub->GetXaxis()->SetLabelSize(0.05);
  dh1_E_cal_pimi_sub->GetYaxis()->SetLabelSize(0.05);
  dh1_E_cal_pimi_sub->GetXaxis()->SetTitleSize(0.05);
  dh1_E_cal_pimi_sub->GetYaxis()->SetTitleSize(0.05);
  dh1_E_cal_pimi_sub->GetXaxis()->SetTitleOffset(0.95);
  dh1_E_cal_pimi_sub->GetYaxis()->SetTitleOffset(1.);
  dh1_E_cal_pimi_sub->GetYaxis()->CenterTitle(true);

  dh1_E_cal_pipl_sub->GetXaxis()->SetTitle("E [GeV]");
  dh1_E_cal_pipl_sub->GetYaxis()->SetTitle("Fraction of Total Events");
  dh1_E_cal_pipl_sub->GetXaxis()->SetLabelSize(0.05);
  dh1_E_cal_pipl_sub->GetYaxis()->SetLabelSize(0.05);
  dh1_E_cal_pipl_sub->GetXaxis()->SetTitleSize(0.05);
  dh1_E_cal_pipl_sub->GetYaxis()->SetTitleSize(0.05);
  dh1_E_cal_pipl_sub->GetXaxis()->SetTitleOffset(0.95);
  dh1_E_cal_pipl_sub->GetYaxis()->SetTitleOffset(1.);
  dh1_E_cal_pipl_sub->GetYaxis()->CenterTitle(true);

  char *label = new char; 
  if(strcmp(argv[1], "C12")==0){label = "^{12}C";}
  else if(strcmp(argv[1], "56Fe")==0){label = "^{56}Fe";}
  else {label = "Pain";}

  if(strcmp(argv[3], "2")== 0) // If it is 2.2 GeV
  {
    cout<< "Drawing...\n";
    //Setup bounds on all of our histograms 
    dh1_E_cal_pimi_sub->SetAxisRange(0, 2.5, "X");
    dh1_E_cal_pipl_sub->SetAxisRange(0, 2.5, "X");
    gh1_E_cal_pimi_sub->SetAxisRange(0, 2.5, "X");
    gh1_E_cal_pipl_sub->SetAxisRange(0, 2.5, "X");
    dh1_E_cal_pimi_sub->SetMinimum(0.0);
    dh1_E_cal_pipl_sub->SetMinimum(0.0);
    gh1_E_cal_pimi_sub->SetMinimum(0.0);
    gh1_E_cal_pipl_sub->SetMinimum(0.0);
    dh1_E_cal_pimi_sub->SetLineWidth(2);
    dh1_E_cal_pipl_sub->SetLineWidth(2);
    gh1_E_cal_pimi_sub->SetLineWidth(2);
    gh1_E_cal_pipl_sub->SetLineWidth(2);

    //Change genie color to put on diff plot 
    gh1_E_cal_pimi_sub->SetLineColor(2);
    gh1_E_cal_pipl_sub->SetLineColor(2);

    double localMax = 0.0;
    char* localLabel = new char;
    localLabel = "";

    //data pimi sub drawing 
    localLabel = "#pi^{-}";
    TCanvas *c1 = new TCanvas("c1","",567,370);
    localMax = dh1_E_cal_pimi_sub->GetMaximum();
    dh1_E_cal_pimi_sub->Draw("HIST");
    TLatex text(.1, localMax*.8, Form("%s %s at 2.2GeV",label, localLabel));
    text.SetTextSize(.08);
    text.DrawClone();
    c1->Print(Form("data_%s_1p1pi_sub_pimi_thesis_2GeV.png", argv[1]));

    localLabel = "#pi^{+}";
    TCanvas *c2 = new TCanvas("c2","",567,370);
    localMax = dh1_E_cal_pipl_sub->GetMaximum();
    dh1_E_cal_pipl_sub->Draw("HIST");
    TLatex text2(.1, localMax*.8, Form("%s %s at 2.2GeV",label, localLabel));
    text2.SetTextSize(.08);
    text2.DrawClone();
    c2->Print(Form("data_%s_1p1pi_sub_pipl_thesis_2GeV.png", argv[1]));

    localLabel = "#pi^{-}";
    TCanvas *c3 = new TCanvas("c3","",567,370);
    localMax = gh1_E_cal_pimi_sub->GetMaximum();
    gh1_E_cal_pimi_sub->Draw("HIST");
    TLatex text3(.1, localMax*.8, Form("%s %s at 2.2GeV",label, localLabel));
    text3.SetTextSize(.08);
    text3.DrawClone();
    c3->Print(Form("genie_%s_1p1pi_sub_pimi_thesis_2GeV.png", argv[1]));

    localLabel = "#pi^{+}";
    TCanvas *c4 = new TCanvas("c4","",567,370);
    localMax = gh1_E_cal_pipl_sub->GetMaximum();
    gh1_E_cal_pipl_sub->Draw("HIST");
    TLatex text4(.1, localMax*.8, Form("%s %s at 2.2GeV",label, localLabel));
    text4.SetTextSize(.08);
    text4.DrawClone();
    c4->Print(Form("genie_%s_1p1pi_sub_pipl_thesis_2GeV.png", argv[1]));

    //Adjust the maximum so that the plot doesn't look wonky
    double pimiMax = findMaximum(gh1_E_cal_pimi_sub, dh1_E_cal_pimi_sub);
    double piplMax = findMaximum(gh1_E_cal_pipl_sub, dh1_E_cal_pipl_sub);
    gh1_E_cal_pimi_sub->SetMaximum(pimiMax);
    gh1_E_cal_pipl_sub->SetMaximum(piplMax);

    localLabel = "#pi^{-}";
    TCanvas *c5 = new TCanvas("c5", "", 567,370);
    gh1_E_cal_pimi_sub->Draw("HIST");
    dh1_E_cal_pimi_sub->Draw("HIST SAME");
    TLatex text5(.1, pimiMax*.8, Form("%s %s at 2.2GeV",label, localLabel));
    text5.SetTextSize(.08);
    text5.DrawClone();
    TLatex textg1(.1, pimiMax*.2, "#splitline{#color[4]{CLAS Data}}{#color[2]{GENIE}}");
    textg1.SetTextSize(.07);
    textg1.DrawClone();
    c5->Print(Form("combined_%s_1p1pi_sub_pimi_thesis_2GeV.png", argv[1]));

    localLabel = "#pi^{+}";
    TCanvas *c6 = new TCanvas("c6", "", 567,370);
    gh1_E_cal_pipl_sub->Draw("HIST");
    dh1_E_cal_pipl_sub->Draw("HIST SAME");
    TLatex text6(.1, piplMax*.8, Form("%s %s at 2.2GeV",label, localLabel));
    text6.SetTextSize(.08);
    text6.DrawClone();
    TLatex textg2(.1, piplMax*.2, "#splitline{#color[4]{CLAS Data}}{#color[2]{GENIE}}");
    textg2.SetTextSize(.07);
    textg2.DrawClone();
    c6->Print(Form("combined_%s_1p1pi_sub_pipl_thesis_2GeV.png", argv[1]));
  }
  else if(strcmp(argv[3], "4")==0)
  {
    //Setup bounds on all of our histograms 
    dh1_E_cal_pimi_sub->SetAxisRange(0, 5.0, "X");
    dh1_E_cal_pipl_sub->SetAxisRange(0, 5.0, "X");
    gh1_E_cal_pimi_sub->SetAxisRange(0, 5.0, "X");
    gh1_E_cal_pipl_sub->SetAxisRange(0, 5.0, "X");
    dh1_E_cal_pimi_sub->SetMinimum(0.0);
    dh1_E_cal_pipl_sub->SetMinimum(0.0);
    gh1_E_cal_pimi_sub->SetMinimum(0.0);
    gh1_E_cal_pipl_sub->SetMinimum(0.0);
    dh1_E_cal_pimi_sub->SetLineWidth(2);
    dh1_E_cal_pipl_sub->SetLineWidth(2);
    gh1_E_cal_pimi_sub->SetLineWidth(2);
    gh1_E_cal_pipl_sub->SetLineWidth(2);

    //Change genie color to put on diff plot 
    gh1_E_cal_pimi_sub->SetLineColor(2);
    gh1_E_cal_pipl_sub->SetLineColor(2);

    double localMax = 0.0;
    char* localLabel = new char; 
    localLabel =  "";
    //data pimi sub drawing 
    localLabel = "#pi^{-}";
    TCanvas *c1 = new TCanvas("c1","",567,370);
    localMax = dh1_E_cal_pimi_sub->GetMaximum();
    dh1_E_cal_pimi_sub->Draw("HIST");
    TLatex text(.2, localMax*.8, Form("%s %s at 4.4GeV",label, localLabel));
    text.SetTextSize(.08);
    text.DrawClone();
    c1->Print(Form("data_%s_1p1pi_sub_pimi_thesis_4GeV.png", argv[1]));

    localLabel = "#pi^{+}";
    TCanvas *c2 = new TCanvas("c2","",567,370);
    localMax = dh1_E_cal_pipl_sub->GetMaximum();
    dh1_E_cal_pipl_sub->Draw("HIST");
    TLatex text2(.2, localMax*.8, Form("%s %s at 4.4GeV",label, localLabel));
    text2.SetTextSize(.08);
    text2.DrawClone();
    c2->Print(Form("data_%s_1p1pi_sub_pipl_thesis_4GeV.png", argv[1]));

    localLabel = "#pi^{-}";
    TCanvas *c3 = new TCanvas("c3","",567,370);
    localMax = gh1_E_cal_pimi_sub->GetMaximum();
    gh1_E_cal_pimi_sub->Draw("HIST");
    TLatex text3(.2, localMax*.8, Form("%s %s at 4.4GeV",label, localLabel));
    text3.SetTextSize(.08);
    text3.DrawClone();
    c3->Print(Form("genie_%s_1p1pi_sub_pimi_thesis_4GeV.png", argv[1]));

    localLabel = "#pi^{+}";
    TCanvas *c4 = new TCanvas("c4","",567,370);
    localMax = gh1_E_cal_pipl_sub->GetMaximum();
    gh1_E_cal_pipl_sub->Draw("HIST");
    TLatex text4(.2, localMax*.8, Form("%s %s at 4.4GeV",label, localLabel));
    text4.SetTextSize(.08);
    text4.DrawClone();
    c4->Print(Form("genie_%s_1p1pi_sub_pipl_thesis_4GeV.png", argv[1]));

    //Adjust the maximum so that the plot doesn't look wonky
    double pimiMax = findMaximum(gh1_E_cal_pimi_sub, dh1_E_cal_pimi_sub);
    double piplMax = findMaximum(gh1_E_cal_pipl_sub, dh1_E_cal_pipl_sub);
    gh1_E_cal_pimi_sub->SetMaximum(pimiMax);
    gh1_E_cal_pipl_sub->SetMaximum(piplMax);

    localLabel = "#pi^{-}";
    TCanvas *c5 = new TCanvas("c5", "", 567,370);
    gh1_E_cal_pimi_sub->Draw("HIST");
    dh1_E_cal_pimi_sub->Draw("HIST SAME");
    TLatex text5(.2, pimiMax*.8, Form("%s %s at 4.4GeV",label, localLabel));
    text5.SetTextSize(.08);
    text5.DrawClone();
    TLatex textg1(.2, pimiMax*.2, "#splitline{#color[4]{CLAS Data}}{#color[2]{GENIE}}");
    textg1.SetTextSize(.07);
    textg1.DrawClone();
    c5->Print(Form("combined_%s_1p1pi_sub_pimi_thesis_4GeV.png", argv[1]));

    localLabel = "#pi^{+}";
    TCanvas *c6 = new TCanvas("c6", "", 567,370);
    gh1_E_cal_pipl_sub->Draw("HIST");
    dh1_E_cal_pipl_sub->Draw("HIST SAME");
    TLatex text6(.2, piplMax*.8, Form("%s %s at 4.4GeV",label, localLabel));
    text6.SetTextSize(.08);
    text6.DrawClone();
    TLatex textg2(.2, piplMax*.2, "#splitline{#color[4]{CLAS Data}}{#color[2]{GENIE}}");
    textg2.SetTextSize(.07);
    textg2.DrawClone();
    c6->Print(Form("combined_%s_1p1pi_sub_pipl_thesis_4GeV.png", argv[1]));
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

//Find Maximum of histos 

double findMaximum(TH1F *h1, TH1F *h2) 
{
  double val1 = h1->GetMaximum(); 
  double val2 = h2->GetMaximum();
  if(val1 > val2) {return val1*1.2;}
  else {return val2 *1.2;}
}