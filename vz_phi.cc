#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>


using namespace std;

void vz_phi()
{
  double theta_peak=28;
  const int n_bins = 40;
  gStyle->SetOptFit(1);
  TF1 *f_Gaus = new TF1("f_Gaus", "[0]*TMath::Gaus(x, [1], [2])", -6., 3.);


  TFile *file_in = new TFile("e2a_ep_3He_2261_neutrino6_united2_test2.root");
  //  TFile *file_in = new TFile("e2a_empty_vert.root");
  // TFile *file_out = new TFile("vz_fits.root", "Recreate");
TFile *file_out = new TFile("vz_3He_2261_2ndgroup.root", "Recreate");

  TH2D *h2_el_phi_vert = (TH2D*)file_in->Get("h2_el_phi_vert_uncorr");
  const int n_gr = h2_el_phi_vert->GetNbinsY()/n_bins;
  cout<<"n_gr = "<<n_gr<<endl;

  TH1D *h_proj_vz_[n_bins];
  TGraphErrors *gr1 = new TGraphErrors();
  TCanvas *c1=new TCanvas("c1","",500,500);
  int gr_ind = 0;

  for( int i = 0; i < n_bins; i++ ){
    h_proj_vz_[i] = (TH1D*)h2_el_phi_vert->ProjectionX(Form("h_proj_vz_%d", i), n_gr*i + 1, n_gr*(i+1) );
    h_proj_vz_[i]->SetAxisRange(-5., -1., "X");
    
    if(h_proj_vz_[i]->GetEntries() > 100){
      
      double phi_val = h2_el_phi_vert->GetYaxis()->GetBinCenter(n_gr*i + n_gr/2);
      double max = h_proj_vz_[i]->GetMaximum();
      double mean = h_proj_vz_[i]->GetBinCenter(h_proj_vz_[i]->GetMaximumBin());
      f_Gaus->SetParameters(max, mean, 0.3);
      h_proj_vz_[i]->Fit(f_Gaus, "MeV", "", mean -0.5, mean + 0.5);
      
      double vz_val = f_Gaus->GetParameter(1);
 double vz_val_error = f_Gaus->GetParError(1);
      gr1->SetPoint(gr_ind, phi_val, vz_val);
   gr1->SetPointError(gr_ind, 0, vz_val_error);
      gr_ind = gr_ind + 1;
    }
    
  }
  
  TF1 *f_vz = new TF1("f_vz", Form("[0] + [1]*cos((x - [2])*TMath::DegToRad())/tan(%1.1f*TMath::DegToRad())",theta_peak), 0., 360.);

  gr1->SetMarkerStyle(23);
  gr1->SetMarkerColor(2);
  gr1->GetXaxis()->SetTitle("#phi [Deg.]");
  gr1->GetYaxis()->SetTitle("Z-vert [cm]");
  gr1->Draw("AP");
  
  gr1->Fit(f_vz, "MeV", "", 0., 360.);
  f_vz->Write();
  gr1->SetName("gr1");

 TLatex *lat1 = new TLatex();
 lat1->SetNDC();
 lat1->DrawLatex(0.2, 0.6, Form("#scale[0.8]{F=a+b*cos(#phi-#phi_{0})/tg(%1.0f^{o})}",theta_peak));
 lat1->DrawLatex(0.2, 0.54, Form("#scale[0.8]{a=%1.2fcm}",f_vz->GetParameter(0)));
 lat1->DrawLatex(0.4, 0.54, Form("#scale[0.8]{b=%1.2fcm}",f_vz->GetParameter(1)));
 lat1->DrawLatex(0.6, 0.54, Form("#scale[0.8]{#phi_{0}=%1.2f^{o}}",f_vz->GetParameter(2)));



 c1->Write();
  c1->Print("Plots/vert_corr.pdf");

  TCanvas *c2= new TCanvas("c2","",500,500);
  h2_el_phi_vert->Draw("colz");
  h2_el_phi_vert->GetXaxis()->SetTitle("Z vert[cm]");
  h2_el_phi_vert->GetYaxis()->SetTitle("#phi [Deg.]");
  h2_el_phi_vert->GetYaxis()->SetTitleOffset(1.3);
  c2->Print("Plots/el_phi_vert.pdf");


  gDirectory->Write();
}



