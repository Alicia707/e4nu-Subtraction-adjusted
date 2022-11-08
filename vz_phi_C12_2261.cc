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

void vz_phi_C12_2261()
{
  const int n_bins = 40;
  const double theta_peak_pos = 28;
  double xmin=0.5,xmax=2; 

  gStyle->SetOptStat(0);
  TF1 *f_Gaus = new TF1("f_Gaus", "[0]*TMath::Gaus(x, [1], [2])", xmin, xmax);


  TFile *file_in = new TFile("e2a_ep_C12_2261_neutrino6_united4_photsub_test.root");
  //TFile *file_in1 = new TFile("vz_fits17868.root");        //run17868 for He4, 18522 for He3
  //  TFile *file_in = new TFile("e2a_empty_vert.root");
  // TFile *file_out = new TFile("vz_fits.root", "Recreate");
  TFile *file_out = new TFile("vz_C12_2261.root", "Recreate");

  // TF1 *old_corr_func=(TF1 *)file_in1->Get("f_vz");
  TH2D *h2_el_phi_vert = (TH2D*)file_in->Get("h2_el_phi_vert_uncorr");
  const int n_gr = h2_el_phi_vert->GetNbinsY()/n_bins;
  cout<<"n_gr = "<<n_gr<<endl;

  int n_points=6;
  double phi_interval=60,phi_init_val=30;
  // double x_min[6]={0.85, 1.2, 1.4,1.2 , 0.6, 0.3};    //3He
  //double x_max[6]={1.3, 2, 2.2, 2, 1.4, 1.4};   //3He
  // double x_min[6]={2.1, 2.6, 3.1, 3 , 2.4, 1.8};    //4He
  //double x_max[6]={2.8, 3.3, 3.9, 3.5, 3, 2.5};   //4He
    double x_min[6]={3.5, 3.5, 3.5, 3.5 ,3.5, 3.5};    //C12
  double x_max[6]={6.5, 6.5, 6.5, 6.5, 6.5, 6.5};   //C12
  TH1D *h_proj_vz_[n_points];
  TGraphErrors *gr1 = new TGraphErrors();
  TCanvas *c1=new TCanvas("c1","",500,500);
  int gr_ind = 0;


  for(int i = 0; i < n_points; i++){
    h_proj_vz_[i] = (TH1D*)h2_el_phi_vert->ProjectionX(Form("h_proj_vz_%d", i),h2_el_phi_vert->GetYaxis()->FindBin(phi_init_val+i*phi_interval)-n_gr , h2_el_phi_vert->GetYaxis()->FindBin(phi_init_val+i*phi_interval)+n_gr);
    // h_proj_vz_[i]->SetAxisRange(xmin, xmax, "X");  //use for iron
     h_proj_vz_[i]->Rebin(2);   //for liquid targets
    h_proj_vz_[i]->SetAxisRange(x_min[i], x_max[i], "X");   //for liquid targets
  

    if(h_proj_vz_[i]->GetEntries() > 100){
      
   
      double phi_val = h2_el_phi_vert->GetYaxis()->GetBinCenter(h2_el_phi_vert->GetYaxis()->FindBin(phi_init_val+i*phi_interval));
      double max = h_proj_vz_[i]->GetMaximum();
      double mean = h_proj_vz_[i]->GetBinCenter(h_proj_vz_[i]->GetMaximumBin());
      f_Gaus->SetParameters(max, mean, 0.25);

      h_proj_vz_[i]->Fit(f_Gaus, "MeV", "", mean -0.2, mean + 0.2);

 
      
      double vz_val = f_Gaus->GetParameter(1);
      double vz_val_error = f_Gaus->GetParError(1);
      gr1->SetPoint(gr_ind, phi_val, vz_val);
      gr1->SetPointError(gr_ind, 0, vz_val_error);
      gr_ind = gr_ind + 1;
    }
    
  }
  
  TF1 *f_vz = new TF1("f_vz",Form("[0] + [1]*cos((x - [2])*TMath::DegToRad())/tan(%1.1f*TMath::DegToRad())",theta_peak_pos), 0., 360.);

  gr1->SetMarkerStyle(23);
  gr1->SetMarkerColor(2);
  gr1->GetYaxis()->SetTitleOffset(1.3);
  gr1->GetXaxis()->SetTitle("#phi [Deg.]");
  gr1->GetYaxis()->SetTitle("Z-vert [cm]");
  gr1->Draw("AP");
  
  gr1->Fit(f_vz, "MeV", "", 0., 360.);
  f_vz->Write();
  // old_corr_func->SetLineColor(kBlue);
  //old_corr_func->Draw("Same");
  gr1->SetName("gr1");

 TLatex *lat1 = new TLatex();
 lat1->SetNDC();
 lat1->DrawLatex(0.2, 0.6, Form("#scale[0.8]{F=a+b*cos(#phi-#phi_{0})/tg(%1.0f^{o})}",theta_peak_pos));
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

