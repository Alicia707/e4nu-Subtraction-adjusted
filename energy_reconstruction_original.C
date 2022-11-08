#define energy_reconstruction_cxx
#include "energy_reconstruction.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH3.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace std;

string fbeam_E, target_name;
map<std::string,double>en_beam;
map<std::string,double>vert_min;
map<std::string,double>vert_max;
map<std::string,double>vertdiff_min;
map<std::string,double>vertdiff_max;
map<std::string,double>EC_photon_beta;
map<std::string,double>LEC_photon_beta;
map<std::pair<std::string, int>, double> EC_time_offset;
map<std::string,double>bind_en;
map<std::string,double>target_mass;
map<std::string,double>residual_target_mass;
map<std::string, double> Ecal_offset;

void SetFiducialCutParameters();
Float_t ProtonMomCorrection_He3_4Cell(std::string atarget, TLorentzVector V4Pr, Float_t vertex_p );
TLorentzVector EMomentumCorrection(TLorentzVector V4el);
TVector3 FindUVW(TVector3 xyz);
Bool_t CutUVW(TVector3 ecxyz);
Bool_t GetEPhiLimits(Float_t momentum, Float_t theta, Int_t sector,Float_t *EPhiMin, Float_t *EPhiMax);
Bool_t PiplFiducialCut(TVector3 momentum,Float_t *philow,Float_t *phiup);
Bool_t PimiFiducialCut(TVector3 momentum, Float_t *pimi_philow, Float_t *pimi_phiup);
Bool_t EFiducialCut(TVector3 momentum);
Bool_t PFiducialCut(TVector3 momentum);
bool Phot_fid(TVector3 V3_phot);
void rot_2pi_1p(TVector3 V3_pi[2], double qpi[2], TVector3 V3_prot, TVector3 V3q, double N1pi1p[2], double *N2pi1p, int N_tot);
void rot_1pi_2p(TVector3 V3_pi, double qpi, TVector3 V3_prot[2], TVector3 V3q, double N1pi1p[2], double *N1pi2p, int N_tot);
void rot_1pi_3p(TVector3 V3_pi, double qpi, TVector3 V3_prot[3], TVector3 V3q, double N1pi1p[3], double N1pi2p[3], double *N1pi3p, int N_tot);
void rot_3pi_1p(TVector3 V3_pi[3], double qpi[3], TVector3 V3_prot, TVector3 V3q, double N1pi1p[3], double N2pi1p[3], double *N3pi1p, int N_tot);
void rot_1phot_1pi_1p(TVector3 V3_phot, TVector3 V3_pi, double qpi, TVector3 V3_prot, TVector3 V3q, bool radstat, double *N1pi1p0phot, double *N1pi1p1phot, int N_tot);
void rot_1phot_1pi_2p(TVector3 V3_phot, TVector3 V3_pi, double qpi, TVector3 V3_prot[2], TVector3 V3q, bool radstat, double N1pi1p0phot[2], double N1pi1p1phot[2], double *N1pi2p0phot, double *N1pi2p1phot, int N_tot);
void rot_2phot_1pi_1p(TVector3 V3_phot[2], TVector3 V3_pi, double qpi, TVector3 V3_prot, TVector3 V3q, bool radstat[2], double *N1pi1p0phot, double *N1pi1p1phot, double *N1pi1p2phot, int N_tot);
void rot_1phot_2pi_1p(TVector3 V3_phot, TVector3 V3_pi[2], double qpi[2], TVector3 V3_prot, TVector3 V3q, bool radstat, double N1pi1p0phot[2], double N1pi1p1phot[2], double *N2pi1p0phot, double *N2pi1p1phot, int N_tot);
void rot_1phot_1pi_3p(TVector3 V3_phot, TVector3 V3_pi, double qpi, TVector3 V3_prot[3], TVector3 V3q, bool radstat, double N1pi1p0phot[3], double N1pi1p1phot[3], double N1pi2p0phot[3], double N1pi2p1phot[3], double *N1pi3p1phot, double *N1pi3p0phot, int N_tot);
double vz_corr(double phi,double theta);

TF1* vz_corr_func;
TF1 *up_lim1_ec, *up_lim2_ec,*up_lim3_ec,*up_lim4_ec, *up_lim5_ec,*up_lim6_ec,*low_lim1_ec,*low_lim2_ec,*low_lim3_ec, *low_lim4_ec,*low_lim5_ec,*low_lim6_ec;
TF1 *leftside_lim1_ec, *leftside_lim2_ec,*leftside_lim3_ec, *leftside_lim4_ec,*leftside_lim5_ec, *leftside_lim6_ec,*rightside_lim1_ec, *rightside_lim2_ec,*rightside_lim3_ec, *rightside_lim4_ec,*rightside_lim5_ec, *rightside_lim6_ec;
TF1 *pipl_deltat_sig,*pipl_deltat_mean,*pimi_deltat_sig,*pimi_deltat_mean, *fsum_pimi,*fsub_pimi, *fsum_pipl,*fsub_pipl, *prot_deltat_sig, *prot_deltat_mean,*fsum_prot,*fsub_prot,*el_Epratio_sig,*el_Epratio_mean,*fsum_e,*fsub_e;


//e- E_tot/p vs p PID cut
Double_t FSum_e(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return el_Epratio_mean->EvalPar(x)+par[0]*el_Epratio_sig->EvalPar(x);
  if(x[0]>=par[1]) return el_Epratio_mean->Eval(par[1])+par[0]*el_Epratio_sig->Eval(par[1]);
}

Double_t FSub_e(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return el_Epratio_mean->EvalPar(x)-par[0]*el_Epratio_sig->EvalPar(x);
  if(x[0]>=par[1]) return el_Epratio_mean->Eval(par[1])-par[0]*el_Epratio_sig->Eval(par[1]);
}

//proton Delta_t vs momentum PID cut

Double_t FSum_prot(Double_t *x, Double_t *par){   //the 2 parameters are the cut range and momentum limit
  if(x[0] < par[1])    return prot_deltat_mean->EvalPar(x)+par[0]*prot_deltat_sig->EvalPar(x);
  if(x[0] >= par[1])   return prot_deltat_mean->Eval(par[1])+par[0]*prot_deltat_sig->Eval(par[1]);
}
Double_t FSub_prot(Double_t *x,Double_t *par){
  if(x[0] < par[1])    return prot_deltat_mean->EvalPar(x)-par[0]*prot_deltat_sig->EvalPar(x);
  if(x[0] >= par[1])   return prot_deltat_mean->Eval(par[1])-par[0]*prot_deltat_sig->Eval(par[1]);
}

//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions

Double_t FSum_pimi(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pimi_deltat_mean->EvalPar(x)+par[0]*pimi_deltat_sig->EvalPar(x);
  if(x[0]>=par[1])return pimi_deltat_mean->Eval(par[1])+par[0]*pimi_deltat_sig->Eval(par[1]);
}
Double_t FSub_pimi(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pimi_deltat_mean->EvalPar(x)-par[0]*pimi_deltat_sig->EvalPar(x);
  if(x[0]>=par[1])return pimi_deltat_mean->Eval(par[1])-par[0]*pimi_deltat_sig->Eval(par[1]);
}

//To Draw two sigma pid cuts lines on Delta t vs p distribution of negative pions

Double_t FSum_pipl(Double_t *x,Double_t *par){
  if(x[0]<par[1])  return pipl_deltat_mean->EvalPar(x)+par[0]*pipl_deltat_sig->EvalPar(x);
  if(x[0]>=par[1]) return pipl_deltat_mean->Eval(par[1])+par[0]*pipl_deltat_sig->Eval(par[1]);
}
Double_t FSub_pipl(Double_t *x,Double_t *par){
  if(x[0]<par[1]) return pipl_deltat_mean->EvalPar(x)-par[0]*pipl_deltat_sig->EvalPar(x);
  if(x[0]>=par[1])return pipl_deltat_mean->Eval(par[1])-par[0]*pipl_deltat_sig->Eval(par[1]);
}

int fTorusCurrent;
bool SCpdcut=true;
Float_t fgPar_4Gev_2250_Phi[6][3],fgPar_4Gev_2250_Theta[6][4],fgPar_4Gev_2250_Efid_t0_p[6][2];// 4GeV e- fiducial cut parameters
Float_t fgPar_4Gev_2250_Efid_t1_p[6][6],fgPar_4Gev_2250_Efid_b_p[6][2][6],fgPar_4Gev_2250_Efid_a_p[6][2][6];
double fgPar_Efid_Theta_S5_extra[8][4],fgPar_Efid_Theta_S4_extra[2][4],fgPar_Efid_Theta_S3_extra[4][4],fgPar_Efid_Theta_S5[8][8],fgPar_Efid_Theta_S4[2][8], fgPar_Efid_Theta_S3[4][8];
double fgPar_Pfid_ScpdS2[2][6],fgPar_Pfid_ScpdS3[8][6], fgPar_Pfid_ScpdS4[4][6],fgPar_Pfid_ScpdS5[8][6], fgPar_Pfid_ScpdS2_extra[2][4],fgPar_Pfid_ScpdS3_extra[8][4], fgPar_Pfid_ScpdS4_extra[4][4],fgPar_Pfid_ScpdS5_extra[8][4], fgPar_4Gev_2250_Pfidft1l[6][6],fgPar_4Gev_2250_Pfidft1r[6][6], fgPar_4Gev_2250_Pfidft2l[6][6], fgPar_4Gev_2250_Pfidft2r[6][6],fgPar_4Gev_2250_Pfidbt1l[6][6],fgPar_4Gev_2250_Pfidbt1r[6][6],fgPar_4Gev_2250_Pfidbt2l[6][6], fgPar_4Gev_2250_Pfidbt2r[6][6], fgPar_4Gev_2250_Pfidbl[6][6], fgPar_4Gev_2250_Pfidbr[6][6];  // 4GeV proton fiducial cut parameters
double   fgPar_1gev_1500_Pfid[6][5][6],fgPar_1gev_750_Pfid[6][5][6];   // 1.1 GeV proton fiducial cut parameters
 double fgPar_1gev_1500_Pfid_ScpdS2[2][6],fgPar_1gev_1500_Pfid_ScpdS3[8][6],fgPar_1gev_1500_Pfid_ScpdS4[4][6],fgPar_1gev_1500_Pfid_ScpdS5[8][6];
 double fgPar_1gev_750_Pfid_ScpdS2[2][6],fgPar_1gev_750_Pfid_ScpdS3[8][6],fgPar_1gev_750_Pfid_ScpdS4[4][6],fgPar_1gev_750_Pfid_ScpdS5[8][6];
double   fgPar_1gev_1500_Piplfid[6][5][6],fgPar_1gev_750_Piplfid[6][5][6];
 double fgPar_1gev_1500_Piplfid_ScpdS2[2][6],fgPar_1gev_1500_Piplfid_ScpdS3[8][6],fgPar_1gev_1500_Piplfid_ScpdS4[4][6],fgPar_1gev_1500_Piplfid_ScpdS5[8][6];
 double fgPar_1gev_750_Piplfid_ScpdS2[2][6],fgPar_1gev_750_Piplfid_ScpdS3[8][6],fgPar_1gev_750_Piplfid_ScpdS4[4][6],fgPar_1gev_750_Piplfid_ScpdS5[8][6];
 double fgPar_1gev_1500_Efid[6][5][6],fgPar_1gev_1500_Efid_Theta_S3[4][8],fgPar_1gev_1500_Efid_Theta_S4[2][8],fgPar_1gev_1500_Efid_Theta_S5[8][8];  // 1.1 GeV e- fiducial cut parameters
 double fgPar_1gev_750_Efid[6][5][6],fgPar_1gev_750_Efid_Theta_S3[4][8],fgPar_1gev_750_Efid_Theta_S4[2][8],fgPar_1gev_750_Efid_Theta_S5[8][8];

 double fgPar_1gev_1500_Pimfid[6][5][6],fgPar_1gev_1500_Pimfid_Theta_S3[4][8],fgPar_1gev_1500_Pimfid_Theta_S4[2][8],fgPar_1gev_1500_Pimfid_Theta_S5[8][8];
double fgPar_1gev_750_Pimfid[6][5][6],fgPar_1gev_750_Pimfid_Theta_S3[4][8],fgPar_1gev_750_Pimfid_Theta_S4[2][8],fgPar_1gev_750_Pimfid_Theta_S5[8][8],fgPar_1gev_750_Pimfid_Theta_S3_extra[4][4],fgPar_1gev_750_Pimfid_Theta_S4_extra[2][4],fgPar_1gev_750_Pimfid_Theta_S5_extra[8][4],fgPar_1gev_1500_Pimfid_Theta_S3_extra[4][4],fgPar_1gev_1500_Pimfid_Theta_S4_extra[2][4],fgPar_1gev_1500_Pimfid_Theta_S5_extra[8][4];
const int n_slice=3,nsect=6;
double prot_accept_mom_lim,chargedpi_accept_mom_lim,phot_accept_en_lim;
double prot_mom_lim;
double min_good_mom;
double max_mom;
double pipl_maxmom, pimi_maxmom,pimi_delt_cutrange,pipl_delt_cutrange;
int N_pperp,N_Ecal;
double *pperp_cut,*Ecal_lowlim,*Ecal_uplim;
double   Q2_mean,Q2_step;
Double_t elmom_corr_fact[nsect];
const double Mpi = 0.139570, e_mass=0.000510998;
const double m_prot= 0.9382720813,m_neut=0.939565,H3_bind_en=0.008481,He4_bind_en=0.0283,C12_bind_en=0.09215, B_bind_en=0.0762,He3_bind_en=0.0077,D2_bind_en=0.00222,Fe_bind_en=0.49226,Mn_bind_en=0.4820764;
const double  mH3 = 2*m_neut+ m_prot-H3_bind_en;
const double  mHe4 = 2*m_neut+ 2*m_prot-He4_bind_en;
double ece, el_sccc_timediff;
int el_segment ,el_cc_sector;
const double eps = .02;
const double m_delta = 1.232;
double m_a_minus = (5*m_prot + 6*m_neut - B_bind_en);
double m_C12 = (6*m_prot + 6*m_neut - C12_bind_en);
int N_tot = 2000;
TLorentzVector V4_target(0,0,0,6*m_prot+6*m_neut-C12_bind_en);
const Double_t c = 2.99792E+10;
Float_t cphil = 0;
Float_t cphir = 0;
Double_t ns_to_s = 1.0E-9;
double bett = 0;
double deltt = 0;
double hist_max = 0;

void energy_reconstruction::Loop()
{
  if (fChain == 0) return;
  en_beam["2261"]=2.261;
  en_beam["4461"]=4.461;
  fbeam_E=fbeam_en;
  TLorentzVector V4_beam(0,0,en_beam[fbeam_E],en_beam[fbeam_E]);

  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.)
    {
      hist_max = 5;
      prot_accept_mom_lim=0.3;
      chargedpi_accept_mom_lim=0.15;
      phot_accept_en_lim=0.3;
      prot_mom_lim=2.15;
      min_good_mom=0.55;
      max_mom=2.1;
      pipl_maxmom=1.4;
      pimi_maxmom=1.3;
      pimi_delt_cutrange=3.;
      pipl_delt_cutrange=3.;
      N_pperp=2,N_Ecal=6;
      pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;  pperp_cut[1]=0.2;
      for (int i=0;i<N_Ecal;i++){
	Ecal_lowlim[i]=0.75+i*0.25;
	Ecal_uplim[i]=1.+i*0.25;
      }
      Ecal_lowlim[5]=0.;
      Ecal_uplim[5]=2.;
      for (int i=0;i<N_Ecal;i++)	cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<endl;
      Q2_mean=0.6;
      Q2_step=0.4;

      vert_min["3He"]=-3.29;
      vert_min["4He"]=-2.53;
      vert_min["C12"]=4.8;
      vert_min["56Fe"]=4.6;
      vert_max["3He"]=-0.23;
      vert_max["4He"]=1.73;
      vert_max["C12"]=5.5;
      vert_max["56Fe"]=5.3;

      vertdiff_min["3He"]=-1.;
      vertdiff_min["4He"]=-1.;
      vertdiff_min["C12"]=-1.;
      vertdiff_min["56Fe"]=-1.;

      vertdiff_max["3He"]=1.;
      vertdiff_max["4He"]=1.;
      vertdiff_max["C12"]=1.;
      vertdiff_max["56Fe"]=1.;


//original cuts
 EC_photon_beta["3He"]=0.93;
 EC_photon_beta["4He"]=0.92;
 EC_photon_beta["C12"]=0.92;
 EC_photon_beta["56Fe"]=0.90;

      /*
//2.25sigma
cout<<"!!!!!!!! Modified photon PID!!!!!!!!!!!!!"<<endl;
 EC_photon_beta["3He"]=0.92;
 EC_photon_beta["4He"]=0.91;
 EC_photon_beta["C12"]=0.91;
 EC_photon_beta["56Fe"]=0.89;
      */
      /*
//1.75sigma
cout<<"!!!!!!!! Modified photon PID!!!!!!!!!!!!!"<<endl;
 EC_photon_beta["3He"]=0.94;
 EC_photon_beta["4He"]=0.93;
 EC_photon_beta["C12"]=0.93;
 EC_photon_beta["56Fe"]=0.91;
      */

 LEC_photon_beta["3He"]=0.96;
 LEC_photon_beta["4He"]=0.94;
 LEC_photon_beta["C12"]=0.94;
 LEC_photon_beta["56Fe"]=0.95;


    //map<std::pair<std::string, int>, double> EC_time_offset;
 EC_time_offset[std::make_pair("3He",1)]=-1.37;  EC_time_offset[std::make_pair("3He",2)]=-1.42; EC_time_offset[std::make_pair("3He",3)]=-1.55;
 EC_time_offset[std::make_pair("3He",4)]=-1.53;  EC_time_offset[std::make_pair("3He",5)]=-1.49; EC_time_offset[std::make_pair("3He",6)]=-1.44;

 EC_time_offset[std::make_pair("4He",1)]=0.72;  EC_time_offset[std::make_pair("4He",2)]=0.27; EC_time_offset[std::make_pair("4He",3)]=0.16;
 EC_time_offset[std::make_pair("4He",4)]=0.21;  EC_time_offset[std::make_pair("4He",5)]=0.22; EC_time_offset[std::make_pair("4He",6)]=0.21;

 EC_time_offset[std::make_pair("C12",1)]=0.50;  EC_time_offset[std::make_pair("C12",2)]=0.39; EC_time_offset[std::make_pair("C12",3)]=0.29;
 EC_time_offset[std::make_pair("C12",4)]=0.29;  EC_time_offset[std::make_pair("C12",5)]=0.32; EC_time_offset[std::make_pair("C12",6)]=0.33;

 EC_time_offset[std::make_pair("56Fe",1)]=0.75;  EC_time_offset[std::make_pair("56Fe",2)]=0.49; EC_time_offset[std::make_pair("56Fe",3)]=0.37;
 EC_time_offset[std::make_pair("56Fe",4)]=0.39;  EC_time_offset[std::make_pair("56Fe",5)]=0.43; EC_time_offset[std::make_pair("56Fe",6)]=0.44;

 elmom_corr_fact[0]=1.001;elmom_corr_fact[1]=0.991;elmom_corr_fact[2]=1.005;elmom_corr_fact[3]=1.004;elmom_corr_fact[4]=1.006;elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
    }


  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5)
    {
  hist_max = 6;
  prot_accept_mom_lim=0.3;
  chargedpi_accept_mom_lim=0.15;
  phot_accept_en_lim=0.3;
  prot_mom_lim=2.7;
  min_good_mom=1.1;
  if(ftarget=="3He") min_good_mom=1.3;//the EC threshold was different for He3
  max_mom=3.7;
  pipl_maxmom=1.9;
  pipl_delt_cutrange=3.;
  pimi_maxmom=1.6;
  pimi_delt_cutrange=3.;
  N_pperp=2,N_Ecal=6;
  pperp_cut=new double[N_pperp];
      Ecal_lowlim=new double[N_Ecal];
      Ecal_uplim=new double[N_Ecal];
      pperp_cut[0]=0.;  pperp_cut[1]=0.2;
      for (int i=0;i<N_Ecal;i++){
	Ecal_lowlim[i]=1.5+i*0.5;
	Ecal_uplim[i]=2.+i*0.5;
      }
      Ecal_lowlim[5]=0.;
      Ecal_uplim[5]=4.;
      for (int i=0;i<N_Ecal;i++)	cout<<Ecal_lowlim[i]<<"  to  "<<Ecal_uplim[i]<<endl;
      Q2_mean=1.5;
      Q2_step=1;

vert_min["3He"]=-3.27;
 vert_min["4He"]=-2.51;
 vert_min["C12"]=4.7;
 vert_min["56Fe"]=4.6;
 vert_max["3He"]=0.07;
 vert_max["4He"]=1.71;
 vert_max["C12"]=5.3;
 vert_max["56Fe"]=5.4;

 vertdiff_min["3He"]=-1.;
 vertdiff_min["4He"]=-1;
 vertdiff_min["C12"]=-1;
 vertdiff_min["56Fe"]=-1;

 vertdiff_max["3He"]=1.;
 vertdiff_max["4He"]=1.;
 vertdiff_max["C12"]=1;
 vertdiff_max["56Fe"]=1;


//original cut
EC_photon_beta["3He"]=0.92;
 EC_photon_beta["4He"]=0.91;
 EC_photon_beta["C12"]=0.92;
 EC_photon_beta["56Fe"]=0.91;

 /*
 //2.25sigma
cout<<"!!!!!!!! Modified photon PID!!!!!!!!!!!!!"<<endl;
EC_photon_beta["3He"]=0.91;
 EC_photon_beta["4He"]=0.9;
 EC_photon_beta["C12"]=0.91;
  EC_photon_beta["56Fe"]=0.9;
 */
 /*
//1.75sigma
cout<<"!!!!!!!! Modified photon PID!!!!!!!!!!!!!"<<endl;
EC_photon_beta["3He"]=0.93;
 EC_photon_beta["4He"]=0.92;
 EC_photon_beta["C12"]=0.93;
 EC_photon_beta["56Fe"]=0.92;
 */

 LEC_photon_beta["3He"]=0.97;
 LEC_photon_beta["4He"]=0.97;
 LEC_photon_beta["C12"]=0.95;
 LEC_photon_beta["56Fe"]=0.96;


 EC_time_offset[std::make_pair("3He",1)]=-0.15;  EC_time_offset[std::make_pair("3He",2)]=-0.26; EC_time_offset[std::make_pair("3He",3)]=-0.41;
 EC_time_offset[std::make_pair("3He",4)]=-0.29;  EC_time_offset[std::make_pair("3He",5)]=-0.25; EC_time_offset[std::make_pair("3He",6)]=-0.23;

 EC_time_offset[std::make_pair("4He",1)]=-0.01;  EC_time_offset[std::make_pair("4He",2)]=-0.11; EC_time_offset[std::make_pair("4He",3)]=-0.23;
 EC_time_offset[std::make_pair("4He",4)]=-0.26;  EC_time_offset[std::make_pair("4He",5)]=-0.21; EC_time_offset[std::make_pair("4He",6)]=-0.09;

 EC_time_offset[std::make_pair("C12",1)]=-0.01;  EC_time_offset[std::make_pair("C12",2)]=-0.11; EC_time_offset[std::make_pair("C12",3)]=-0.23;
 EC_time_offset[std::make_pair("C12",4)]=-0.27;  EC_time_offset[std::make_pair("C12",5)]=-0.21; EC_time_offset[std::make_pair("C12",6)]=-0.08;

 EC_time_offset[std::make_pair("56Fe",1)]=-0.49;  EC_time_offset[std::make_pair("56Fe",2)]=-0.14; EC_time_offset[std::make_pair("56Fe",3)]=-0.32;
 EC_time_offset[std::make_pair("56Fe",4)]=-0.25;  EC_time_offset[std::make_pair("56Fe",5)]=-0.17; EC_time_offset[std::make_pair("56Fe",6)]=-0.35;

 elmom_corr_fact[0]=1.001;elmom_corr_fact[1]=0.991;elmom_corr_fact[2]=1.005;elmom_corr_fact[3]=1.004;elmom_corr_fact[4]=1.006;elmom_corr_fact[5]=1.005;//a constant to multiply e- momentum with to correct the location of n peak in MM(3He(e,e'pp)n)
    }
  Ecal_offset["3He"]=0.004;
  Ecal_offset["4He"]=0.005;
  Ecal_offset["C12"]=0.005;
  Ecal_offset["56Fe"]=0.011;

  bind_en["3He"]=He3_bind_en-D2_bind_en+Ecal_offset["3He"]; //the offset is used to shift the peak to be at 0
  bind_en["C12"]=C12_bind_en-B_bind_en+Ecal_offset["C12"];
  bind_en["56Fe"]=Fe_bind_en-Mn_bind_en+Ecal_offset["56Fe"];

  target_mass["3He"]=2*m_prot+m_neut-He3_bind_en;
  target_mass["C12"]=6*m_prot+6*m_neut-C12_bind_en;
  target_mass["56Fe"]=26*m_prot+30*m_neut-Fe_bind_en;

  residual_target_mass["3He"]=m_prot+m_neut-D2_bind_en;
  residual_target_mass["C12"]=5*m_prot+6*m_neut-B_bind_en;
  residual_target_mass["56Fe"]=25*m_prot+30*m_neut-Mn_bind_en;

  up_lim1_ec=new TF1("up_lim1_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim2_ec=new TF1("up_lim2_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim3_ec=new TF1("up_lim3_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim4_ec=new TF1("up_lim4_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim5_ec=new TF1("up_lim5_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  up_lim6_ec=new TF1("up_lim6_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim1_ec=new TF1("low_lim1_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim2_ec=new TF1("low_lim2_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim3_ec=new TF1("low_lim3_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim4_ec=new TF1("low_lim4_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim5_ec=new TF1("low_lim5_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  low_lim6_ec=new TF1("low_lim6_ec","[0]+(x-[1])*(x-[1])*[2]",0,360);
  rightside_lim1_ec=new TF1("rightside_lim1_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim1_ec=new TF1("leftside_lim1_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim2_ec=new TF1("rightside_lim2_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim2_ec=new TF1("leftside_lim2_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim3_ec=new TF1("rightside_lim3_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim3_ec=new TF1("leftside_lim3_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim4_ec=new TF1("rightside_lim4_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim4_ec=new TF1("leftside_lim4_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim5_ec=new TF1("rightside_lim5_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim5_ec=new TF1("leftside_lim5_ec","[0]*(x+[1])+[2]",0,360);
  rightside_lim6_ec=new TF1("rightside_lim6_ec","[0]*(x+[1])+[2]",0,360);
  leftside_lim6_ec=new TF1("leftside_lim6_ec","[0]*(x+[1])+[2]",0,360);


  up_lim1_ec->SetParameters(0.995,30,-0.0001);
  up_lim2_ec->SetParameters(0.995,90,-0.0001);
  up_lim3_ec->SetParameters(0.995,150,-0.0001);
  up_lim4_ec->SetParameters(0.995,210,-0.0001);
  up_lim5_ec->SetParameters(0.995,270,-0.0001);
  up_lim6_ec->SetParameters(0.995,330,-0.0001);
  low_lim1_ec->SetParameters(0.7,30,-0.00005);
  low_lim2_ec->SetParameters(0.7,90,-0.00005);
  low_lim3_ec->SetParameters(0.7,150,-0.00005);
  low_lim4_ec->SetParameters(0.7,210,-0.00005);
  low_lim5_ec->SetParameters(0.7,270,-0.00005);
  low_lim6_ec->SetParameters(0.7,330,-0.00005);
  leftside_lim1_ec->SetParameters(0.11,0,0.03);
  rightside_lim1_ec->SetParameters(-0.11,-60,0.03);
  leftside_lim2_ec->SetParameters(0.11,-60,0.03);
  rightside_lim2_ec->SetParameters(-0.11,-120,0.03);
  leftside_lim3_ec->SetParameters(0.11,-120,0.03);
  rightside_lim3_ec->SetParameters(-0.11,-180,0.03);
  leftside_lim4_ec->SetParameters(0.11,-180,0.03);
  rightside_lim4_ec->SetParameters(-0.11,-240,0.03);
  leftside_lim5_ec->SetParameters(0.11,-240,0.03);
  rightside_lim5_ec->SetParameters(-0.11,-300,0.03);
  leftside_lim6_ec->SetParameters(0.11,-300,0.03);
  rightside_lim6_ec->SetParameters(-0.11,-360,0.03);

  Long64_t nentries = fChain->GetEntriesFast();
  nentries = 10000000;

  TFile *file_in=new TFile(Form("el_Epratio_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in1=new TFile(Form("protdeltat_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in2=new TFile(Form("pimideltat_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in4=new TFile(Form("pipldeltat_mom_%s.root",fbeam_en.c_str()));
  TFile *file_in3 = new TFile(Form("vz_%s_%s.root",ftarget.c_str(),fbeam_en.c_str()));
  TFile *file_in6 = new TFile("vz_3He_2261_2ndrungroup.root");
  TFile *file_in7 = new TFile("vz_56Fe_2261_badruns.root");
  TFile *file_out = new TFile(Form("%s_%s.root", ftarget.c_str(),fbeam_en.c_str()), "Recreate");
  TF1* f_gaus=new TF1("f_gaus","TMath::Landau(x,[0],[1])", 0,20);//smearing with Landau
  vz_corr_func = (TF1 *)file_in3->Get("f_vz");
  el_Epratio_mean=(TF1*)file_in->Get("f_mean");
  el_Epratio_sig=(TF1*)file_in->Get("f_sig");
  prot_deltat_sig=(TF1*)file_in1->Get("sig_pol9");
  prot_deltat_mean=(TF1*)file_in1->Get("mean_pol9");
  pipl_deltat_sig=(TF1*)file_in4->Get("sig_pol9");
  pipl_deltat_mean=(TF1*)file_in4->Get("mean_pol9");
  pimi_deltat_sig=(TF1*)file_in2->Get("sig_pol9");
  pimi_deltat_mean=(TF1*)file_in2->Get("mean_pol9");

  fsum_pimi=new TF1("fsum_pimi",FSum_pimi,0.,5.,2);
  fsub_pimi=new TF1("fsub_pimi",FSub_pimi,0.,5.,2);
  fsum_pipl=new TF1("fsum_pipl",FSum_pipl,0.,5.,2);
  fsub_pipl=new TF1("fsub_pipl",FSub_pipl,0.,5.,2);
  fsum_prot=new TF1("fsum_prot",FSum_prot,0.,5.,2);
  fsub_prot=new TF1("fsub_pprot",FSub_prot,0.,5.,2);
  fsum_e=new TF1("fsum_e",FSum_e,0.,5.,2);
  fsub_e=new TF1("fsub_e",FSub_e,0.,5.,2);

  int N_qe;
    double *x_qe;
  // Setting parameters for each beam energy and target
  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
   N_qe=109;
   x_qe=new double[N_qe+1];
   for (int i=0;i<=64;i++)x_qe[i]=-1+i*0.015;
   for (int i=0;i<=44;i++)x_qe[i+65]=-0.04+(i+1)*0.01;
 }


if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){

  if(ftarget=="56Fe"){

    N_qe=111;
    x_qe=new double[N_qe+1];
    for (int i=0;i<=58;i++)x_qe[i]=-1+i*0.015;
    for (int i=0;i<=52;i++)x_qe[i+59]=-0.13+(i+1)*0.01;
  }
  else{
   N_qe=82;
   x_qe=new double[N_qe+1];
   for (int i=0;i<=29;i++)x_qe[i]=-1+i*0.03;
   for (int i=0;i<=52;i++)x_qe[i+30]=-0.13+(i+1)*0.01;
  }
 }
 //Variable binning
 int n_bins;
 double *x_values;

 if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
  n_bins=38;
  x_values=new double[n_bins+1];

    for (int i=0;i<=17;i++)x_values[i]=0.4+i*0.04;
    for (int i=0;i<=20;i++)x_values[i+18]=1.08+(i+1)*0.02;
  }


  if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
  n_bins=54;
  x_values=new double[n_bins+1];

    for (int i=0;i<=23;i++)x_values[i]=i*0.09;
    for (int i=0;i<=30;i++)x_values[i+24]=2.07+(i+1)*0.03;
  }


  if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
  n_bins=38;
  x_values=new double[n_bins+1];

    for (int i=0;i<=21;i++) x_values[i]=i*0.2;
    for (int i=0;i<=16;i++)  x_values[i+22]=4.2+(i+1)*0.05;
    }
  // declaring histograms
  TH1F *h1_en_recon1=new TH1F("calorimetric","", n_bins, x_values);
  h1_en_recon1->Sumw2();
  TH1F *h1_en_recon1_pimi=new TH1F("cal_pimi", "", n_bins, x_values);
  h1_en_recon1_pimi->Sumw2();
  TH1F *h1_en_recon1_pipl=new TH1F("cal_pipl", "",n_bins, x_values);
  h1_en_recon1_pipl->Sumw2();
  TH1F *h1_en_recon2=new TH1F("kin_e","",n_bins, x_values);
  h1_en_recon2->Sumw2();
  TH1F *h1_en_recon2_pimi=new TH1F("kin_e_pimi","",n_bins, x_values);
  h1_en_recon2_pimi->Sumw2();
  TH1F *h1_en_recon2_pipl=new TH1F("kin_e_pipl","",n_bins, x_values);
  h1_en_recon2_pipl->Sumw2();
  TH1F *h1_en_recon3=new TH1F("kin_e_pi","",n_bins, x_values);
  h1_en_recon3->Sumw2();
  TH1F *h1_en_recon3_pimi=new TH1F("kin_e_pi_pimi","",n_bins, x_values);
  h1_en_recon3_pimi->Sumw2();
  TH1F *h1_en_recon3_pipl=new TH1F("kin_e_pi_pipl","",n_bins, x_values);
  h1_en_recon3_pipl->Sumw2();
  TH1F *h1_rot1_2pi_1p=new TH1F("2pi_1p_rot_cal","",n_bins, x_values);
  h1_rot1_2pi_1p->Sumw2();
  TH1F *h1_rot1_2pi_1p_pimi=new TH1F("2pi_1p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_2pi_1p_pimi->Sumw2();
  TH1F *h1_rot1_2pi_1p_pipl=new TH1F("2pi_1p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_2pi_1p_pipl->Sumw2();
  TH1F *h1_rot2_2pi_1p=new TH1F("2pi_1p_rot_kin_e","",n_bins, x_values);
  h1_rot2_2pi_1p->Sumw2();
  TH1F *h1_rot2_2pi_1p_pimi=new TH1F("2pi_1p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_2pi_1p_pimi->Sumw2();
  TH1F *h1_rot2_2pi_1p_pipl=new TH1F("2pi_1p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_2pi_1p_pipl->Sumw2();
  TH1F *h1_rot3_2pi_1p=new TH1F("2pi_1p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_2pi_1p->Sumw2();
  TH1F *h1_rot3_2pi_1p_pimi=new TH1F("2pi_1p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_2pi_1p_pimi->Sumw2();
  TH1F *h1_rot3_2pi_1p_pipl=new TH1F("2pi_1p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_2pi_1p_pipl->Sumw2();
  TH1F *h1_rot1_1pi_2p=new TH1F("1pi_2p_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_2p->Sumw2();
  TH1F *h1_rot1_1pi_2p_pimi=new TH1F("1pi_2p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_2p_pimi->Sumw2();
  TH1F *h1_rot1_1pi_2p_pipl=new TH1F("1pi_2p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_2p_pipl->Sumw2();
  TH1F *h1_rot2_1pi_2p=new TH1F("1pi_2p_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_2p->Sumw2();
  TH1F *h1_rot2_1pi_2p_pimi=new TH1F("1pi_2p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_2p_pimi->Sumw2();
  TH1F *h1_rot2_1pi_2p_pipl=new TH1F("1pi_2p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_2p_pipl->Sumw2();
  TH1F *h1_rot3_1pi_2p=new TH1F("1pi_2p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_2p->Sumw2();
  TH1F *h1_rot3_1pi_2p_pimi=new TH1F("1pi_2p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_2p_pimi->Sumw2();
  TH1F *h1_rot3_1pi_2p_pipl=new TH1F("1pi_2p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_2p_pipl->Sumw2();
  TH1F *h1_rot1_2pi_2p=new TH1F("2pi_2p_rot_cal","",n_bins, x_values);
  h1_rot1_2pi_2p->Sumw2();
  TH1F *h1_rot1_2pi_2p_pimi=new TH1F("2pi_2p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_2pi_2p_pimi->Sumw2();
  TH1F *h1_rot1_2pi_2p_pipl=new TH1F("2pi_2p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_2pi_2p_pipl->Sumw2();
  TH1F *h1_rot2_2pi_2p=new TH1F("2pi_2p_rot_kin_e","",n_bins, x_values);
  h1_rot2_2pi_2p->Sumw2();
  TH1F *h1_rot2_2pi_2p_pimi=new TH1F("2pi_2p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_2pi_2p_pimi->Sumw2();
  TH1F *h1_rot2_2pi_2p_pipl=new TH1F("2pi_2p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_2pi_2p_pipl->Sumw2();
  TH1F *h1_rot3_2pi_2p=new TH1F("2pi_2p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_2pi_2p->Sumw2();
  TH1F *h1_rot3_2pi_2p_pimi=new TH1F("2pi_2p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_2pi_2p_pimi->Sumw2();
  TH1F *h1_rot3_2pi_2p_pipl=new TH1F("2pi_2p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_2pi_2p_pipl->Sumw2();
  TH1F *h1_rot1_1pi_3p=new TH1F("1pi_3p_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_3p->Sumw2();
  TH1F *h1_rot1_1pi_3p_pimi=new TH1F("1pi_3p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_3p_pimi->Sumw2();
  TH1F *h1_rot1_1pi_3p_pipl=new TH1F("1pi_3p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_3p_pipl->Sumw2();
  TH1F *h1_rot2_1pi_3p=new TH1F("1pi_3p_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_3p->Sumw2();
  TH1F *h1_rot2_1pi_3p_pimi=new TH1F("1pi_3p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_3p_pimi->Sumw2();
  TH1F *h1_rot2_1pi_3p_pipl=new TH1F("1pi_3p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_3p_pipl->Sumw2();
  TH1F *h1_rot3_1pi_3p=new TH1F("1pi_3p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_3p->Sumw2();
  TH1F *h1_rot3_1pi_3p_pimi=new TH1F("1pi_3p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_3p_pimi->Sumw2();
  TH1F *h1_rot3_1pi_3p_pipl=new TH1F("1pi_3p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_3p_pipl->Sumw2();
  TH1F *h1_rot1_3pi_1p=new TH1F("3pi_1p_rot_cal","",n_bins, x_values);
  h1_rot1_3pi_1p->Sumw2();
  TH1F *h1_rot1_3pi_1p_pimi=new TH1F("3pi_1p_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_3pi_1p_pimi->Sumw2();
  TH1F *h1_rot1_3pi_1p_pipl=new TH1F("3pi_1p_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_3pi_1p_pipl->Sumw2();
  TH1F *h1_rot2_3pi_1p=new TH1F("3pi_1p_rot_kin_e","",n_bins, x_values);
  h1_rot2_3pi_1p->Sumw2();
  TH1F *h1_rot2_3pi_1p_pimi=new TH1F("3pi_1p_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_3pi_1p_pimi->Sumw2();
  TH1F *h1_rot2_3pi_1p_pipl=new TH1F("3pi_1p_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_3pi_1p_pipl->Sumw2();
  TH1F *h1_rot3_3pi_1p=new TH1F("3pi_1p_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_3pi_1p->Sumw2();
  TH1F *h1_rot3_3pi_1p_pimi=new TH1F("3pi_1p_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_3pi_1p_pimi->Sumw2();
  TH1F *h1_rot3_3pi_1p_pipl=new TH1F("3pi_1p_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_3pi_1p_pipl->Sumw2();
  TH1F *h1_rot1_1pi_1p_1phot=new TH1F("1pi_1p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_1p_1phot->Sumw2();
  TH1F *h1_rot1_1pi_1p_1phot_pimi=new TH1F("1pi_1p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_1p_1phot_pimi->Sumw2();
  TH1F *h1_rot1_1pi_1p_1phot_pipl=new TH1F("1pi_1p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_1p_1phot_pipl->Sumw2();
  TH1F *h1_rot2_1pi_1p_1phot=new TH1F("1pi_1p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_1p_1phot->Sumw2();
  TH1F *h1_rot2_1pi_1p_1phot_pimi=new TH1F("1pi_1p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_1p_1phot_pimi->Sumw2();
  TH1F *h1_rot2_1pi_1p_1phot_pipl=new TH1F("1pi_1p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_1p_1phot_pipl->Sumw2();
  TH1F *h1_rot3_1pi_1p_1phot=new TH1F("1pi_1p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_1p_1phot->Sumw2();
  TH1F *h1_rot3_1pi_1p_1phot_pimi=new TH1F("1pi_1p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_1p_1phot_pimi->Sumw2();
  TH1F *h1_rot3_1pi_1p_1phot_pipl=new TH1F("1pi_1p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_1p_1phot_pipl->Sumw2();
  TH1F *h1_rot1_1pi_2p_1phot=new TH1F("1pi_2p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_2p_1phot->Sumw2();
  TH1F *h1_rot1_1pi_2p_1phot_pimi=new TH1F("1pi_2p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_2p_1phot_pimi->Sumw2();
  TH1F *h1_rot1_1pi_2p_1phot_pipl=new TH1F("1pi_2p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_2p_1phot_pipl->Sumw2();
  TH1F *h1_rot2_1pi_2p_1phot=new TH1F("1pi_2p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_2p_1phot->Sumw2();
  TH1F *h1_rot2_1pi_2p_1phot_pimi=new TH1F("1pi_2p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_2p_1phot_pimi->Sumw2();
  TH1F *h1_rot2_1pi_2p_1phot_pipl=new TH1F("1pi_2p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_2p_1phot_pipl->Sumw2();
  TH1F *h1_rot3_1pi_2p_1phot=new TH1F("1pi_2p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_2p_1phot->Sumw2();
  TH1F *h1_rot3_1pi_2p_1phot_pimi=new TH1F("1pi_2p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_2p_1phot_pimi->Sumw2();
  TH1F *h1_rot3_1pi_2p_1phot_pipl=new TH1F("1pi_2p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_2p_1phot_pipl->Sumw2();
  TH1F *h1_rot1_1pi_1p_2phot=new TH1F("1pi_1p_2phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_1p_2phot->Sumw2();
  TH1F *h1_rot1_1pi_1p_2phot_pimi=new TH1F("1pi_1p_2phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_1p_2phot_pimi->Sumw2();
  TH1F *h1_rot1_1pi_1p_2phot_pipl=new TH1F("1pi_1p_2phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_1p_2phot_pipl->Sumw2();
  TH1F *h1_rot2_1pi_1p_2phot=new TH1F("1pi_1p_2phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_1p_2phot->Sumw2();
  TH1F *h1_rot2_1pi_1p_2phot_pimi=new TH1F("1pi_1p_2phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_1p_2phot_pimi->Sumw2();
  TH1F *h1_rot2_1pi_1p_2phot_pipl=new TH1F("1pi_1p_2phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_1p_2phot_pipl->Sumw2();
  TH1F *h1_rot3_1pi_1p_2phot=new TH1F("1pi_1p_2phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_1p_2phot->Sumw2();
  TH1F *h1_rot3_1pi_1p_2phot_pimi=new TH1F("1pi_1p_2phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_1p_2phot_pimi->Sumw2();
  TH1F *h1_rot3_1pi_1p_2phot_pipl=new TH1F("1pi_1p_2phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_1p_2phot_pipl->Sumw2();
  TH1F *h1_rot1_2pi_1p_1phot=new TH1F("2pi_1p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_2pi_1p_1phot->Sumw2();
  TH1F *h1_rot1_2pi_1p_1phot_pimi=new TH1F("2pi_1p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_2pi_1p_1phot_pimi->Sumw2();
  TH1F *h1_rot1_2pi_1p_1phot_pipl=new TH1F("2pi_1p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_2pi_1p_1phot_pipl->Sumw2();
  TH1F *h1_rot2_2pi_1p_1phot=new TH1F("2pi_1p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_2pi_1p_1phot->Sumw2();
  TH1F *h1_rot2_2pi_1p_1phot_pimi=new TH1F("2pi_1p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_2pi_1p_1phot_pimi->Sumw2();
  TH1F *h1_rot2_2pi_1p_1phot_pipl=new TH1F("2pi_1p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_2pi_1p_1phot_pipl->Sumw2();
  TH1F *h1_rot3_2pi_1p_1phot=new TH1F("2pi_1p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_2pi_1p_1phot->Sumw2();
  TH1F *h1_rot3_2pi_1p_1phot_pimi=new TH1F("2pi_1p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_2pi_1p_1phot_pimi->Sumw2();
  TH1F *h1_rot3_2pi_1p_1phot_pipl=new TH1F("2pi_1p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_2pi_1p_1phot_pipl->Sumw2();
  TH1F *h1_rot1_1pi_3p_1phot=new TH1F("1pi_3p_1phot_rot_cal","",n_bins, x_values);
  h1_rot1_1pi_3p_1phot->Sumw2();
  TH1F *h1_rot1_1pi_3p_1phot_pimi=new TH1F("1pi_3p_1phot_rot_cal_pimi","",n_bins, x_values);
  h1_rot1_1pi_3p_1phot_pimi->Sumw2();
  TH1F *h1_rot1_1pi_3p_1phot_pipl=new TH1F("1pi_3p_1phot_rot_cal_pipl","",n_bins, x_values);
  h1_rot1_1pi_3p_1phot_pipl->Sumw2();
  TH1F *h1_rot2_1pi_3p_1phot=new TH1F("1pi_3p_1phot_rot_kin_e","",n_bins, x_values);
  h1_rot2_1pi_3p_1phot->Sumw2();
  TH1F *h1_rot2_1pi_3p_1phot_pimi=new TH1F("1pi_3p_1phot_rot_kin_e_pimi","",n_bins, x_values);
  h1_rot2_1pi_3p_1phot_pimi->Sumw2();
  TH1F *h1_rot2_1pi_3p_1phot_pipl=new TH1F("1pi_3p_1phot_rot_kin_e_pipl","",n_bins, x_values);
  h1_rot2_1pi_3p_1phot_pipl->Sumw2();
  TH1F *h1_rot3_1pi_3p_1phot=new TH1F("1pi_3p_1phot_rot_kin_e_pi","",n_bins, x_values);
  h1_rot3_1pi_3p_1phot->Sumw2();
  TH1F *h1_rot3_1pi_3p_1phot_pimi=new TH1F("1pi_3p_1phot_rot_kin_e_pi_pimi","",n_bins, x_values);
  h1_rot3_1pi_3p_1phot_pimi->Sumw2();
  TH1F *h1_rot3_1pi_3p_1phot_pipl=new TH1F("1pi_3p_1phot_rot_kin_e_pi_pipl","",n_bins, x_values);
  h1_rot3_1pi_3p_1phot_pipl->Sumw2();
  TH1F *h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
  TH1F *h1_p_vert_corr=new TH1F("h1_p_vert_corr", "", 300, -10, 10);
  TH1F *h1_pimi_vert_corr=new TH1F("h1_pimi_vert_corr", "", 300, -10, 10);
  TH1F *h1_pipl_vert_corr=new TH1F("h1_pipl_vert_corr", "", 300, -10, 10);
  TH1F *h1_p_vert_corr_cut=new TH1F("h1_p_vert_corr_cut", "", 300, -10, 10);
  TH1F *h1_pimi_vert_corr_cut=new TH1F("h1_pimi_vert_corr_cut", "", 300, -10, 10);
  TH1F *h1_pipl_vert_corr_cut=new TH1F("h1_pipl_vert_corr_cut", "", 300, -10, 10);
  TH1F *h1_el_vertuncorr=new TH1F("h1_el_vertuncorr","",200,-10,10);
  TH1F *h1_el_vertcorr=new TH1F("h1_el_vertcorr","",200,-10,10);
  TH1F *h1_el_vertcorr_cut=new TH1F("h1_el_vertcorr_cut","",200,-10,10);
  TH1F *h1_ec_beta_corr=new TH1F("h1_ec_beta_corr","",300,0,2);
  TH1F *h1_ec_beta_corr_cut=new TH1F("h1_ec_beta_corr_cut","",300,0,2);
  TH1F *h1_Q2 = new TH1F("h1_Q2","",200,0,6);
  h1_Q2->Sumw2();
  TH1F *h1_omega = new TH1F("h1_omega","",200,0,5);
  h1_omega->Sumw2();
  TH1F *h1_Wvar = new TH1F("h1_Wvar","",200,0,3);
  h1_Wvar->Sumw2();
  TH1F *h1_Q2_sub = new TH1F("h1_Q2_sub","",200,0,6);
  h1_Q2_sub->Sumw2();
  TH1F *h1_omega_sub = new TH1F("h1_omega_sub","",200,0,5);
  h1_omega_sub->Sumw2();
  TH1F *h1_Wvar_sub = new TH1F("h1_Wvar_sub","",200,0,3);
  h1_Wvar_sub->Sumw2();
  TH1F *h1_Q2_cut = new TH1F("h1_Q2_cut","",200,0,6);
  h1_Q2_cut->Sumw2();
  TH1F *h1_omega_cut = new TH1F("h1_omega_cut","",200,0,5);
  h1_omega_cut->Sumw2();
  TH1F *h1_Wvar_cut = new TH1F("h1_Wvar_cut","",200,0,3);
  h1_Wvar_cut->Sumw2();
  TH1F *h1_feed_down_cal = new TH1F("h1_feed_down_cal","",200,0,1);
  h1_feed_down_cal->Sumw2();
  TH1F *h1_feed_down_kin = new TH1F("h1_feed_down_kin","",200,0,1);
  h1_feed_down_kin->Sumw2();
  TH1F *h1_p_perp = new TH1F("h1_p_perp","",400,0,1);
  h1_p_perp->Sumw2();
  TH1F *h1_p_perp_sub = new TH1F("h1_p_perp_sub","",400,0,1);
  h1_p_perp_sub->Sumw2();
  TH1F *h1_p_perp_cut = new TH1F("h1_p_perp_cut","",400,0,1);
  h1_p_perp_cut->Sumw2();
  TH1F *h1_cal_p_slice1 = new TH1F("h1_cal_p_slice1","",n_bins, x_values);
  h1_cal_p_slice1->Sumw2();
  TH1F *h1_cal_p_slice2 = new TH1F("h1_cal_p_slice2","",n_bins, x_values);
  h1_cal_p_slice2->Sumw2();
  TH1F *h1_cal_p_slice3 = new TH1F("h1_cal_p_slice3","",n_bins, x_values);
  h1_cal_p_slice3->Sumw2();
  TH1F *h1_cal_p_slice1_sub = new TH1F("h1_cal_p_slice1_sub","",n_bins, x_values);
  h1_cal_p_slice1_sub->Sumw2();
  TH1F *h1_cal_p_slice2_sub = new TH1F("h1_cal_p_slice2_sub","",n_bins, x_values);
  h1_cal_p_slice2_sub->Sumw2();
  TH1F *h1_cal_p_slice3_sub = new TH1F("h1_cal_p_slice3_sub","",n_bins, x_values);
  h1_cal_p_slice3_sub->Sumw2();
  TH2F *h2_e_ec_xy = new TH2F("h2_e_ec_xy","",100,-600,600,100,-600,600);
  TH2F *h2_e_ec_xy_fidcut = new TH2F("h2_e_ec_xy_fidcut","",100,-600,600,100,-600,600);
  TH2F *h2_neut_costheta_phi=new TH2F("h2_neut_costheta_phi","",200,0,360,200,0,1.1);
  TH2F *h2_neut_costheta_phi_cut=new TH2F("h2_neut_costheta_phi_cut","",200,0,360,200,0,1.1);
  TH2F *h2_e_phi_theta=new TH2F("h2_e_phi_theta", "", 200,0,360,200,0,180);
  TH2F *h2_e_phi_theta_cut=new TH2F("h2_e_phi_theta_cut", "", 200,0,360,200,0,180);
  TH2F *h2_e_phi_costheta=new TH2F("h2_e_phi_costheta", "", 200,0,360,200,-1,1);
  TH2F *h2_e_phi_costheta_cut=new TH2F("h2_e_phi_costheta_cut", "", 200,0,360,200,-1,1);
  TH2F *h2_p_phi_theta=new TH2F("h2_p_phi_theta", "", 200,0,360,200,0,180);
  TH2F *h2_p_phi_theta_cut=new TH2F("h2_p_phi_theta_cut", "", 200,0,360,200,0,180);
  TH2F *h2_pimi_phi_theta=new TH2F("h2_pimi_phi_theta", "", 200,0,360,200,0,180);
  TH2F *h2_pimi_phi_theta_cut=new TH2F("h2_pimi_phi_theta_cut", "", 200,0,360,200,0,180);
  TH2F *h2_pipl_phi_theta=new TH2F("h2_pipl_phi_theta", "", 200,0,360,200,0,180);
  TH2F *h2_pipl_phi_theta_cut=new TH2F("h2_pipl_phi_theta_cut", "", 200,0,360,200,0,180);
  TH2F *h2_Np_Npi=new TH2F("h2_Np_Npi","",11,-0.5,4.5,11,-0.5, 4.5);
  TH2F *h2_phot_pi_1p=new TH2F("h2_phot_pi_1p", "",11,-0.5,4.5,11,-0.5,4.5);
  TH2F *h2_phot_pi_2p=new TH2F("h2_phot_pi_2p", "",11,-0.5,4.5,11,-0.5,4.5);
  TH2F *h2_phot_pi_3p=new TH2F("h2_phot_pi_3p", "",11,-0.5,4.5,11,-0.5,4.5);
  TH2F *prot_Deltat_p= new TH2F("prot_Deltat_p","",300,0.,2.7,300,-15,15);
  TH2F *prot_Deltat_p_cut= new TH2F("prot_Deltat_p_cut","",300,0.,2.7,300,-15,15);
  TH2F *pimi_delt_p= new TH2F("pimi_delt_p","",300,0.,2.5,300,-15,15);
  TH2F *pimi_delt_p_cut= new TH2F("pimi_delt_p_cut","",300,0.,2.5,300,-15,15);
  TH2F *pipl_delt_p= new TH2F("pipl_delt_p","",300,0.,2.5,300,-15,15);
  TH2F *pipl_delt_p_cut= new TH2F("pipl_delt_p_cut","",300,0.,2.5,300,-15,15);
  TH2F *h2_el_E_p_ratio = new TH2F("h2_el_E_p_ratio","",200,0,4.5,200,0,0.5);
  TH2F *h2_el_E_p_ratio_cut = new TH2F("h2_el_E_p_ratio_cut","",200,0,4.5,200,0,0.5);
  TH2F *h2_Wvar_Q2 = new TH2F("h2_Wvar_Q2","", 200,0,3,200,0,5);
  h2_Wvar_Q2->Sumw2();
  TH2F *h2_Wvar_Q2_sect1 = new TH2F("h2_Wvar_Q2_sect1","", 200,0,3,200,0,5);
  h2_Wvar_Q2_sect1->Sumw2();
  TH2F *h2_omega_Q2 = new TH2F("h2_omega_Q2","",200,0,3.5,200,0,5);
  h2_omega_Q2->Sumw2();
  TH2F *h2_omega_Q2_sect1 = new TH2F("h2_omega_Q2_sect1","",200,0,3.5,200,0,5);
  h2_omega_Q2_sect1->Sumw2();
  TH2F *h2_Wvar_Q2_sub = new TH2F("h2_Wvar_Q2_sub","", 200,0,3,200,0,5);
  h2_Wvar_Q2_sub->Sumw2();
  TH2F *h2_omega_Q2_sub = new TH2F("h2_omega_Q2_sub","",200,0,3.5,200,0,5);
  h2_omega_Q2_sub->Sumw2();
  TH2F *h2_kin_e_Wvar = new TH2F("h2_kin_e_Wvar","",200,0,5,200,0,5);
  h2_kin_e_Wvar->Sumw2();
  TH2F *h2_kin_e_pi_Wvar = new TH2F("h2_kin_e_pi_Wvar","",200,0,5,200,0,5);
  h2_kin_e_pi_Wvar->Sumw2();
  TH2F *h2_cal_Wvar = new TH2F("h2_cal_Wvar","",200,0,5,200,0,5);
  h2_cal_Wvar->Sumw2();
  TH2F *h2_pperp_cal = new TH2F("h2_pperp_cal","",200,0,1,200,0,6);
  h2_pperp_cal->Sumw2();
  TH2F *h2_pperp_kin = new TH2F("h2_pperp_kin","",200,0,1,200,0,6);
  h2_pperp_kin->Sumw2();
  TH2F *h2_p_perp_Wvar = new TH2F("h2_p_perp_Wvar","",200,0,1,200,0,5);
  h2_p_perp_Wvar->Sumw2();
  TH3F *h3_Npi_Np_Nphoton=new TH3F("h3_Npi_Np_Nphoton", "", 11, -0.5, 4.5, 11, -0.5, 4.5, 11, -0.5, 4.5);

  // Pulling parameters from f_vz file
  double pars[3];
  if(ftarget=="56Fe"){
    pars[0]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(0);
    pars[1]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(1);
    pars[2]=((TF1 *)file_in7->Get("f_vz"))->GetParameter(2);
  }
  if(ftarget=="3He")
  {
    pars[0]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(0);
    pars[1]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(1);
    pars[2]=((TF1 *)file_in6->Get("f_vz"))->GetParameter(2);
  }
  int n_evnt=1;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if( jentry%200000 == 0 )
    {
	     gDirectory->Write("hist_Files", TObject::kOverwrite);
	     cout<<"Processing event "<<jentry<<endl;
    }
    if (runnb==18258 || runnb==18259 || (runnb>18382 && runnb<18438) || (runnb>18220 && runnb<18253)) vz_corr_func->SetParameters(pars);

    if(runnb==18258 || runnb==18259)   {   //setting appropriate e- vertex cut range for the runs with the same target and bema energy, but different vertex correction
      vert_max["56Fe"]=6.;
      vert_min["56Fe"]=5.2;//runs with exploaded cell
    }
    if(runnb>18382 && runnb<18438){
      vert_max["3He"]=0.01;
      vert_min["3He"]=-3.31; //runs with thin exit window
    }

    if((runnb>18283 && runnb<18289) || (runnb>18300 && runnb<18304) || (runnb>18317 && runnb<18329))fTorusCurrent=750;    //setting appropriate torrus magnet current
    else if ((runnb>18293 && runnb<18301) || (runnb>18305 && runnb<18317) || (runnb>18328 && runnb<18336))fTorusCurrent=1500;
    else fTorusCurrent=2250;
    fTorusCurrent=2250;



    if(n_evnt==1){SetFiducialCutParameters();}
    n_evnt++;


    int n_elec = 0;
    const int ind_em=0;

    double prot_delt_cutrange=3.;
    double min_good_mom=0.55;
    int nsect = 6;
    Double_t sc_cc_delt_cut_sect[6]={-2,-5,-8,-8,-2,2};
    TVector3 el_mom1(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em] ,p[ind_em]*cz[ind_em]);
    double sc_time = sc_t[sc[ind_em] - 1];    // time it took the electron to reach the time of fligh (TOF) detectors from the production vertex
    double sc_path = sc_r[sc[ind_em] - 1]; //path length to TOF detectors
    int sc_paddle = sc_pd[sc[ind_em] - 1];
    int sc_sector = sc_sect[sc[ind_em] - 1];  //The TOF sector that the electron has been detected in (1-6)
    float el_vert=vz[ind_em];
    double ec_x=ech_x[ind_em];
    double ec_y=ech_y[ind_em];
    double ec_z=ech_z[ind_em];
    TVector3 e_ec_xyz1(ech_x[ec[ind_em]-1],ech_y[ec[ind_em]-1],ech_z[ec[ind_em]-1]);
    double el_theta=TMath::ACos(cz[ind_em])*TMath::RadToDeg();
    double el_phi_mod=TMath::ATan2(cy[ind_em],cx[ind_em])*TMath::RadToDeg()+30;
    if(el_phi_mod<0)el_phi_mod=el_phi_mod+360;
    int el_ec_sector = ec_sect[ec[ind_em] - 1];
    double el_vert_corr=el_vert+vz_corr(el_phi_mod,el_theta);

    ece = TMath::Max(ec_ei[ec[ind_em] - 1] + ec_eo[ec[ind_em] - 1],
		      etot[ec[ind_em] - 1]);
    el_segment=int((cc_segm[cc[ind_em]-1]-int(cc_segm[cc[ind_em]-1]/1000)*1000)/10);
    el_cc_sector=cc_sect[cc[ind_em]-1];
    el_sccc_timediff=sc_t[cc[ind_em]-1]-cc_t[cc[ind_em]-1]-(sc_r[cc[ind_em]-1]-cc_r[cc[ind_em]-1])/(c*ns_to_s);


    TLorentzVector V4_el_uncorr(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em],p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass));
    TLorentzVector V4_el(p[ind_em]*cx[ind_em],p[ind_em]*cy[ind_em],p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass));
    if(el_ec_sector>=1 && el_ec_sector<=6) V4_el.SetPxPyPzE(elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cx[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cy[ind_em],elmom_corr_fact[el_ec_sector-1]*p[ind_em]*cz[ind_em],TMath::Sqrt(p[ind_em]*p[ind_em]*elmom_corr_fact[el_ec_sector-1]*elmom_corr_fact[el_ec_sector-1]+e_mass*e_mass));
    TVector3 V3_el = V4_el.Vect();
    h1_el_mom_ratio->Fill(V4_el.Rho()/V4_el_uncorr.Rho());
    double en_recon2 = (m_delta*m_delta-(m_prot-eps)*(m_prot-eps)+2*(m_prot-eps)*V4_el.E())/(2*(m_prot-eps-V4_el.E()+V4_el.Rho()*cz[ind_em]));
    TVector3 V3_pipl;
    TVector3 V3_pimi;
    int index_p[20],ind_p,index_pi[20],ind_pi, ind_pi_phot[20];
    int num_p = 0,num_pi=0,num_pimi=0,num_pipl=0, num_pi_phot=0;
    int index_n[20]={},ec_index_n[20]={},index_pipl[20]={},index_pimi[20]={};
    int num_n = 0,ec_num_n = 0;
    double pi_phi,pi_phi_mod, pi_theta,pimi_phi,pimi_phi_mod,pimi_theta,pipl_phi,pipl_phi_mod, pipl_theta;
    double fine_struc_const=0.007297;
    double Mott_cross_sec=(fine_struc_const*fine_struc_const*(cz[ind_em]+1))/(2*V4_el.E()*V4_el.E()*(1-cz[ind_em])*(1-cz[ind_em]));
    const double pimi_vertcut=2.5,pipl_vertcut=2.5;
    double pimi_vert, pipl_vert, prot_vert, pimi_vert_corr, pipl_vert_corr, prot_vert_corr, prot_theta, prot_phi, prot_phi_mod;
    double neut_zvert, neut_yvert, neut_xvert, neut_ecx, neut_ecy, neut_ecz;
    TVector3 phot_ec_xyz, V3_phot, V3_el_angles_SC;
    double ec_deltt, neut_ecpath_corr, neut_ectime_corr, neut_beta_corr, neutr_phi_mod, photon_ece;
    double phot_rad_cut = 40, phot_e_phidiffcut = 30;
    bool ec_radstat_n[20]={false};
    double epratio_sig_cutrange=3.;
    double max_mom=2.1;

    fsum_e->SetParameters(epratio_sig_cutrange, max_mom);
    fsub_e->SetParameters(epratio_sig_cutrange, max_mom);

    h2_el_E_p_ratio->Fill(p[ind_em], ece/p[ind_em]);

    if(ec[ind_em] < 0.5 || cc[ind_em] < 0.5 ||  sc[ind_em] < 0.5 || q[ind_em] >=0 || ec_ei[ec[ind_em] - 1] < 0.06 || ece/p[ind_em]<fsub_e->Eval(p[ind_em]) || ece/p[ind_em] >fsum_e->Eval(p[ind_em]) || p[ind_em]<min_good_mom || cc_c2[cc[ind_em]-1]>=0.1 || el_sccc_timediff<sc_cc_delt_cut_sect[el_cc_sector-1]  || TMath::Sqrt(p[ind_em]*p[ind_em]+e_mass*e_mass)>en_beam[fbeam_en]) //electron pid cuts
        {continue;}
    h2_el_E_p_ratio_cut->Fill(p[ind_em], ece/p[ind_em]);

    h2_e_phi_theta->Fill(el_phi_mod, el_theta);
    h2_e_phi_costheta->Fill(el_phi_mod, cz[ind_em]);
    h2_e_ec_xy->Fill(ec_x,ec_y);
    if(!EFiducialCut(el_mom1))continue;//theta, phi cuts
    if(!CutUVW(e_ec_xyz1))continue;
    h2_e_ec_xy_fidcut->Fill(ec_x,ec_y);
    h2_e_phi_theta_cut->Fill(el_phi_mod, el_theta);
    h2_e_phi_costheta_cut->Fill(el_phi_mod, cz[ind_em]);

  for( int i = 1; i < TMath::Min(gpart, 20); i++ )   //i is the index of the particle in the event
      {
	if( sc[i] > 0 && stat[i] > 0 &&  id[i] == 2212 )// looking for protons, id ==2212
	  {
      bett=p[i]/TMath::Sqrt(p[i]*p[i]+m_prot*m_prot);
      deltt=sc_t[sc[i]-1]-sc_r[sc[i]-1]/(bett*c*ns_to_s) - tr_time;
      fsum_prot->SetParameters(prot_delt_cutrange,prot_mom_lim);
	    fsub_prot->SetParameters(prot_delt_cutrange,prot_mom_lim);

	    prot_Deltat_p->Fill(p[i],deltt);

	    if(deltt<fsum_prot->Eval(p[i]) && deltt>fsub_prot->Eval(p[i]) && p[i]>=prot_accept_mom_lim){//Proton PID cut
      prot_Deltat_p_cut->Fill(p[i],deltt);
      TVector3 V3_p;
      V3_p.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
      prot_vert = vz[i];
      prot_phi=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
      prot_phi_mod=prot_phi+30;
      if (prot_phi_mod<0)prot_phi_mod=prot_phi_mod+360;
      prot_theta=TMath::ACos(cz[i])*TMath::RadToDeg();
      prot_vert_corr = prot_vert + vz_corr(prot_phi_mod,prot_theta);
      h2_p_phi_theta->Fill(prot_phi_mod, prot_theta);


      if(PFiducialCut(V3_p)){
        h2_p_phi_theta_cut->Fill(prot_phi_mod, prot_theta);
        h1_p_vert_corr->Fill(el_vert_corr-prot_vert_corr);
      if((el_vert_corr-prot_vert_corr)>vertdiff_min[ftarget] && (el_vert_corr-prot_vert_corr)<vertdiff_max[ftarget]){
        h1_p_vert_corr_cut->Fill(el_vert_corr-prot_vert_corr);
      num_p = num_p + 1;
	    index_p[num_p-1]=i;
      }
      }
	   }
    }
	if(q[i]!=0 &&  sc[i] > 0 && dc[i]>0 && stat[i] > 0 &&  id[i] == -211) //looking for negative pions, id=-211
	  {
      V3_pimi.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
      bett=p[i]/TMath::Sqrt(p[i]*p[i]+Mpi*Mpi);
	    deltt=sc_t[sc[i]-1]-sc_r[sc[i]-1]/(bett*c*ns_to_s) - tr_time;

      fsub_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);
	    fsum_pimi->SetParameters(pimi_delt_cutrange,pimi_maxmom);

	    pimi_delt_p->Fill(p[i],deltt);

	    if(deltt<fsum_pimi->Eval(p[i]) && deltt>fsub_pimi->Eval(p[i])){//Negative pion PID cut

      pimi_delt_p_cut->Fill(p[i],deltt);

      pimi_vert = vz[i];
      pimi_phi=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
      pimi_phi_mod=pimi_phi+30;
      if (pimi_phi_mod<0)pimi_phi_mod=pimi_phi_mod+360;
      pimi_theta=TMath::ACos(cz[i])*TMath::RadToDeg();
      pimi_vert_corr = pimi_vert + vz_corr(pimi_phi_mod,pimi_theta);
      h2_pimi_phi_theta->Fill(pimi_phi_mod, pimi_theta);
      if(PimiFiducialCut(V3_pimi, &cphil, &cphir)){
        h2_pimi_phi_theta_cut->Fill(pimi_phi_mod, pimi_theta);
        h1_pimi_vert_corr->Fill(el_vert_corr-pimi_vert_corr);
      if(abs(el_vert_corr-pimi_vert_corr)<pimi_vertcut){
        h1_pimi_vert_corr_cut->Fill(el_vert_corr-pimi_vert_corr);
      num_pimi = num_pimi+1;
      num_pi = num_pi+1;
	    index_pimi[num_pimi - 1]=i;
      index_pi[num_pi - 1]=i;
    }
    }
    }
	  }
	if(q[i]!=0 &&  sc[i] > 0 && dc[i]>0 && stat[i] > 0 &&  id[i] == 211) //looking for positive pions, id=211
	  {
      V3_pipl.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);

      bett=p[i]/TMath::Sqrt(p[i]*p[i]+Mpi*Mpi);
	    deltt=sc_t[sc[i]-1]-sc_r[sc[i]-1]/(bett*c*ns_to_s) - tr_time;

	    fsub_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);
	    fsum_pipl->SetParameters(pipl_delt_cutrange,pipl_maxmom);

      pipl_delt_p->Fill(p[i],deltt);

      if(deltt<fsum_pipl->Eval(p[i]) && deltt>fsub_pipl->Eval(p[i])){//Positive pion PID cut
      pipl_delt_p_cut->Fill(p[i],deltt);
      pipl_vert = vz[i];
      pipl_phi=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg();
      pipl_phi_mod=pipl_phi+30;
      if (pipl_phi_mod<0)pipl_phi_mod=pipl_phi_mod+360;
      pipl_theta=TMath::ACos(cz[i])*TMath::RadToDeg();
      pipl_vert_corr = pipl_vert + vz_corr(pipl_phi_mod,pipl_theta);
      h2_pipl_phi_theta->Fill(pipl_phi_mod, pipl_theta);
      if(PiplFiducialCut(V3_pipl, &cphil, &cphir)){
        h1_pipl_vert_corr->Fill(el_vert_corr-pipl_vert_corr);
        h2_pipl_phi_theta_cut->Fill(pipl_phi_mod, pipl_theta);
      if(abs(el_vert_corr-pipl_vert_corr)<pipl_vertcut){
        h1_pipl_vert_corr_cut->Fill(el_vert_corr-pipl_vert_corr);
      num_pipl = num_pipl+1;
      num_pi = num_pi + 1;
	    index_pipl[num_pipl - 1]=i;
      index_pi[num_pi - 1]=i;
    }
    }
	  }
}
  if(ec[i] > 0 && dc[i]<=0  && sc[i]<=0  && stat[i] > 0 &&  q[i] == 0) //looking for neutral particles
  {
    neut_zvert=vz[i];
    neut_yvert=vy[i];
    neut_xvert=vx[i];
    neut_ecx=ech_x[ec[i]-1];
    neut_ecy=ech_y[ec[i]-1];
    neut_ecz=ech_z[ec[i]-1];
    phot_ec_xyz.SetXYZ(ech_x[ec[i]-1],ech_y[ec[i]-1],ech_z[ec[i]-1]);

    neut_ecpath_corr=TMath::Sqrt((neut_ecx-neut_xvert)*(neut_ecx-neut_xvert)+(neut_ecy-neut_yvert)*(neut_ecy-neut_yvert)+(neut_ecz-neut_zvert)*(neut_ecz-neut_zvert));
    neut_ectime_corr=neut_ecpath_corr/(b[i]*c*ns_to_s)-EC_time_offset[make_pair(ftarget,ec_sect[ec[i]-1])];
    neut_beta_corr=neut_ecpath_corr/(neut_ectime_corr*c*ns_to_s);
    neutr_phi_mod=TMath::ATan2(cy[i],cx[i])*TMath::RadToDeg()+30;
    if (neutr_phi_mod<0)neutr_phi_mod=neutr_phi_mod+360;
    V3_phot.SetXYZ(p[i]*cx[i],p[i]*cy[i],p[i]*cz[i]);
    V3_el_angles_SC.SetXYZ(dc_cxsc[dc[ind_em]-1],dc_cysc[dc[ind_em]-1],dc_czsc[dc[ind_em]-1]);
    ec_deltt=ec_t[ec[i]-1]-neut_ecpath_corr/(c*ns_to_s)+EC_time_offset[make_pair(ftarget,ec_sect[ec[i]-1])] - tr_time;
    h1_ec_beta_corr->Fill(neut_beta_corr);
    if(neut_beta_corr >EC_photon_beta[ftarget]) //looking for photons
    {
      h1_ec_beta_corr_cut->Fill(neut_beta_corr);
      h2_neut_costheta_phi->Fill(neutr_phi_mod, cz[i]);
      if(Phot_fid(V3_phot))
      {
        h2_neut_costheta_phi_cut->Fill(neutr_phi_mod, cz[i]);
        ec_num_n = ec_num_n + 1;
        ec_index_n[ec_num_n - 1]=i;
        if(V3_phot.Angle(V3_el)*TMath::RadToDeg()<phot_rad_cut && abs(neutr_phi_mod-el_phi_mod)<phot_e_phidiffcut) ec_radstat_n[ec_num_n - 1]=true;
        photon_ece = TMath::Max(ec_ei[ec[i] - 1] + ec_eo[ec[i] - 1], etot[ec[i] - 1]);
      }
    }
  }
}
  h2_Np_Npi->Fill(num_pi,num_p);
  h3_Npi_Np_Nphoton->Fill(num_pi,num_p,ec_num_n);
  h1_el_vertcorr->Fill(el_vert_corr);
  h1_el_vertuncorr->Fill(el_vert);

  TVector3 V3_q=(V4_beam-V4_el).Vect();
  double q2 = -(V4_beam-V4_el).Mag2();
  double omega= (V4_beam-V4_el).E();
  double Wvar=TMath::Sqrt((m_prot+omega)*(m_prot+omega)-V3_q*V3_q);

  h1_Q2->Fill(q2, 1/Mott_cross_sec);
  h1_omega->Fill(omega, 1/Mott_cross_sec);
  h1_Wvar->Fill(Wvar, 1/Mott_cross_sec);
  h2_Wvar_Q2->Fill(Wvar, q2, 1/Mott_cross_sec);
  if(el_phi_mod<30 && el_phi_mod>-30)
  {
      h2_Wvar_Q2_sect1->Fill(Wvar, q2, 1/Mott_cross_sec);
      h2_omega_Q2_sect1->Fill(omega, q2, 1/Mott_cross_sec);
  }
  h2_omega_Q2->Fill(omega,q2, 1/Mott_cross_sec);
if(el_vert_corr<vert_max[ftarget] && el_vert_corr>vert_min[ftarget]){

  h1_el_vertcorr_cut->Fill(el_vert_corr);

  double N_all = 0;
  double qpi[2]={0};
  bool pi2_stat[2]={false};
  double N_2pi_1p=0,N_1pi_1p[2]={0};
  Float_t pimi_phimin, pimi_phimax;
  double p_kin;
    //Beginning of rotations
    //Requiring one proton
    if (num_p == 1)
      {
      h2_phot_pi_1p->Fill(num_pi, ec_num_n);
      TLorentzVector V4_p;
      TVector3 V3_p;
      V4_p.SetPxPyPzE(p[index_p[0]]*cx[index_p[0]],p[index_p[0]]*cy[index_p[0]],p[index_p[0]]*cz[index_p[0]],TMath::Sqrt(m_prot*m_prot+p[index_p[0]]*p[index_p[0]]));
      V3_p = V4_p.Vect();
      TLorentzVector V4_p_corr;
      double prot_mom_corr = ProtonMomCorrection_He3_4Cell(ftarget,V4_p,prot_vert_corr);
      V4_p_corr.SetPxPyPzE(prot_mom_corr*cx[index_p[0]],prot_mom_corr*cy[index_p[0]],prot_mom_corr*cz[index_p[0]],TMath::Sqrt(m_prot*m_prot+prot_mom_corr*prot_mom_corr));
      p_kin = V4_p_corr.E() - m_prot;
      //Requiring 3 pions, 0 photons
      if(num_pi==3 && ec_num_n==0 && num_n==0)
      {
        TLorentzVector V4_pi[3], V4_total[3];
        TVector3 V3_pi[3];
        double q_pi[3] = {0};
        double p_perp[3] = {0};
        for(int i=0; i<3;i++)
        {
          V4_pi[i].SetPxPyPzE(p[index_pi[i]]*cx[index_pi[i]],p[index_pi[i]]*cy[index_pi[i]],p[index_pi[i]]*cz[index_pi[i]],TMath::Sqrt(p[index_pi[i]]*p[index_pi[i]]+Mpi*Mpi));
          V3_pi[i] = V4_pi[i].Vect();
          q_pi[i] = q[index_pi[i]];
          V4_total[i] = V4_pi[i] + V4_p_corr + V4_el;
          p_perp[i] = TMath::Sqrt(V4_total[i].Px()*V4_total[i].Px()+V4_total[i].Py()*V4_total[i].Py());
        }

        double en_recon1[3];
        double en_recon3[3];

        for(int g=0; g<3; g++)
        {
          en_recon1[g] = V4_el.E() + p_kin + V4_pi[g].E();
          en_recon3[g] = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi[g].E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi[g].E() + 2*V4_el.Rho()*V4_pi[g].Rho()*cos(V3_pi[g].Angle(V3_el)) + Mpi*Mpi)/
                                                (2*(m_neut - eps - V4_el.E() - V4_pi[g].E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi[g].Rho()*cz[index_pi[g]]));
        }
        double N1pi1p[3] = {0};
        double N2pi1p[3] = {0};
        double N3pi1p = 0;

        rot_3pi_1p(V3_pi, q_pi, V3_p, V3_q, N1pi1p, N2pi1p, &N3pi1p, N_tot);

        //fill the histograms here

        if(N3pi1p!=0)
        {
          for(int j=0;j<3;j++)
          {
            h1_rot1_3pi_1p->Fill(en_recon1[j], (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            if(q_pi[j]>0)
            {
              h1_rot1_3pi_1p_pipl->Fill(en_recon1[j], (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot1_3pi_1p_pimi->Fill(en_recon1[j], (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h2_cal_Wvar->Fill(en_recon1[j], Wvar, -(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h1_Q2_sub->Fill(q2,-(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h1_omega_sub->Fill(omega,-(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h1_Wvar_sub->Fill(Wvar,-(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h2_omega_Q2_sub->Fill(omega,q2,-(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              if(p_perp[j] > 0 && p_perp[j] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[j], -(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              if(p_perp[j] > 0.2 && p_perp[j] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[j], -(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              if(p_perp[j] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[j], -(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));

            }

            h1_rot2_3pi_1p->Fill(en_recon2, (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            if(q_pi[j]>0)
            {
              h1_rot2_3pi_1p_pipl->Fill(en_recon2, (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot2_3pi_1p_pimi->Fill(en_recon2, (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            }

            h1_rot3_3pi_1p->Fill(en_recon3[j], (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            if(q_pi[j]>0)
            {
              h1_rot3_3pi_1p_pipl->Fill(en_recon3[j], (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot3_3pi_1p_pimi->Fill(en_recon3[j], (N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
              h2_kin_e_pi_Wvar->Fill(en_recon3[j], Wvar, -(N1pi1p[j]/N3pi1p)*(1/Mott_cross_sec));
            }
          }
        }

        double N_2pion_1prot = 0;
        double N_1pion_1prot[2] = {0};
        int count = 0;
        TVector3 V3_pi1[2];
        double q_pi1[2] = {0};
        for(int i=0; i<3; i++)
        {
          for(int j=0; j<3; j++)
          {
            if(i<j)
            {
              N_1pion_1prot[0] = N_1pion_1prot[1] = 0;
              N_2pion_1prot = 0;
              V3_pi1[0] = V3_pi[i];
              V3_pi1[1] = V3_pi[j];
              q_pi1[0] = q_pi[i];
              q_pi1[1] = q_pi[j];
              V4_total[0] = V4_pi[i] + V4_p_corr + V4_el;
              p_perp[0] = TMath::Sqrt(V4_total[0].Px()*V4_total[0].Px()+V4_total[0].Py()*V4_total[0].Py());
              V4_total[1] = V4_pi[j] + V4_p_corr + V4_el;
              p_perp[1] = TMath::Sqrt(V4_total[1].Px()*V4_total[1].Px()+V4_total[1].Py()*V4_total[1].Py());
              rot_2pi_1p (V3_pi1, q_pi1, V3_p, V3_q, N_1pion_1prot, &N_2pion_1prot, N_tot);
              en_recon1[0] = V4_el.E() + p_kin + V4_pi[i].E();
              en_recon1[1] = V4_el.E() + p_kin + V4_pi[j].E();
              en_recon3[0] = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi[i].E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi[i].E() + 2*V4_el.Rho()*V4_pi[i].Rho()*cos(V3_pi[i].Angle(V3_el)) + Mpi*Mpi)/
                                                    (2*(m_neut - eps - V4_el.E() - V4_pi[i].E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi[i].Rho()*cz[index_pi[i]]));
              en_recon3[1] = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi[j].E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi[j].E() + 2*V4_el.Rho()*V4_pi[j].Rho()*cos(V3_pi[j].Angle(V3_el)) + Mpi*Mpi)/
                                                                                          (2*(m_neut - eps - V4_el.E() - V4_pi[j].E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi[j].Rho()*cz[index_pi[j]]));
              if(N_2pion_1prot!=0 && N3pi1p !=0)
              {
                h1_rot1_3pi_1p->Fill(en_recon1[0], -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));

              if(q_pi1[0]>0)
              {
                h1_rot1_3pi_1p_pipl->Fill(en_recon1[0], -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot1_3pi_1p_pimi->Fill(en_recon1[0], -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h1_Q2_sub->Fill(q2,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h1_omega_sub->Fill(omega,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h1_Wvar_sub->Fill(Wvar,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_omega_Q2_sub->Fill(omega,q2,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              h1_rot1_3pi_1p->Fill(en_recon1[1], -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              if(q_pi1[1]>0)
              {
                h1_rot1_3pi_1p_pipl->Fill(en_recon1[1], -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot1_3pi_1p_pimi->Fill(en_recon1[1], -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h1_Q2_sub->Fill(q2,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h1_omega_sub->Fill(omega,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h1_Wvar_sub->Fill(Wvar,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_omega_Q2_sub->Fill(omega,q2,(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              h1_rot2_3pi_1p->Fill(en_recon2, -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              if(q_pi1[0]>0)
              {
                h1_rot2_3pi_1p_pipl->Fill(en_recon2, -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot2_3pi_1p_pimi->Fill(en_recon2, -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              h1_rot2_3pi_1p->Fill(en_recon2, -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              if(q_pi1[1]>0)
              {
                h1_rot2_3pi_1p_pipl->Fill(en_recon2, -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot2_3pi_1p_pimi->Fill(en_recon2, -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              h1_rot3_3pi_1p->Fill(en_recon3[0], -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              if(q_pi1[0]>0)
              {
                h1_rot3_3pi_1p_pipl->Fill(en_recon3[0], -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot3_3pi_1p_pimi->Fill(en_recon3[0], -(N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, (N_1pion_1prot[0]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              h1_rot3_3pi_1p->Fill(en_recon3[1], -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              if(q_pi1[1]>0)
              {
                h1_rot3_3pi_1p_pipl->Fill(en_recon3[1], -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot3_3pi_1p_pimi->Fill(en_recon3[1], -(N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
                h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, (N_1pion_1prot[1]/N_2pion_1prot)*(N2pi1p[count]/N3pi1p)*(1/Mott_cross_sec));
              }
            }
              count=count+1;
            }
          }
        }
      }//end of 3pi 0 photon statetment
      //Requiring 2 pions, 1 photon
      if (num_pi==2 && ec_num_n==1)
      {
        TLorentzVector V4_pi[2], V4_total[2];
        TVector3 V3_pi[2];
        double p_perp[2] = {0};
        V4_pi[0].SetPxPyPzE(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]],TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
        V4_pi[1].SetPxPyPzE(p[index_pi[1]]*cx[index_pi[1]],p[index_pi[1]]*cy[index_pi[1]],p[index_pi[1]]*cz[index_pi[1]],TMath::Sqrt(p[index_pi[1]]*p[index_pi[1]]+Mpi*Mpi));
        V3_pi[0] = V4_pi[0].Vect();
        V3_pi[1] = V4_pi[1].Vect();
        qpi[0] = q[index_pi[0]];
        qpi[1] = q[index_pi[1]];
        V4_total[0] = V4_pi[0] + V4_p_corr + V4_el;
        p_perp[0] = TMath::Sqrt(V4_total[0].Px()*V4_total[0].Px()+V4_total[0].Py()*V4_total[0].Py());
        V4_total[1] = V4_pi[1] + V4_p_corr + V4_el;
        p_perp[1] = TMath::Sqrt(V4_total[1].Px()*V4_total[1].Px()+V4_total[1].Py()*V4_total[1].Py());
        TVector3 V3_phot(p[ec_index_n[0]]*cx[ec_index_n[0]],p[ec_index_n[0]]*cy[ec_index_n[0]],p[ec_index_n[0]]*cz[ec_index_n[0]]);

        double en_recon1[2];
        double en_recon3[2];

        for(int g=0; g<2; g++)
        {
          en_recon1[g] = V4_el.E() + p_kin + V4_pi[g].E();
          en_recon3[g] = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi[g].E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi[g].E() + 2*V4_el.Rho()*V4_pi[g].Rho()*cos(V3_pi[g].Angle(V3_el)) + Mpi*Mpi)/
                                                (2*(m_neut - eps - V4_el.E() - V4_pi[g].E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi[g].Rho()*cz[index_pi[g]]));
        }
        double N_1pi_1p_0phot[2] = {0};
        double N_1pi_1p_1phot[2] = {0};
        double N_2pi_1p_0phot = 0;
        double N_2pi_1p_1phot = 0;

        rot_1phot_2pi_1p(V3_phot, V3_pi, qpi, V3_p, V3_q, ec_radstat_n[0], N_1pi_1p_0phot, N_1pi_1p_1phot, &N_2pi_1p_0phot, &N_2pi_1p_1phot, N_tot);

        if(N_2pi_1p_1phot!=0){
            // First reconstruction method
            h1_rot1_2pi_1p_1phot->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(qpi[0]>0)
            {
              h1_rot1_2pi_1p_1phot_pipl->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot1_2pi_1p_1phot_pimi->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h1_Q2_sub->Fill(q2,-(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h1_omega_sub->Fill(omega,-(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            h1_rot1_2pi_1p_1phot->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(qpi[1]>0)
            {
              h1_rot1_2pi_1p_1phot_pipl->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot1_2pi_1p_1phot_pimi->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h1_Q2_sub->Fill(q2,-(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h1_omega_sub->Fill(omega,-(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            //Second reconstruction method
            h1_rot2_2pi_1p_1phot->Fill(en_recon2, (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(qpi[0]>0)
            {
              h1_rot2_2pi_1p_1phot_pipl->Fill(en_recon2, (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot2_2pi_1p_1phot_pimi->Fill(en_recon2, (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            h1_rot2_2pi_1p_1phot->Fill(en_recon2, (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(qpi[1]>0)
            {
              h1_rot2_2pi_1p_1phot_pipl->Fill(en_recon2, (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot2_2pi_1p_1phot_pimi->Fill(en_recon2, (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            //Third reconstruction method
            h1_rot3_2pi_1p_1phot->Fill(en_recon3[0], (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(qpi[0]>0)
            {
              h1_rot3_2pi_1p_1phot_pipl->Fill(en_recon3[0], (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot3_2pi_1p_1phot_pimi->Fill(en_recon3[0], (N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, -(N_1pi_1p_0phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            h1_rot3_2pi_1p_1phot->Fill(en_recon3[1], (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(qpi[1]>0)
            {
              h1_rot3_2pi_1p_1phot_pipl->Fill(en_recon3[1], (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
            else
            {
              h1_rot3_2pi_1p_1phot_pimi->Fill(en_recon3[1], (N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
              h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, -(N_1pi_1p_0phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            }
      }
      double N1pi1p0phot[2] = {0};
      double N1pi1p1phot[2] = {0};
      rot_1phot_1pi_1p(V3_phot, V3_pi[0], qpi[0], V3_p, V3_q, ec_radstat_n[0], &N1pi1p0phot[0], &N1pi1p1phot[0], N_tot);
      rot_1phot_1pi_1p(V3_phot, V3_pi[1], qpi[1], V3_p, V3_q, ec_radstat_n[0], &N1pi1p0phot[1], &N1pi1p1phot[1], N_tot);

      if(N_2pi_1p_1phot!=0 && N1pi1p1phot[0]!=0 && N1pi1p1phot[1]!=0){
          // First reconstruction method
          h1_rot1_2pi_1p_1phot->Fill(en_recon1[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(qpi[0]>0)
          {
            h1_rot1_2pi_1p_1phot_pipl->Fill(en_recon1[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot1_2pi_1p_1phot_pimi->Fill(en_recon1[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h1_Q2_sub->Fill(q2,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h1_omega_sub->Fill(omega,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h1_Wvar_sub->Fill(Wvar,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_Wvar_Q2_sub->Fill(Wvar, q2,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_omega_Q2_sub->Fill(omega,q2,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          h1_rot1_2pi_1p_1phot->Fill(en_recon1[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(qpi[1]>0)
          {
            h1_rot1_2pi_1p_1phot_pipl->Fill(en_recon1[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot1_2pi_1p_1phot_pimi->Fill(en_recon1[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h1_Q2_sub->Fill(q2,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h1_omega_sub->Fill(omega,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h1_Wvar_sub->Fill(Wvar,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_Wvar_Q2_sub->Fill(Wvar, q2,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_omega_Q2_sub->Fill(omega,q2,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          //Second reconstruction method
          h1_rot2_2pi_1p_1phot->Fill(en_recon2, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(qpi[0]>0)
          {
            h1_rot2_2pi_1p_1phot_pipl->Fill(en_recon2, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot2_2pi_1p_1phot_pimi->Fill(en_recon2, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          h1_rot2_2pi_1p_1phot->Fill(en_recon2, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(qpi[1]>0)
          {
            h1_rot2_2pi_1p_1phot_pipl->Fill(en_recon2, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot2_2pi_1p_1phot_pimi->Fill(en_recon2, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          //Third reconstruction method
          h1_rot3_2pi_1p_1phot->Fill(en_recon3[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(qpi[0]>0)
          {
            h1_rot3_2pi_1p_1phot_pipl->Fill(en_recon3[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot3_2pi_1p_1phot_pimi->Fill(en_recon3[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          h1_rot3_2pi_1p_1phot->Fill(en_recon3[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(qpi[1]>0)
          {
            h1_rot3_2pi_1p_1phot_pipl->Fill(en_recon3[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot3_2pi_1p_1phot_pimi->Fill(en_recon3[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
            h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          }
    }
    double N_1pi_1p[2] = {0};
    double N_2pi_1p = 0;
    rot_2pi_1p (V3_pi, qpi, V3_p, V3_q, N_1pi_1p, &N_2pi_1p, N_tot);

    if(N_2pi_1p_1phot!=0 && N_2pi_1p!=0){
        // First reconstruction method
        h1_rot1_2pi_1p_1phot->Fill(en_recon1[0], -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot1_2pi_1p_1phot_pipl->Fill(en_recon1[0], -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_1p_1phot_pimi->Fill(en_recon1[0], -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_1p_1phot->Fill(en_recon1[1], -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot1_2pi_1p_1phot_pipl->Fill(en_recon1[1], -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_1p_1phot_pimi->Fill(en_recon1[1], -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        //Second reconstruction method
        h1_rot2_2pi_1p_1phot->Fill(en_recon2, -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot2_2pi_1p_1phot_pipl->Fill(en_recon2, -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_1p_1phot_pimi->Fill(en_recon2, -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        h1_rot2_2pi_1p_1phot->Fill(en_recon2, -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot2_2pi_1p_1phot_pipl->Fill(en_recon2, -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_1p_1phot_pimi->Fill(en_recon2, -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        //Third reconstruction method
        h1_rot3_2pi_1p_1phot->Fill(en_recon3[0], -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot3_2pi_1p_1phot_pipl->Fill(en_recon3[0], -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_1p_1phot_pimi->Fill(en_recon3[0], -(N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, (N_1pi_1p[0]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        h1_rot3_2pi_1p_1phot->Fill(en_recon3[1], -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot3_2pi_1p_1phot_pipl->Fill(en_recon3[1], -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_1p_1phot_pimi->Fill(en_recon3[1], -(N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, (N_1pi_1p[1]/N_2pi_1p)*(N_2pi_1p_0phot/N_2pi_1p_1phot)*(1/Mott_cross_sec));
        }
  }
    }//end of 2 pi 1 photon statement
    //Requiring 2 pions, 0 photons
      if (num_pi==2 && ec_num_n==0 && num_n==0 ) {
        TLorentzVector V4_pi[2], V4_total[2];
        TVector3 V3_pi[2];
        double p_perp[2] = {0};
        V4_pi[0].SetPxPyPzE(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]],TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
        V4_pi[1].SetPxPyPzE(p[index_pi[1]]*cx[index_pi[1]],p[index_pi[1]]*cy[index_pi[1]],p[index_pi[1]]*cz[index_pi[1]],TMath::Sqrt(p[index_pi[1]]*p[index_pi[1]]+Mpi*Mpi));
        V3_pi[0] = V4_pi[0].Vect();
        V3_pi[1] = V4_pi[1].Vect();
        qpi[0] = q[index_pi[0]];
        qpi[1] = q[index_pi[1]];
        V4_total[0] = V4_pi[0] + V4_p_corr + V4_el;
        p_perp[0] = TMath::Sqrt(V4_total[0].Px()*V4_total[0].Px()+V4_total[0].Py()*V4_total[0].Py());
        V4_total[1] = V4_pi[1] + V4_p_corr + V4_el;
        p_perp[1] = TMath::Sqrt(V4_total[1].Px()*V4_total[1].Px()+V4_total[1].Py()*V4_total[1].Py());

        double en_recon1[2];
        double en_recon3[2];

        for(int g=0; g<2; g++)
        {
          en_recon1[g] = V4_el.E() + p_kin + V4_pi[g].E();
          en_recon3[g] = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi[g].E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi[g].E() + 2*V4_el.Rho()*V4_pi[g].Rho()*cos(V3_pi[g].Angle(V3_el)) + Mpi*Mpi)/
                                                (2*(m_neut - eps - V4_el.E() - V4_pi[g].E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi[g].Rho()*cz[index_pi[g]]));
        }
        rot_2pi_1p (V3_pi, qpi, V3_p, V3_q, N_1pi_1p, &N_2pi_1p, N_tot);

        //fill the histograms here

          if(N_2pi_1p!=0){
              // First reconstruction method
              h1_rot1_2pi_1p->Fill(en_recon1[0], (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              if(qpi[0]>0)
              {
                h1_rot1_2pi_1p_pipl->Fill(en_recon1[0], (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot1_2pi_1p_pimi->Fill(en_recon1[0], (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h1_Q2_sub->Fill(q2,-(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h1_omega_sub->Fill(omega,-(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              h1_rot1_2pi_1p->Fill(en_recon1[1], (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              if(qpi[1]>0)
              {
                h1_rot1_2pi_1p_pipl->Fill(en_recon1[1], (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot1_2pi_1p_pimi->Fill(en_recon1[1], (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h1_Q2_sub->Fill(q2,-(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h1_omega_sub->Fill(omega,-(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              //Second reconstruction method
              h1_rot2_2pi_1p->Fill(en_recon2, (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              if(qpi[0]>0)
              {
                h1_rot2_2pi_1p_pipl->Fill(en_recon2, (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot2_2pi_1p_pimi->Fill(en_recon2, (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              h1_rot2_2pi_1p->Fill(en_recon2, (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              if(qpi[1]>0)
              {
                h1_rot2_2pi_1p_pipl->Fill(en_recon2, (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot2_2pi_1p_pimi->Fill(en_recon2, (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              //Third reconstruction method
              h1_rot3_2pi_1p->Fill(en_recon3[0], (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              if(qpi[0]>0)
              {
                h1_rot3_2pi_1p_pipl->Fill(en_recon3[0], (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot3_2pi_1p_pimi->Fill(en_recon3[0], (N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, -(N_1pi_1p[0]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              h1_rot3_2pi_1p->Fill(en_recon3[1], (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              if(qpi[1]>0)
              {
                h1_rot3_2pi_1p_pipl->Fill(en_recon3[1], (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot3_2pi_1p_pimi->Fill(en_recon3[1], (N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
                h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, -(N_1pi_1p[1]/N_2pi_1p)*(1/Mott_cross_sec));
              }
        }

}//end of 2pi 0 photon statement
    //Requiring 1 pion, 0 photons
    if(num_pi==1 && ec_num_n==0 && num_n==0)
	  {
      TLorentzVector V4_total;
      double p_perp = 0;
      if(num_pipl==1 )
      {
        TLorentzVector V4_pi(p[index_pipl[0]]*cx[index_pipl[0]],p[index_pipl[0]]*cy[index_pipl[0]],p[index_pipl[0]]*cz[index_pipl[0]], TMath::Sqrt(p[index_pipl[0]]*p[index_pipl[0]]+Mpi*Mpi));
        TVector3 V3_pi = V4_pi.Vect();
        // First reconstruction method
        double en_recon1 = V4_el.E() + p_kin + V4_pi.E();
        h1_en_recon1->Fill(en_recon1, 1/Mott_cross_sec);
        h1_en_recon1_pipl->Fill(en_recon1, 1/Mott_cross_sec);

        //Second reconstruction method
        h1_en_recon2->Fill(en_recon2, 1/Mott_cross_sec);
        h1_en_recon2_pipl->Fill(en_recon2, 1/Mott_cross_sec);

        //Third reconstruction method
        double en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                              (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pipl[0]]));
        h1_en_recon3->Fill(en_recon3, 1/Mott_cross_sec);
        h1_en_recon3_pipl->Fill(en_recon3, 1/Mott_cross_sec);
      }
      if(num_pimi==1)
      {
        TLorentzVector V4_pi(p[index_pimi[0]]*cx[index_pimi[0]],p[index_pimi[0]]*cy[index_pimi[0]],p[index_pimi[0]]*cz[index_pimi[0]], TMath::Sqrt(p[index_pimi[0]]*p[index_pimi[0]]+Mpi*Mpi));
        TVector3 V3_pi = V4_pi.Vect();
        V4_total = V4_pi + V4_p_corr + V4_el;
        p_perp = TMath::Sqrt(V4_total.Px()*V4_total.Px()+V4_total.Py()*V4_total.Py());

        h1_p_perp_cut->Fill(p_perp, 1/Mott_cross_sec);
        TVector3 V3_q=(V4_beam-V4_el).Vect();
        double q2 = -(V4_beam-V4_el).Mag2();
        double omega= (V4_beam-V4_el).E();
        double Wvar=TMath::Sqrt((m_prot+omega)*(m_prot+omega)-V3_q*V3_q);

        h1_Q2_cut->Fill(q2, 1/Mott_cross_sec);
        h1_omega_cut->Fill(omega, 1/Mott_cross_sec);
        h1_Wvar_cut->Fill(Wvar, 1/Mott_cross_sec);
        // First reconstruction method
        double en_recon1 = V4_el.E() + p_kin + V4_pi.E();
        h1_en_recon1->Fill(en_recon1, 1/Mott_cross_sec);
        h1_en_recon1_pimi->Fill(en_recon1, 1/Mott_cross_sec);
        h2_cal_Wvar->Fill(en_recon1, Wvar, 1/Mott_cross_sec);
        //Second reconstruction method
        h1_en_recon2->Fill(en_recon2, 1/Mott_cross_sec);
        h1_en_recon2_pimi->Fill(en_recon2, 1/Mott_cross_sec);
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, 1/Mott_cross_sec);
        //Third reconstruction method
        double en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                              (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pimi[0]]));
        h1_en_recon3->Fill(en_recon3, 1/Mott_cross_sec);
        h1_en_recon3_pimi->Fill(en_recon3, 1/Mott_cross_sec);
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, 1/Mott_cross_sec);
        h1_Q2_sub->Fill(q2,1/Mott_cross_sec);
        h1_omega_sub->Fill(omega,1/Mott_cross_sec);
        h1_Wvar_sub->Fill(Wvar,1/Mott_cross_sec);
        h2_Wvar_Q2_sub->Fill(Wvar, q2,1/Mott_cross_sec);
        h2_omega_Q2_sub->Fill(omega,q2,1/Mott_cross_sec);
        double p_perp = TMath::Sqrt(V4_total.Px()*V4_total.Px()+V4_total.Py()*V4_total.Py());
        h1_p_perp->Fill(p_perp,1/Mott_cross_sec);
        h1_p_perp_sub->Fill(p_perp,1/Mott_cross_sec);
        h2_p_perp_Wvar->Fill(p_perp,Wvar,1/Mott_cross_sec);
        h2_pperp_cal->Fill(p_perp,en_recon1,1/Mott_cross_sec);
        h2_pperp_kin->Fill(p_perp,en_recon3,1/Mott_cross_sec);
        if(p_perp > 0 && p_perp < 0.2) h1_cal_p_slice1->Fill(en_recon1, 1/Mott_cross_sec);
        if(p_perp > 0.2 && p_perp < 0.4) h1_cal_p_slice2->Fill(en_recon1, 1/Mott_cross_sec);
        if(p_perp > 0.4) h1_cal_p_slice3->Fill(en_recon1, 1/Mott_cross_sec);
        if(p_perp > 0 && p_perp < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1, 1/Mott_cross_sec);
        if(p_perp > 0.2 && p_perp < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1, 1/Mott_cross_sec);
        if(p_perp > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1, 1/Mott_cross_sec);
    }
  }//end of 1pi 0 photon statement
  //Requiring 1 pion, 1 photon
  if(num_pi==1 && ec_num_n==1)
  {
    TLorentzVector V4_total;
    TLorentzVector V4_pi(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]], TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
    TVector3 V3_pi = V4_pi.Vect();
    V4_total = V4_pi + V4_p_corr + V4_el;
    double p_perp = TMath::Sqrt(V4_total.Px()*V4_total.Px()+V4_total.Py()*V4_total.Py());
    TVector3 V3_phot(p[ec_index_n[0]]*cx[ec_index_n[0]],p[ec_index_n[0]]*cy[ec_index_n[0]],p[ec_index_n[0]]*cz[ec_index_n[0]]);
    double N_1pi_1p_0phot = 0;
    double N_1pi_1p_1phot = 0;

    double en_recon1 = V4_el.E() + p_kin + V4_pi.E();
    double en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                          (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pipl[0]]));

    rot_1phot_1pi_1p(V3_phot, V3_pi, q[index_pi[0]], V3_p, V3_q, ec_radstat_n[0], &N_1pi_1p_0phot, &N_1pi_1p_1phot, N_tot);

    //fill histograms here
    if(N_1pi_1p_1phot!=0)
    {
      h1_rot1_1pi_1p_1phot->Fill(en_recon1, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot1_1pi_1p_1phot_pipl->Fill(en_recon1, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot1_1pi_1p_1phot_pimi->Fill(en_recon1, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,-(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,-(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1, Wvar, -(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp,-(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        if(p_perp > 0 && p_perp < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1, -(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        if(p_perp > 0.2 && p_perp < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1, -(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        if(p_perp > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1, -(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot2_1pi_1p_1phot->Fill(en_recon2, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot2_1pi_1p_1phot_pipl->Fill(en_recon2, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot2_1pi_1p_1phot_pimi->Fill(en_recon2, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot3_1pi_1p_1phot->Fill(en_recon3, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot3_1pi_1p_1phot_pipl->Fill(en_recon3, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot3_1pi_1p_1phot_pimi->Fill(en_recon3, (N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_1p_0phot/N_1pi_1p_1phot)*(1/Mott_cross_sec));
      }
    }
  }//end of 1 pi 1 photon statement
  //Requiring 1 pion, 2 photons
  if(num_pi==1 && ec_num_n==2)
  {
    TLorentzVector V4_total;
    TLorentzVector V4_pi(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]], TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
    TVector3 V3_pi = V4_pi.Vect();
    V4_total = V4_pi + V4_p_corr + V4_el;
    double p_perp = TMath::Sqrt(V4_total.Px()*V4_total.Px()+V4_total.Py()*V4_total.Py());
    TVector3 V3_phot[2];
    V3_phot[0].SetXYZ(p[ec_index_n[0]]*cx[ec_index_n[0]],p[ec_index_n[0]]*cy[ec_index_n[0]],p[ec_index_n[0]]*cz[ec_index_n[0]]);
    V3_phot[1].SetXYZ(p[ec_index_n[1]]*cx[ec_index_n[1]],p[ec_index_n[1]]*cy[ec_index_n[1]],p[ec_index_n[1]]*cz[ec_index_n[1]]);
    double N_1pi_1p_0phot = 0;
    double N_1pi_1p_1phot = 0;
    double N_1pi_1p_2phot = 0;
    bool radstat[2];
    radstat[0] = ec_radstat_n[0];
    radstat[1] = ec_radstat_n[1];

    double en_recon1 = V4_el.E() + p_kin + V4_pi.E();
    double en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                          (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pipl[0]]));

    rot_2phot_1pi_1p(V3_phot, V3_pi, q[index_pi[0]], V3_p, V3_q, radstat, &N_1pi_1p_0phot, &N_1pi_1p_1phot, &N_1pi_1p_2phot, N_tot);

    if(N_1pi_1p_2phot!=0)
    {
      h1_rot1_1pi_1p_2phot->Fill(en_recon1, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot1_1pi_1p_2phot_pipl->Fill(en_recon1, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot1_1pi_1p_2phot_pimi->Fill(en_recon1, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,-(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,-(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1, Wvar, -(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp,-(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        if(p_perp > 0 && p_perp < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1, -(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        if(p_perp > 0.2 && p_perp < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1, -(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        if(p_perp > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1, -(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));

      }

      h1_rot2_1pi_1p_2phot->Fill(en_recon2, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot2_1pi_1p_2phot_pipl->Fill(en_recon2, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot2_1pi_1p_2phot_pimi->Fill(en_recon2, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }

      h1_rot3_1pi_1p_2phot->Fill(en_recon3, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot3_1pi_1p_2phot_pipl->Fill(en_recon3, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot3_1pi_1p_2phot_pimi->Fill(en_recon3, (N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_1p_0phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
    }
    double N1pi1p0phot = 0;
    double N1pi1p1phot = 0;
    rot_1phot_1pi_1p(V3_phot[0], V3_pi, q[index_pi[0]], V3_p, V3_q, ec_radstat_n[0], &N_1pi_1p_0phot, &N1pi1p1phot, N_tot);
    rot_1phot_1pi_1p(V3_phot[1], V3_pi, q[index_pi[0]], V3_p, V3_q, ec_radstat_n[1], &N_1pi_1p_0phot, &N1pi1p1phot, N_tot);
    if(N_1pi_1p_2phot!=0 && N1pi1p1phot!=0)
    {
      h1_rot1_1pi_1p_2phot->Fill(en_recon1, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot1_1pi_1p_2phot_pipl->Fill(en_recon1, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot1_1pi_1p_2phot_pimi->Fill(en_recon1, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1, Wvar, (N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp,(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        if(p_perp > 0 && p_perp < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1, (N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        if(p_perp > 0.2 && p_perp < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1, (N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        if(p_perp > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1, (N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }

      h1_rot2_1pi_1p_2phot->Fill(en_recon2, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot2_1pi_1p_2phot_pipl->Fill(en_recon2, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot2_1pi_1p_2phot_pimi->Fill(en_recon2, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }

      h1_rot3_1pi_1p_2phot->Fill(en_recon3, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot3_1pi_1p_2phot_pipl->Fill(en_recon3, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot3_1pi_1p_2phot_pimi->Fill(en_recon3, -(N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N1pi1p0phot/N1pi1p1phot)*(N_1pi_1p_1phot/N_1pi_1p_2phot)*(1/Mott_cross_sec));
      }
    }
  }//end of 1 pi 2 photon statement
}//end of 1proton statement

TVector3 V3_1pi_rot, V3_p_rot[2], V3_pi, V3_p[2];
TLorentzVector V4_pi, V4_p[2];
double p2_kin[2];

//Outer loop requiring 2 protons
if(num_p == 2)
{
  h2_phot_pi_2p->Fill(num_pi, ec_num_n);
  double prot_phi[2];
  double prot_phi_mod[2];
  double prot_theta[2];
  double prot_vert[2];
  double prot_vert_corr[2];
  double prot_mom_corr[2];
  TLorentzVector V4_p_corr[2];
  V4_p[0].SetPxPyPzE(p[index_p[0]]*cx[index_p[0]],p[index_p[0]]*cy[index_p[0]],p[index_p[0]]*cz[index_p[0]],TMath::Sqrt(m_prot*m_prot+p[index_p[0]]*p[index_p[0]]));
  V4_p[1].SetPxPyPzE(p[index_p[1]]*cx[index_p[1]],p[index_p[1]]*cy[index_p[1]],p[index_p[1]]*cz[index_p[1]],TMath::Sqrt(m_prot*m_prot+p[index_p[1]]*p[index_p[1]]));
  V3_p[0] = V4_p[0].Vect();
  V3_p[1] = V4_p[1].Vect();
  for(int i=0; i<2; i++)
  {
    prot_vert[i] = vz[index_p[i]];
    prot_phi[i]=TMath::ATan2(cy[index_p[i]],cx[index_p[i]])*TMath::RadToDeg();
    prot_phi_mod[i]=prot_phi[i]+30;
    if (prot_phi_mod[i]<0)prot_phi_mod[i]=prot_phi_mod[i]+360;
    prot_theta[i]=TMath::ACos(cz[index_p[i]])*TMath::RadToDeg();
    prot_vert_corr[i] = prot_vert[i] + vz_corr(prot_phi_mod[i],prot_theta[i]);
    prot_mom_corr[i] = ProtonMomCorrection_He3_4Cell(ftarget,V4_p[i],prot_vert_corr[i]);
    V4_p_corr[i].SetPxPyPzE(prot_mom_corr[i]*cx[index_p[i]],prot_mom_corr[i]*cy[index_p[i]],prot_mom_corr[i]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+prot_mom_corr[i]*prot_mom_corr[i]));
    p2_kin[i] = V4_p_corr[i].E() - m_prot;
  }
  //Requiring 1 pion, 1 photon
  if(num_pi==1 && ec_num_n == 1)
  {
    double p_perp[2] = {0};
    TLorentzVector V4_total[2];
    TLorentzVector V4_pi(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]], TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
    TVector3 V3_pi = V4_pi.Vect();
    V4_total[0] = V4_pi + V4_p_corr[0] + V4_el;
    V4_total[1] = V4_pi + V4_p_corr[1] + V4_el;
    p_perp[0] = TMath::Sqrt(V4_total[0].Px()*V4_total[0].Px()+V4_total[0].Py()*V4_total[0].Py());
    p_perp[1] = TMath::Sqrt(V4_total[1].Px()*V4_total[1].Px()+V4_total[1].Py()*V4_total[1].Py());
    TVector3 V3_phot(p[ec_index_n[0]]*cx[ec_index_n[0]],p[ec_index_n[0]]*cy[ec_index_n[0]],p[ec_index_n[0]]*cz[ec_index_n[0]]);
    double N_1pi_1p_0phot[2] = {0};
    double N_1pi_1p_1phot[2] = {0};
    double N_1pi_2p_1phot = 0;
    double N_1pi_2p_0phot = 0;

    double en_recon1[2];
    en_recon1[0] = V4_el.E() + p2_kin[0] + V4_pi.E();
    en_recon1[1] = V4_el.E() + p2_kin[1] + V4_pi.E();
    double en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                          (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pipl[0]]));

    rot_1phot_1pi_2p(V3_phot, V3_pi, q[index_pi[0]], V3_p, V3_q, ec_radstat_n[0], N_1pi_1p_0phot, N_1pi_1p_1phot, &N_1pi_2p_0phot, &N_1pi_2p_1phot, N_tot);

    //fill histograms here
    if(N_1pi_2p_1phot!=0)
    {
      h1_rot1_1pi_2p_1phot->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      h1_rot1_1pi_2p_1phot->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot1_1pi_2p_1phot_pipl->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_2p_1phot_pipl->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot1_1pi_2p_1phot_pimi->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_2p_1phot_pimi->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,-(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,-(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,-(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,-(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[0],-(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_1p_0phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[1],-(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot2_1pi_2p_1phot->Fill(en_recon2, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot2_1pi_2p_1phot_pipl->Fill(en_recon2, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot2_1pi_2p_1phot_pimi->Fill(en_recon2, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, -((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot3_1pi_2p_1phot->Fill(en_recon3, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot3_1pi_2p_1phot_pipl->Fill(en_recon3, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot3_1pi_2p_1phot_pimi->Fill(en_recon3, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1])/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
    }
    double N1pi1p0phot[2] = {0};
    double N1pi1p1phot[2] = {0};
    rot_1phot_1pi_1p(V3_phot, V3_pi, q[index_pi[0]], V3_p[0], V3_q, ec_radstat_n[0], &N1pi1p0phot[0], &N1pi1p1phot[0], N_tot);
    rot_1phot_1pi_1p(V3_phot, V3_pi, q[index_pi[0]], V3_p[1], V3_q, ec_radstat_n[0], &N1pi1p0phot[1], &N1pi1p1phot[1], N_tot);
    if(N_1pi_2p_1phot!=0 && N1pi1p1phot[0]!=0 && N1pi1p1phot[1]!=0)
    {
      h1_rot1_1pi_2p_1phot->Fill(en_recon1[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      h1_rot1_1pi_2p_1phot->Fill(en_recon1[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot1_1pi_2p_1phot_pipl->Fill(en_recon1[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_2p_1phot_pipl->Fill(en_recon1[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot1_1pi_2p_1phot_pimi->Fill(en_recon1[0], -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_2p_1phot_pimi->Fill(en_recon1[1], -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[0],(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[1],(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot2_1pi_2p_1phot->Fill(en_recon2, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      h1_rot2_1pi_2p_1phot->Fill(en_recon2, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot2_1pi_2p_1phot_pipl->Fill(en_recon2, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot2_1pi_2p_1phot_pipl->Fill(en_recon2, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot2_1pi_2p_1phot_pimi->Fill(en_recon2, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot2_1pi_2p_1phot_pimi->Fill(en_recon2, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot3_1pi_2p_1phot->Fill(en_recon3, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      h1_rot3_1pi_2p_1phot->Fill(en_recon3, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot3_1pi_2p_1phot_pipl->Fill(en_recon3, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot3_1pi_2p_1phot_pipl->Fill(en_recon3, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot3_1pi_2p_1phot_pimi->Fill(en_recon3, -(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot3_1pi_2p_1phot_pimi->Fill(en_recon3, -(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot[0]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot[1]/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
    }
    double N_1pi_1p[2] = {0};
    double N_1pi_2p = 0;
    rot_1pi_2p(V3_pi, q[index_pi[0]], V3_p, V3_q, N_1pi_1p, &N_1pi_2p, N_tot);

    if(N_1pi_2p_1phot!=0 && N_1pi_2p!=0)
    {
      h1_rot1_1pi_2p_1phot->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      h1_rot1_1pi_2p_1phot->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot1_1pi_2p_1phot_pipl->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_2p_1phot_pipl->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot1_1pi_2p_1phot_pimi->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_2p_1phot_pimi->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[0],(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[1],(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot2_1pi_2p_1phot->Fill(en_recon2, -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      h1_rot2_1pi_2p_1phot->Fill(en_recon2, -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot2_1pi_2p_1phot_pipl->Fill(en_recon2, -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot2_1pi_2p_1phot_pipl->Fill(en_recon2, -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot2_1pi_2p_1phot_pimi->Fill(en_recon2, -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot2_1pi_2p_1phot_pimi->Fill(en_recon2, -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot3_1pi_2p_1phot->Fill(en_recon3, -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      h1_rot3_1pi_2p_1phot->Fill(en_recon3, -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot3_1pi_2p_1phot_pipl->Fill(en_recon3, -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot3_1pi_2p_1phot_pipl->Fill(en_recon3, -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot3_1pi_2p_1phot_pimi->Fill(en_recon3, -(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h1_rot3_1pi_2p_1phot_pimi->Fill(en_recon3, -(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot/N_1pi_2p_1phot)*(1/Mott_cross_sec));
      }
    }
  }
  //end of 1pi 1phot statement
  //Requiring 1 pion, 0 photon
  if (num_pi==1 && ec_num_n==0 && num_n==0) {
    double p_perp[2] = {0};
    TLorentzVector V4_total[2];
    V4_pi.SetPxPyPzE(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]],TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
    V3_pi = V4_pi.Vect();
    V4_total[0] = V4_pi + V4_p_corr[0] + V4_el;
    p_perp[0] = TMath::Sqrt(V4_total[0].Px()*V4_total[0].Px()+V4_total[0].Py()*V4_total[0].Py());
    V4_total[1] = V4_pi + V4_p_corr[1] + V4_el;
    p_perp[1] = TMath::Sqrt(V4_total[1].Px()*V4_total[1].Px()+V4_total[1].Py()*V4_total[1].Py());
    double qpi = q[index_pi[0]];

    N_1pi_1p[0]=N_1pi_1p[1]=0;
    double N_1pi_2p = 0;
    double en_recon1[2];
    double en_recon3;
    double rot_angle;

    en_recon1[0] = V4_el.E() + p2_kin[0] + V4_pi.E();
    en_recon1[1] = V4_el.E() + p2_kin[1] + V4_pi.E();
    en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                            (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pi[0]]));


rot_1pi_2p(V3_pi, qpi, V3_p, V3_q, N_1pi_1p, &N_1pi_2p, N_tot);
    //fill the histograms here
      if(N_1pi_2p!=0){
          // First reconstruction method
          h1_rot1_1pi_2p->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot1_1pi_2p_pipl->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot1_1pi_2p_pimi->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_Q2_sub->Fill(q2,-(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_omega_sub->Fill(omega,-(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_p_perp_sub->Fill(p_perp[0],-(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          h1_rot1_1pi_2p->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot1_1pi_2p_pipl->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot1_1pi_2p_pimi->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_Q2_sub->Fill(q2,-(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_omega_sub->Fill(omega,-(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h1_p_perp_sub->Fill(p_perp[1],-(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          //Second reconstruction method
          h1_rot2_1pi_2p->Fill(en_recon2, (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot2_1pi_2p_pipl->Fill(en_recon2, (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot2_1pi_2p_pimi->Fill(en_recon2, (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          h1_rot2_1pi_2p->Fill(en_recon2, (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot2_1pi_2p_pipl->Fill(en_recon2, (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot2_1pi_2p_pimi->Fill(en_recon2, (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          //Third reconstruction method
          h1_rot3_1pi_2p->Fill(en_recon3, (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot3_1pi_2p_pipl->Fill(en_recon3, (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot3_1pi_2p_pimi->Fill(en_recon3, (N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          h1_rot3_1pi_2p->Fill(en_recon3, (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot3_1pi_2p_pipl->Fill(en_recon3, (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot3_1pi_2p_pimi->Fill(en_recon3, (N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
            h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
          }
    }
}//end of 1pi 0 photon statement

TLorentzVector V4_pi[2];
TVector3 V3_pi[2], V3_pi_rot[2];
double N_all2 = 0;
double N_p1_pi1 = 0, N_p1_pi2 = 0, N_p2_pi2 = 0, N_p2_pi1 = 0;
//Requiring 2 pions, 0 photons
if(num_pi == 2 && ec_num_n==0 && num_n==0)
{
  double p_perp[4];
  TLorentzVector V4_total[4];
  V4_pi[0].SetPxPyPzE(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]],TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
  V4_pi[1].SetPxPyPzE(p[index_pi[1]]*cx[index_pi[1]],p[index_pi[1]]*cy[index_pi[1]],p[index_pi[1]]*cz[index_pi[1]],TMath::Sqrt(p[index_pi[1]]*p[index_pi[1]]+Mpi*Mpi));
  V3_pi[0] = V4_pi[0].Vect();
  V3_pi[1] = V4_pi[1].Vect();
  V4_total[0] = V4_pi[0] + V4_p_corr[0] + V4_el;
  p_perp[0] = TMath::Sqrt(V4_total[0].Px()*V4_total[0].Px()+V4_total[0].Py()*V4_total[0].Py());
  V4_total[1] = V4_pi[1] + V4_p_corr[0] + V4_el;
  p_perp[1] = TMath::Sqrt(V4_total[1].Px()*V4_total[1].Px()+V4_total[1].Py()*V4_total[1].Py());
  V4_total[2] = V4_pi[0] + V4_p_corr[1] + V4_el;
  p_perp[2] = TMath::Sqrt(V4_total[2].Px()*V4_total[2].Px()+V4_total[2].Py()*V4_total[2].Py());
  V4_total[3] = V4_pi[1] + V4_p_corr[1] + V4_el;
  p_perp[3] = TMath::Sqrt(V4_total[3].Px()*V4_total[3].Px()+V4_total[3].Py()*V4_total[3].Py());
  double N_1pi_1p_[4]={0}, N_1pi_2p_[2] = {0}, N_2pi_1p_[2] = {0};
  double en_recon1[4];
  double en_recon3[2];
  double rot_angle;
  double qpi[2];
  qpi[0] = q[index_pi[0]];
  qpi[1] = q[index_pi[1]];

  en_recon1[0] = V4_el.E() + p2_kin[0] + V4_pi[0].E();
  en_recon1[1] = V4_el.E() + p2_kin[0] + V4_pi[1].E();
  en_recon1[2] = V4_el.E() + p2_kin[1] + V4_pi[0].E();
  en_recon1[3] = V4_el.E() + p2_kin[1] + V4_pi[1].E();
for(int g = 0; g<2;g++){
    en_recon3[g] = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi[g].E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi[g].E() + 2*V4_el.Rho()*V4_pi[g].Rho()*cos(V3_pi[g].Angle(V3_el)) + Mpi*Mpi)/
                                          (2*(m_neut - eps - V4_el.E() - V4_pi[g].E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi[g].Rho()*cz[index_pi[g]]));
}


for(int g=0; g<N_tot; g++){

rot_angle=gRandom->Uniform(0,2*TMath::Pi());

    V3_pi_rot[0]=V3_pi[0];
    V3_pi_rot[1]=V3_pi[1];
    V3_p_rot[0]=V3_p[0];
    V3_p_rot[1]=V3_p[1];
    V3_pi_rot[0].Rotate(rot_angle,V3_q);
    V3_pi_rot[1].Rotate(rot_angle,V3_q);
    V3_p_rot[0].Rotate(rot_angle,V3_q);
    V3_p_rot[1].Rotate(rot_angle,V3_q);

    for(int z=0;z<2;z++){
      if(q[index_pi[z]]>0) pi2_stat[z]=PiplFiducialCut(V3_pi_rot[z], &cphil, &cphir);
      else  pi2_stat[z]=PimiFiducialCut(V3_pi_rot[z], &pimi_phimin, &pimi_phimax);
    }
    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  pi2_stat[0] && pi2_stat[1])  N_all2=N_all2+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) &&  pi2_stat[0] && !pi2_stat[1])  N_1pi_1p_[0]=N_1pi_1p_[0]+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) &&  !pi2_stat[0] && pi2_stat[1])  N_1pi_1p_[1]=N_1pi_1p_[1]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  pi2_stat[0] && !pi2_stat[1])  N_1pi_1p_[2]=N_1pi_1p_[2]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  !pi2_stat[0] && pi2_stat[1])  N_1pi_1p_[3]=N_1pi_1p_[3]+1;
    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  pi2_stat[0] && !pi2_stat[1])  N_1pi_2p_[0]=N_1pi_2p_[0]+1;
    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  !pi2_stat[0] && pi2_stat[1])  N_1pi_2p_[1]=N_1pi_2p_[1]+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) &&  pi2_stat[0] && pi2_stat[1])  N_2pi_1p_[0]=N_2pi_1p_[0]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  pi2_stat[0] && pi2_stat[1])  N_2pi_1p_[1]=N_2pi_1p_[1]+1;
  }//end of N_tot loop

  if(N_all2 !=0){

//---------------------------------------------------2p2pi->1p2pi-------------------------------------------------------
      double N_2pi_1prot=0;
      double P_2p2pi_1p2pi[2][2]={0};
      double N_1pi_1prot[2] = {0};
      for(int i=0;i<2;i++){

      rot_2pi_1p (V3_pi, qpi, V3_p[i], V3_q, N_1pi_1prot, &N_2pi_1prot, N_tot);

      if(N_2pi_1prot!=0){
        P_2p2pi_1p2pi[i][0]=(N_2pi_1p_[i]/N_all2)*(N_1pi_1prot[0]/N_2pi_1prot);
        P_2p2pi_1p2pi[i][1]=(N_2pi_1p_[i]/N_all2)*(N_1pi_1prot[1]/N_2pi_1prot);
      }
    }

//---------------------------------------------------2p2pi->2p1pi-------------------------------------------------------
    double N_1pion_2prot = 0;
    double N_1pion_1prot[2] = {0};
    double P_2p2pi_2p1pi[2][2] = {0};
    for(int i=0;i<2;i++){

      rot_1pi_2p (V3_pi[i], qpi[i], V3_p, V3_q, N_1pion_1prot, &N_1pion_2prot, N_tot);

      if(N_1pion_2prot!=0){
      P_2p2pi_2p1pi[i][0] = (N_1pi_2p_[i]/N_all2)*(N_1pion_1prot[0]/N_1pion_2prot);
      P_2p2pi_2p1pi[i][1] = (N_1pi_2p_[i]/N_all2)*(N_1pion_1prot[1]/N_1pion_2prot);

    }
    }

//---------------------------------------------------2p2pi->1p1pi-------------------------------------------------------

//fill the histograms here

        // First reconstruction method
        h1_rot1_2pi_2p->Fill(en_recon1[0], (N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[0], (N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[0], (N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,-(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,-(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[0],-(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_1p_[0]/N_all2)*(1/Mott_cross_sec));

        }
        h1_rot1_2pi_2p->Fill(en_recon1[1], (N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[1], (N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[1], (N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,-(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,-(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[1],-(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[1]> 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_1p_[1]/N_all2)*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[2], (N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[2], (N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[2], (N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,-(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,-(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[2], Wvar, -(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[2],-(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[2] > 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], -(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], -(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[2], -(N_1pi_1p_[2]/N_all2)*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[3], (N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[3], (N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[3], (N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,-(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,-(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[3], Wvar, -(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[3],-(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[3] > 0 && p_perp[3] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[3], -(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[3] > 0.2 && p_perp[3] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[3], -(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
          if(p_perp[3] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[3], -(N_1pi_1p_[3]/N_all2)*(1/Mott_cross_sec));
        }

        h1_rot1_2pi_2p->Fill(en_recon1[0], -(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[0], -(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[0], -(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[0], Wvar, (P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[0],(P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
          if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (P_2p2pi_1p2pi[0][0])*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[1], -(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[1], -(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[1], -(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[1], Wvar, (P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[1],(P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
          if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (P_2p2pi_1p2pi[0][1])*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[2], -(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[2], -(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[2], -(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[2], Wvar, (P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[2],(P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          if(p_perp[2] > 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], (P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], (P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[2], (P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[3], -(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[3], -(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[3], -(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[3], Wvar, (P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[3],(P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          if(p_perp[3] > 0 && p_perp[3] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[3], (P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          if(p_perp[3] > 0.2 && p_perp[3] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[3], (P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          if(p_perp[3] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[3], (P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        }

        h1_rot1_2pi_2p->Fill(en_recon1[0], -(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[0], -(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[0], -(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[0], Wvar, (P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[0],(P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
          if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (P_2p2pi_2p1pi[0][0])*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[1], -(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[1], -(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[1], -(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[1], Wvar, (P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[1],(P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
          if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (P_2p2pi_2p1pi[0][1])*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[2], -(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[2], -(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[2], -(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[2], Wvar, (P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[2],(P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          if(p_perp[2] > 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], (P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], (P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[2], (P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        }
        h1_rot1_2pi_2p->Fill(en_recon1[3], -(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot1_2pi_2p_pipl->Fill(en_recon1[3], -(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot1_2pi_2p_pimi->Fill(en_recon1[3], -(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h1_Q2_sub->Fill(q2,(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h1_omega_sub->Fill(omega,(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h1_Wvar_sub->Fill(Wvar,(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h2_Wvar_Q2_sub->Fill(Wvar, q2,(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h2_omega_Q2_sub->Fill(omega,q2,(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h2_cal_Wvar->Fill(en_recon1[3], Wvar, (P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h1_p_perp_sub->Fill(p_perp[3],(P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          if(p_perp[3] > 0 && p_perp[3] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[3], (P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          if(p_perp[3] > 0.2 && p_perp[3] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[3], (P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          if(p_perp[3] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[3], (P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        }
        //Second reconstruction method
        h1_rot2_2pi_2p->Fill(en_recon2, ((N_1pi_1p_[0]+N_1pi_1p_[1]+N_1pi_1p_[2]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot2_2pi_2p_pipl->Fill(en_recon2, ((N_1pi_1p_[0]+N_1pi_1p_[2])/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_2p_pimi->Fill(en_recon2, ((N_1pi_1p_[0]+N_1pi_1p_[2])/N_all2)*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, -((N_1pi_1p_[0]+N_1pi_1p_[2])/N_all2)*(1/Mott_cross_sec));
        }
        if(qpi[1]>0)
        {
          h1_rot2_2pi_2p_pipl->Fill(en_recon2, ((N_1pi_1p_[1]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_2p_pimi->Fill(en_recon2, ((N_1pi_1p_[1]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, -((N_1pi_1p_[1]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
        }
        h1_rot2_2pi_2p->Fill(en_recon2, -(P_2p2pi_1p2pi[0][0] +  P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][0] +  P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot2_2pi_2p_pipl->Fill(en_recon2, -(P_2p2pi_1p2pi[0][0] + P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_2p_pimi->Fill(en_recon2, -(P_2p2pi_1p2pi[0][0] + P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, (P_2p2pi_1p2pi[0][0] + P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        }
        if(qpi[1]>0)
        {
          h1_rot2_2pi_2p_pipl->Fill(en_recon2, -(P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_2p_pimi->Fill(en_recon2, -(P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, (P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        }
        h1_rot2_2pi_2p->Fill(en_recon2, -(P_2p2pi_2p1pi[0][0] +  P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][0] +  P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot2_2pi_2p_pipl->Fill(en_recon2, -(P_2p2pi_2p1pi[0][0] + P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_2p_pimi->Fill(en_recon2, -(P_2p2pi_2p1pi[0][0] + P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, (P_2p2pi_2p1pi[0][0] + P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        }
        if(qpi[1]>0)
        {
          h1_rot2_2pi_2p_pipl->Fill(en_recon2, -(P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot2_2pi_2p_pimi->Fill(en_recon2, -(P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h2_kin_e_Wvar->Fill(en_recon2, Wvar, (P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        }

        //Third reconstruction method
        //1st pi
        h1_rot3_2pi_2p->Fill(en_recon3[0], ((N_1pi_1p_[0]+N_1pi_1p_[2])/N_all2)*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot3_2pi_2p_pipl->Fill(en_recon3[0], ((N_1pi_1p_[0]+N_1pi_1p_[2])/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_2p_pimi->Fill(en_recon3[0], ((N_1pi_1p_[0]+N_1pi_1p_[2])/N_all2)*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, -((N_1pi_1p_[0]+N_1pi_1p_[2])/N_all2)*(1/Mott_cross_sec));
        }
        h1_rot3_2pi_2p->Fill(en_recon3[0], -(P_2p2pi_1p2pi[0][0] + P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot3_2pi_2p_pipl->Fill(en_recon3[0], -(P_2p2pi_1p2pi[0][0] + P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_2p_pimi->Fill(en_recon3[0], -(P_2p2pi_1p2pi[0][0] + P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, (P_2p2pi_1p2pi[0][0] + P_2p2pi_1p2pi[1][0])*(1/Mott_cross_sec));
        }
        h1_rot3_2pi_2p->Fill(en_recon3[0], -(P_2p2pi_2p1pi[0][0] + P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        if(qpi[0]>0)
        {
          h1_rot3_2pi_2p_pipl->Fill(en_recon3[0], -(P_2p2pi_2p1pi[0][0] + P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_2p_pimi->Fill(en_recon3[0], -(P_2p2pi_2p1pi[0][0] + P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[0], Wvar, (P_2p2pi_2p1pi[0][0] + P_2p2pi_2p1pi[1][0])*(1/Mott_cross_sec));
        }
        //2nd pi
        h1_rot3_2pi_2p->Fill(en_recon3[1], ((N_1pi_1p_[1]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot3_2pi_2p_pipl->Fill(en_recon3[1], ((N_1pi_1p_[1]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_2p_pimi->Fill(en_recon3[1], ((N_1pi_1p_[1]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, -((N_1pi_1p_[1]+N_1pi_1p_[3])/N_all2)*(1/Mott_cross_sec));
        }
        h1_rot3_2pi_2p->Fill(en_recon3[1], -(P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot3_2pi_2p_pipl->Fill(en_recon3[1], -(P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_2p_pimi->Fill(en_recon3[1], -(P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, (P_2p2pi_1p2pi[0][1] + P_2p2pi_1p2pi[1][1])*(1/Mott_cross_sec));
        }
        h1_rot3_2pi_2p->Fill(en_recon3[1], -(P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        if(qpi[1]>0)
        {
          h1_rot3_2pi_2p_pipl->Fill(en_recon3[1], -(P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        }
        else
        {
          h1_rot3_2pi_2p_pimi->Fill(en_recon3[1], -(P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
          h2_kin_e_pi_Wvar->Fill(en_recon3[1], Wvar, (P_2p2pi_2p1pi[0][1] + P_2p2pi_2p1pi[1][1])*(1/Mott_cross_sec));
        }

}//N_all=0 requirement
}//end of 2pi statement
}//end of 2p statetment

TLorentzVector V4_prot[3];
double p3_kin[3];
TVector3 V3_prot[3];
//Requiring 3 protons
if(num_p == 3){
  h2_phot_pi_3p->Fill(num_pi, ec_num_n);
  V4_prot[0].SetPxPyPzE(p[index_p[0]]*cx[index_p[0]],p[index_p[0]]*cy[index_p[0]],p[index_p[0]]*cz[index_p[0]],TMath::Sqrt(m_prot*m_prot+p[index_p[0]]*p[index_p[0]]));
  V4_prot[1].SetPxPyPzE(p[index_p[1]]*cx[index_p[1]],p[index_p[1]]*cy[index_p[1]],p[index_p[1]]*cz[index_p[1]],TMath::Sqrt(m_prot*m_prot+p[index_p[1]]*p[index_p[1]]));
  V4_prot[2].SetPxPyPzE(p[index_p[2]]*cx[index_p[2]],p[index_p[2]]*cy[index_p[2]],p[index_p[2]]*cz[index_p[2]],TMath::Sqrt(m_prot*m_prot+p[index_p[2]]*p[index_p[2]]));
  V3_prot[0] = V4_prot[0].Vect();
  V3_prot[1] = V4_prot[1].Vect();
  V3_prot[2] = V4_prot[2].Vect();

  double prot_phi[3];
  double prot_phi_mod[3];
  double prot_theta[3];
  double prot_vert[3];
  double prot_vert_corr[3];
  double prot_mom_corr[3];
  TLorentzVector V4_p_corr[3];

  for(int i=0; i<3; i++)
  {
    prot_vert[i] = vz[index_p[i]];
    prot_phi[i]=TMath::ATan2(cy[index_p[i]],cx[index_p[i]])*TMath::RadToDeg();
    prot_phi_mod[i]=prot_phi[i]+30;
    if (prot_phi_mod[i]<0)prot_phi_mod[i]=prot_phi_mod[i]+360;
    prot_theta[i]=TMath::ACos(cz[index_p[i]])*TMath::RadToDeg();
    prot_vert_corr[i] = prot_vert[i] + vz_corr(prot_phi_mod[i],prot_theta[i]);
    prot_mom_corr[i] = ProtonMomCorrection_He3_4Cell(ftarget,V4_p[i],prot_vert_corr[i]);
    V4_p_corr[i].SetPxPyPzE(prot_mom_corr[i]*cx[index_p[i]],prot_mom_corr[i]*cy[index_p[i]],prot_mom_corr[i]*cz[index_p[i]],TMath::Sqrt(m_prot*m_prot+prot_mom_corr[i]*prot_mom_corr[i]));
    p3_kin[i] = V4_p_corr[i].E() - m_prot;
  }

  TLorentzVector V4_pi;
  TVector3 V3_pi;
  //Requiring 1 pion, 1 photon
  if(num_pi == 1 && ec_num_n == 1)
  {
    TLorentzVector V4_total[3];
    double p_perp[3] = {0};
    V4_pi.SetPxPyPzE(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]],TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
    V3_pi = V4_pi.Vect();
    for(int j=0; j<3; j++)
    {
      V4_total[j] = V4_pi + V4_p_corr[j] + V4_el;
      p_perp[j] = TMath::Sqrt(V4_total[j].Px()*V4_total[j].Px()+V4_total[j].Py()*V4_total[j].Py());
    }
    TVector3 V3_phot(p[ec_index_n[0]]*cx[ec_index_n[0]],p[ec_index_n[0]]*cy[ec_index_n[0]],p[ec_index_n[0]]*cz[ec_index_n[0]]);

    double en_recon1[3];
    double en_recon3;
    double qpi = q[index_pi[0]];
    for(int i = 0;i<3;i++){
      en_recon1[i] = V4_el.E() + p3_kin[i] + V4_pi.E();
    }
    en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                            (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pi[0]]));

    double N_1pi_1p_0phot[3] = {0};
    double N_1pi_1p_1phot[3] = {0};
    double N_1pi_2p_0phot[3] = {0};
    double N_1pi_2p_1phot[3] = {0};
    double N_1pi_3p_1phot = 0;
    double N_1pi_3p_0phot = 0;

    rot_1phot_1pi_3p(V3_phot, V3_pi, qpi, V3_prot, V3_q, ec_radstat_n[0], N_1pi_1p_0phot, N_1pi_1p_1phot, N_1pi_2p_0phot, N_1pi_2p_1phot, &N_1pi_3p_1phot, &N_1pi_3p_0phot, N_tot);

    if(N_1pi_3p_1phot!=0)
    {
      h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      h1_rot1_1pi_3p_1phot->Fill(en_recon1[2], (N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[2], (N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], (N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], (N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[2], (N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_Q2_sub->Fill(q2,-((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_omega_sub->Fill(omega,-((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_Wvar_sub->Fill(Wvar,-((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h2_Wvar_Q2_sub->Fill(Wvar, q2,-((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h2_omega_Q2_sub->Fill(omega,q2,-((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h2_cal_Wvar->Fill(en_recon1[2], Wvar, -(N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[0], -(N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[1], -(N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h1_p_perp_sub->Fill(p_perp[2], -(N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_1p_0phot[0]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_1p_0phot[1]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[2] > 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], -(N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], -(N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[2], -(N_1pi_1p_0phot[2]/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot2_1pi_3p_1phot->Fill(en_recon2, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_Wvar->Fill(en_recon2, Wvar, -((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      }

      h1_rot3_1pi_3p_1phot->Fill(en_recon3, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      if(q[index_pi[0]]>0)
      {
        h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      }
      else
      {
        h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, ((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -((N_1pi_1p_0phot[0]+N_1pi_1p_0phot[1]+N_1pi_1p_0phot[2])/N_1pi_3p_1phot)*(1/Mott_cross_sec));
      }
    }
    TVector3 V3_prot[2];
    double N_1pi_1p[2] = {0};
    double N_1pi_2p = 0;
    V3_prot[0] = V3_p[0];
    V3_prot[1] = V3_p[1];
    rot_1pi_2p(V3_pi, qpi, V3_prot, V3_q, N_1pi_1p, &N_1pi_2p, N_tot);
        //fill the histograms here
          if(N_1pi_2p!=0 && N_1pi_3p_1phot!=0){
              // First reconstruction method
              h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
              h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              if(qpi>0)
              {
                h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_Q2_sub->Fill(q2,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_omega_sub->Fill(omega,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_Wvar_sub->Fill(Wvar,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_Q2_sub->Fill(q2,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_omega_sub->Fill(omega,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_Wvar_sub->Fill(Wvar,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_p_perp_sub->Fill(p_perp[0], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_p_perp_sub->Fill(p_perp[1], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              }

              //Second reconstruction method
              h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
              h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              if(qpi>0)
              {
                h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              }

              //Third reconstruction method
              h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
              h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              if(qpi>0)
              {
                h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              }
              else
              {
                h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
              }
        }
        N_1pi_1p[0] = 0;
        N_1pi_1p[1] = 0;
        N_1pi_2p = 0;
        V3_prot[0] = V3_p[0];
        V3_prot[1] = V3_p[2];
        rot_1pi_2p(V3_pi, qpi, V3_prot, V3_q, N_1pi_1p, &N_1pi_2p, N_tot);
            //fill the histograms here
              if(N_1pi_2p!=0 && N_1pi_3p_1phot!=0){
                  // First reconstruction method
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[2], -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  if(qpi>0)
                  {
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[2], -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[2], -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[2], Wvar, (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[0], (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[2], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[0]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    if(p_perp[2] > 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  }

                  //Second reconstruction method
                  h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                  h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  if(qpi>0)
                  {
                    h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  }

                  //Third reconstruction method
                  h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                  h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  if(qpi>0)
                  {
                    h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(1/Mott_cross_sec));
                    h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_2p_0phot[1]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                  }
            }
            N_1pi_1p[0] = 0;
            N_1pi_1p[1] = 0;
            N_1pi_2p = 0;
            V3_prot[0] = V3_p[1];
            V3_prot[1] = V3_p[2];
            rot_1pi_2p(V3_pi, qpi, V3_prot, V3_q, N_1pi_1p, &N_1pi_2p, N_tot);
                //fill the histograms here
                  if(N_1pi_2p!=0 && N_1pi_3p_1phot!=0){
                      // First reconstruction method
                      h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                      h1_rot1_1pi_3p_1phot->Fill(en_recon1[2], -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[2], -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[2], -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_Q2_sub->Fill(q2,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_omega_sub->Fill(omega,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_Wvar_sub->Fill(Wvar,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_Q2_sub->Fill(q2,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_omega_sub->Fill(omega,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_Wvar_sub->Fill(Wvar,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_cal_Wvar->Fill(en_recon1[2], Wvar, (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_p_perp_sub->Fill(p_perp[1], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_p_perp_sub->Fill(p_perp[2], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        if(p_perp[2]> 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[2], (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      }

                      //Second reconstruction method
                      h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                      h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      }

                      //Third reconstruction method
                      h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                      h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(1/Mott_cross_sec));
                        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_2p_0phot[2]/N_1pi_3p_1phot)*(N_1pi_1p[2]/N_1pi_2p)*(1/Mott_cross_sec));
                      }
                }

                double N_1pi_2p_[3] = {0};
                double N_1pi_1p_[3] = {0};
                double N_1pi_3p_ = 0;
                rot_1pi_3p(V3_pi, qpi, V3_prot, V3_q, N_1pi_1p_, N_1pi_2p_, &N_1pi_3p_, N_tot);

                if(N_1pi_3p_!=0 && N_1pi_3p_1phot!=0){
                h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_Q2_sub->Fill(q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_omega_sub->Fill(omega,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_Wvar_sub->Fill(Wvar,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_p_perp_sub->Fill(p_perp[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_Q2_sub->Fill(q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_omega_sub->Fill(omega,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_Wvar_sub->Fill(Wvar,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_p_perp_sub->Fill(p_perp[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                h1_rot1_1pi_3p_1phot->Fill(en_recon1[2], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[2], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[2], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_Q2_sub->Fill(q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_omega_sub->Fill(omega,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_Wvar_sub->Fill(Wvar,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_cal_Wvar->Fill(en_recon1[2], Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h1_p_perp_sub->Fill(p_perp[2], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[2] > 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[2], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                }

                h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                }

                h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[0]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[1]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                if(qpi>0)
                {
                  h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
                else
                {
                  h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                  h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_[2]/N_1pi_3p_)*(1/Mott_cross_sec));
                }
              }

                double N_1pion_2prot = 0;
                double N_1pion_1prot[2] = {0};
                int count = 0;

                TVector3 V3_prot1[2];
                for(int i=0; i<3; i++)
                {
                  for(int j=0; j<3; j++)
                  {
                    if(i<j)
                    {
                      N_1pion_2prot = 0;
                      N_1pion_1prot[0] = N_1pion_1prot[1] = 0;
                      V3_prot1[0] = V3_prot[i];
                      V3_prot1[1] = V3_prot[j];
                      rot_1pi_2p (V3_pi, qpi, V3_prot1, V3_q, N_1pion_1prot, &N_1pion_2prot, N_tot);
                      en_recon1[0] = V4_el.E() + p3_kin[i] + V4_pi.E();
                      en_recon1[1] = V4_el.E() + p3_kin[j] + V4_pi.E();
                      if(N_1pi_3p_1phot!=0 && N_1pion_2prot!=0 && N_1pi_3p_!=0)
                      {

                      h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h1_Q2_sub->Fill(q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h1_omega_sub->Fill(omega,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h1_Wvar_sub->Fill(Wvar,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h1_p_perp_sub->Fill(p_perp[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }

                      h1_rot2_1pi_3p_1phot->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      h1_rot2_1pi_3p_1phot->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }

                      h1_rot3_1pi_3p_1phot->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      h1_rot3_1pi_3p_1phot->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      if(qpi>0)
                      {
                        h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      else
                      {
                        h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                        h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p_[count]/N_1pi_3p_)*(1/Mott_cross_sec));
                      }
                      count=count+1;
                    }
                  }
                }
              }
                TVector3 V3_p1[2];
                for(int i=0; i<3; i++)
                {
                  for(int j=0; j<3; j++)
                  {
                    if(i<j)
                    {
                      V3_p1[0] = V3_p[i];
                      V3_p1[1] = V3_p[j];
                double N_1pi_1p_0phot_[2] = {0};
                double N_1pi_1p_1phot_[2] = {0};
                double N_1pi_2p_1phot_ = 0;
                double N_1pi_2p_0phot_ = 0;
                en_recon1[0] = V4_el.E() + p3_kin[i] + V4_pi.E();
                en_recon1[1] = V4_el.E() + p3_kin[j] + V4_pi.E();
                rot_1phot_1pi_2p(V3_phot, V3_pi, q[index_pi[0]], V3_p1, V3_q, ec_radstat_n[0], N_1pi_1p_0phot_, N_1pi_1p_1phot_, &N_1pi_2p_0phot_, &N_1pi_2p_1phot_, N_tot);

                //fill histograms here
                if(N_1pi_2p_1phot_!=0 && N_1pi_3p_1phot!=0)
                {
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p_0phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }

                  h1_rot2_1pi_3p_1phot->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }

                  h1_rot3_1pi_3p_1phot->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                      h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*((N_1pi_1p_0phot_[0]+N_1pi_1p_0phot_[1])/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                }
                double N1pi1p0phot[2] = {0};
                double N1pi1p1phot[2] = {0};
                rot_1phot_1pi_1p(V3_phot, V3_pi, q[index_pi[0]], V3_p[0], V3_q, ec_radstat_n[0], &N1pi1p0phot[0], &N1pi1p1phot[0], N_tot);
                rot_1phot_1pi_1p(V3_phot, V3_pi, q[index_pi[0]], V3_p[1], V3_q, ec_radstat_n[0], &N1pi1p0phot[1], &N1pi1p1phot[1], N_tot);
                if(N_1pi_2p_1phot_!=0 && N_1pi_3p_1phot!=0 && N1pi1p1phot[0]!=0 && N1pi1p1phot[1]!=0 && N_1pi_2p_1phot_!=0)
                {
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }

                  h1_rot2_1pi_3p_1phot->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  h1_rot2_1pi_3p_1phot->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }

                  h1_rot3_1pi_3p_1phot->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  h1_rot3_1pi_3p_1phot->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[0]/N1pi1p1phot[0])*(N_1pi_1p_1phot_[0]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N1pi1p0phot[1]/N1pi1p1phot[1])*(N_1pi_1p_1phot_[1]/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                }
                N_1pi_1p[0] = 0;
                N_1pi_1p[1] = 0;
                N_1pi_2p = 0;
                rot_1pi_2p(V3_pi, q[index_pi[0]], V3_p, V3_q, N_1pi_1p, &N_1pi_2p, N_tot);

                if(N_1pi_2p_1phot_!=0 && N_1pi_2p!=0 && N_1pi_3p_1phot!=0)
                {
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  h1_rot1_1pi_3p_1phot->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pipl->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[0], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot1_1pi_3p_1phot_pimi->Fill(en_recon1[1], (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Q2_sub->Fill(q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_omega_sub->Fill(omega,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_Wvar_sub->Fill(Wvar,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_p_perp_sub->Fill(p_perp[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }

                  h1_rot2_1pi_3p_1phot->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  h1_rot2_1pi_3p_1phot->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot2_1pi_3p_1phot_pipl->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot2_1pi_3p_1phot_pimi->Fill(en_recon2, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }

                  h1_rot3_1pi_3p_1phot->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  h1_rot3_1pi_3p_1phot->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  if(q[index_pi[0]]>0)
                  {
                    h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot3_1pi_3p_1phot_pipl->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                  else
                  {
                    h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h1_rot3_1pi_3p_1phot_pimi->Fill(en_recon3, (N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[0]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                    h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_3p_0phot/N_1pi_3p_1phot)*(N_1pi_1p[1]/N_1pi_2p)*(N_1pi_2p_0phot_/N_1pi_2p_1phot_)*(1/Mott_cross_sec));
                  }
                }
              }
            }
          }
}
//End of 1 pi 0 phot statement
  //Requiring 1 pion, 0 photon
  if(num_pi == 1 && ec_num_n==0 && num_n==0)
  {
    double p_perp[3] = {0};
    TLorentzVector V4_total[3];
    V4_pi.SetPxPyPzE(p[index_pi[0]]*cx[index_pi[0]],p[index_pi[0]]*cy[index_pi[0]],p[index_pi[0]]*cz[index_pi[0]],TMath::Sqrt(p[index_pi[0]]*p[index_pi[0]]+Mpi*Mpi));
    V3_pi = V4_pi.Vect();
    for(int j=0;j<3;j++)
    {
      V4_total[j] = V4_pi + V4_p_corr[j] + V4_el;
      p_perp[j] = TMath::Sqrt(V4_total[j].Px()*V4_total[j].Px()+V4_total[j].Py()*V4_total[j].Py());
    }

    double en_recon1[3];
    double en_recon3;
    double qpi = 0;
    qpi = q[index_pi[0]];
    for(int i = 0;i<3;i++){
      en_recon1[i] = V4_el.E() + p3_kin[i] + V4_pi.E();
    }
    en_recon3 = (m_prot*m_prot - (m_prot - eps)*(m_prot - eps) + 2*(m_neut - eps)*(V4_el.E() + V4_pi.E()) - e_mass*e_mass - 2*V4_el.E()*V4_pi.E() + 2*V4_el.Rho()*V4_pi.Rho()*cos(V3_pi.Angle(V3_el)) + Mpi*Mpi)/
                                            (2*(m_neut - eps - V4_el.E() - V4_pi.E()) + 2*(V4_el.Rho()*cz[ind_em] + V4_pi.Rho()*cz[index_pi[0]]));

    double N_1pi_1p[3] = {0};
    double N_1pi_2p[3] = {0};
    double N_1pi_3p = 0;
    rot_1pi_3p(V3_pi, qpi, V3_prot, V3_q, N_1pi_1p, N_1pi_2p, &N_1pi_3p, N_tot);

    if(N_1pi_3p!=0){
    h1_rot1_1pi_3p->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot1_1pi_3p_pipl->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot1_1pi_3p_pimi->Fill(en_recon1[0], (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_Q2_sub->Fill(q2,-(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_omega_sub->Fill(omega,-(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_cal_Wvar->Fill(en_recon1[0], Wvar, -(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_p_perp_sub->Fill(p_perp[0],-(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], -(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    h1_rot1_1pi_3p->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot1_1pi_3p_pipl->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot1_1pi_3p_pimi->Fill(en_recon1[1], (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_Q2_sub->Fill(q2,-(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_omega_sub->Fill(omega,-(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_cal_Wvar->Fill(en_recon1[1], Wvar, -(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_p_perp_sub->Fill(p_perp[1],-(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], -(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    h1_rot1_1pi_3p->Fill(en_recon1[2], (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot1_1pi_3p_pipl->Fill(en_recon1[2], (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot1_1pi_3p_pimi->Fill(en_recon1[2], (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_Q2_sub->Fill(q2,-(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_omega_sub->Fill(omega,-(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_Wvar_sub->Fill(Wvar,-(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_Wvar_Q2_sub->Fill(Wvar, q2,-(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_omega_Q2_sub->Fill(omega,q2,-(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_cal_Wvar->Fill(en_recon1[2], Wvar, -(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h1_p_perp_sub->Fill(p_perp[2],-(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[2] > 0 && p_perp[2] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[2], -(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[2] > 0.2 && p_perp[2] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[2], -(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      if(p_perp[2] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[2], -(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    }

    h1_rot2_1pi_3p->Fill(en_recon2, (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot2_1pi_3p_pipl->Fill(en_recon2, (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot2_1pi_3p_pimi->Fill(en_recon2, (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    h1_rot2_1pi_3p->Fill(en_recon2, (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot2_1pi_3p_pipl->Fill(en_recon2, (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot2_1pi_3p_pimi->Fill(en_recon2, (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    h1_rot2_1pi_3p->Fill(en_recon2, (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot2_1pi_3p_pipl->Fill(en_recon2, (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot2_1pi_3p_pimi->Fill(en_recon2, (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_kin_e_Wvar->Fill(en_recon2, Wvar, -(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    h1_rot3_1pi_3p->Fill(en_recon3, (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot3_1pi_3p_pipl->Fill(en_recon3, (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot3_1pi_3p_pimi->Fill(en_recon3, (N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_1p[0]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    h1_rot3_1pi_3p->Fill(en_recon3, (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot3_1pi_3p_pipl->Fill(en_recon3, (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot3_1pi_3p_pimi->Fill(en_recon3, (N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_1p[1]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    h1_rot3_1pi_3p->Fill(en_recon3, (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    if(qpi>0)
    {
      h1_rot3_1pi_3p_pipl->Fill(en_recon3, (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    }
    else
    {
      h1_rot3_1pi_3p_pimi->Fill(en_recon3, (N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
      h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, -(N_1pi_1p[2]/N_1pi_3p)*(1/Mott_cross_sec));
    }
  }

    double N_1pion_2prot = 0;
    double N_1pion_1prot[2] = {0};
    int count = 0;

    TVector3 V3_prot1[2];
    for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
      {
        if(i<j)
        {
          N_1pion_2prot = 0;
          N_1pion_1prot[0] = N_1pion_1prot[1] = 0;
          V3_prot1[0] = V3_prot[i];
          V3_prot1[1] = V3_prot[j];
          rot_1pi_2p (V3_pi, qpi, V3_prot1, V3_q, N_1pion_1prot, &N_1pion_2prot, N_tot);
          en_recon1[0] = V4_el.E() + p3_kin[i] + V4_pi.E();
          en_recon1[1] = V4_el.E() + p3_kin[j] + V4_pi.E();
          V4_total[0] = V4_pi + V4_p_corr[i] + V4_el;
          p_perp[0] = TMath::Sqrt(V4_total[0].Px()*V4_total[0].Px()+V4_total[0].Py()*V4_total[0].Py());
          V4_total[1] = V4_pi + V4_p_corr[j] + V4_el;
          p_perp[1] = TMath::Sqrt(V4_total[1].Px()*V4_total[1].Px()+V4_total[1].Py()*V4_total[1].Py());
          if(N_1pion_2prot!=0 && N_1pi_3p!=0)
          {
          h1_rot1_1pi_3p->Fill(en_recon1[0], -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot1_1pi_3p_pipl->Fill(en_recon1[0], -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot1_1pi_3p_pimi->Fill(en_recon1[0], -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_Q2_sub->Fill(q2,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_omega_sub->Fill(omega,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_Wvar_sub->Fill(Wvar,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_omega_Q2_sub->Fill(omega,q2,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_p_perp_sub->Fill(p_perp[0],(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          h1_rot1_1pi_3p->Fill(en_recon1[1], -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot1_1pi_3p_pipl->Fill(en_recon1[1], -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot1_1pi_3p_pimi->Fill(en_recon1[1], -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_Q2_sub->Fill(q2,(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_omega_sub->Fill(omega,(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_Wvar_sub->Fill(Wvar,(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_omega_Q2_sub->Fill(omega,q2,(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_cal_Wvar->Fill(en_recon1[1], Wvar, (N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h1_p_perp_sub->Fill(p_perp[1],(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            if(p_perp[1] > 0 && p_perp[1] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[1], (N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            if(p_perp[1] > 0.2 && p_perp[1] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[1], (N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            if(p_perp[1] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[1], (N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          h1_rot2_1pi_3p->Fill(en_recon2, -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot2_1pi_3p_pipl->Fill(en_recon2, -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot2_1pi_3p_pimi->Fill(en_recon2, -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          h1_rot2_1pi_3p->Fill(en_recon2, -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot2_1pi_3p_pipl->Fill(en_recon2, -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot2_1pi_3p_pimi->Fill(en_recon2, -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_kin_e_Wvar->Fill(en_recon2, Wvar, (N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          h1_rot3_1pi_3p->Fill(en_recon3, -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot3_1pi_3p_pipl->Fill(en_recon3, -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot3_1pi_3p_pimi->Fill(en_recon3, -(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          h1_rot3_1pi_3p->Fill(en_recon3, -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          if(qpi>0)
          {
            h1_rot3_1pi_3p_pipl->Fill(en_recon3, -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          else
          {
            h1_rot3_1pi_3p_pimi->Fill(en_recon3, -(N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
            h2_kin_e_pi_Wvar->Fill(en_recon3, Wvar, (N_1pion_1prot[1]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
          }
          count=count+1;
        }
      }
      }
    }
  }//end of 1pi statement
}//end of 3p statement

}//electron vertex cut
}
//Here I clone the original histograms and then add (subtract) the subtraction histograms
TH1F* h1_sub_cal_1pi_2p = (TH1F*) h1_en_recon1->Clone("sub_cal_1pi_2p");
h1_sub_cal_1pi_2p->Add(h1_rot1_1pi_2p, -1);

TH1F* h1_sub_cal_1pi_2p_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_1pi_2p_pimi");
h1_sub_cal_1pi_2p_pimi->Add(h1_rot1_1pi_2p_pimi, -1);

TH1F* h1_sub_cal_1pi_2p_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_1pi_2p_pipl");
h1_sub_cal_1pi_2p_pipl->Add(h1_rot1_1pi_2p_pipl, -1);

TH1F* h1_sub_kin_e_1pi_2p = (TH1F*) h1_en_recon2->Clone("sub_kin_e_1pi_2p");
h1_sub_kin_e_1pi_2p->Add(h1_rot2_1pi_2p, -1);

TH1F* h1_sub_kin_e_1pi_2p_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_1pi_2p_pimi");
h1_sub_kin_e_1pi_2p_pimi->Add(h1_rot2_1pi_2p_pimi, -1);

TH1F* h1_sub_kin_e_1pi_2p_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_1pi_2p_pipl");
h1_sub_kin_e_1pi_2p_pipl->Add(h1_rot2_1pi_2p_pipl, -1);

TH1F* h1_sub_kin_e_pi_1pi_2p = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_1pi_2p");
h1_sub_kin_e_pi_1pi_2p->Add(h1_rot3_1pi_2p, -1);

TH1F* h1_sub_kin_e_pi_1pi_2p_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_1pi_2p_pimi");
h1_sub_kin_e_pi_1pi_2p_pimi->Add(h1_rot3_1pi_2p_pimi, -1);

TH1F* h1_sub_kin_e_pi_1pi_2p_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_1pi_2p_pipl");
h1_sub_kin_e_pi_1pi_2p_pipl->Add(h1_rot3_1pi_2p_pipl, -1);

TH1F* h1_sub_cal_2pi_1p = (TH1F*) h1_en_recon1->Clone("sub_cal_2pi_1p");
h1_sub_cal_2pi_1p->Add(h1_rot1_2pi_1p, -1);

TH1F* h1_sub_cal_2pi_1p_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_2pi_1p_pimi");
h1_sub_cal_2pi_1p_pimi->Add(h1_rot1_2pi_1p_pimi, -1);

TH1F* h1_sub_cal_2pi_1p_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_2pi_1p_pipl");
h1_sub_cal_2pi_1p_pipl->Add(h1_rot1_2pi_1p_pipl, -1);

TH1F* h1_sub_kin_e_2pi_1p = (TH1F*) h1_en_recon2->Clone("sub_kin_e_2pi_1p");
h1_sub_kin_e_2pi_1p->Add(h1_rot2_2pi_1p, -1);

TH1F* h1_sub_kin_e_2pi_1p_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_2pi_1p_pimi");
h1_sub_kin_e_2pi_1p_pimi->Add(h1_rot2_2pi_1p_pimi, -1);

TH1F* h1_sub_kin_e_2pi_1p_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_2pi_1p_pipl");
h1_sub_kin_e_2pi_1p_pipl->Add(h1_rot2_2pi_1p_pipl, -1);

TH1F* h1_sub_kin_e_pi_2pi_1p = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_2pi_1p");
h1_sub_kin_e_pi_2pi_1p->Add(h1_rot3_2pi_1p, -1);

TH1F* h1_sub_kin_e_pi_2pi_1p_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_2pi_1p_pimi");
h1_sub_kin_e_pi_2pi_1p_pimi->Add(h1_rot3_2pi_1p_pimi, -1);

TH1F* h1_sub_kin_e_pi_2pi_1p_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_2pi_1p_pipl");
h1_sub_kin_e_pi_2pi_1p_pipl->Add(h1_rot3_2pi_1p_pipl, -1);

TH1F* h1_sub_cal_2pi_2p = (TH1F*) h1_en_recon1->Clone("sub_cal_2pi_2p");
h1_sub_cal_2pi_2p->Add(h1_rot1_2pi_2p, -1);

TH1F* h1_sub_cal_2pi_2p_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_2pi_2p_pimi");
h1_sub_cal_2pi_2p_pimi->Add(h1_rot1_2pi_2p_pimi, -1);

TH1F* h1_sub_cal_2pi_2p_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_2pi_2p_pipl");
h1_sub_cal_2pi_2p_pipl->Add(h1_rot1_2pi_2p_pipl, -1);

TH1F* h1_sub_kin_e_2pi_2p = (TH1F*) h1_en_recon2->Clone("sub_kin_e_2pi_2p");
h1_sub_kin_e_2pi_2p->Add(h1_rot2_2pi_2p, -1);

TH1F* h1_sub_kin_e_2pi_2p_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_2pi_2p_pimi");
h1_sub_kin_e_2pi_2p_pimi->Add(h1_rot2_2pi_2p_pimi, -1);

TH1F* h1_sub_kin_e_2pi_2p_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_2pi_2p_pipl");
h1_sub_kin_e_2pi_2p_pipl->Add(h1_rot2_2pi_2p_pipl, -1);

TH1F* h1_sub_kin_e_pi_2pi_2p = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_2pi_2p");
h1_sub_kin_e_pi_2pi_2p->Add(h1_rot3_2pi_2p, -1);

TH1F* h1_sub_kin_e_pi_2pi_2p_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_2pi_2p_pimi");
h1_sub_kin_e_pi_2pi_2p_pimi->Add(h1_rot3_2pi_2p_pimi, -1);

TH1F* h1_sub_kin_e_pi_2pi_2p_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_2pi_2p_pipl");
h1_sub_kin_e_pi_2pi_2p_pipl->Add(h1_rot3_2pi_2p_pipl, -1);

TH1F* h1_sub_cal_1pi_3p = (TH1F*) h1_en_recon1->Clone("sub_cal_1pi_3p");
h1_sub_cal_1pi_3p->Add(h1_rot1_1pi_3p, -1);

TH1F* h1_sub_cal_1pi_3p_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_1pi_3p_pimi");
h1_sub_cal_1pi_3p_pimi->Add(h1_rot1_1pi_3p_pimi, -1);

TH1F* h1_sub_cal_1pi_3p_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_1pi_3p_pipl");
h1_sub_cal_1pi_3p_pipl->Add(h1_rot1_1pi_3p_pipl, -1);

TH1F* h1_sub_kin_e_1pi_3p = (TH1F*) h1_en_recon2->Clone("sub_kin_e_1pi_3p");
h1_sub_kin_e_1pi_3p->Add(h1_rot2_1pi_3p, -1);

TH1F* h1_sub_kin_e_1pi_3p_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_1pi_3p_pimi");
h1_sub_kin_e_1pi_3p_pimi->Add(h1_rot2_1pi_3p_pimi, -1);

TH1F* h1_sub_kin_e_1pi_3p_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_1pi_3p_pipl");
h1_sub_kin_e_1pi_3p_pipl->Add(h1_rot2_1pi_3p_pipl, -1);

TH1F* h1_sub_kin_e_pi_1pi_3p = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_1pi_3p");
h1_sub_kin_e_pi_1pi_3p->Add(h1_rot3_1pi_3p, -1);

TH1F* h1_sub_kin_e_pi_1pi_3p_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_1pi_3p_pimi");
h1_sub_kin_e_pi_1pi_3p_pimi->Add(h1_rot3_1pi_3p_pimi, -1);

TH1F* h1_sub_kin_e_pi_1pi_3p_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_1pi_3p_pipl");
h1_sub_kin_e_pi_1pi_3p_pipl->Add(h1_rot3_1pi_3p_pipl, -1);

TH1F* h1_sub_cal_3pi_1p = (TH1F*) h1_en_recon1->Clone("sub_cal_3pi_1p");
h1_sub_cal_3pi_1p->Add(h1_rot1_3pi_1p, -1);

TH1F* h1_sub_cal_3pi_1p_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_3pi_1p_pimi");
h1_sub_cal_3pi_1p_pimi->Add(h1_rot1_3pi_1p_pimi, -1);

TH1F* h1_sub_cal_3pi_1p_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_3pi_1p_pipl");
h1_sub_cal_3pi_1p_pipl->Add(h1_rot1_3pi_1p_pipl, -1);

TH1F* h1_sub_kin_e_3pi_1p = (TH1F*) h1_en_recon2->Clone("sub_kin_e_3pi_1p");
h1_sub_kin_e_3pi_1p->Add(h1_rot2_3pi_1p, -1);

TH1F* h1_sub_kin_e_3pi_1p_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_3pi_1p_pimi");
h1_sub_kin_e_3pi_1p_pimi->Add(h1_rot2_3pi_1p_pimi, -1);

TH1F* h1_sub_kin_e_3pi_1p_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_3pi_1p_pipl");
h1_sub_kin_e_3pi_1p_pipl->Add(h1_rot2_3pi_1p_pipl, -1);

TH1F* h1_sub_kin_e_pi_3pi_1p = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_3pi_1p");
h1_sub_kin_e_pi_3pi_1p->Add(h1_rot3_3pi_1p, -1);

TH1F* h1_sub_kin_e_pi_3pi_1p_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_3pi_1p_pimi");
h1_sub_kin_e_pi_3pi_1p_pimi->Add(h1_rot3_3pi_1p_pimi, -1);

TH1F* h1_sub_kin_e_pi_3pi_1p_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_3pi_1p_pipl");
h1_sub_kin_e_pi_3pi_1p_pipl->Add(h1_rot3_3pi_1p_pipl, -1);

TH1F* h1_sub_cal_1pi_1p_1phot = (TH1F*) h1_en_recon1->Clone("sub_cal_1pi_1p_1phot");
h1_sub_cal_1pi_1p_1phot->Add(h1_rot1_1pi_1p_1phot, -1);

TH1F* h1_sub_cal_1pi_1p_1phot_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_1pi_1p_1phot_pimi");
h1_sub_cal_1pi_1p_1phot_pimi->Add(h1_rot1_1pi_1p_1phot_pimi, -1);

TH1F* h1_sub_cal_1pi_1p_1phot_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_1pi_1p_1phot_pipl");
h1_sub_cal_1pi_1p_1phot_pipl->Add(h1_rot1_1pi_1p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_1pi_1p_1phot = (TH1F*) h1_en_recon2->Clone("sub_kin_e_1pi_1p_1phot");
h1_sub_kin_e_1pi_1p_1phot->Add(h1_rot2_1pi_1p_1phot, -1);

TH1F* h1_sub_kin_e_1pi_1p_1phot_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_1pi_1p_1phot_pimi");
h1_sub_kin_e_1pi_1p_1phot_pimi->Add(h1_rot2_1pi_1p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_1pi_1p_1phot_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_1pi_1p_1phot_pipl");
h1_sub_kin_e_1pi_1p_1phot_pipl->Add(h1_rot2_1pi_1p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_pi_1pi_1p_1phot = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_1pi_1p_1phot");
h1_sub_kin_e_pi_1pi_1p_1phot->Add(h1_rot3_1pi_1p_1phot, -1);

TH1F* h1_sub_kin_e_pi_1pi_1p_1phot_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_1pi_1p_1phot_pimi");
h1_sub_kin_e_pi_1pi_1p_1phot_pimi->Add(h1_rot3_1pi_1p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_pi_1pi_1p_1phot_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_1pi_1p_1phot_pipl");
h1_sub_kin_e_pi_1pi_1p_1phot_pipl->Add(h1_rot3_1pi_1p_1phot_pipl, -1);

TH1F* h1_sub_cal_1pi_2p_1phot = (TH1F*) h1_en_recon1->Clone("sub_cal_1pi_2p_1phot");
h1_sub_cal_1pi_2p_1phot->Add(h1_rot1_1pi_2p_1phot, -1);

TH1F* h1_sub_cal_1pi_2p_1phot_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_1pi_2p_1phot_pimi");
h1_sub_cal_1pi_2p_1phot_pimi->Add(h1_rot1_1pi_2p_1phot_pimi, -1);

TH1F* h1_sub_cal_1pi_2p_1phot_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_1pi_2p_1phot_pipl");
h1_sub_cal_1pi_2p_1phot_pipl->Add(h1_rot1_1pi_2p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_1pi_2p_1phot = (TH1F*) h1_en_recon2->Clone("sub_kin_e_1pi_2p_1phot");
h1_sub_kin_e_1pi_2p_1phot->Add(h1_rot2_1pi_2p_1phot, -1);

TH1F* h1_sub_kin_e_1pi_2p_1phot_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_1pi_2p_1phot_pimi");
h1_sub_kin_e_1pi_2p_1phot_pimi->Add(h1_rot2_1pi_2p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_1pi_2p_1phot_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_1pi_2p_1phot_pipl");
h1_sub_kin_e_1pi_2p_1phot_pipl->Add(h1_rot2_1pi_2p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_pi_1pi_2p_1phot = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_1pi_2p_1phot");
h1_sub_kin_e_pi_1pi_2p_1phot->Add(h1_rot3_1pi_2p_1phot, -1);

TH1F* h1_sub_kin_e_pi_1pi_2p_1phot_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_1pi_2p_1phot_pimi");
h1_sub_kin_e_pi_1pi_2p_1phot_pimi->Add(h1_rot3_1pi_2p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_pi_1pi_2p_1phot_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_1pi_2p_1phot_pipl");
h1_sub_kin_e_pi_1pi_2p_1phot_pipl->Add(h1_rot3_1pi_2p_1phot_pipl, -1);

TH1F* h1_sub_cal_1pi_1p_2phot = (TH1F*) h1_en_recon1->Clone("sub_cal_1pi_1p_2phot");
h1_sub_cal_1pi_1p_2phot->Add(h1_rot1_1pi_1p_2phot, -1);

TH1F* h1_sub_cal_1pi_1p_2phot_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_1pi_1p_2phot_pimi");
h1_sub_cal_1pi_1p_2phot_pimi->Add(h1_rot1_1pi_1p_2phot_pimi, -1);

TH1F* h1_sub_cal_1pi_1p_2phot_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_1pi_1p_2phot_pipl");
h1_sub_cal_1pi_1p_2phot_pipl->Add(h1_rot1_1pi_1p_2phot_pipl, -1);

TH1F* h1_sub_kin_e_1pi_1p_2phot = (TH1F*) h1_en_recon2->Clone("sub_kin_e_1pi_1p_2phot");
h1_sub_kin_e_1pi_1p_2phot->Add(h1_rot2_1pi_1p_2phot, -1);

TH1F* h1_sub_kin_e_1pi_1p_2phot_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_1pi_1p_2phot_pimi");
h1_sub_kin_e_1pi_1p_2phot_pimi->Add(h1_rot2_1pi_1p_2phot_pimi, -1);

TH1F* h1_sub_kin_e_1pi_1p_2phot_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_1pi_1p_2phot_pipl");
h1_sub_kin_e_1pi_1p_2phot_pipl->Add(h1_rot2_1pi_1p_2phot_pipl, -1);

TH1F* h1_sub_kin_e_pi_1pi_1p_2phot = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_1pi_1p_2phot");
h1_sub_kin_e_pi_1pi_1p_2phot->Add(h1_rot3_1pi_1p_2phot, -1);

TH1F* h1_sub_kin_e_pi_1pi_1p_2phot_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_1pi_1p_2phot_pimi");
h1_sub_kin_e_pi_1pi_1p_2phot_pimi->Add(h1_rot3_1pi_1p_2phot_pimi, -1);

TH1F* h1_sub_kin_e_pi_1pi_1p_2phot_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_1pi_1p_2phot_pipl");
h1_sub_kin_e_pi_1pi_1p_2phot_pipl->Add(h1_rot3_1pi_1p_2phot_pipl, -1);

TH1F* h1_sub_cal_2pi_1p_1phot = (TH1F*) h1_en_recon1->Clone("sub_cal_2pi_1p_1phot");
h1_sub_cal_2pi_1p_1phot->Add(h1_rot1_2pi_1p_1phot, -1);

TH1F* h1_sub_cal_2pi_1p_1phot_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_2pi_1p_1phot_pimi");
h1_sub_cal_2pi_1p_1phot_pimi->Add(h1_rot1_2pi_1p_1phot_pimi, -1);

TH1F* h1_sub_cal_2pi_1p_1phot_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_2pi_1p_1phot_pipl");
h1_sub_cal_2pi_1p_1phot_pipl->Add(h1_rot1_2pi_1p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_2pi_1p_1phot = (TH1F*) h1_en_recon2->Clone("sub_kin_e_2pi_1p_1phot");
h1_sub_kin_e_2pi_1p_1phot->Add(h1_rot2_2pi_1p_1phot, -1);

TH1F* h1_sub_kin_e_2pi_1p_1phot_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_2pi_1p_1phot_pimi");
h1_sub_kin_e_2pi_1p_1phot_pimi->Add(h1_rot2_2pi_1p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_2pi_1p_1phot_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_2pi_1p_1phot_pipl");
h1_sub_kin_e_2pi_1p_1phot_pipl->Add(h1_rot2_2pi_1p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_pi_2pi_1p_1phot = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_2pi_1p_1phot");
h1_sub_kin_e_pi_2pi_1p_1phot->Add(h1_rot3_2pi_1p_1phot, -1);

TH1F* h1_sub_kin_e_pi_2pi_1p_1phot_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_2pi_1p_1phot_pimi");
h1_sub_kin_e_pi_2pi_1p_1phot_pimi->Add(h1_rot3_2pi_1p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_pi_2pi_1p_1phot_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_2pi_1p_1phot_pipl");
h1_sub_kin_e_pi_2pi_1p_1phot_pipl->Add(h1_rot3_2pi_1p_1phot_pipl, -1);

TH1F* h1_sub_cal_1pi_3p_1phot = (TH1F*) h1_en_recon1->Clone("sub_cal_1pi_3p_1phot");
h1_sub_cal_1pi_3p_1phot->Add(h1_rot1_1pi_3p_1phot, -1);

TH1F* h1_sub_cal_1pi_3p_1phot_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_1pi_3p_1phot_pimi");
h1_sub_cal_1pi_3p_1phot_pimi->Add(h1_rot1_1pi_3p_1phot_pimi, -1);

TH1F* h1_sub_cal_1pi_3p_1phot_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_1pi_3p_1phot_pipl");
h1_sub_cal_1pi_3p_1phot_pipl->Add(h1_rot1_1pi_3p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_1pi_3p_1phot = (TH1F*) h1_en_recon2->Clone("sub_kin_e_1pi_3p_1phot");
h1_sub_kin_e_1pi_3p_1phot->Add(h1_rot2_1pi_3p_1phot, -1);

TH1F* h1_sub_kin_e_1pi_3p_1phot_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_1pi_3p_1phot_pimi");
h1_sub_kin_e_1pi_3p_1phot_pimi->Add(h1_rot2_1pi_3p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_1pi_3p_1phot_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_1pi_3p_1phot_pipl");
h1_sub_kin_e_1pi_3p_1phot_pipl->Add(h1_rot2_1pi_3p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_pi_1pi_3p_1phot = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_1pi_3p_1phot");
h1_sub_kin_e_pi_1pi_3p_1phot->Add(h1_rot3_1pi_3p_1phot, -1);

TH1F* h1_sub_kin_e_pi_1pi_3p_1phot_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_1pi_3p_1phot_pimi");
h1_sub_kin_e_pi_1pi_3p_1phot_pimi->Add(h1_rot3_1pi_3p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_pi_1pi_3p_1phot_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_1pi_3p_1phot_pipl");
h1_sub_kin_e_pi_1pi_3p_1phot_pipl->Add(h1_rot3_1pi_3p_1phot_pipl, -1);

//Fully subtracted histograms
TH1F* h1_sub_cal_all = (TH1F*) h1_en_recon1->Clone("sub_cal_all");
h1_sub_cal_all->Add(h1_rot1_1pi_2p, -1);
h1_sub_cal_all->Add(h1_rot1_2pi_1p, -1);
h1_sub_cal_all->Add(h1_rot1_2pi_2p, -1);
h1_sub_cal_all->Add(h1_rot1_1pi_3p, -1);
h1_sub_cal_all->Add(h1_rot1_3pi_1p, -1);
h1_sub_cal_all->Add(h1_rot1_1pi_1p_1phot, -1);
h1_sub_cal_all->Add(h1_rot1_1pi_2p_1phot, -1);
h1_sub_cal_all->Add(h1_rot1_1pi_1p_2phot, -1);
h1_sub_cal_all->Add(h1_rot1_2pi_1p_1phot, -1);
h1_sub_cal_all->Add(h1_rot1_1pi_3p_1phot, -1);


TH1F* h1_sub_cal_all_pimi = (TH1F*) h1_en_recon1_pimi->Clone("sub_cal_all_pimi");
h1_sub_cal_all_pimi->Add(h1_rot1_1pi_2p_pimi, -1);
h1_sub_cal_all_pimi->Write("process1");
h1_sub_cal_all_pimi->Add(h1_rot1_2pi_1p_pimi, -1);
h1_sub_cal_all_pimi->Write("process2");
h1_sub_cal_all_pimi->Add(h1_rot1_2pi_2p_pimi, -1);
h1_sub_cal_all_pimi->Write("process3");
h1_sub_cal_all_pimi->Add(h1_rot1_1pi_3p_pimi, -1);
h1_sub_cal_all_pimi->Write("process4");
h1_sub_cal_all_pimi->Add(h1_rot1_3pi_1p_pimi, -1);
h1_sub_cal_all_pimi->Write("process5");
h1_sub_cal_all_pimi->Add(h1_rot1_1pi_1p_1phot_pimi, -1);
h1_sub_cal_all_pimi->Add(h1_rot1_1pi_2p_1phot_pimi, -1);
h1_sub_cal_all_pimi->Add(h1_rot1_1pi_1p_2phot_pimi, -1);
h1_sub_cal_all_pimi->Add(h1_rot1_2pi_1p_1phot_pimi, -1);
h1_sub_cal_all_pimi->Add(h1_rot1_1pi_3p_1phot_pimi, -1);

TH1F* h1_sub_cal_all_pipl = (TH1F*) h1_en_recon1_pipl->Clone("sub_cal_all_pipl");
h1_sub_cal_all_pipl->Add(h1_rot1_1pi_2p_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_2pi_1p_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_2pi_2p_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_1pi_3p_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_3pi_1p_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_1pi_1p_1phot_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_1pi_2p_1phot_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_1pi_1p_2phot_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_2pi_1p_1phot_pipl, -1);
h1_sub_cal_all_pipl->Add(h1_rot1_1pi_3p_1phot_pipl, -1);

TH1F* h1_sub_kin_e_all = (TH1F*) h1_en_recon2->Clone("sub_kin_e_all");
h1_sub_kin_e_all->Add(h1_rot2_1pi_2p, -1);
h1_sub_kin_e_all->Add(h1_rot2_2pi_1p, -1);
h1_sub_kin_e_all->Add(h1_rot2_2pi_2p, -1);
h1_sub_kin_e_all->Add(h1_rot2_1pi_3p, -1);
h1_sub_kin_e_all->Add(h1_rot2_3pi_1p, -1);
h1_sub_kin_e_all->Add(h1_rot2_1pi_1p_1phot, -1);
h1_sub_kin_e_all->Add(h1_rot2_1pi_2p_1phot, -1);
h1_sub_kin_e_all->Add(h1_rot2_1pi_1p_2phot, -1);
h1_sub_kin_e_all->Add(h1_rot2_2pi_1p_1phot, -1);
h1_sub_kin_e_all->Add(h1_rot2_1pi_3p_1phot, -1);

TH1F* h1_sub_kin_e_all_pimi = (TH1F*) h1_en_recon2_pimi->Clone("sub_kin_e_all_pimi");
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_2p_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_2pi_1p_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_2pi_2p_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_3p_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_3pi_1p_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_1p_1phot_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_2p_1phot_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_1p_2phot_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_1p_1phot_pimi, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_3p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_all_pipl = (TH1F*) h1_en_recon2_pipl->Clone("sub_kin_e_all_pipl");
h1_sub_kin_e_all_pipl->Add(h1_rot2_1pi_2p_pipl, -1);
h1_sub_kin_e_all_pipl->Add(h1_rot2_2pi_1p_pipl, -1);
h1_sub_kin_e_all_pipl->Add(h1_rot2_2pi_2p_pipl, -1);
h1_sub_kin_e_all_pipl->Add(h1_rot2_1pi_3p_pipl, -1);
h1_sub_kin_e_all_pipl->Add(h1_rot2_3pi_1p_pipl, -1);
h1_sub_kin_e_all_pipl->Add(h1_rot2_1pi_1p_1phot_pipl, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_2p_1phot_pimi, -1);
h1_sub_kin_e_all_pipl->Add(h1_rot2_1pi_1p_2phot_pipl, -1);
h1_sub_kin_e_all_pipl->Add(h1_rot2_2pi_1p_1phot_pipl, -1);
h1_sub_kin_e_all_pimi->Add(h1_rot2_1pi_3p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_pi_all = (TH1F*) h1_en_recon3->Clone("sub_kin_e_pi_all");
h1_sub_kin_e_pi_all->Add(h1_rot3_1pi_2p, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_2pi_1p, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_2pi_2p, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_1pi_3p, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_3pi_1p, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_1pi_1p_1phot, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_1pi_2p_1phot, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_1pi_1p_2phot, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_2pi_1p_1phot, -1);
h1_sub_kin_e_pi_all->Add(h1_rot3_1pi_3p_1phot, -1);

TH1F* h1_sub_kin_e_pi_all_pimi = (TH1F*) h1_en_recon3_pimi->Clone("sub_kin_e_pi_all_pimi");
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_1pi_2p_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_2pi_1p_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_2pi_2p_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_1pi_3p_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_3pi_1p_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_1pi_1p_1phot_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_1pi_2p_1phot_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_1pi_1p_2phot_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_2pi_1p_1phot_pimi, -1);
h1_sub_kin_e_pi_all_pimi->Add(h1_rot3_1pi_3p_1phot_pimi, -1);

TH1F* h1_sub_kin_e_pi_all_pipl = (TH1F*) h1_en_recon3_pipl->Clone("sub_kin_e_pi_all_pipl");
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_1pi_2p_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_2pi_1p_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_2pi_2p_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_1pi_3p_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_3pi_1p_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_1pi_1p_1phot_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_1pi_2p_1phot_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_1pi_1p_2phot_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_2pi_1p_1phot_pipl, -1);
h1_sub_kin_e_pi_all_pipl->Add(h1_rot3_1pi_3p_1phot_pipl, -1);

gDirectory->Write("hist_Files", TObject::kOverwrite);
}

Bool_t GetEPhiLimits(Float_t momentum, Float_t theta, Int_t sector,Float_t *EPhiMin, Float_t *EPhiMax){
//Begin_Html
/*</pre>
 Information for electron fiducial cut,
    returns the minimum and maximum phi accepted for a given momentum, theta and sector
    momentum is in GeV/c, theta is in degrees, 0 <= sector <= 5
    EPhiMin and EPhiMax are in degrees
    Function returns False if inputs are out of bounds
    1.1 GeV not implemented yet
 tested against EFiducialCut to make sure the limits are identical
    2.2 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.1 < Ef < 2.261
             2 inconsistent events out of 10^6
    4.4 GeV: tested for 10 < theta < 65, -30 < phi < 360, 0.3 < Ef < 4.461
             0 inconsistent events out of 10^6
 Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
 For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
<pre>
*/
//End_Html
  if (sector < 0 || sector > 5) return kFALSE;    // bad input

  if(en_beam[fbeam_E]>4. && en_beam[fbeam_E]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){// 4.4GeV fiducial cuts by protopop@jlab.org
    if ((theta < 15.) || (momentum < 0.9)) return kFALSE;         // out of range
    Float_t t0, t1, b[2], a[2];

    if (momentum > 3.7) momentum = 3.7; // don't extrapolate past the data


    // uncomment this if you want 100MeV energy bins
    //Enrgy = 0.100*int(Enrgy/0.100);


    // calculates parameters of cut functions for this energy
    t0 = fgPar_4Gev_2250_Efid_t0_p[sector][0]/pow(momentum, fgPar_4Gev_2250_Efid_t0_p[sector][1]);
    t1 = 0.; for(int k=0; k<6; k++) t1 += (fgPar_4Gev_2250_Efid_t1_p[sector][k]*pow(momentum, k));
    for(int l=0; l<2; l++){
      b[l] = 0.; for(int k=0; k<6; k++) b[l] += (fgPar_4Gev_2250_Efid_b_p[sector][l][k]*pow(momentum, k));
      a[l] = 0.; for(int k=0; k<6; k++) a[l] += (fgPar_4Gev_2250_Efid_a_p[sector][l][k]*pow(momentum, k));
    }



    // adjust upper limit according to hardware
    if(t1 < 45.) t1 = 45.;
    if(t0 < theta && theta < t1){

      *EPhiMin = 60.*sector - b[0]*(1. - 1/((theta - t0)/(b[0]/a[0]) + 1.));
      *EPhiMax = 60.*sector + b[1]*(1. - 1/((theta - t0)/(b[1]/a[1]) + 1.));
      // if(momentum<1.65 && momentum>1.60)cout<<sector<<"  "<<a[0]<<"    "<<a[1]<<"    "<<a[2]<<endl;
    }
    else {
      *EPhiMin = 60.*sector;
      *EPhiMax = 60.*sector;
    }


  }   // 4.4 GeV e2a
  else {
    return kFALSE;     // wrong beam energy/torus
  }
  return kTRUE;
}

Bool_t PimiFiducialCut(TVector3 momentum, Float_t *pimi_philow, Float_t *pimi_phiup)
{
  // Electron fiducial cut, return kTRUE if pass or kFALSE if not
  Bool_t status = kTRUE;

 if(en_beam[fbeam_E]>1. &&  en_beam[fbeam_E]<2.) {


   TVector3 mom = momentum;
double phi = mom.Phi();
  if (phi < -M_PI/6.) phi+= 2.*M_PI;
  int sector = (phi+M_PI/6.)/(M_PI/3.);
  sector = sector%6;
  double phi_deg = phi * 180./M_PI;
  phi_deg -= sector*60;

  double theta = mom.Theta();
  double theta_deg = theta * 180./M_PI;
  double mom_e = mom.Mag();


  if( fTorusCurrent>1490 && fTorusCurrent<1510){

    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    if (mom_e < .15)
      mom_e = .15;
    if (mom_e > 1.1)
      mom_e = 1.1;

    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_1500_Pimfid[sector][phipar][mompar]*pow(mom_e,mompar);
        //std::cout << mom_e << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi_deg<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta_deg-phipars[4])/phipars[3]+1.)));
      status = ((phi_deg>phicutoff) && (theta_deg>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta_deg-phipars[4])/phipars[2]+1.)));
      status = ((phi_deg<phicutoff) && (theta_deg>phipars[4]));
    }
    if (mom_e >= .3)
      {
        bool SCpdcut = true;
        if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
          if (status){
            int tsector = sector + 1;
            // sector 3 has two bad paddles
            if (tsector == 3){
              float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
              for (int i=0; i<4; i++){
                badpar3[i] = 0;
                // calculate the parameters using pol7
                for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom_e + fgPar_1gev_1500_Pimfid_Theta_S3[i][d];}
              }
              for(int ipar=0;ipar<2;ipar++)
                status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
            }
            // sector 4 has one bad paddle
            else if (tsector == 4){
              float badpar4[2];     // 2 parameters to determine the position of the theta gap
              for (int i=0; i<2; i++){
                badpar4[i] = 0;
                // calculate the parameters using pol7
                for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom_e + fgPar_1gev_1500_Pimfid_Theta_S4[i][d];}
              }
              status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
            }
            // sector 5 has four bad paddles
            else if (tsector == 5){
              Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
              for (Int_t i=0; i<8; i++){
                badpar5[i] = 0;
                // calculate the parameters using pol7
                for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom_e + fgPar_1gev_1500_Pimfid_Theta_S5[i][d];}
              }
              if (mom_e<1.25) badpar5[0] = 23.4*1500/2250;
              if (mom_e<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.
              for(Int_t ipar=0;ipar<4;ipar++)
                status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
            }
          }
        }
        return status;
      }
    else{
      bool SCpdcut = true;
      if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
        if (status){
          int tsector = sector + 1;
          // sector 3 has two bad paddles
          if (tsector == 3){
            float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
            for (int i=0; i<4; i++){
              badpar3[i] = 0;
              // calculate the parameters using 1/p
              badpar3[i] = fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][1]/mom_e + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            for(int ipar=0;ipar<2;ipar++)
              status = status && !(theta_deg>badpar3[2*ipar] && theta_deg<badpar3[2*ipar+1]);
          }
          // sector 4 has one bad paddle
          else if (tsector == 4){
            float badpar4[2];     // 2 parameters to determine the position of the theta gap
            for (int i=0; i<2; i++){
              badpar4[i] = 0;
              // calculate the parameters using 1/p
              badpar4[i] = fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][1]/mom_e + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            status = !(theta_deg>badpar4[0] && theta_deg<badpar4[1]);
          }
          // sector 5 has four bad paddles
          else if (tsector == 5){
            Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
            for (Int_t i=0; i<8; i++){
              badpar5[i] = 0;
              // calculate the parameters using 1/p
              badpar5[i] = fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][0] + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][1]/mom_e + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][2]/(mom_e*mom_e) + fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][3]/(mom_e*mom_e*mom_e);
            }
            if (mom_e<1.25) badpar5[0] = 23.4*1500/2250;
            if (mom_e<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.
            for(Int_t ipar=0;ipar<4;ipar++)
              status = status && !(theta_deg>badpar5[2*ipar] && theta_deg<badpar5[2*ipar+1]);
          }
        }
      }
      return (status);
    }
  }
  if (fTorusCurrent>740 && fTorusCurrent<760){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    if (mom_e < .15)
      mom_e = .15;
    if (mom_e > 1.1)
      mom_e = 1.1;
    //    std::cout << mom_e << std::endl;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_750_Pimfid[sector][phipar][mompar]*pow(mom_e,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }


    Float_t p_theta=mom_e, thetamax=0;
   if (p_theta>0.7)  p_theta=0.7;
   else if(p_theta<0.1)   p_theta=0.1;
   for(int i=4;i>=0;i--)thetamax=thetamax*p_theta+pimi_thetamax1[i];

    Int_t uplow;
    Double_t phicutoff;
    if(phi_deg<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta_deg-phipars[4])/phipars[3]+1.)));
      status = ((phi_deg>phicutoff) && (theta_deg>phipars[4]) && theta_deg<=thetamax);
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta_deg-phipars[4])/phipars[2]+1.)));
      status = ((phi_deg<phicutoff) && (theta_deg>phipars[4]) && theta_deg<=thetamax);
    }
    ///pasted here


 bool SCpdcut = true;
	if (SCpdcut ){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

	  if (status){

	    int tsector = sector + 1;
	    mom_e = mom.Mag();
	    //if(mom_e<0.2)mom_e=0.2;


	    //sector 2 has one gap
	    if(tsector == 2){
	      double parsec2_l,parsec2_h;
	      if(mom_e<0.4)mom_e=0.4;
		parsec2_l= fid_1gev_750_efid_S2[0][0]+fid_1gev_750_efid_S2[0][1]/mom_e +fid_1gev_750_efid_S2[0][2]/(mom_e*mom_e) +fid_1gev_750_efid_S2[0][3]/(mom_e*mom_e*mom_e);
		parsec2_h= fid_1gev_750_efid_S2[1][0]+fid_1gev_750_efid_S2[1][1]/mom_e +fid_1gev_750_efid_S2[1][2]/(mom_e*mom_e) +fid_1gev_750_efid_S2[1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta_deg>parsec2_l && theta_deg<parsec2_h);

	    }
	    //sector 3 has four gaps, first gap is due to CC so should be applied only on pimi
	    if(tsector == 3){
	      double parsec3_l[4],parsec3_h[4];
	      for(int d=1;d<4;d++){    //first gap is due to CC so should be applied only on pimi
		mom_e = mom.Mag();
		if(d>1 && mom_e>0.3 )mom_e=0.3;
		else if(d<2 && mom_e<0.45 )mom_e=0.45;
		parsec3_l[d]= fid_1gev_750_efid_S3[d][0][0]+fid_1gev_750_efid_S3[d][0][1]/mom_e +fid_1gev_750_efid_S3[d][0][2]/(mom_e*mom_e) +fid_1gev_750_efid_S3[d][0][3]/(mom_e*mom_e*mom_e);
		parsec3_h[d]= fid_1gev_750_efid_S3[d][1][0]+fid_1gev_750_efid_S3[d][1][1]/mom_e +fid_1gev_750_efid_S3[d][1][2]/(mom_e*mom_e) +fid_1gev_750_efid_S3[d][1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta_deg>parsec3_l[d] && theta_deg<parsec3_h[d]);
	      }
	    }
	    //sector 4 has two gaps
	    else if(tsector == 4){
	      double parsec4_l[2],parsec4_h[2];
	      for(int d=0;d<2;d++){
		mom_e = mom.Mag();
		if(d==1 && mom_e>0.25 )mom_e=0.25;
		else if(d==0 && mom_e<0.775 )mom_e=0.775;
		parsec4_l[d]= fid_1gev_750_efid_S4[d][0][0]+fid_1gev_750_efid_S4[d][0][1]/mom_e +fid_1gev_750_efid_S4[d][0][2]/(mom_e*mom_e) +fid_1gev_750_efid_S4[d][0][3]/(mom_e*mom_e*mom_e);
		parsec4_h[d]= fid_1gev_750_efid_S4[d][1][0]+fid_1gev_750_efid_S4[d][1][1]/mom_e +fid_1gev_750_efid_S4[d][1][2]/(mom_e*mom_e) +fid_1gev_750_efid_S4[d][1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta_deg>parsec4_l[d] && theta_deg<parsec4_h[d]);
	      }
	    }
	    //sector 5 has three gaps
	    else if(tsector == 5){
	      double parsec5_l[3],parsec5_h[3];
	      for(int d=0;d<3;d++){
		mom_e = mom.Mag();
		if(d==0 && mom_e>0.3)mom_e=0.3;
		if(d>0 && mom_e<0.5)mom_e=0.5;
		parsec5_l[d]= fid_1gev_750_efid_S5[d][0][0]+fid_1gev_750_efid_S5[d][0][1]/mom_e +fid_1gev_750_efid_S5[d][0][2]/(mom_e*mom_e) +fid_1gev_750_efid_S5[d][0][3]/(mom_e*mom_e*mom_e);
		parsec5_h[d]= fid_1gev_750_efid_S5[d][1][0]+fid_1gev_750_efid_S5[d][1][1]/mom_e +fid_1gev_750_efid_S5[d][1][2]/(mom_e*mom_e) +fid_1gev_750_efid_S5[d][1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta_deg>parsec5_l[d] && theta_deg<parsec5_h[d]);
	      }
	    }

	  }
	}

 return (status);
  }
 }


   if ( en_beam[fbeam_E]>2. &&  en_beam[fbeam_E]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){

 Float_t phi=momentum.Phi()*180./TMath::Pi();
    if(phi<-30.) phi+=360.;
    Int_t sector = (Int_t)((phi+30.)/60.);
    if(sector<0)sector=0;
    if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180./TMath::Pi();
   Float_t mom = momentum.Mag(), phimin, phimax;

    Float_t p_theta=mom, thetamax=0;
   if (p_theta>2.075)  p_theta=2.075;
   else if(p_theta<0.1)   p_theta=0.1;
   for(int i=4;i>=0;i--)thetamax=thetamax*p_theta+pimi_thetamax2and4[i]; //upper theta limit for pi- at different p


   if(mom>0.35){   //theta vs phi outline for high p region obtained by Bin

   if(mom>2.)mom=2.; //to extrapolate the cut to higher momenta for pimi

    Float_t par[6];               // six parameters to determine the outline of Theta vs Phi
    for (Int_t i=0; i<6; i++){
      par[i] = 0;
      for (Int_t d=8; d>=0; d--){
	par[i] = par[i]*mom + fgPar_2GeV_2250_Efid[sector][i][d];
      }                          // calculate the parameters using pol8
    }
    if (phi < 0) {
      Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
      phimin =  par[3]/((theta-par[0])+par[3]/par[2])-par[2];
      phimax =  par[2]-par[3]/((theta-par[0])+par[3]/par[2]);
      *pimi_philow = phimin;
      *pimi_phiup = phimax;
      status = (theta>tmptheta && tmptheta>=par[0] && theta<=thetamax);
    }
    else {
      Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
      phimin =  par[5]/((theta-par[0])+par[5]/par[4])-par[4];
      phimax =  par[4]-par[5]/((theta-par[0])+par[5]/par[4]);
      *pimi_philow = phimin;
      *pimi_phiup = phimax;
      status = (theta>tmptheta && tmptheta>=par[0] && theta<=thetamax);
    }


   }//end of high momentum cut



   if(mom<=0.35){     //theta vs phi outline for low p obtained by me

 if(mom>0.325)mom=0.325;
 else if (mom<0.125)mom=0.125;

  Float_t params[6];               // six parameters to determine the outline of Theta vs Phi
    for (Int_t i=0; i<6; i++){
      params[i] = 0;
      for (Int_t d=4; d>=0; d--){
	params[i] = params[i]*mom + fid_2gev_2250_pimifid_outline[sector][i][d];
      }                          // calculate the parameters using pol4
    }

  if (phi < 0) {
      phimin =  params[3]/((theta-params[0])+params[3]/params[1])-params[1];
      status = (phi>phimin && theta>params[0] && theta<=thetamax);
    }
    else {
      phimax = params[2]-params[4]/((theta-params[0])+params[4]/params[2]);
      status = (phi<phimax && theta>params[0]  && theta<=thetamax);
    }


    }//end of low momentum cut



 ////////////////////////////////Remove bad TOF paddles //////////////////////////////////

   bool SCpdcut=true;
   Int_t tsector = sector + 1;


  // by now, we have checked if the electron is within the outline of theta vs phi plot
    if (SCpdcut){  // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
      if (status){

	//gaps obtained with e-

	mom=momentum.Mag();
	if(mom>2.)mom=2.; //to extrapolate the cut to higher momenta for pimi
	if(mom<0.35)mom=0.35; //to extrapolate the cut to higher momenta for pimi

	if (tsector == 3){               // sector 3 has two bad paddles
	  Float_t badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	  for (Int_t i=0; i<4; i++){
	    badpar3[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar3[i] = badpar3[i]*mom + fgPar_2GeV_2250_EfidTheta_S3[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  Int_t ipar;
	  for(Int_t ipar=0;ipar<2;ipar++){
	    if(!(ipar==1 && momentum.Mag()<0.35))  status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
}
	}
	else if (tsector == 4){         // sector 4 has one bad paddle
	  Float_t badpar4[2];           // 2 parameters to determine the position of the theta gap
	  for (Int_t i=0; i<2; i++){
	    badpar4[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar4[i] = badpar4[i]*mom + fgPar_2GeV_2250_EfidTheta_S4[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  status = !(theta>badpar4[0] && theta<badpar4[1]);
	}
	else if (tsector == 5){         // sector 5 has four bad paddles
	  Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar5[i] = badpar5[i]*mom + fgPar_2GeV_2250_EfidTheta_S5[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  if (mom<1.25) badpar5[0] = 23.4;
	  if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	  for(Int_t ipar=0;ipar<4;ipar++)
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}

	//gaps obtained with pi-
	    double mom_e = momentum.Mag();



	    //sector 4 has one gap
	    if(tsector == 4){
	      double parsec4_l,parsec4_h;
	      if(mom_e<0.15)mom_e=0.15;
	      else if(mom_e>0.5)mom_e=0.5;
	      parsec4_l= fid_2gev_2250_efid_S4[0][0]+fid_2gev_2250_efid_S4[0][1]/mom_e +fid_2gev_2250_efid_S4[0][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S4[0][3]/(mom_e*mom_e*mom_e);
		parsec4_h= fid_2gev_2250_efid_S4[1][0]+fid_2gev_2250_efid_S4[1][1]/mom_e +fid_2gev_2250_efid_S4[1][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S4[1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta>parsec4_l && theta<parsec4_h);

	    }
	    //sector 3 has 2 gaps
	    else if(tsector == 3){
	      double parsec3_l[2],parsec3_h[2];
	      for(int d=0;d<2;d++){    //first gap is due to CC so should be applied only on pimi
		mom_e = momentum.Mag();
		if(mom_e<0.15 )mom_e=0.15;
		else if(d==0 && mom_e>0.55 )mom_e=0.55;
		else if(d==1 && mom_e<0.25 )mom_e=0.25;
		parsec3_l[d]= fid_2gev_2250_efid_S3[d][0][0]+fid_2gev_2250_efid_S3[d][0][1]/mom_e +fid_2gev_2250_efid_S3[d][0][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S3[d][0][3]/(mom_e*mom_e*mom_e);
		parsec3_h[d]= fid_2gev_2250_efid_S3[d][1][0]+fid_2gev_2250_efid_S3[d][1][1]/mom_e +fid_2gev_2250_efid_S3[d][1][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S3[d][1][3]/(mom_e*mom_e*mom_e);
		if(!(d==1 && mom_e>0.35))status=status && !(theta>parsec3_l[d] && theta<parsec3_h[d]);

	      }
	    }
	    //sector 5 has two gaps
	    else if(tsector == 5){
	      double parsec5_l[2],parsec5_h[2];
	      for(int d=0;d<2;d++){
		mom_e = momentum.Mag();
		if( mom_e<0.275)mom_e=0.275;
		else if(d==0 && mom_e>0.425)mom_e=0.425;
		else if(d==1 && mom_e>0.475)mom_e=0.475;
		parsec5_l[d]= fid_2gev_2250_efid_S5[d][0][0]+fid_2gev_2250_efid_S5[d][0][1]/mom_e +fid_2gev_2250_efid_S5[d][0][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S5[d][0][3]/(mom_e*mom_e*mom_e);
		parsec5_h[d]= fid_2gev_2250_efid_S5[d][1][0]+fid_2gev_2250_efid_S5[d][1][1]/mom_e +fid_2gev_2250_efid_S5[d][1][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S5[d][1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta>parsec5_l[d] && theta<parsec5_h[d]);
	      }
	    }
	    //sector 1 has two gaps
	    else if(tsector == 1){
	      double parsec1_l[2],parsec1_h[2];
	      for(int d=0;d<2;d++){
		mom_e = momentum.Mag();
		if(mom_e<0.1 )mom_e=0.1;
		else if(mom_e>0.3)mom_e=0.3;
		parsec1_l[d]= fid_2gev_2250_efid_S1[d][0][0]+fid_2gev_2250_efid_S1[d][0][1]/mom_e +fid_2gev_2250_efid_S1[d][0][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S1[d][0][3]/(mom_e*mom_e*mom_e);
		parsec1_h[d]= fid_2gev_2250_efid_S1[d][1][0]+fid_2gev_2250_efid_S1[d][1][1]/mom_e +fid_2gev_2250_efid_S1[d][1][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S1[d][1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta>parsec1_l[d] && theta<parsec1_h[d]);
	      }
	    }


	    //sector 6 has one gap
	    if(tsector == 6){
	      double parsec6_l,parsec6_h;
	      mom_e = momentum.Mag();
	      if(mom_e<0.175)mom_e=0.175;
	      else if(mom_e>0.275)mom_e=0.275;
	      parsec6_l= fid_2gev_2250_efid_S6[0][0]+fid_2gev_2250_efid_S6[0][1]/mom_e +fid_2gev_2250_efid_S6[0][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S6[0][3]/(mom_e*mom_e*mom_e);
		parsec6_h= fid_2gev_2250_efid_S6[1][0]+fid_2gev_2250_efid_S6[1][1]/mom_e +fid_2gev_2250_efid_S6[1][2]/(mom_e*mom_e) +fid_2gev_2250_efid_S6[1][3]/(mom_e*mom_e*mom_e);
		status=status && !(theta>parsec6_l && theta<parsec6_h);
	    }


	    ///end of gaps with pi-


      }
    }


    ////////////////////end of bad TOF paddle removal




   }   //end of 2 and 4GEV fiducial
  return status;
}

Bool_t PiplFiducialCut(TVector3 momentum, Float_t *philow, Float_t *phiup){
  //Positive Hadron Fiducial Cut
  //Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).

  Bool_t status = kTRUE;

 if(en_beam[fbeam_E]>1. && en_beam[fbeam_E]<2.){

	Float_t theta = momentum.Theta()*180/M_PI;
	Float_t phi   = momentum.Phi()  *180/M_PI;
	if(phi<-30) phi+=360;
	Int_t sector = Int_t ((phi+30)/60);
	if(sector<0) sector=0;
	if(sector>5) sector=5;
	phi -= sector*60;
	Float_t p = momentum.Mag();


if (fTorusCurrent < 1510 && fTorusCurrent > 1490){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < 0.15)p=0.15;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_1500_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = p;          // Momentum for bad sc paddles cuts
			if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
			if(tsector==2){      // sector 2 has one bad paddle
				Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
				for (Int_t i=0; i<2; i++){
					badpar2[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS2[i][d];
					}                // calculate the parameters using pol5
				}
				status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
			}
			else if(tsector==3){ // sector 3 has four bad paddles
				Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar3[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS3[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
				}
			}
			else if(tsector==4){ // sector 4 has two bad paddles
				Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<4; i++){
          if (i==0 || i==1)
            if (mom_scpd > .65)
              mom_scpd = .65;
					badpar4[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS4[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<2;ipar++){
					status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
				}
			}
			else if(tsector==5){ // sector 5 has four bad paddles
				Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar5[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS5[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
				}
			}


		}

    return status;
  }


  if (fTorusCurrent < 760 && fTorusCurrent > 740){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < 0.15)p=0.15;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_750_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    Float_t p_theta=p, thetamax_p=0;
   if (p_theta>0.95)  p_theta=0.95;
   else if(p_theta<0.1)   p_theta=0.1;
   for(int i=4;i>=0;i--)thetamax_p=thetamax_p*p_theta+pipl_thetamax1[i];


    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
    }

    //have pasted here
  if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd =momentum.Mag();         // momentum for bad sc paddles cuts
		       	if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
			//  if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c





  //sector 1 has two gaps , 1st gap appears only at smaller momenta and theta>110 and affects only pipl
	     if(tsector == 1){
	      double parsec1_l[2],parsec1_h[2];
	      for(int d=0;d<2;d++){

		mom_scpd =momentum.Mag();
		if(d==0 && mom_scpd>0.45 )mom_scpd=0.45;
		else if(d==0 && mom_scpd<0.15 )mom_scpd=0.15;
		else if(d==1 && mom_scpd>0.55) mom_scpd=0.55;
		else if(d==1 && mom_scpd<0.3) mom_scpd=0.3;
		parsec1_l[d]= fid_1gev_750_pfid_S1[d][0][0]+fid_1gev_750_pfid_S1[d][0][1]/mom_scpd +fid_1gev_750_pfid_S1[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
		parsec1_h[d]= fid_1gev_750_pfid_S1[d][1][0]+fid_1gev_750_pfid_S1[d][1][1]/mom_scpd +fid_1gev_750_pfid_S1[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
		status=status && !(theta>parsec1_l[d] && theta<parsec1_h[d]);
	      }
	    }
	  //sector 2 has 1 gap1
			else  if(tsector == 2){
			    double parsec2_l,parsec2_h;

			      parsec2_l= fid_1gev_750_pfid_S2[0][0]+fid_1gev_750_pfid_S2[0][1]/mom_scpd +fid_1gev_750_pfid_S2[0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec2_h= fid_1gev_750_pfid_S2[1][0]+fid_1gev_750_pfid_S2[1][1]/mom_scpd +fid_1gev_750_pfid_S2[1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec2_l && theta<parsec2_h);
			  }

			  //sector 3 has four gaps
			  else if(tsector == 3){
			    double parsec3_l[4],parsec3_h[4];
			    for(int d=0;d<4;d++){
			      parsec3_l[d]= fid_1gev_750_pfid_S3[d][0][0]+fid_1gev_750_pfid_S3[d][0][1]/mom_scpd +fid_1gev_750_pfid_S3[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec3_h[d]= fid_1gev_750_pfid_S3[d][1][0]+fid_1gev_750_pfid_S3[d][1][1]/mom_scpd +fid_1gev_750_pfid_S3[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec3_l[d] && theta<parsec3_h[d]);
			    }
			  }

	  //sector 4 has two gaps
			  else if(tsector == 4){
			    double parsec4_l[2],parsec4_h[2];
			    for(int d=0;d<2;d++){
			      parsec4_l[d]= fid_1gev_750_pfid_S4[d][0][0]+fid_1gev_750_pfid_S4[d][0][1]/mom_scpd +fid_1gev_750_pfid_S4[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec4_h[d]= fid_1gev_750_pfid_S4[d][1][0]+fid_1gev_750_pfid_S4[d][1][1]/mom_scpd +fid_1gev_750_pfid_S4[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec4_l[d] && theta<parsec4_h[d]);
			    }
			  }

  //sector 5 has four gaps
			  else if(tsector == 5){
			    double parsec5_l[4],parsec5_h[4];
			    for(int d=0;d<4;d++){
			      mom_scpd=momentum.Mag();
			      //  if(d==0 && d==1 && mom_scpd>0.6)mom_scpd=0.6;
			      if(d==2 && mom_scpd<0.5)mom_scpd=0.5;
			      if(d==3 && mom_scpd>0.3)mom_scpd=0.3;
			      parsec5_l[d]= fid_1gev_750_pfid_S5[d][0][0]+fid_1gev_750_pfid_S5[d][0][1]/mom_scpd +fid_1gev_750_pfid_S5[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec5_h[d]= fid_1gev_750_pfid_S5[d][1][0]+fid_1gev_750_pfid_S5[d][1][1]/mom_scpd +fid_1gev_750_pfid_S5[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec5_l[d] && theta<parsec5_h[d]);
			    }
			  }

  //sector 6 has two gaps
			 else if(tsector == 6){
			   double parsec6_l[2],parsec6_h[2];
			   for(int d=0;d<2;d++){
			     mom_scpd = momentum.Mag();
			     if(mom_scpd>0.6 )mom_scpd=0.6;
			     else if(mom_scpd<0.3)mom_scpd=0.3;
			     parsec6_l[d]= fid_1gev_750_pfid_S6[d][0][0]+fid_1gev_750_pfid_S6[d][0][1]/mom_scpd +fid_1gev_750_pfid_S6[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			     parsec6_h[d]= fid_1gev_750_pfid_S6[d][1][0]+fid_1gev_750_pfid_S6[d][1][1]/mom_scpd +fid_1gev_750_pfid_S6[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			     status=status && !(theta>parsec6_l[d] && theta<parsec6_h[d]);
			   }
			 }


		}

    return status;
  }

}




  if (en_beam[fbeam_E]>2. && en_beam[fbeam_E]<3. && fTorusCurrent>2240 && fTorusCurrent<2260){

    Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();
    Float_t mom_for = p;              // momentum for forward constraints
    if (mom_for<0.3) mom_for = 0.3;   // momentum smaller than 300 MeV/c, use 300 MeV/c
    if (mom_for>1.6) mom_for = 1.6;   // momentum greater than 1.6 GeV/c, use 1.6 GeV/c
    Float_t mom_bak = p;              // momentum for backward constraints
    if (mom_bak<0.2) mom_bak = 0.2;   // momentum smaller than 200 MeV/c, use 200 MeV/c
    if (mom_bak>1.0) mom_bak = 1.0;   // momentum greater than 1.0 GeV/c, use 1.0 GeV/c
    Float_t theta0 = 8.5;
    Float_t phi_lower = -24.0;
    Float_t phi_upper = 24.0;
    Float_t phimin, phimax;
    Float_t par_for[4], par_bak[4];
    for (Int_t i=0; i<4; i++){
      par_for[i] = 0; par_bak[i] = 0;
      for (Int_t d=6; d>=0; d--){
	par_for[i] = par_for[i]*mom_for + fgPar_2GeV_2250_Pfid_For[sector][i][d];
	par_bak[i] = par_bak[i]*mom_bak + fgPar_2GeV_2250_Pfid_Bak[sector][i][d];
      }
    }
    if (phi < 0) {
      Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
       phimin = par_for[1]/((theta-theta0)+par_for[1]/par_for[0])-par_for[0];
       phimax = par_for[0]-par_for[1]/((theta-theta0)+par_for[1]/par_for[0]);
       *philow = phimin;
       *phiup = phimax;
      status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
    }
    else {
      Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
      phimin = par_for[3]/(theta-theta0+par_for[3]/par_for[2])-par_for[2];
      phimax = par_for[2]-par_for[3]/(theta-theta0+par_for[3]/par_for[2]);
      *phiup = phimax;
      *philow = phimin;
      status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);
    }                     // now the forward constrains are checked
    if ( status ) {       // now check the backward constrains
      if(theta>par_bak[0]) status = kFALSE;
      else if(theta>par_bak[1]) status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]);
    }

    if(status && SCpdcut){ // cut bad scintillator paddles

      Int_t tsector = sector + 1;
      Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
      if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if(tsector==2){      // sector 2 has one bad paddle
	Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
		for (Int_t i=0; i<2; i++){
	  badpar2[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar2[i] = badpar2[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS2[i][d];
	  }                // calculate the parameters using pol5
	}
	status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
      }
      else if(tsector==3){ // sector 3 has four bad paddles
	Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar3[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar3[i] = badpar3[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS3[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
      }
      else if(tsector==4){ // sector 4 has two bad paddles
	Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<4; i++){
	  badpar4[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar4[i] = badpar4[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS4[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<2;ipar++){
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	}
      }
      else if(tsector==5){ // sector 5 has four bad paddles
	Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar5[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar5[i] = badpar5[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS5[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }

  }


  if (en_beam[fbeam_E]>4. && en_beam[fbeam_E]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){//4 GeV Fiducial Cut Rustam Niyazov

   Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = Int_t ((phi+30)/60); if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();

    Float_t parfidl[3];for(Int_t i=0; i<3; i++){parfidl[i]=0;}
    Float_t parfidr[3];for(Int_t i=0; i<3; i++){parfidr[i]=0;}
    Float_t parfidbl[2];for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
    Float_t parfidbr[2];for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
    Float_t cphil=0;Float_t cphir=0;
    Float_t phi45l=0; Float_t phi45r=0;
    Float_t phi60l=0; Float_t phi60r=0;
    Float_t theta_min=11;

    bool Forward=kFALSE; //defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
    Int_t thetab=45; //this variable defines the edge point for Forward<->Backward regions
    Float_t p1=0.575; //last bin momentum for region p<0.6 GeV/c
    Float_t theta_max=140;
    if(p<0.2)p=0.2; //momentum less than 0.2 GeV/c, use 0.2 GeV/c
    if(p>4.4)p=4.4; //momentum greater than 4.4 GeV/c, use 4.4 GeV/c

    //get parametrized values of theta_max for p<0.6 GeV/c region
    if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
    //get parametrized values of theta_max for p>0.6 GeV/c region
    else{theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p1;}

    //Get the momentum dependent parameters for Forward Region (theta <45 deg)
    Forward=kTRUE;
    if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
      //parameters for hyperbolic function
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;
      }
    }
    else{//forward2 defines  regions of momenta and p>0.6 GeV/c
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;
      }
    }
    phi45l=parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
    phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
    if(theta>thetab){//backward region defined by theta >45 deg.
      if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
      if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

      //Get the momentum dependent parameters for Backward Region

      Forward=kFALSE;
      if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;
        }
      }
      else{//backward2 defines  regions of momenta p>0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;
        }
      }
    }

    if(Forward){//Forward region
      if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.
        cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
        cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
    }
    else{//Backward region
      phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
      phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

      if(theta<60){
        cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
        cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
      }
      Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0;
      //dr and er are theta_flat and phi_edge parameters for phi>0;
      dl=parfidbl[0];el=parfidbl[1];
      dr=parfidbr[0];er=parfidbr[1];

      if(theta>45&&theta<60){ //BackwardA region
        //try to match parametrized values from Forward region to Backward region parameters
        if(cphil>phi45l)cphil=phi45l;
        if(cphir<phi45r)cphir=phi45r;
      }
      //BackwardB region & phi<0
      else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant
      else if(theta>dl&&theta<=theta_max){
        cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line
      else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
      //BackwardB region & phi>0
      if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant
      else if(theta>dr&&theta<=theta_max){
        cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line
      else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
    }//Backward Region


    if(phi<0) status=(phi>cphil); //check the constrains
    else if(phi>=0) {status=(phi<cphir);
  }

    if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min

    if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

   //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for
   //some range of momentum, where edge does not look good.
    bool s1s4=(theta<11.7&&(sector==0||sector==3));
    bool s5=(theta<12.2&&sector==4);
    bool s6=(theta<11.4&&sector==5);
    if(p>=0.6&&p<1.5&&(s1s4||s5||s6)) status=kFALSE;

    *philow = cphil;
    *phiup = cphir;



bool SCpdcut = true;
    if(status && SCpdcut){ // cut bad scintillator paddles
      if(p < 1.0){
        Int_t tsector = sector + 1;
        Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
        if (mom_scpd<0.3)mom_scpd=0.3; // momentum smaller than 200 MeV/c, use 200 MeV/c
        if(tsector==2){      // sector 2 has one bad paddle
          Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
          for (Int_t i=0; i<2; i++){
            badpar2[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar2[i] = badpar2[i]*mom_scpd + fgPar_Pfid_ScpdS2[i][d];
            }                // calculate the parameters using pol5
          }
          status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
        }
        else if(tsector==3){ // sector 3 has four bad paddles
          Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar3[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar3[i] = badpar3[i]*mom_scpd + fgPar_Pfid_ScpdS3[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
        }
        else if(tsector==4){ // sector 4 has two bad paddles
          Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<4; i++){
            badpar4[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar4[i] = badpar4[i]*mom_scpd + fgPar_Pfid_ScpdS4[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<2;ipar++){
            status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
        }
        else if(tsector==5){ // sector 5 has four bad paddles
          Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar5[i] = badpar5[i]*mom_scpd + fgPar_Pfid_ScpdS5[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
          }
        }
      }
      else{
        int tsector = sector + 1;
        double mom_scpd =p;
        // sector 2 has one bad paddles
        if (tsector == 2){
          float badpar2[2];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<2; i++){
            badpar2[i] = 0;
            // calculate the parameters using 1/p
            badpar2[i] = fgPar_Pfid_ScpdS2_extra[i][0] + fgPar_Pfid_ScpdS2_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS2_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS2_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<1;ipar++)
            status = status && !(theta>badpar2[2*ipar] && theta<badpar2[2*ipar+1]);
        }
        if (tsector == 3){
          float badpar3[8];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<8; i++){
            badpar3[i] = 0;
            // calculate the parameters using 1/p
            badpar3[i] = fgPar_Pfid_ScpdS3_extra[i][0] + fgPar_Pfid_ScpdS3_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS3_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS3_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
        }
        // sector 4 has two bad paddle
        else if (tsector == 4){
          float badpar4[4];     // 2 parameters to determine the position of the theta gap
          for (int i=0; i<4; i++){
            badpar4[i] = 0;
            // calculate the parameters using 1/p
            badpar4[i] = fgPar_Pfid_ScpdS4_extra[i][0] + fgPar_Pfid_ScpdS4_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS4_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS4_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<2;ipar++)
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
        // sector 5 has four bad paddles
        else if (tsector == 5){
          Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            // calculate the parameters using 1/p
            badpar5[i] = fgPar_Pfid_ScpdS5_extra[i][0] + fgPar_Pfid_ScpdS5_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS5_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS5_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(Int_t ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
        }
      }
    }



  }


  return status;
}

Bool_t PFiducialCut(TVector3 momentum){
  //Positive Hadron Fiducial Cut
  //Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/pfiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).

   Bool_t status = kTRUE;



  if(en_beam[fbeam_E]>1. && en_beam[fbeam_E]<2.){

	Float_t theta = momentum.Theta()*180/M_PI;
	Float_t phi   = momentum.Phi()  *180/M_PI;
	if(phi<-30) phi+=360;
	Int_t sector = Int_t ((phi+30)/60);
	if(sector<0) sector=0;
	if(sector>5) sector=5;
	phi -= sector*60;
	Float_t p = momentum.Mag();


if ( fTorusCurrent < 1510 && fTorusCurrent > 1490){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < .3)
      return false;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_1500_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]));
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]));
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = p;          // Momentum for bad sc paddles cuts
			if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c
			if(tsector==2){      // sector 2 has one bad paddle
				Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
				for (Int_t i=0; i<2; i++){
					badpar2[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar2[i] = badpar2[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS2[i][d];
					}                // calculate the parameters using pol5
				}
				status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
			}
			else if(tsector==3){ // sector 3 has four bad paddles
				Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar3[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar3[i] = badpar3[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS3[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
				}
			}
			else if(tsector==4){ // sector 4 has two bad paddles
				Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<4; i++){
          if (i==0 || i==1)
            if (mom_scpd > .65)
              mom_scpd = .65;
					badpar4[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar4[i] = badpar4[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS4[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<2;ipar++){
					status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
				}
			}
			else if(tsector==5){ // sector 5 has four bad paddles
				Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
				for (Int_t i=0; i<8; i++){
					badpar5[i] = 0;
					for (Int_t d=5; d>=0; d--){
						badpar5[i] = badpar5[i]*mom_scpd + fgPar_1gev_1500_Pfid_ScpdS5[i][d];
					}                // calculate the parameters using pol5
				}
				for (Int_t ipar=0;ipar<4;ipar++){
					status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
				}
			}
		}
    return status;
  }
  if (fTorusCurrent < 760 && fTorusCurrent > 740){
    Double_t phipars[5]={0,0,0,0,0};
    status = true;
    bool SCpdcut = true;
    if (p < .3)
      return false;
    if (p > 1)
      p = 1;
    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t phipar=0;phipar<5;phipar++) {
        phipars[phipar]+=fgPar_1gev_750_Pfid[sector][phipar][mompar]*pow(p,mompar);
        //std::cout << p << " " << mompar << " " << phipar << " " << phipars[1] << " " << phipars[2] << " " << phipars[3] << " " << phipars[4] << " " << phipars[5] << std::endl;
      }
    }
    Int_t uplow;
    Double_t phicutoff;
    Float_t p_theta=p, thetamax_p=0;
   if (p_theta>0.95)  p_theta=0.95;
   else if(p_theta<0.1)   p_theta=0.1;
   for(int i=4;i>=0;i--)thetamax_p=thetamax_p*p_theta+pipl_thetamax1[i];


    if(phi<=0) {
      phicutoff = phipars[1]*(1.-(1./((theta-phipars[4])/phipars[3]+1.)));
      //std::cout << "bottom " << theta << std::endl;
      status = ((phi>phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
    }
    else {
      phicutoff = phipars[0]*(1.-(1./((theta-phipars[4])/phipars[2]+1.)));
      //std::cout << "top " << phicutoff << std::endl;
      status = ((phi<phicutoff) && (theta>phipars[4]) && theta<=thetamax_p);
    }
    if(status && SCpdcut){ // cut bad scintillator paddles
			Int_t tsector = sector + 1;
			Float_t mom_scpd = momentum.Mag();          // momentum for bad sc paddles cuts
		       	if (mom_scpd<0.15)mom_scpd=0.15; // momentum smaller than 200 MeV/c, use 200 MeV/c
			//  if (mom_scpd>1.0)mom_scpd=1.0; // momentum greater than 1000 MeV/c, use 1000 MeV/c




			//sector 1 has two gaps , 1st gap appears only at smaller momenta and theta>110 and affects only pipl
			if(tsector == 1){
			  double parsec1_l[2],parsec1_h[2];
			  for(int d=0;d<2;d++){

			    mom_scpd =momentum.Mag();
			    if(d==0 && mom_scpd>0.45 )mom_scpd=0.45;
			    else if(d==0 && mom_scpd<0.15 )mom_scpd=0.15;
			    else if(d==1 && mom_scpd>0.55) mom_scpd=0.55;
			    else if(d==1 && mom_scpd<0.3) mom_scpd=0.3;
			    parsec1_l[d]= fid_1gev_750_pfid_S1[d][0][0]+fid_1gev_750_pfid_S1[d][0][1]/mom_scpd +fid_1gev_750_pfid_S1[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			    parsec1_h[d]= fid_1gev_750_pfid_S1[d][1][0]+fid_1gev_750_pfid_S1[d][1][1]/mom_scpd +fid_1gev_750_pfid_S1[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S1[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			    status=status && !(theta>parsec1_l[d] && theta<parsec1_h[d]);
			  }
			}

	  //sector 2 has 1 gap1
			 else if(tsector == 2){
			    double parsec2_l,parsec2_h;

			      parsec2_l= fid_1gev_750_pfid_S2[0][0]+fid_1gev_750_pfid_S2[0][1]/mom_scpd +fid_1gev_750_pfid_S2[0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec2_h= fid_1gev_750_pfid_S2[1][0]+fid_1gev_750_pfid_S2[1][1]/mom_scpd +fid_1gev_750_pfid_S2[1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S2[1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec2_l && theta<parsec2_h);
			  }

			  //sector 3 has four gaps
			 else if(tsector == 3){
			    double parsec3_l[4],parsec3_h[4];
			    for(int d=0;d<4;d++){
			      parsec3_l[d]= fid_1gev_750_pfid_S3[d][0][0]+fid_1gev_750_pfid_S3[d][0][1]/mom_scpd +fid_1gev_750_pfid_S3[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec3_h[d]= fid_1gev_750_pfid_S3[d][1][0]+fid_1gev_750_pfid_S3[d][1][1]/mom_scpd +fid_1gev_750_pfid_S3[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S3[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec3_l[d] && theta<parsec3_h[d]);
			    }
			  }

	  //sector 4 has two gaps
			else  if(tsector == 4){
			    double parsec4_l[2],parsec4_h[2];
			    for(int d=0;d<2;d++){
			      parsec4_l[d]= fid_1gev_750_pfid_S4[d][0][0]+fid_1gev_750_pfid_S4[d][0][1]/mom_scpd +fid_1gev_750_pfid_S4[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec4_h[d]= fid_1gev_750_pfid_S4[d][1][0]+fid_1gev_750_pfid_S4[d][1][1]/mom_scpd +fid_1gev_750_pfid_S4[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S4[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec4_l[d] && theta<parsec4_h[d]);
			    }
			  }

  //sector 5 has four gaps
			 else if(tsector == 5){//the fourth bad TOF pd. can be seen only below b=0.3 and so there are just three bad TOFs for p
			    double parsec5_l[3],parsec5_h[3];
			    for(int d=0;d<3;d++){
			      //  if(d==0 && d==1 && mom_scpd>0.6)mom_scpd=0.6;
			      mom_scpd = momentum.Mag();
			      if(d==2 && mom_scpd<0.5)mom_scpd=0.5;
			      parsec5_l[d]= fid_1gev_750_pfid_S5[d][0][0]+fid_1gev_750_pfid_S5[d][0][1]/mom_scpd +fid_1gev_750_pfid_S5[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			      parsec5_h[d]= fid_1gev_750_pfid_S5[d][1][0]+fid_1gev_750_pfid_S5[d][1][1]/mom_scpd +fid_1gev_750_pfid_S5[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S5[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			      status=status && !(theta>parsec5_l[d] && theta<parsec5_h[d]);
			    }
			  }

			  //sector 6 has two gaps
			 else if(tsector == 6){
			   double parsec6_l[2],parsec6_h[2];
			   for(int d=0;d<2;d++){
			     mom_scpd = momentum.Mag();
			     if(mom_scpd>0.6 )mom_scpd=0.6;
			     else if(mom_scpd<0.3)mom_scpd=0.3;
			     parsec6_l[d]= fid_1gev_750_pfid_S6[d][0][0]+fid_1gev_750_pfid_S6[d][0][1]/mom_scpd +fid_1gev_750_pfid_S6[d][0][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][0][3]/(mom_scpd*mom_scpd*mom_scpd);
			     parsec6_h[d]= fid_1gev_750_pfid_S6[d][1][0]+fid_1gev_750_pfid_S6[d][1][1]/mom_scpd +fid_1gev_750_pfid_S6[d][1][2]/(mom_scpd*mom_scpd) +fid_1gev_750_pfid_S6[d][1][3]/(mom_scpd*mom_scpd*mom_scpd);
			     status=status && !(theta>parsec6_l[d] && theta<parsec6_h[d]);
			   }
			 }


		}
    return status;
  }

}




  if (en_beam[fbeam_E]>2. && en_beam[fbeam_E]<3. && fTorusCurrent>2240 && fTorusCurrent<2260){
    Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();
    Float_t mom_for = p;              // momentum for forward constraints
    if (mom_for<0.3) mom_for = 0.3;   // momentum smaller than 300 MeV/c, use 300 MeV/c
    if (mom_for>1.6) mom_for = 1.6;   // momentum greater than 1.6 GeV/c, use 1.6 GeV/c
    Float_t mom_bak = p;              // momentum for backward constraints
    if (mom_bak<0.2) mom_bak = 0.2;   // momentum smaller than 200 MeV/c, use 200 MeV/c
    if (mom_bak>1.0) mom_bak = 1.0;   // momentum greater than 1.0 GeV/c, use 1.0 GeV/c
    Float_t theta0 = 8.5;
    Float_t phi_lower = -24.0;
    Float_t phi_upper = 24.0;
    Float_t par_for[4], par_bak[4];
    for (Int_t i=0; i<4; i++){
      par_for[i] = 0; par_bak[i] = 0;
      for (Int_t d=6; d>=0; d--){
	par_for[i] = par_for[i]*mom_for + fgPar_2GeV_2250_Pfid_For[sector][i][d];
	par_bak[i] = par_bak[i]*mom_bak + fgPar_2GeV_2250_Pfid_Bak[sector][i][d];
      }
    }
    if (phi < 0) {
      Float_t tmptheta = theta0 - par_for[1]/par_for[0] + par_for[1]/(par_for[0]+phi);
      status = (theta>tmptheta && tmptheta>=theta0 && phi>=phi_lower);
    }
    else {
      Float_t tmptheta = theta0 - par_for[3]/par_for[2] + par_for[3]/(par_for[2]-phi);
      status = (theta>tmptheta && tmptheta>=theta0 && phi<=phi_upper);
    }                     // now the forward constrains are checked
    if ( status ) {       // now check the backward constrains
      if(theta>par_bak[0]) status = kFALSE;
      else if(theta>par_bak[1]) status = (phi-phi_lower)/(theta-par_bak[1])>=(par_bak[2]-phi_lower)/(par_bak[0]-par_bak[1]) && (phi-phi_upper)/(theta-par_bak[1])<=(par_bak[3]-phi_upper)/(par_bak[0]-par_bak[1]);
    }

    if(status && SCpdcut){ // cut bad scintillator paddles

      Int_t tsector = sector + 1;
      Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
      if (mom_scpd<0.2)mom_scpd=0.2; // momentum smaller than 200 MeV/c, use 200 MeV/c
      if(tsector==2){      // sector 2 has one bad paddle
	Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
	for (Int_t i=0; i<2; i++){
	  badpar2[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar2[i] = badpar2[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS2[i][d];
	  }                // calculate the parameters using pol5
	}
	status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
      }
      else if(tsector==3){ // sector 3 has four bad paddles
	Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar3[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar3[i] = badpar3[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS3[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
      }
      else if(tsector==4){ // sector 4 has two bad paddles
	Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<4; i++){
	  badpar4[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar4[i] = badpar4[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS4[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<2;ipar++){
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
	}
      }
      else if(tsector==5){ // sector 5 has four bad paddles
	Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
	for (Int_t i=0; i<8; i++){
	  badpar5[i] = 0;
	  for (Int_t d=5; d>=0; d--){
	    badpar5[i] = badpar5[i]*mom_scpd + fgPar_2GeV_2250_Pfid_ScpdS5[i][d];
	  }                // calculate the parameters using pol5
	}
	for (Int_t ipar=0;ipar<4;ipar++){
	  status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }
  }



  if (en_beam[fbeam_E]>4. && en_beam[fbeam_E]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){//4 GeV Fiducial Cut Rustam Niyazov

    Float_t phi=momentum.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
    Int_t sector = Int_t ((phi+30)/60); if(sector<0)sector=0; if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180/TMath::Pi();
    Float_t p = momentum.Mag();

    Float_t parfidl[3];for(Int_t i=0; i<3; i++){parfidl[i]=0;}
    Float_t parfidr[3];for(Int_t i=0; i<3; i++){parfidr[i]=0;}
    Float_t parfidbl[2];for(Int_t i=0; i<2; i++){parfidbl[i]=0;}
    Float_t parfidbr[2];for(Int_t i=0; i<2; i++){parfidbr[i]=0;}
    Float_t cphil=0;Float_t cphir=0;
    Float_t phi45l=0; Float_t phi45r=0;
    Float_t phi60l=0; Float_t phi60r=0;
    Float_t theta_min=11;

    bool Forward=kFALSE; //defines if particle in Forward (Forward=kTRUE) or Backward (Forward=kFALSE) region.
    Int_t thetab=45; //this variable defines the edge point for Forward<->Backward regions
    Float_t p1=0.575; //last bin momentum for region p<0.6 GeV/c
    Float_t theta_max=140;
    if(p<0.2)p=0.2; //momentum less than 0.2 GeV/c, use 0.2 GeV/c
    if(p>4.4)p=4.4; //momentum greater than 4.4 GeV/c, use 4.4 GeV/c

    //get parametrized values of theta_max for p<0.6 GeV/c region
    if(p<0.6){theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p;}
    //get parametrized values of theta_max for p>0.6 GeV/c region
    else{theta_max=fgPar_4Gev_2250_Pfidbl[sector][4]+fgPar_4Gev_2250_Pfidbl[sector][5]*p1;}

    //Get the momentum dependent parameters for Forward Region (theta <45 deg)
    Forward=kTRUE;
    if(p<0.6){//forward1 defines  regions of momenta p<0.6 GeV/c
      //parameters for hyperbolic function
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft1l[sector][j]+fgPar_4Gev_2250_Pfidft1l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft1r[sector][j]+fgPar_4Gev_2250_Pfidft1r[sector][j+1]/p;
      }
    }
    else{//forward2 defines  regions of momenta and p>0.6 GeV/c
      for (Int_t i=0; i<3; i++){
        Int_t j=2*i;
        parfidl[i]=fgPar_4Gev_2250_Pfidft2l[sector][j]+fgPar_4Gev_2250_Pfidft2l[sector][j+1]/p;
        parfidr[i]=fgPar_4Gev_2250_Pfidft2r[sector][j]+fgPar_4Gev_2250_Pfidft2r[sector][j+1]/p;
      }
    }
    phi45l=parfidl[0]*(parfidl[2]-45)/(45-parfidl[2]+(parfidl[1]/parfidl[0])); //parametrized value of phi at theta=45 deg.
    phi45r=-parfidr[0]*(parfidr[2]-45)/(45-parfidr[2]+(parfidr[1]/parfidr[0]));
    if(theta>thetab){//backward region defined by theta >45 deg.
      if(theta>140) theta =140; //theta greater than 140 degrees, use 140 degrees
      if(p>1)p=1.; //momentum greater than 1.0 GeV/c, use 1.0 GeV/c

      //Get the momentum dependent parameters for Backward Region

      Forward=kFALSE;
      if(p<0.6){//backward1 defines  regions of momenta p<0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt1l[sector][j]+fgPar_4Gev_2250_Pfidbt1l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt1r[sector][j]+fgPar_4Gev_2250_Pfidbt1r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p<0.6 GeV/c
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p;
        }
      }
      else{//backward2 defines  regions of momenta p>0.6 GeV/c
        //parameters for quadratic function
        for (Int_t i=0; i<3; i++){
          Int_t j=2*i;
          parfidl[i]=fgPar_4Gev_2250_Pfidbt2l[sector][j]+fgPar_4Gev_2250_Pfidbt2l[sector][j+1]/p;
          parfidr[i]=fgPar_4Gev_2250_Pfidbt2r[sector][j]+fgPar_4Gev_2250_Pfidbt2r[sector][j+1]/p;
        }
        //these parameters determine theta_flat and phi_edge at p=0.575 GeV/c momentum
        for (Int_t i=0; i<2; i++){
          Int_t j=2*i;
          parfidbl[i]=fgPar_4Gev_2250_Pfidbl[sector][j]+fgPar_4Gev_2250_Pfidbl[sector][j+1]/p1;
          parfidbr[i]=fgPar_4Gev_2250_Pfidbr[sector][j]+fgPar_4Gev_2250_Pfidbr[sector][j+1]/p1;
        }
      }
    }

    if(Forward){//Forward region
      if(p<0.6) theta_min=14; else theta_min=11;//for p<0.6 GeV/c Region theta starts from 14 deg., otherwise 11 deg.
        cphil=parfidl[0]*(parfidl[2]-theta)/(theta-parfidl[2]+(parfidl[1]/parfidl[0]));//hyperbolic function
        cphir=-parfidr[0]*(parfidr[2]-theta)/(theta-parfidr[2]+(parfidr[1]/parfidr[0]));
    }
    else{//Backward region
      phi60l=parfidl[0]+ parfidl[1]*60.+ parfidl[2]*3600.;//parametrized value of phi at theta=60 deg.
      phi60r=-(parfidr[0]+ parfidr[1]*60.+ parfidr[2]*3600.);

      if(theta<60){
        cphil=parfidl[0]+ parfidl[1]*theta+ parfidl[2]*theta*theta; //quadratic function
        cphir=-(parfidr[0]+ parfidr[1]*theta+ parfidr[2]*theta*theta);
      }
      Float_t dl,el,dr,er; //dl and el are theta_flat and phi_edge parameters for phi<0;
      //dr and er are theta_flat and phi_edge parameters for phi>0;
      dl=parfidbl[0];el=parfidbl[1];
      dr=parfidbr[0];er=parfidbr[1];

      if(theta>45&&theta<60){ //BackwardA region
        //try to match parametrized values from Forward region to Backward region parameters
        if(cphil>phi45l)cphil=phi45l;
        if(cphir<phi45r)cphir=phi45r;
      }
      //BackwardB region & phi<0
      else if(theta>=60&&theta<=dl){cphil=phi60l;} //phi=constant
      else if(theta>dl&&theta<=theta_max){
        cphil=(140-theta)*(phi60l-el)/(140-dl) +el;}//phi=stright line
      else if(theta>theta_max){cphil=0;} //cut out if theta>theta_max
      //BackwardB region & phi>0
      if(theta>=60&&theta<=dr){cphir=phi60r;} //phi=constant
      else if(theta>dr&&theta<=theta_max){
        cphir=(140-theta)*(phi60r-er)/(140-dr) +er;}//phi=stright line
      else if(theta>theta_max){cphir=0;} //cut out if theta>theta_max
    }//Backward Region

    if(phi<0) status=(phi>cphil); //check the constrains
    else if(phi>=0) {status=(phi<cphir);
    }

    if(theta<theta_min) status=kFALSE; //Cutting out events below theta_min

    if(Forward && p<0.6 && theta<20.6-11.4*p)status=kFALSE; //function defines cut of the edge at low theta for p<0.6 GeV/c

   //p>0.6 GeV/c. Cut of the edge at low theta  for some sectors and for
   //some range of momentum, where edge does not look good.
    bool s1s4=(theta<11.7&&(sector==0||sector==3));
    bool s5=(theta<12.2&&sector==4);
    bool s6=(theta<11.4&&sector==5);
    if(p>=0.6&&p<1.5&&(s1s4||s5||s6)) status=kFALSE;



bool SCpdcut = true;
    if(status && SCpdcut){ // cut bad scintillator paddles
      if(p < 1.0){
        Int_t tsector = sector + 1;
        Float_t mom_scpd = p;          // momentum for bad sc paddles cuts
        if (mom_scpd<0.3)mom_scpd=0.3; // momentum smaller than 200 MeV/c, use 200 MeV/c
        if(tsector==2){      // sector 2 has one bad paddle
          Float_t badpar2[2];// 2 parameters to determine the position of the theta gap
          for (Int_t i=0; i<2; i++){
            badpar2[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar2[i] = badpar2[i]*mom_scpd + fgPar_Pfid_ScpdS2[i][d];
            }                // calculate the parameters using pol5
          }
          status = status && !(theta>badpar2[0]&&theta<badpar2[1]);
        }
        else if(tsector==3){ // sector 3 has four bad paddles
          Float_t badpar3[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar3[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar3[i] = badpar3[i]*mom_scpd + fgPar_Pfid_ScpdS3[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
          }
        }
        else if(tsector==4){ // sector 4 has two bad paddles
          Float_t badpar4[4];// 4 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<4; i++){
            badpar4[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar4[i] = badpar4[i]*mom_scpd + fgPar_Pfid_ScpdS4[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<2;ipar++){
            status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
          }
        }
        else if(tsector==5){ // sector 5 has four bad paddles
          Float_t badpar5[8];// 8 parameters to determine the positions of the theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            for (Int_t d=5; d>=0; d--){
              badpar5[i] = badpar5[i]*mom_scpd + fgPar_Pfid_ScpdS5[i][d];
            }                // calculate the parameters using pol5
          }
          for (Int_t ipar=0;ipar<4;ipar++){
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
          }
        }
      }
      else{
        int tsector = sector + 1;
        double mom_scpd =p;
        // sector 2 has one bad paddles
        if (tsector == 2){
          float badpar2[2];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<2; i++){
            badpar2[i] = 0;
            // calculate the parameters using 1/p
            badpar2[i] = fgPar_Pfid_ScpdS2_extra[i][0] + fgPar_Pfid_ScpdS2_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS2_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS2_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<1;ipar++)
	    status = status && !(theta>badpar2[2*ipar] && theta<badpar2[2*ipar+1]);
	  // status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
        if (tsector == 3){
          float badpar3[8];            // 4 parameters to determine the positions of the two theta gaps
          for (int i=0; i<8; i++){
            badpar3[i] = 0;
            // calculate the parameters using 1/p
            badpar3[i] = fgPar_Pfid_ScpdS3_extra[i][0] + fgPar_Pfid_ScpdS3_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS3_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS3_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
        }
        // sector 4 has two bad paddle
        else if (tsector == 4){
          float badpar4[4];     // 2 parameters to determine the position of the theta gap
          for (int i=0; i<4; i++){
            badpar4[i] = 0;
            // calculate the parameters using 1/p
            badpar4[i] = fgPar_Pfid_ScpdS4_extra[i][0] + fgPar_Pfid_ScpdS4_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS4_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS4_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(int ipar=0;ipar<2;ipar++)
	  status = status && !(theta>badpar4[2*ipar] && theta<badpar4[2*ipar+1]);
        }
        // sector 5 has four bad paddles
        else if (tsector == 5){
          Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
          for (Int_t i=0; i<8; i++){
            badpar5[i] = 0;
            // calculate the parameters using 1/p
            badpar5[i] = fgPar_Pfid_ScpdS5_extra[i][0] + fgPar_Pfid_ScpdS5_extra[i][1]/mom_scpd + fgPar_Pfid_ScpdS5_extra[i][2]/(mom_scpd*mom_scpd) + fgPar_Pfid_ScpdS5_extra[i][3]/(mom_scpd*mom_scpd*mom_scpd);
          }
          for(Int_t ipar=0;ipar<4;ipar++)
            status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
        }
      }
    }



  }
  return status;
}

Float_t ProtonMomCorrection_He3_4Cell(std::string atarget, TLorentzVector V4Pr, Float_t vertex_p ){

  // Low energy proton momentum correction function
  // to be used with He3 target (4th target cell) (RUN # 18338-18438)
  // Input: Proton momentum 4 vector, and Z coord of proton vertex.
  // Returns the corrected MAGNITUDE of the proton momentum,

  Float_t  up_parm[6]   = {2.001,  -14.94,  47.2,   -77.59,  65.73,  -22.85};
  Float_t  down_parm[6] = {1.4165, -13.004, 48.897, -92.443, 86.984, -32.424};

  Float_t  proton_p     = V4Pr.Vect().Mag();
  Float_t  thetta_p     = V4Pr.Vect().Theta()*57.3;
  Float_t  polinom_up   = (((((up_parm[5]*proton_p+up_parm[4])*proton_p+up_parm[3])
			  *proton_p+up_parm[2])*proton_p+up_parm[1])*proton_p+up_parm[0]);

  Float_t polinom_down = (((((down_parm[5]*proton_p+down_parm[4])*proton_p+down_parm[3])
			  *proton_p+down_parm[2])*proton_p+down_parm[1])*proton_p+down_parm[0]);


  if(polinom_up<0.  ) polinom_up   = 0;
  if(polinom_down<0.) polinom_down = 0;

  Float_t  p_corr_up   = proton_p + proton_p*polinom_up;
  Float_t  p_corr_down = proton_p + proton_p*polinom_down;

  p_corr_down=p_corr_down*4./3-proton_p/3;//artificial cut to match with Bins distribution
  p_corr_up=p_corr_down;//artificial cut to match with Bins distribution

  if((thetta_p>=70.)) return p_corr_up;

  if((thetta_p < 30.)||(vertex_p>=(1/20.*thetta_p-5/2))||
     (thetta_p<=(-200*proton_p+86))){
    return p_corr_down;
  }



  if((thetta_p<=70.)&&(thetta_p>=30)&&
     (thetta_p>(20*vertex_p+50))){
    return p_corr_up;
  }else if(proton_p<0.57){
    return p_corr_down;
  } else { return p_corr_up;}

  return -1.;
}

TLorentzVector EMomentumCorrection(TLorentzVector V4el)
{


  // Electron Momentum correction, Pass the electron 4 vector, return corrected 4 Vector pointer.

  TLorentzVector V4ecor= V4el;
  TVector3       V3el(V4el.Vect());
  Float_t p = V3el.Mag();
  Float_t cz =V3el.CosTheta();
  Float_t phi = 180*V3el.Phi()/TMath::Pi(); // (MWH 5/6/2000) Convert to degrees.
 if (phi<-30.) phi += 360;
Float_t theta =  57.29578*(V3el.Theta()); // (DP 9/12/2000) Needed for my function

  Int_t sectInd = (Int_t)(phi+30)/60;
  if(sectInd>5) sectInd = 5; if(sectInd<0) sectInd = 0;


  if(en_beam[fbeam_E]>2. && en_beam[fbeam_E]<3. && fTorusCurrent>2240 && fTorusCurrent<2260){
    //
    //the correction for 2.2Gev and 2250A data
    //
    phi -= 60*sectInd;
    p = p*(fgPar_2Gev_2250_Phi[sectInd][0] + fgPar_2Gev_2250_Phi[sectInd][1]*phi + fgPar_2Gev_2250_Phi[sectInd][2]*phi*phi);
    if(cz>0.8 && cz<0.885){
      p=p*(fgPar_2Gev_2250_Theta[0][0]*sin(fgPar_2Gev_2250_Theta[0][1]*(cz+fgPar_2Gev_2250_Theta[0][2])) + fgPar_2Gev_2250_Theta[0][3]);
    }
    if(cz>0.885 && cz<0.935){
      p=p*(fgPar_2Gev_2250_Theta[1][0] + fgPar_2Gev_2250_Theta[1][1]*cz + fgPar_2Gev_2250_Theta[1][2]*cz*cz);
    }
    if(cz>0.935 && cz<0.97){
      p=p*(fgPar_2Gev_2250_Theta[2][0] + fgPar_2Gev_2250_Theta[2][1]*cz + fgPar_2Gev_2250_Theta[2][2]*cz*cz);
    }
  }

 if(p!=0.) V3el.SetMag(p);
  V4ecor.SetVectMag(V3el,e_mass);
  return V4ecor;
}



bool Phot_fid(TVector3 V3_phot){

  bool status=true;
  //  double costheta=V3_phot.Pz();
  double costheta=V3_phot.CosTheta();
  double theta_deg=TMath::ACos(V3_phot.Pz())*TMath::RadToDeg();
  double phi_deg=TMath::ATan2(V3_phot.Py(),V3_phot.Px())*TMath::RadToDeg()+30;
  if(phi_deg<0)phi_deg=phi_deg+360;
  bool hot_spot=(phi_deg>185 && phi_deg<191 && costheta<0.71 && costheta>0.67) || (phi_deg>221 && phi_deg<236 && costheta<0.73 && costheta>0.67); // used to kill the two hot spots in sector 4


  if((costheta>low_lim1_ec->Eval(phi_deg) && costheta<up_lim1_ec->Eval(phi_deg) && costheta<rightside_lim1_ec->Eval(phi_deg) && costheta<leftside_lim1_ec->Eval(phi_deg))  ||
     (costheta>low_lim2_ec->Eval(phi_deg) && costheta<up_lim2_ec->Eval(phi_deg) && costheta<rightside_lim2_ec->Eval(phi_deg) && costheta<leftside_lim2_ec->Eval(phi_deg))       ||
     (costheta>low_lim3_ec->Eval(phi_deg) && costheta<up_lim3_ec->Eval(phi_deg) && costheta<rightside_lim3_ec->Eval(phi_deg) && costheta<leftside_lim3_ec->Eval(phi_deg))       ||
     (costheta>low_lim4_ec->Eval(phi_deg) && costheta<up_lim4_ec->Eval(phi_deg) && costheta<rightside_lim4_ec->Eval(phi_deg) && costheta<leftside_lim4_ec->Eval(phi_deg) && !hot_spot)       ||
     (costheta>low_lim5_ec->Eval(phi_deg) && costheta<up_lim5_ec->Eval(phi_deg) && costheta<rightside_lim5_ec->Eval(phi_deg) && costheta<leftside_lim5_ec->Eval(phi_deg))       ||
     (costheta>low_lim6_ec->Eval(phi_deg) && costheta<up_lim6_ec->Eval(phi_deg) && costheta<rightside_lim6_ec->Eval(phi_deg) && costheta<leftside_lim6_ec->Eval(phi_deg))
 ){ //EC only
      status=true;}
    else {
      status=false;}

     return status;
}

TVector3 FindUVW(TVector3 xyz)
{
  // get the U V W distance to EC edge for the purpose of geometry cut
  // ported from Stepan's function ec_xyz_duvw. the input is lab coordinates
  // of the EC hit.
  Float_t x = xyz.X(); Float_t y = xyz.Y(); Float_t z = xyz.Z();
  Float_t xi,yi,zi,u,v,w;
  Float_t ec_the = 0.4363323;
  Float_t ylow = -182.974; Float_t yhi = 189.956;
  Float_t tgrho=1.95325; Float_t sinrho=0.8901256; Float_t cosrho=0.455715;
  Float_t phi=xyz.Phi()*180./TMath::Pi(); if(phi<-30) phi+=360;
  Int_t ec_sect = (phi+30)/60.; if(ec_sect<0)ec_sect=0; if(ec_sect>5)ec_sect=5;
  Float_t ec_phi = ec_sect*TMath::Pi()/3.;
  xi = -x*sin(ec_phi) + y*cos(ec_phi);
  yi = x*cos(ec_the)*cos(ec_phi) + y*cos(ec_the)*sin(ec_phi) - z*sin(ec_the);
  zi = x*sin(ec_the)*cos(ec_phi) + y*sin(ec_the)*sin(ec_phi) + z*cos(ec_the);
  zi -= 510.32;
  u = (yi-ylow)/sinrho;
  v = (yhi-ylow)/tgrho - xi + (yhi-yi)/tgrho;
  w = ((yhi-ylow)/tgrho + xi + (yhi-yi)/tgrho)/2./cosrho;
  TVector3 uvw(u,v,w);
  return uvw;
}

Bool_t CutUVW(TVector3 ecxyz)
{
  // Cut the edges of EC according to UVW distance threshold defined by par_EcUVW array.
  // If it passes the cut, return true, if not return false
  TVector3 ecuvw = FindUVW(ecxyz);
  Float_t phi=ecxyz.Phi()*180/TMath::Pi(); if(phi<-30) phi+=360;
  Int_t sector = (phi+30)/60; if(sector<0)sector=0; if(sector>5) sector=5;
  return (ecuvw.X()>par_EcUVW[sector][0] && ecuvw.Y()<par_EcUVW[sector][1] && ecuvw.Z()<par_EcUVW[sector][2]);
}

Bool_t EFiducialCut(TVector3 momentum)
{

// Electron fiducial cut, return kTRUE if pass or kFALSE if not
  Bool_t status = kTRUE;


 if(en_beam[fbeam_E]>1. &&  en_beam[fbeam_E]<2. && fTorusCurrent>740 && fTorusCurrent<1510) {

  Float_t phiMin, phiMax;
  Float_t mom = momentum.Mag();
  Float_t phi = momentum.Phi()*180./TMath::Pi();
  if(phi<-30.) phi += 360.;
  Float_t theta = momentum.Theta()*180./TMath::Pi();
  Int_t  sector = (Int_t)((phi+30.)/60.);
  if(sector < 0) sector = 0;
  if(sector > 5) sector = 5; // to match array index


    phi -= sector*60;
    Double_t elmom = (momentum.Mag())*1000;
    Double_t thetapars[5]={0,0,0,0,0};

    for(Int_t mompar=0;mompar<6;mompar++) {
      for(Int_t thetapar=0;thetapar<5;thetapar++) {
	if((fTorusCurrent>1490) && (fTorusCurrent<1510)) {
	  // 1500A torus current
	  thetapars[thetapar]+=fgPar_1gev_1500_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
	}
	if((fTorusCurrent>740) && (fTorusCurrent<760)) {
	  // 750A torus current
	  thetapars[thetapar]+=fgPar_1gev_750_Efid[sector][thetapar][mompar]*pow(elmom,mompar);
	}
      }
    }

    Int_t uplow;
    Double_t thetacutoff;
    Float_t p_thetae=mom, thetamax_e=0;
    if (p_thetae>1.05)  p_thetae=1.05;
    else if(p_thetae<0.4)   p_thetae=0.4;
    for(int i=4;i>=0;i--)thetamax_e=thetamax_e*p_thetae+el_thetamax1[i];

    if(phi<=0) {
      uplow=1;
      thetacutoff=((phi*(thetapars[0]-(thetapars[1]/thetapars[2])))+
		   (double(uplow)*thetapars[2]*thetapars[0]))/(phi+(double(uplow)*thetapars[2]));
    }
    else {
      uplow=-1;
      thetacutoff=( (phi*(thetapars[0]-(thetapars[3]/thetapars[4]))) +
		   (double(uplow)*thetapars[4]*thetapars[0]))/(phi+(double(uplow)*thetapars[4]) );
    }

    status = (theta>thetacutoff) && (thetacutoff>=thetapars[0]) && (elmom>300) && (elmom<=1100)  && theta<=thetamax_e;



		bool SCpdcut = true;
		if (SCpdcut && (fTorusCurrent>1490) && (fTorusCurrent<1510) ){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
		  if (status){
		    int tsector = sector + 1;
		    // sector 3 has two bad paddles
		    if (tsector == 3){
		      float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
		      for (int i=0; i<4; i++){
			badpar3[i] = 0;
			// calculate the parameters using pol7
			for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom + fgPar_1gev_1500_Efid_Theta_S3[i][d];}
		      }
		      for(int ipar=0;ipar<2;ipar++)
			status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
		    }
		    // sector 4 has one bad paddle
		    else if (tsector == 4){
		      float badpar4[2];     // 2 parameters to determine the position of the theta gap
		      for (int i=0; i<2; i++){
			badpar4[i] = 0;
			// calculate the parameters using pol7
			for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom + fgPar_1gev_1500_Efid_Theta_S4[i][d];}
		      }
		      status = !(theta>badpar4[0] && theta<badpar4[1]);
		    }
		    // sector 5 has four bad paddles
		    else if (tsector == 5){
		      Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
		      for (Int_t i=0; i<8; i++){
			badpar5[i] = 0;
			// calculate the parameters using pol7
			for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom + fgPar_1gev_1500_Efid_Theta_S5[i][d];}
		      }
		      if (mom<1.25) badpar5[0] = 23.4*1500/2250;
		      if (mom<1.27) badpar5[1] = 24.0*1500/2250; // some dummy constants. see fiducial cuts webpage.

		      for(Int_t ip=0;ip<4;ip++)status = status && !(theta>badpar5[2*ip] && theta<badpar5[2*ip+1]);
		    }
		  }
		}


	if (SCpdcut && (fTorusCurrent>740) && (fTorusCurrent<760) ){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.

	  if (status){

	    int tsector = sector + 1;
	    mom = momentum.Mag();

 //sector 2 has one gap
	    if(tsector == 2){
	      double parsec2_l,parsec2_h;
	      if(mom<0.4)mom=0.4;
		parsec2_l= fid_1gev_750_efid_S2[0][0]+fid_1gev_750_efid_S2[0][1]/mom +fid_1gev_750_efid_S2[0][2]/(mom*mom) +fid_1gev_750_efid_S2[0][3]/(mom*mom*mom);
		parsec2_h= fid_1gev_750_efid_S2[1][0]+fid_1gev_750_efid_S2[1][1]/mom +fid_1gev_750_efid_S2[1][2]/(mom*mom) +fid_1gev_750_efid_S2[1][3]/(mom*mom*mom);
		status=status && !(theta>parsec2_l && theta<parsec2_h);

	    }
	    //sector 3 has four gaps, the last two appear only at low momenta (p<0.3) and affect only pimi
	    if(tsector == 3){
	      double parsec3_l[4],parsec3_h[4];
	      for(int d=0;d<4;d++){
		mom = momentum.Mag();
		if((d==2 || d==3) && mom>0.3 )mom=0.3;
		else if(d<2 && mom<0.4)mom=0.4;
		parsec3_l[d]= fid_1gev_750_efid_S3[d][0][0]+fid_1gev_750_efid_S3[d][0][1]/mom +fid_1gev_750_efid_S3[d][0][2]/(mom*mom) +fid_1gev_750_efid_S3[d][0][3]/(mom*mom*mom);
		parsec3_h[d]= fid_1gev_750_efid_S3[d][1][0]+fid_1gev_750_efid_S3[d][1][1]/mom +fid_1gev_750_efid_S3[d][1][2]/(mom*mom) +fid_1gev_750_efid_S3[d][1][3]/(mom*mom*mom);
		status=status && !(theta>parsec3_l[d] && theta<parsec3_h[d]);
	      }
	    }
	    //sector 4 has two gaps , second gap appears only at p<0.3 and theta>105 and affects only pimi
	    else if(tsector == 4){
	      double parsec4_l[2],parsec4_h[2];
	      for(int d=0;d<2;d++){
		mom = momentum.Mag();
		if(d==0 && mom<0.775 )mom=0.775;
		else if(d==1 && mom>0.3) mom=0.3;
		parsec4_l[d]= fid_1gev_750_efid_S4[d][0][0]+fid_1gev_750_efid_S4[d][0][1]/mom +fid_1gev_750_efid_S4[d][0][2]/(mom*mom) +fid_1gev_750_efid_S4[d][0][3]/(mom*mom*mom);
		parsec4_h[d]= fid_1gev_750_efid_S4[d][1][0]+fid_1gev_750_efid_S4[d][1][1]/mom +fid_1gev_750_efid_S4[d][1][2]/(mom*mom) +fid_1gev_750_efid_S4[d][1][3]/(mom*mom*mom);
		status=status && !(theta>parsec4_l[d] && theta<parsec4_h[d]);
	      }


	    }
	    //sector 5 has three gaps,
	    else if(tsector == 5){
	      double parsec5_l[3],parsec5_h[3];
	      for(int d=0;d<3;d++){

		mom = momentum.Mag();
		if(d==0 && mom>0.3)mom=0.3; //first one shows up only for pimi at p<0.3 and theta~128
		else if(d>0 && mom<0.5)mom=0.5;
		parsec5_l[d]= fid_1gev_750_efid_S5[d][0][0]+fid_1gev_750_efid_S5[d][0][1]/mom +fid_1gev_750_efid_S5[d][0][2]/(mom*mom) +fid_1gev_750_efid_S5[d][0][3]/(mom*mom*mom);
		parsec5_h[d]= fid_1gev_750_efid_S5[d][1][0]+fid_1gev_750_efid_S5[d][1][1]/mom +fid_1gev_750_efid_S5[d][1][2]/(mom*mom) +fid_1gev_750_efid_S5[d][1][3]/(mom*mom*mom);
		status=status && !(theta>parsec5_l[d] && theta<parsec5_h[d]);
	      }
	    }



	  }
	}


    return status;
  }

  if ( en_beam[fbeam_E]>2. &&  en_beam[fbeam_E]<3. && fTorusCurrent>2240 && fTorusCurrent<2260){
    Float_t phi=momentum.Phi()*180./TMath::Pi();
    if(phi<-30.) phi+=360.;
    Int_t sector = (Int_t)((phi+30.)/60.);
    if(sector<0)sector=0;
    if(sector>5) sector=5;
    phi -= sector*60;
    Float_t theta = momentum.Theta()*180./TMath::Pi();
    Float_t mom = momentum.Mag();
    Float_t par[6];               // six parameters to determine the outline of Theta vs Phi
    for (Int_t i=0; i<6; i++){
      par[i] = 0;
      for (Int_t d=8; d>=0; d--){
	par[i] = par[i]*mom + fgPar_2GeV_2250_Efid[sector][i][d];
      }                          // calculate the parameters using pol8
    }
    if (phi < 0) {
      Float_t tmptheta = par[0] - par[3]/par[2] + par[3]/(par[2]+phi);
      status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);
    }
    else {
      Float_t tmptheta = par[0] - par[5]/par[4] + par[5]/(par[4]-phi);
      status = (theta>tmptheta && tmptheta>=par[0] && theta<par[1]);
    }


    // by now, we have checked if the electron is within the outline of theta vs phi plot
    if (SCpdcut){  // if the kESCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
      if (status){
	Int_t tsector = sector + 1;
	if (tsector == 3){               // sector 3 has two bad paddles
	  Float_t badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	  for (Int_t i=0; i<4; i++){
	    badpar3[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar3[i] = badpar3[i]*mom + fgPar_2GeV_2250_EfidTheta_S3[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  for(Int_t ipar=0;ipar<2;ipar++)
	    status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
	else if (tsector == 4){         // sector 4 has one bad paddle
	  Float_t badpar4[2];           // 2 parameters to determine the position of the theta gap
	  for (Int_t i=0; i<2; i++){
	    badpar4[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar4[i] = badpar4[i]*mom + fgPar_2GeV_2250_EfidTheta_S4[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  status = !(theta>badpar4[0] && theta<badpar4[1]);
	}
	else if (tsector == 5){         // sector 5 has four bad paddles
	  Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    for (Int_t d=7; d>=0; d--){
	      badpar5[i] = badpar5[i]*mom + fgPar_2GeV_2250_EfidTheta_S5[i][d];
	    }                           // calculate the parameters using pol7
	  }
	  if (mom<1.25) badpar5[0] = 23.4;
	  if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	  for(Int_t ipar=0;ipar<4;ipar++)
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }
  }


  if ( en_beam[fbeam_E]>4. &&  en_beam[fbeam_E]<5. && fTorusCurrent>2240 && fTorusCurrent<2260){


//Begin_Html
/*</pre>
  Electron fiducial cut, return kTRUE if the electron is in the fiducial volume
  modified 14 May 2001 lbw
  Now calls GetEPhiLimits for 2.2 and 4.4 GeV
  tested against EFiducialCut for both 2.2 (with and without bad scintillator cuts) and 4.4 GeV
  discrepancy less than 2 in 10^6 events
  Please refer to <A HREF="http://www.jlab.org/Hall-B/secure/e2/bzh/efiducialcut.html">Electron Fiducial Cuts</A> -- Bin Zhang (MIT).
  For 4.4GeV please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu (UNH)
  Please refer to <a href="http://www.jlab.org/Hall-B/secure/e2/stevenmc/FiducialCuts/index.html">1.1 GeV fiducial cuts</a> -- Steven McLauchlan (GU).
<pre>
*/
//End_Html

  Float_t phiMin, phiMax;
  Float_t mom = momentum.Mag();
  Float_t phi = momentum.Phi()*180./TMath::Pi();
  if(phi<-30.) phi += 360.;
  Float_t theta = momentum.Theta()*180./TMath::Pi();
  Int_t  sector = (Int_t)((phi+30.)/60.);
  if(sector < 0) sector = 0;
  if(sector > 5) sector = 5; // to match array index
  // all the work is now done in GetEPhiLimits

  status = GetEPhiLimits(mom, theta, sector, &phiMin, &phiMax);

  if (status) {
    status = status && (phi > phiMin) && (phi < phiMax);
  }


  if(mom <= 2.0)
    {
      bool SCpdcut = true;
      if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
	if (status){
	  int tsector = sector + 1;
	  // sector 3 has two bad paddles
	  if (tsector == 3){
	    float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	    for (int i=0; i<4; i++){
	      badpar3[i] = 0;
	      // calculate the parameters using pol7
	      for (int d=7; d>=0; d--){badpar3[i] = badpar3[i]*mom + fgPar_Efid_Theta_S3[i][d];}
	    }
	    for(int ipar=0;ipar<2;ipar++)
	      status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	  }
	  // sector 4 has one bad paddle
	  else if (tsector == 4){
	    float badpar4[2];     // 2 parameters to determine the position of the theta gap
	    for (int i=0; i<2; i++){
	      badpar4[i] = 0;
	      // calculate the parameters using pol7
	      for (int d=7; d>=0; d--){badpar4[i] = badpar4[i]*mom + fgPar_Efid_Theta_S4[i][d];}
	    }
	    status = !(theta>badpar4[0] && theta<badpar4[1]);
	  }
	  // sector 5 has four bad paddles
	  else if (tsector == 5){
	    Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	    for (Int_t i=0; i<8; i++){
	      badpar5[i] = 0;
	      // calculate the parameters using pol7
	      for (Int_t d=7; d>=0; d--){badpar5[i] = badpar5[i]*mom + fgPar_Efid_Theta_S5[i][d];}
	    }
	    if (mom<1.25) badpar5[0] = 23.4;
	    if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	    for(Int_t ipar=0;ipar<4;ipar++)
	      status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	  }
	}
      }
      return (status && (phi < phiMax) && (phi>phiMin));
    }
  else{
    bool SCpdcut = true;
    if (SCpdcut){  // if the SCpdCut bit is set, take off the bad SC paddle by strictly cutting off a theta gap.
      if (status){
	int tsector = sector + 1;
	// sector 3 has two bad paddles
	if (tsector == 3){
	  float badpar3[4];            // 4 parameters to determine the positions of the two theta gaps
	  for (int i=0; i<4; i++){
	    badpar3[i] = 0;
	    // calculate the parameters using 1/p
	    badpar3[i] = fgPar_Efid_Theta_S3_extra[i][0] + fgPar_Efid_Theta_S3_extra[i][1]/mom + fgPar_Efid_Theta_S3_extra[i][2]/(mom*mom) + fgPar_Efid_Theta_S3_extra[i][3]/(mom*mom*mom);
	  }
	  for(int ipar=0;ipar<2;ipar++)
	    status = status && !(theta>badpar3[2*ipar] && theta<badpar3[2*ipar+1]);
	}
	// sector 4 has one bad paddle
	else if (tsector == 4){
	  float badpar4[2];     // 2 parameters to determine the position of the theta gap
	  for (int i=0; i<2; i++){
	    badpar4[i] = 0;
	    // calculate the parameters using 1/p
	    badpar4[i] = fgPar_Efid_Theta_S4_extra[i][0] + fgPar_Efid_Theta_S4_extra[i][1]/mom + fgPar_Efid_Theta_S4_extra[i][2]/(mom*mom) + fgPar_Efid_Theta_S4_extra[i][3]/(mom*mom*mom);
	  }
	  status = !(theta>badpar4[0] && theta<badpar4[1]);
	}
	// sector 5 has four bad paddles
	else if (tsector == 5){
	  Float_t badpar5[8];           // 8 parameters to determine the positions of the four theta gaps
	  for (Int_t i=0; i<8; i++){
	    badpar5[i] = 0;
	    // calculate the parameters using 1/p
	    badpar5[i] = fgPar_Efid_Theta_S5_extra[i][0] + fgPar_Efid_Theta_S5_extra[i][1]/mom + fgPar_Efid_Theta_S5_extra[i][2]/(mom*mom) + fgPar_Efid_Theta_S5_extra[i][3]/(mom*mom*mom);
	  }
	  if (mom<1.25) badpar5[0] = 23.4;
	  if (mom<1.27) badpar5[1] = 24.0; // some dummy constants. see fiducial cuts webpage.
	  for(Int_t ipar=0;ipar<4;ipar++)
	    status = status && !(theta>badpar5[2*ipar] && theta<badpar5[2*ipar+1]);
	}
      }
    }
    return (status && (phi < phiMax) && (phi>phiMin));
  }





  }

  return status;
}
//Histogram filling functions
void HistFill(double *q2, double *omega, double *Wvar, double *en_recon1, double *weight, double *p_perp)
{
  h1_Q2_sub->Fill(q2,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  h1_omega_sub->Fill(omega,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  h1_Wvar_sub->Fill(Wvar,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  h2_Wvar_Q2_sub->Fill(Wvar, q2,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  h2_omega_Q2_sub->Fill(omega,q2,(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  h2_cal_Wvar->Fill(en_recon1[0], Wvar, (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  h1_p_perp_sub->Fill(p_perp[0],(N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  if(p_perp[0] > 0 && p_perp[0] < 0.2) h1_cal_p_slice1_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  if(p_perp[0] > 0.2 && p_perp[0] < 0.4) h1_cal_p_slice2_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
  if(p_perp[0] > 0.4) h1_cal_p_slice3_sub->Fill(en_recon1[0], (N_1pion_1prot[0]/N_1pion_2prot)*(N_1pi_2p[count]/N_1pi_3p)*(1/Mott_cross_sec));
}

//Rotation functions start here
void rot_2pi_1p (TVector3 V3_pi[2], double qpi[2], TVector3 V3_prot, TVector3 V3q, double N1pi1p[2], double *N2pi1p, int N_tot)
{
  double N_all = 0;
  N1pi1p[0] = 0;
  N1pi1p[1] = 0;
  TVector3 V3_2pi_rot[2];
  TVector3 V3_p_rot;
  bool pi2_stat[2]={0};
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle;
  for(int g=0; g<N_tot; g++)
    {

    rot_angle=gRandom->Uniform(0,2*TMath::Pi());

        V3_2pi_rot[0]=V3_pi[0];
        V3_2pi_rot[1]=V3_pi[1];
        V3_p_rot=V3_prot;
        V3_2pi_rot[0].Rotate(rot_angle,V3q);
        V3_2pi_rot[1].Rotate(rot_angle,V3q);
        V3_p_rot.Rotate(rot_angle,V3q);

        for(int z=0;z<2;z++)
          {
            if(qpi[z]>0) pi2_stat[z]=PiplFiducialCut(V3_2pi_rot[z], &cphil, &cphir);
            else  pi2_stat[z]=PimiFiducialCut(V3_2pi_rot[z], &pimi_phimin, &pimi_phimax);
          }
        if(PFiducialCut(V3_p_rot) &&  pi2_stat[0] && !pi2_stat[1])  N1pi1p[0]=N1pi1p[0]+1;
        if(PFiducialCut(V3_p_rot) && !pi2_stat[0] &&  pi2_stat[1])  N1pi1p[1]=N1pi1p[1]+1;
        if(PFiducialCut(V3_p_rot) &&  pi2_stat[0] &&  pi2_stat[1])  N_all=N_all+1;

    }
  if(N_all!=0)  *N2pi1p = N_all;
  else *N2pi1p=0;
}

void rot_1pi_2p (TVector3 V3_pi, double qpi, TVector3 V3_prot[2], TVector3 V3q, double N1pi1p[2], double *N1pi2p, int N_tot)
{
  double N_all = 0;
  N1pi1p[0] = 0;
  N1pi1p[1] = 0;
  TVector3 V3_pi_rot;
  TVector3 V3_p_rot[2];
  bool pi_stat = 0;
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;

  for(int g=0; g<N_tot; g++)
    {
    rot_angle=gRandom->Uniform(0,2*TMath::Pi());
    V3_pi_rot = V3_pi;
    V3_p_rot[0] = V3_prot[0];
    V3_p_rot[1] = V3_prot[1];
    V3_pi_rot.Rotate(rot_angle,V3q);
    V3_p_rot[0].Rotate(rot_angle,V3q);
    V3_p_rot[1].Rotate(rot_angle,V3q);

    if(qpi>0) pi_stat=PiplFiducialCut(V3_pi_rot, &cphil, &cphir);
    else pi_stat=PimiFiducialCut(V3_pi_rot, &pimi_phimin, &pimi_phimax);

    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  pi_stat)  N_all=N_all+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) &&  pi_stat)  N1pi1p[0]=N1pi1p[0]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) &&  pi_stat)  N1pi1p[1]=N1pi1p[1]+1;
  }
if(N_all!=0)  *N1pi2p = N_all;
else *N1pi2p=0;
}
void rot_1pi_3p(TVector3 V3_pi, double qpi, TVector3 V3_prot[3], TVector3 V3q, double N1pi1p[3], double N1pi2p[3], double *N1pi3p, int N_tot)
{
  double N_all = 0;
  TVector3 V3_pi_rot;
  TVector3 V3_p_rot[3];
  bool pi_stat = 0;
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;
  for(int g=0; g<N_tot; g++)
  {
    rot_angle = gRandom->Uniform(0,2*TMath::Pi());
    V3_pi_rot = V3_pi;
    V3_pi_rot.Rotate(rot_angle, V3q);
    V3_p_rot[0] = V3_prot[0];
    V3_p_rot[0].Rotate(rot_angle, V3q);
    V3_p_rot[1] = V3_prot[1];
    V3_p_rot[1].Rotate(rot_angle, V3q);
    V3_p_rot[2] = V3_prot[2];
    V3_p_rot[2].Rotate(rot_angle, V3q);
      if(qpi>0) pi_stat = PiplFiducialCut(V3_pi_rot, &cphil, &cphir);
      else pi_stat = PimiFiducialCut(V3_pi_rot, &pimi_phimin, &pimi_phimax);

      if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat) N_all = N_all+1;
      if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat) N1pi1p[0]=N1pi1p[0]+1;
      if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat) N1pi1p[1]=N1pi1p[1]+1;
      if(!PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat) N1pi1p[2]=N1pi1p[2]+1;
      if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat) N1pi2p[0]=N1pi2p[0]+1;
      if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat) N1pi2p[1]=N1pi2p[1]+1;
      if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat) N1pi2p[2]=N1pi2p[2]+1;
  }
  *N1pi3p = N_all;
}

void rot_3pi_1p(TVector3 V3_pi[3], double qpi[3], TVector3 V3_prot, TVector3 V3q, double N1pi1p[3], double N2pi1p[3], double *N3pi1p, int N_tot)
{
  double N_all = 0;
  TVector3 V3_pi_rot[3];
  TVector3 V3_p_rot;
  bool pi_stat[3] = {0};
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;
  for(int g=0; g<N_tot; g++)
  {
    rot_angle = gRandom->Uniform(0,2*TMath::Pi());
    V3_p_rot = V3_prot;
    V3_p_rot.Rotate(rot_angle, V3q);
    V3_pi_rot[0] = V3_pi[0];
    V3_pi_rot[0].Rotate(rot_angle, V3q);
    V3_pi_rot[1] = V3_pi[1];
    V3_pi_rot[1].Rotate(rot_angle, V3q);
    V3_pi_rot[2] = V3_pi[2];
    V3_pi_rot[2].Rotate(rot_angle, V3q);
    for(int z=0; z<3; z++){
      if(qpi[z]>0) pi_stat[z] = PiplFiducialCut(V3_pi_rot[z], &cphil, &cphir);
      else pi_stat[z] = PimiFiducialCut(V3_pi_rot[z], &pimi_phimin, &pimi_phimax);
    }


    if(PFiducialCut(V3_p_rot) && pi_stat[0] && pi_stat[1] && pi_stat[2]) N_all = N_all+1;
    if(PFiducialCut(V3_p_rot) && pi_stat[0] && !pi_stat[1] && !pi_stat[2]) N1pi1p[0]=N1pi1p[0]+1;
    if(PFiducialCut(V3_p_rot) && !pi_stat[0] && pi_stat[1] && !pi_stat[2]) N1pi1p[1]=N1pi1p[1]+1;
    if(PFiducialCut(V3_p_rot) && !pi_stat[0] && !pi_stat[1] && pi_stat[2]) N1pi1p[2]=N1pi1p[2]+1;
    if(PFiducialCut(V3_p_rot) && pi_stat[0] && pi_stat[1] && !pi_stat[2]) N2pi1p[0]=N2pi1p[0]+1;
    if(PFiducialCut(V3_p_rot) && pi_stat[0] && !pi_stat[1] && pi_stat[2]) N2pi1p[1]=N2pi1p[1]+1;
    if(PFiducialCut(V3_p_rot) && !pi_stat[0] && pi_stat[1] && pi_stat[2]) N2pi1p[2]=N2pi1p[2]+1;
  }
  *N3pi1p = N_all;
}

void rot_1phot_1pi_1p(TVector3 V3_phot, TVector3 V3_pi, double qpi, TVector3 V3_prot, TVector3 V3q, bool radstat, double *N1pi1p0phot, double *N1pi1p1phot, int N_tot)
{
  double N_all = 0;
  double N_1pi_1p_0phot = 0;
  TVector3 V3_pi_rot;
  TVector3 V3_p_rot;
  TVector3 V3_phot_rot;
  bool pi_stat = true;
  bool phot_stat = true;
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;

  for(int g=0; g<N_tot; g++)
  {
    rot_angle = gRandom->Uniform(0,2*TMath::Pi());
    V3_p_rot = V3_prot;
    V3_p_rot.Rotate(rot_angle, V3q);
    V3_pi_rot = V3_pi;
    V3_pi_rot.Rotate(rot_angle, V3q);

    if(qpi>0) pi_stat = PiplFiducialCut(V3_pi_rot, &cphil, &cphir);
    else pi_stat = PimiFiducialCut(V3_pi_rot, &pimi_phimin, &pimi_phimax);

    if(!radstat)
    {
      V3_phot_rot = V3_phot;
      V3_phot_rot.Rotate(rot_angle, V3q);
      phot_stat = Phot_fid(V3_phot_rot);
    }

    if(PFiducialCut(V3_p_rot) && pi_stat && phot_stat) N_all = N_all+1;
    if(PFiducialCut(V3_p_rot) && pi_stat && !phot_stat) N_1pi_1p_0phot = N_1pi_1p_0phot+1;
  }
  *N1pi1p0phot = N_1pi_1p_0phot;
  *N1pi1p1phot = N_all;
}

void rot_1phot_1pi_2p(TVector3 V3_phot, TVector3 V3_pi, double qpi, TVector3 V3_prot[2], TVector3 V3q, bool radstat, double N1pi1p0phot[2], double N1pi1p1phot[2], double *N1pi2p0phot, double *N1pi2p1phot, int N_tot)
{
  double N_all = 0;
  double N_1pi_2p_0phot = 0;
  TVector3 V3_pi_rot;
  TVector3 V3_p_rot[2];
  TVector3 V3_phot_rot;
  bool pi_stat = true;
  bool phot_stat = true;
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;

  for(int g=0; g<N_tot; g++)
  {
    rot_angle = gRandom->Uniform(0,2*TMath::Pi());
    V3_p_rot[0] = V3_prot[0];
    V3_p_rot[0].Rotate(rot_angle, V3q);
    V3_p_rot[1] = V3_prot[1];
    V3_p_rot[1].Rotate(rot_angle, V3q);
    V3_pi_rot = V3_pi;
    V3_pi_rot.Rotate(rot_angle, V3q);

    if(qpi>0) pi_stat = PiplFiducialCut(V3_pi_rot, &cphil, &cphir);
    else pi_stat = PimiFiducialCut(V3_pi_rot, &pimi_phimin, &pimi_phimax);

    if(!radstat)
    {
      V3_phot_rot = V3_phot;
      V3_phot_rot.Rotate(rot_angle, V3q);
      phot_stat = Phot_fid(V3_phot_rot);
    }

    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && pi_stat && phot_stat) N_all = N_all+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && pi_stat && !phot_stat) N1pi1p0phot[0] = N1pi1p0phot[0]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && pi_stat && !phot_stat) N1pi1p0phot[1] = N1pi1p0phot[1]+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && pi_stat && phot_stat) N1pi1p1phot[0] = N1pi1p1phot[0]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && pi_stat && phot_stat) N1pi1p1phot[1] = N1pi1p1phot[1]+1;
    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && pi_stat && !phot_stat) N_1pi_2p_0phot = N_1pi_2p_0phot+1;
  }
  *N1pi2p1phot = N_all;
  *N1pi2p0phot = N_1pi_2p_0phot;
}

void rot_2phot_1pi_1p(TVector3 V3_phot[2], TVector3 V3_pi, double qpi, TVector3 V3_prot, TVector3 V3q, bool radstat[2], double *N1pi1p0phot, double *N1pi1p1phot, double *N1pi1p2phot, int N_tot)
{
  double N_all = 0;
  double N_1pi_1p_0phot = 0;
  double N_1pi_1p_1phot = 0;
  TVector3 V3_pi_rot;
  TVector3 V3_p_rot;
  TVector3 V3_phot_rot[2];
  bool pi_stat = true;
  bool phot_stat[2] = {true};
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;

  for(int g=0; g<N_tot; g++)
  {
    rot_angle = gRandom->Uniform(0,2*TMath::Pi());
    V3_p_rot = V3_prot;
    V3_p_rot.Rotate(rot_angle, V3q);
    V3_pi_rot = V3_pi;
    V3_pi_rot.Rotate(rot_angle, V3q);

    if(qpi>0) pi_stat = PiplFiducialCut(V3_pi_rot, &cphil, &cphir);
    else pi_stat = PimiFiducialCut(V3_pi_rot, &pimi_phimin, &pimi_phimax);

    for(int i=0; i<2; i++)
    {
      if(!radstat[i])
      {
        V3_phot_rot[i] = V3_phot[i];
        V3_phot_rot[i].Rotate(rot_angle, V3q);
        phot_stat[i] = Phot_fid(V3_phot_rot[i]);
      }
    }
    if(PFiducialCut(V3_p_rot) && pi_stat && phot_stat[0] && phot_stat[1]) N_all = N_all+1;
    if(PFiducialCut(V3_p_rot) && pi_stat && !phot_stat[0] && !phot_stat[1]) N_1pi_1p_0phot = N_1pi_1p_0phot+1;
    if(PFiducialCut(V3_p_rot) && pi_stat && phot_stat[0] && !phot_stat[1]) N_1pi_1p_1phot = N_1pi_1p_1phot+1;
    if(PFiducialCut(V3_p_rot) && pi_stat && !phot_stat[0] && phot_stat[1]) N_1pi_1p_1phot = N_1pi_1p_1phot+1;
  }
  *N1pi1p0phot = N_1pi_1p_0phot;
  *N1pi1p1phot = N_1pi_1p_1phot;
  *N1pi1p2phot = N_all;
}

void rot_1phot_2pi_1p(TVector3 V3_phot, TVector3 V3_pi[2], double qpi[2], TVector3 V3_prot, TVector3 V3q, bool radstat, double N1pi1p0phot[2], double N1pi1p1phot[2], double *N2pi1p0phot, double *N2pi1p1phot, int N_tot)
{
  double N_all = 0;
  double N_2pi_1p_0phot = 0;
  TVector3 V3_pi_rot[2];
  TVector3 V3_p_rot;
  TVector3 V3_phot_rot;
  bool pi_stat[2] = {true};
  bool phot_stat = true;
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;

  for(int g=0; g<N_tot; g++)
  {
    rot_angle = gRandom->Uniform(0,2*TMath::Pi());
    V3_pi_rot[0] = V3_pi[0];
    V3_pi_rot[0].Rotate(rot_angle, V3q);
    V3_pi_rot[1] = V3_pi[1];
    V3_pi_rot[1].Rotate(rot_angle, V3q);
    V3_p_rot = V3_prot;
    V3_p_rot.Rotate(rot_angle, V3q);

    for(int i=0; i<2; i++)
    {
      if(qpi[i]>0) pi_stat[i] = PiplFiducialCut(V3_pi_rot[i], &cphil, &cphir);
      else pi_stat[i] = PimiFiducialCut(V3_pi_rot[i], &pimi_phimin, &pimi_phimax);
    }

    if(!radstat)
    {
      V3_phot_rot = V3_phot;
      V3_phot_rot.Rotate(rot_angle, V3q);
      phot_stat = Phot_fid(V3_phot_rot);
    }

    if(PFiducialCut(V3_p_rot) && pi_stat[0] && pi_stat[1] && phot_stat) N_all = N_all+1;
    if(PFiducialCut(V3_p_rot) && pi_stat[0] && !pi_stat[1] && phot_stat) N1pi1p1phot[0] = N1pi1p1phot[0]+1;
    if(PFiducialCut(V3_p_rot) && !pi_stat[0] && pi_stat[1] && phot_stat) N1pi1p1phot[1] = N1pi1p1phot[1]+1;
    if(PFiducialCut(V3_p_rot) && pi_stat[0] && pi_stat[1] && !phot_stat) N_2pi_1p_0phot = N_2pi_1p_0phot+1;
    if(PFiducialCut(V3_p_rot) && pi_stat[0] && !pi_stat[1] && !phot_stat) N1pi1p0phot[0] = N1pi1p0phot[0]+1;
    if(PFiducialCut(V3_p_rot) && pi_stat[0] && !pi_stat[1] && !phot_stat) N1pi1p0phot[1] = N1pi1p0phot[1]+1;
  }
  *N2pi1p1phot = N_all;
  *N2pi1p0phot = N_2pi_1p_0phot;
}

void rot_1phot_1pi_3p(TVector3 V3_phot, TVector3 V3_pi, double qpi, TVector3 V3_prot[3], TVector3 V3q, bool radstat, double N1pi1p0phot[3], double N1pi1p1phot[3], double N1pi2p0phot[3], double N1pi2p1phot[3], double *N1pi3p1phot, double *N1pi3p0phot, int N_tot)
{
  double N_all = 0;
  double N_1pi_3p_0phot = 0;
  TVector3 V3_pi_rot;
  TVector3 V3_p_rot[3];
  TVector3 V3_phot_rot;
  bool pi_stat = true;
  bool phot_stat = true;
  Float_t pimi_phimin, pimi_phimax, cphil, cphir;
  double rot_angle = 0;

  for(int g=0; g<N_tot; g++)
  {
    rot_angle = gRandom->Uniform(0,2*TMath::Pi());
    V3_p_rot[0] = V3_prot[0];
    V3_p_rot[0].Rotate(rot_angle, V3q);
    V3_p_rot[1] = V3_prot[1];
    V3_p_rot[1].Rotate(rot_angle, V3q);
    V3_p_rot[2] = V3_prot[2];
    V3_p_rot[2].Rotate(rot_angle, V3q);
    V3_pi_rot = V3_pi;
    V3_pi_rot.Rotate(rot_angle, V3q);

    if(qpi>0) pi_stat = PiplFiducialCut(V3_pi_rot, &cphil, &cphir);
    else pi_stat = PimiFiducialCut(V3_pi_rot, &pimi_phimin, &pimi_phimax);

    if(!radstat)
    {
      V3_phot_rot = V3_phot;
      V3_phot_rot.Rotate(rot_angle, V3q);
      phot_stat = Phot_fid(V3_phot_rot);
    }

    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && phot_stat) N_all = N_all+1;
    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat && phot_stat) N1pi2p1phot[0] = N1pi2p1phot[0]+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && phot_stat) N1pi2p1phot[1] = N1pi2p1phot[1]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && phot_stat) N1pi2p1phot[2] = N1pi2p1phot[2]+1;
    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && !phot_stat) N_1pi_3p_0phot = N_1pi_3p_0phot+1;
    if(PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat && !phot_stat) N1pi2p0phot[0] = N1pi2p0phot[0]+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && !phot_stat) N1pi2p0phot[1] = N1pi2p0phot[1]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && !phot_stat) N1pi2p0phot[2] = N1pi2p0phot[2]+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat && phot_stat) N1pi1p1phot[0] = N1pi1p1phot[0]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat && phot_stat) N1pi1p1phot[1] = N1pi1p1phot[1]+1;
    if(!PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && phot_stat) N1pi1p1phot[2] = N1pi1p1phot[2]+1;
    if(PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat && !phot_stat) N1pi1p0phot[0] = N1pi1p0phot[0]+1;
    if(!PFiducialCut(V3_p_rot[0]) && PFiducialCut(V3_p_rot[1]) && !PFiducialCut(V3_p_rot[2]) && pi_stat && !phot_stat) N1pi1p0phot[1] = N1pi1p0phot[1]+1;
    if(!PFiducialCut(V3_p_rot[0]) && !PFiducialCut(V3_p_rot[1]) && PFiducialCut(V3_p_rot[2]) && pi_stat && !phot_stat) N1pi1p0phot[2] = N1pi1p0phot[2]+1;
  }
  *N1pi3p1phot = N_all;
  *N1pi3p0phot = N_1pi_3p_0phot;
}

double vz_corr(double phi,double theta)            //correction function for vertex , takes the arguments in Deg.
{
  //  return (0.2)*cos((phi+47.9)*TMath::DegToRad())/tan(theta*TMath::DegToRad()); // vertex correction function obtained for the empty runs 18393 and 18394, works fine for 3He runs at 2.261[GeV] beam energy
  return (-(vz_corr_func->GetParameter(1)))*cos((phi-(vz_corr_func->GetParameter(2)))*TMath::DegToRad())/tan(theta*TMath::DegToRad()); //vertex correction function for 4He runs at 2.261[GeV] beam energy obtained for the empty run18283

}

void SetFiducialCutParameters(){
// reads from a file the parameters of the fiducial cut functions
// Please refer to <A HREF="http://einstein.unh.edu/protopop/FiducialCuts/fc4E2.html">Fiducial Cuts</A> -- D.Protopopescu(UNH)



 if(en_beam[fbeam_E]>4. && en_beam[fbeam_E]<5.){    //
   // reads FC parameters for 4.4GeV , e- and p fiducial cut parameters at 4GeV
   //

  std::ifstream param_file2(Form("./PFID_%s_%d.dat",fbeam_E.c_str(),fTorusCurrent));//reading the proton fiducial cut parameters at 4GeV
  std::ifstream param_file(Form("./FCP_%s_%d.dat",fbeam_E.c_str(),fTorusCurrent));


   //	std::ifstream param_file("./FCP_4461_2250.dat");
   int param_type, sector;
   double data[6];
   while ( (sector!=6 || param_type!=21))
     {
       param_file >> param_type;
       param_file >> sector >> data[0] >> data[1] >> data[2] >> data[3] >> data[4] >> data[5];
       // Test the type of parameter and assign it to the proper data array
       //  std::cout << param_type << " " << sector << std::endl;
       switch (param_type)
	 {
	 case  0:
	   for(int k=0; k<2; k++) fgPar_4Gev_2250_Efid_t0_p[sector-1][k] = data[k];
	   break;
	 case  1:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_t1_p[sector-1][k] = data[k];
	   break;
	 case 10:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[sector-1][0][k] = data[k];
	   break;
	 case 11:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_b_p[sector-1][1][k] = data[k];
	   break;
	 case 20:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[sector-1][0][k] = data[k];
	   break;
	 case 21:
	   for(int k=0; k<6; k++) fgPar_4Gev_2250_Efid_a_p[sector-1][1][k] = data[k];
	   break;
	 default:
	   printf("Error in Efid parameter file!\nReceived parameter type %d, which is not found.\nAborting!\n\n\n",param_type);
	   break;
	 }
     } // Done reading in Fiducial Region Parameters

// ---
 for(int i = 0 ; i < 4 ; i++){
   for(int j = 0 ; j < 8 ; j++){
     param_file >> fgPar_Efid_Theta_S3[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 2 ; i++){
   for(int j = 0 ; j < 8 ; j++){
     param_file >> fgPar_Efid_Theta_S4[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 8 ; i++){
   for(int j = 0 ; j < 8 ; j++){
     param_file >> fgPar_Efid_Theta_S5[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 4 ; i++){
   for(int j = 0 ; j < 4 ; j++){
     param_file >> fgPar_Efid_Theta_S3_extra[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 2 ; i++){
   for(int j = 0 ; j < 4 ; j++){
     param_file >> fgPar_Efid_Theta_S4_extra[i][j];
   }
 }
 // ---
 for(int i = 0 ; i < 8 ; i++){
   for(int j = 0 ; j < 4 ; j++){
     param_file >> fgPar_Efid_Theta_S5_extra[i][j];
   }
 }
	param_file.close();



   for(int i = 0 ; i < 6 ; i++){
     for(int j = 0 ; j < 6 ; j++){
       param_file2 >> fgPar_4Gev_2250_Pfidft1l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidft1r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidft2l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidft2r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt1l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt1r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt2l[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbt2r[i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbl  [i][j];
       param_file2 >> fgPar_4Gev_2250_Pfidbr  [i][j];
     }
   }

for(int i = 0 ; i < 2 ; i++){//reading the proton bad TOF cuts at 4GeV
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS2[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS3[i][j];
  }
 }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS4[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 6 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS5[i][j];
  }
 }
for(int i = 0 ; i < 2 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS2_extra[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS3_extra[i][j];
  }
 }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS4_extra[i][j];
  }
 }
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 4 ; j++){
    param_file2 >> fgPar_Pfid_ScpdS5_extra[i][j];
  }
 }



 }




 else if(en_beam[fbeam_E]>1. && en_beam[fbeam_E]<2.){

  std::ifstream param_file2(Form("./PFID_%s_%d.dat",fbeam_E.c_str(),fTorusCurrent));//reading the proton fiducial cut parameters at 4GeV
  std::ifstream param_file(Form("./FCP_%s_%d.dat",fbeam_E.c_str(),fTorusCurrent));
  std::ifstream param_file3(Form("./PIPFID_%s_%d.dat",fbeam_E.c_str(),fTorusCurrent));
  std::ifstream param_file4(Form("./PIMFID_%s_%d.dat",fbeam_E.c_str(),fTorusCurrent));
  //
   // reads FC parameters for 1.1GeV , e- fiducial cut parameters at 1GeV
   //

 if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {

for(Int_t sector=0;sector<6;sector++)
  {
    for(Int_t thetapar=0;thetapar<5;thetapar++)
      {
        for(Int_t mompar=0;mompar<6;mompar++)
          {
            param_file >> fgPar_1gev_1500_Efid[sector][thetapar][mompar];
          }
      }
  }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_1500_Efid_Theta_S3[i][j];
  }
 }
// ---
for(int i = 0 ; i < 2 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_1500_Efid_Theta_S4[i][j];
  }
 }
// ---
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_1500_Efid_Theta_S5[i][j];
  }
 }

    }


 if ( fTorusCurrent < 760 && fTorusCurrent > 740){

for(Int_t sector=0;sector<6;sector++)
  {
    for(Int_t thetapar=0;thetapar<5;thetapar++)
      {
        for(Int_t mompar=0;mompar<6;mompar++)
          {
            param_file >> fgPar_1gev_750_Efid[sector][thetapar][mompar];
          }
      }
  }
for(int i = 0 ; i < 4 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_750_Efid_Theta_S3[i][j];
  }
 }
// ---
for(int i = 0 ; i < 2 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_750_Efid_Theta_S4[i][j];
  }
 }
// ---
for(int i = 0 ; i < 8 ; i++){
  for(int j = 0 ; j < 8 ; j++){
    param_file >> fgPar_1gev_750_Efid_Theta_S5[i][j];
  }
 }
 }

	param_file.close();



  //
   // reads FC parameters for 1.1GeV , p fiducial cut parameters at 1GeV
   //

 if ( fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file2 >> fgPar_1gev_1500_Pfid[sector][phipar][mompar];
                  //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
                  //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
                }
            }
        }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_1500_Pfid_ScpdS5[i][j];
        }
      }
    }
  if (fTorusCurrent< 760 && fTorusCurrent > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file2 >> fgPar_1gev_750_Pfid[sector][phipar][mompar];
                  //std::cout << "PFID " << fgPar_1gev_Pfid[sector][phipar][mompar] << std::endl;
                  //std::cout << "EFID " << fgPar_1gev_Efid[sector][phipar][mompar] << std::endl;
                }
            }
        }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file2 >> fgPar_1gev_750_Pfid_ScpdS5[i][j];
        }
      }
    }


  param_file2.close();



  //reads pimi fiducial cut parameters at 1GeV




  if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t thetapar=0;thetapar<5;thetapar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file4 >> fgPar_1gev_1500_Pimfid[sector][thetapar][mompar];
                }
            }
        }
      // ---
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S3[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S4[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S5[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4>> fgPar_1gev_1500_Pimfid_Theta_S3_extra[i][j];
        }
      }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S4_extra[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_1500_Pimfid_Theta_S5_extra[i][j];
        }
      }
    }
  if (fTorusCurrent< 760 && fTorusCurrent > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t thetapar=0;thetapar<5;thetapar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file4 >> fgPar_1gev_750_Pimfid[sector][thetapar][mompar];
                }
            }
        }
      // ---
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S3[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S4[i][j];
        }
      }
      // ---
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 8 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S5[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S3_extra[i][j];
        }
      }
      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S4_extra[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 4 ; j++){
          param_file4 >> fgPar_1gev_750_Pimfid_Theta_S5_extra[i][j];
	  //	  std::cout << fgPar_1gev_750_Pimfid_Theta_S5_extra[i][j] << std::endl;
	}
      }
    }
	param_file4.close();


	//reads fiducial cut parameters for pi+

 if (fTorusCurrent< 1510 && fTorusCurrent > 1490)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file3 >> fgPar_1gev_1500_Piplfid[sector][phipar][mompar];
		  //  std::cout << "PFID " << fgPar_1gev_1500_Pfid[sector][phipar][mompar]  << std::endl;
                }
            }
        }


      for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_1500_Piplfid_ScpdS5[i][j];
        }
      }


    }


  if ( fTorusCurrent < 760 && fTorusCurrent > 740)
    {
      for(Int_t sector=0;sector<6;sector++)
        {
          for(Int_t phipar=0;phipar<5;phipar++)
            {
              for(Int_t mompar=0;mompar<6;mompar++)
                {
                  param_file3 >> fgPar_1gev_750_Piplfid[sector][phipar][mompar];
		  //  std::cout << "PFID " << fgPar_1gev_750_Pfid[sector][phipar][mompar]  << std::endl;
                }
            }
        }

   for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS2[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS3[i][j];
        }
      }
      for(int i = 0 ; i < 4 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS4[i][j];
        }
      }
      for(int i = 0 ; i < 8 ; i++){
        for(int j = 0 ; j < 6 ; j++){
          param_file3 >> fgPar_1gev_750_Piplfid_ScpdS5[i][j];
        }
      }




    }
  param_file3.close();



 }
else printf("There are no fiducial cut parameters to be read at %3.1f GeV!\n", en_beam[fbeam_E]);
 //	param_file2.close();



}
