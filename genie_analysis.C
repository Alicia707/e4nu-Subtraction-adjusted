#define GENIE_ANALYSIS_C

#include "genie_analysis.h"
#include "Constants.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TH1D.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <TMath.h>
#include <exception>
#include <iostream>
#include <fstream>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TH3.h>
#include <TGraph.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;

// ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

vector<double> CalculateCalKineVars(double ECal,TLorentzVector FSElectron) {

	vector<double> CalKineVars; CalKineVars.clear();

	TLorentzVector V4_beam_Cal(0,0,ECal,ECal);
	double nu_Cal = -(FSElectron-V4_beam_Cal).E();
	double Q2_Cal = -(FSElectron-V4_beam_Cal).Mag2();
	double x_bjk_Cal = Q2_Cal/(2*m_prot*nu_Cal);
	TVector3 V3_q_Cal = (V4_beam_Cal-FSElectron).Vect();
	double W_var_Cal = TMath::Sqrt((m_prot+nu_Cal)*(m_prot+nu_Cal)-V3_q_Cal*V3_q_Cal);


	CalKineVars.push_back(nu_Cal); // 0-th element: energy transfer using Ecal
	CalKineVars.push_back(Q2_Cal); // 1st element: Q2 using Ecal
	CalKineVars.push_back(x_bjk_Cal); // 2nd element: xB using Ecal
	CalKineVars.push_back(W_var_Cal); // 3rd element: invariant mass using Ecal

	return CalKineVars;

}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Loading all the constants from Constant.h (e_mass, m_prot, m_pimi, m_pipl, m_pion, m_neut = 0.939565,
// H3_bind_en, He4_bind_en, C12_bind_en, B_bind_en, He3_bind_en, D2_bind_en, Fe_bind_en, Mn_bind_en

void genie_analysis::Loop(Int_t choice) {

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	double histoweight = 1.;
	//Choice = 0 is for analysis of CLAS data while choice = 1 is for the analysis of GENIE Simulation
	if (choice != 1 && choice != 0)
	{
		std::cout << "This parameter value is not implemented in genie_analysis::Loop(). It should be either 0 or 1. The given value is " << choice << std::endl;
		std::exit(0);
	}

	std::map<std::string, double> bind_en;
	std::map<std::string, double> target_mass;
	std::map<std::string, double> residual_target_mass;
	std::map<std::string, double> Ecal_offset; //that might not be necessary for simulation data

	target_name = ftarget; //std string for target name
	en_beam["1161"] = 1.161;
	en_beam["2261"] = 2.261;
	en_beam["4461"] = 4.461;

	en_beam_Ecal["1161"] = 1.161;
	en_beam_Ecal["2261"] = 2.261;
	en_beam_Ecal["4461"] = 4.461;

	en_beam_Eqe["1161"] = 1.161;
	en_beam_Eqe["2261"] = 2.261;
	en_beam_Eqe["4461"] = 4.461;

	if (fChain == 0) return;

	//8.31.22 -> Changed GetEntriesFast() to GetEntries 
	Long64_t nentries = fChain->GetEntries();
	Long64_t runEntries = ((double)per_Events/100) * nentries; 
	std::cout << "Entries: " << runEntries << std::endl;
	
	//Resolutions for Smearing for GENIE simulation data
	double reso_p = 0.01; // smearing for the proton
	double reso_e = 0.005; // smearing for the electrons
	//double reso_pipl = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
	//double reso_pimi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
	double reso_pi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)

	// Resolution defined above seems to be insufficient at 1.1 GeV -> tripled it for all particles
	if(fbeam_en == "1161") { reso_p = 3*reso_p; reso_e = 3*reso_e; reso_pi = 3*reso_pi; }

	double Wcut = 2; //cut for all beam energies < 2
	double Q2cut = 0; // cut for 1.1 GeV > 0.1, for 2.2 GeV > 0.4 and 4.4 GeV > 0.8

	const int n_slice=3; // Stick to the 3 slices
	const double pperp_min[n_slice]={0.,0.2,0.4};
	const double pperp_max[n_slice]={0.2,0.4,10.};

	//	const int n_slice=3; // Stick to the 3 slices
	//	const double pperp_min[n_slice]={0.,0.3,10.};
	//	const double pperp_max[n_slice]={0.3,10.,100.};

	TVector3 V3_rotprot1,V3_rotprot2,V3_rotprot3,V3_rot_pi,V3_rotprot;

	TString E_acc_file;

	if(en_beam[fbeam_en] > 1.0 && en_beam[fbeam_en] < 2.0) //1.1 GeV Configuration parameters and cuts
	{
		E_acc_file="1_161";
		Q2cut = 0.1;
	}


	else if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) //2.2 GeV Configuration parameters and cuts
	{
		E_acc_file="2_261";
		Q2cut = 0.4;
	}

	else if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) //4.4 GeV Configuration parameters and cuts
	{
		E_acc_file="4_461";
		Q2cut = 0.8;
	}

	/**
	 * 12.1.2022
	 * @brief If there is no beam energy within our target range, exit the program
	 */
	else
	{
		std::cout << "Please specify an energy from the following: 1161, 2261, 4461" << std::endl;
		return;
	}

	//Further constants for binding energies and target masses
	Ecal_offset["3He"]  = 0.004;
	Ecal_offset["4He"]  = 0.005;
	Ecal_offset["C12"]  = 0.005;
	Ecal_offset["56Fe"] = 0.011;

	bind_en["3He"] = He3_bind_en-D2_bind_en + Ecal_offset["3He"]; //the offset is used to shift the peak to be at 0
	bind_en["4He"] = He4_bind_en-H3_bind_en + Ecal_offset["4He"];
	bind_en["C12"] = C12_bind_en-B_bind_en	+ Ecal_offset["C12"];
	bind_en["56Fe"]= Fe_bind_en-Mn_bind_en	+ Ecal_offset["56Fe"];
	bind_en["CH2"] = C12_bind_en-B_bind_en;

	target_mass["3He"] = 2*m_prot+m_neut-He3_bind_en;
	target_mass["4He"] = 2*m_prot+2*m_neut-He4_bind_en;
	target_mass["C12"] = 6*m_prot+6*m_neut-C12_bind_en;
	target_mass["56Fe"]= 26*m_prot+30*m_neut-Fe_bind_en;
	target_mass["CH2"] = 6*m_prot+6*m_neut-C12_bind_en;

	residual_target_mass["3He"] = m_prot+m_neut-D2_bind_en;
	residual_target_mass["4He"] = m_prot+2*m_neut-H3_bind_en;
	residual_target_mass["C12"] = 5*m_prot+6*m_neut-B_bind_en;
	residual_target_mass["56Fe"]= 25*m_prot+30*m_neut-Mn_bind_en;
	residual_target_mass["CH2"] = 25*m_prot+30*m_neut-Mn_bind_en;

	//Definition of Histograms
	TH1F *h1_Etot_p_bkgd_slice[n_slice], *h1_Erec_p_bkgd_slice[n_slice];
	TH1F *h1_Etot_Npi0[n_slice],*h1_Erec_Npi0[n_slice];
	TH1F *h1_Etot_Npi1[n_slice],*h1_Erec_Npi1[n_slice];

	TH1F *h1_Erec_bkgd_pipl_pimi_new_fact[n_slice], *h1_Etot_bkgd_pipl_pimi_fact[n_slice];
	TH1F *h1_Etot_bkgd_pipl_pimi_fact_pipl[n_slice], *h1_Etot_bkgd_pipl_pimi_fact_pimi[n_slice];

	TH1F *h1_Etot_bkgd_1p2pi_pipl[n_slice], *h1_Erec_bkgd_1p2pi_pipl[n_slice];
	TH1F *h1_Etot_bkgd_1p2pi_pimi[n_slice], *h1_Erec_bkgd_1p2pi_pimi[n_slice];
	TH1F *h1_Etot_bkgd_1p3pi_pipl[n_slice], *h1_Erec_bkgd_1p3pi_pipl[n_slice];
	TH1F *h1_Etot_bkgd_1p3pi_pimi[n_slice], *h1_Erec_bkgd_1p3pi_pimi[n_slice];

	//TH1F *h1_Etot_p_bkgd_slice_2p2pi[n_slice], *h1_Erec_p_bkgd_slice_2p2pi[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_2p2pi_pimi[n_slice], *h1_Etot_p_bkgd_slice_2p2pi_pipl[n_slice];
	TH1F *h1_Erec_p_bkgd_slice_2p2pi_pimi[n_slice], *h1_Erec_p_bkgd_slice_2p2pi_pipl[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pipl[n_slice], *h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pimi[n_slice], *h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pipl[n_slice], *h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pimi[n_slice];

	TH1F *h1_Etot_piplpimi_subtruct_fact[n_slice],*h1_Erec_piplpimi_subtruct_fact[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub[n_slice],*h1_Erec_p_bkgd_slice_sub[n_slice];
	TH1F *h1_Etot_3pto1p_slice[n_slice],*h1_Erec_3pto1p_slice[n_slice];
	TH1F *h1_Etot_3pto2p_slice[n_slice],*h1_Erec_3pto2p_slice[n_slice];
	TH1F *h1_Etot_3p1pi_slice[n_slice], *h1_Erec_3p1pi_slice_pipl[n_slice], *h1_Erec_3p1pi_slice_pimi[n_slice];
	TH1F *h1_Etot_3p1pi_slice_pipl[n_slice], *h1_Etot_3p1pi_slice_pimi[n_slice];
	TH1F *h1_Etot_4pto3p_slice[n_slice],*h1_Erec_4pto3p_slice[n_slice];
	TH1F *h1_Etot_4pto1p_slice[n_slice],*h1_Erec_4pto1p_slice[n_slice];
	TH1F *h1_Etot_4pto2p_slice[n_slice],*h1_Erec_4pto2p_slice[n_slice];
	TH1F *h1_Etot_43pto1p_slice[n_slice],*h1_Erec_43pto1p_slice[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub1p2pi[n_slice],*h1_Erec_p_bkgd_slice_sub1p2pi[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_1p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_1p[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub2p1pi_2p[n_slice],*h1_Erec_p_bkgd_slice_sub2p1pi_2p[n_slice];


	TH1F *h1_Etot_p_bkgd_slice_sub32[n_slice],*h1_Erec_p_bkgd_slice_sub32[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub31[n_slice],*h1_Erec_p_bkgd_slice_sub31[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub43[n_slice],*h1_Erec_p_bkgd_slice_sub43[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub41[n_slice],*h1_Erec_p_bkgd_slice_sub41[n_slice];
	TH1F *h1_Erec_p_bkgd_slice_sub42[n_slice],*h1_Etot_p_bkgd_slice_sub42[n_slice];
	TH1F *h1_Etot_p_bkgd_slice_sub431[n_slice],*h1_Erec_p_bkgd_slice_sub431[n_slice];
	TH2F *h2_N_pi_phot[20];

	gRandom = new TRandom3();
	gRandom->SetSeed(10);

	TLorentzVector V4_beam(0,0,en_beam[fbeam_en],en_beam[fbeam_en]);
	TLorentzVector V4_target(0,0,0,target_mass[ftarget]);

	//Acceptance Maps

	TString WhichMap = "e2a_maps";
	TFile* file_acceptance;
	TFile* file_acceptance_p;
	TFile* file_acceptance_pip;

	TString Target = "12C";
	if (ftarget.c_str() == "3He") { Target = "3He"; }
	if (ftarget.c_str() == "4He") { Target = "4He"; }

	if (choice == 1)
	{ //Only need acceptance maps for GENIE simulation data
		file_acceptance = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+".root");
		file_acceptance_p = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_p.root");
		file_acceptance_pip = TFile::Open(WhichMap+"/"+WhichMap+"_"+Target+"_E_"+E_acc_file+"_pip.root");
	}

	// ---------------------------------------------------------------------------------------------------------------

	//Output file definition

	TFile *file_out;
	if (choice == 0) { file_out = new TFile(Form("data_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");}
	else { file_out = new TFile(Form("genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str()), "Recreate");}

	// ---------------------------------------------------------------------------------------------------------------

	//initialize Fiducial functions for EC limits
	fiducialcut->InitEClimits();
	std::cout << " Test InitEClimits Loop " << fiducialcut->up_lim1_ec->Eval(60) << std::endl;

	//Definition and initialization of Histograms
	TH1F *h1_el_Mott_crosssec = new TH1F("h1_el_Mott_crosssec","",200,0.,0.01);
	TH1F *h1_Wvar = new TH1F("h1_Wvar","",400,0,3);
	TH1F *h1_xbjk = new TH1F("h1_xbjk","",400,0,3);
	TH1F *h1_Q2 = new TH1F("h1_Q2","",400,0,6);
	TH1F *h1_el_theta = new TH1F("h1_el_theta","",200,0,180);
	TH1F *h1_Nprot=new TH1F("h1_Nprot","",10,-0.5,4.5);
	TH1F *h1_Nprot_NonZeroProt=new TH1F("h1_Nprot_NonZeroProt","",8,0.5,4.5);
	TH1F *h1_Nphot=new TH1F("h1_Nphot","",10,-0.5,4.5);
	TH1F *h1_Npiphot=new TH1F("h1_Npiphot","",10,-0.5,4.5);
	TH1F *h1_Npiphot_norad=new TH1F("h1_Npiphot_norad","",10,-0.5,4.5);
	TH1F *h1_Npi=new TH1F("h1_Npi","",10,-0.5,4.5);
	TH1F *h1_Npi_NonZeroProt=new TH1F("h1_Npi_NonZeroProt","",10,-0.5,4.5);
	TH1F *h1_Npipl=new TH1F("h1_Npipl","",10,-0.5,4.5);
	TH1F *h1_Npimi=new TH1F("h1_Npimi","",10,-0.5,4.5);
	TH1F *h1_MissMomentum_pipl = new TH1F("MissMomentum_pipl","",80,0.,1.);
	TH1F *h1_MissMomentum_pimi = new TH1F("MissMomentum_pimi","",80,0.,1.);
	TH1F *h1_el_mom = new TH1F("h1_el_mom","",100,0.2,6);
	TH1F *h1_el_mom_corr_pipl = new TH1F("h1_el_mom_corr_pipl","",100,0.,5.);
	TH1F *h1_el_mom_corr_pimi = new TH1F("h1_el_mom_corr_pimi","",100,0.,5.);
	TH1F *h1_el_mom_ratio = new TH1F("h1_el_mom_ratio","",50,0.97,1.01);
	TH1F *h1_prot_mom_pipl = new TH1F("h1_prot_mom_pipl","",300,0,3);
	TH1F *h1_prot_mom_pimi = new TH1F("h1_prot_mom_pimi","",300,0,3);
	TH1F *h1_prot_mom_ratio = new TH1F("h1_prot_mom_ratio","",50,0.97,1.2);
	TH1F *h1_Wvar_weight_pipl = new TH1F("h1_Wvar_weight_pipl","",400,0,3);
	TH1F *h1_Wvar_weight_pimi = new TH1F("h1_Wvar_weight_pimi","",400,0,3);
	TH1F *h1_xbjk_weight_pipl = new TH1F("h1_xbjk_weight_pipl","",400,0,3);
	TH1F *h1_xbjk_weight_pimi = new TH1F("h1_xbjk_weight_pimi","",400,0,3);
	TH1F *h1_Q2_weight_pipl = new TH1F("h1_Q2_weight_pipl","",400,0,6);
	TH1F *h1_Q2_weight_pimi = new TH1F("h1_Q2_weight_pimi","",400,0,6);
	TH1F *h1_nu_weight_pipl = new TH1F("h1_nu_weight_pipl","",400,0,4);
	TH1F *h1_nu_weight_pimi = new TH1F("h1_nu_weight_pimi","",400,0,4);
	TH1F *h1_WvarCal_weight_pipl = new TH1F("h1_WvarCal_weight_pipl","",400,0,3);
	TH1F *h1_WvarCal_weight_pimi = new TH1F("h1_WvarCal_weight_pimi","",400,0,3);
	TH1F *h1_xbjkCal_weight_pipl = new TH1F("h1_xbjkCal_weight_pipl","",400,0,3);
	TH1F *h1_xbjkCal_weight_pimi = new TH1F("h1_xbjkCal_weight_pimi","",400,0,3);
	TH1F *h1_Q2Cal_weight_pipl = new TH1F("h1_Q2Cal_weight_pipl","",400,0,6);
	TH1F *h1_Q2Cal_weight_pimi = new TH1F("h1_Q2Cal_weight_pimi","",400,0,6);
	TH1F *h1_nuCal_weight_pipl = new TH1F("h1_nuCal_weight_pipl","",400,0,3);
	TH1F *h1_nuCal_weight_pimi = new TH1F("h1_nuCal_weight_pimi","",400,0,3);

	// -------------------------------------------------------------------------------------------------------

	//Binning for energy reconstruction histograms
	int n_bins;
	double *x_values;
	double *x_qe;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.)
	{
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=17;i++) { x_values[i]=0.4+i*0.04; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=20;i++) { x_values[i+18]=1.08+(i+1)*0.02; x_qe[i+18] = (x_values[i+18] - en_beam[fbeam_en]) / en_beam[fbeam_en]; }
	}

	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.)
	{
		n_bins=54;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=23;i++) { x_values[i]=i*0.09; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=30;i++) { x_values[i+24]=2.07+(i+1)*0.03; x_qe[i+24] = (x_values[i+24] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.)
	{
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=21;i++)	{ x_values[i]=i*0.2; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=16;i++)	{ x_values[i+22]=4.2+(i+1)*0.05; x_qe[i+22] = (x_values[i+22] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	// -------------------------------------------------------------------------------------------------------

	//Definitions of further Histograms

	TH1F *CosDeltaThetaElectronPhotonAboveThreshold=new TH1F("CosDeltaThetaElectronPhotonAboveThreshold","",100,-1.,1.);
	TH1F *CosDeltaPhiElectronPhotonAboveThreshold=new TH1F("CosDeltaPhiElectronPhotonAboveThreshold","",100,-1.,1.);

	TH2F *RadCosThetaGammaEgamma = new TH2F("RadCosThetaGammaEgamma","",100,-1.,1.,600,0.,6.);
	TH2F *RadCosDeltaThetaGammaEgamma = new TH2F("RadCosDeltaThetaGammaEgamma","",100,-1.,1.,600,0.,6.);
	TH2F *NonRadThetaVsPhiGamma = new TH2F("NonRadThetaVsPhiGamma","",360,0.,360.,180,0.,180.);

	TH1F *h1_E_rec_1pi_weight_frac_feed=new TH1F("h1_E_rec_1pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_2pi_weight_frac_feed=new TH1F("h1_E_rec_2pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_3pi_weight_frac_feed=new TH1F("h1_E_rec_3pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pi_weight_frac_feed=new TH1F("h1_E_rec_4pi_weight_frac_feed","",n_bins,x_qe);
	TH1F *h1_E_tot_cut2_fracfeed = new TH1F("h1_E_tot_cut2_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_cut2_new_fracfeed = new TH1F("h1_E_rec_cut2_new_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_p_bkgd_fracfeed = new TH1F("h1_E_tot_p_bkgd_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_p_bkgd_fracfeed = new TH1F("h1_E_rec_p_bkgd_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_2p1pi_1p1pi_fracfeed_pipl = new TH1F("h1_E_tot_2p1pi_1p1pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_tot_2p1pi_1p1pi_fracfeed_pimi = new TH1F("h1_E_tot_2p1pi_1p1pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_rec_2p1pi_1p1pi_fracfeed_pipl = new TH1F("h1_E_rec_2p1pi_1p1pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_rec_2p1pi_1p1pi_fracfeed_pimi = new TH1F("h1_E_rec_2p1pi_1p1pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_tot_3pto2p_fracfeed = new TH1F("h1_E_tot_3pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_3pto2p_fracfeed = new TH1F("h1_E_rec_3pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_3pto1p_fracfeed = new TH1F("h1_E_tot_3pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_3pto1p_fracfeed = new TH1F("h1_E_rec_3pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_4pto3p_fracfeed = new TH1F("h1_E_tot_4pto3p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pto3p_fracfeed = new TH1F("h1_E_rec_4pto3p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_43pto1p_fracfeed = new TH1F("h1_E_tot_43pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_43pto1p_fracfeed = new TH1F("h1_E_rec_43pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_4pto2p_fracfeed = new TH1F("h1_E_tot_4pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pto2p_fracfeed = new TH1F("h1_E_rec_4pto2p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_4pto1p_fracfeed = new TH1F("h1_E_tot_4pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_4pto1p_fracfeed = new TH1F("h1_E_rec_4pto1p_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_1p2pi_fracfeed_pipl = new TH1F("h1_E_tot_1p2pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_tot_1p2pi_fracfeed_pimi = new TH1F("h1_E_tot_1p2pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_rec_1p2pi_fracfeed_pipl = new TH1F("h1_E_rec_1p2pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_rec_1p2pi_fracfeed_pimi = new TH1F("h1_E_rec_1p2pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_tot_1p3pi_fracfeed_pipl = new TH1F("h1_E_tot_1p3pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_tot_1p3pi_fracfeed_pimi = new TH1F("h1_E_tot_1p3pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_rec_1p3pi_fracfeed_pipl = new TH1F("h1_E_rec_1p3pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_rec_1p3pi_fracfeed_pimi = new TH1F("h1_E_rec_1p3pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_tot_2p2pi_fracfeed_pipl = new TH1F("h1_E_tot_2p2pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_tot_2p2pi_fracfeed_pimi = new TH1F("h1_E_tot_2p2pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_rec_2p2pi_fracfeed_pipl = new TH1F("h1_E_rec_2p2pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_rec_2p2pi_fracfeed_pimi = new TH1F("h1_E_rec_2p2pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_tot_3p1pi_fracfeed_pipl = new TH1F("h1_E_tot_3p1pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_tot_3p1pi_fracfeed_pimi = new TH1F("h1_E_tot_3p1pi_fracfeed_pimi","",n_bins,x_qe);
	//TH1F *h1_E_rec_3p1pi_fracfeed = new TH1F("h1_E_rec_3p1pi_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_rec_3p1pi_fracfeed_pipl = new TH1F("h1_E_rec_3p1pi_fracfeed_pipl","",n_bins,x_qe);
	TH1F *h1_E_rec_3p1pi_fracfeed_pimi = new TH1F("h1_E_rec_3p1pi_fracfeed_pimi","",n_bins,x_qe);
	TH1F *h1_E_rec_undetfactor_fracfeed = new TH1F("h1_E_rec_undetfactor_fracfeed","",n_bins,x_qe);
	TH1F *h1_E_tot_undetfactor_fracfeed = new TH1F("h1_E_tot_undetfactor_fracfeed","",n_bins,x_qe);

	TH1F *h1_theta0_pipl=new TH1F("h1_theta0_pipl","",300,0,180);
	TH1F *h1_theta0_pimi=new TH1F("h1_theta0_pimi","",300,0,180);
	TH2F *h2_Ecal_Eres_pipl=new TH2F("h2_Ecal_Eres_pipl","",800,0,8.,800,0,8.);
	TH2F *h2_Ecal_Eres_pimi=new TH2F("h2_Ecal_Eres_pimi","",800,0,8.,800,0,8.);
	TH1F *h1_Ecal_pipl=new TH1F("h1_Ecal_pipl","",n_bins, x_values);
	TH1F *h1_Ecal_pimi=new TH1F("h1_Ecal_pimi","",n_bins, x_values);
	TH1F *h1_Ecal_Reso_pipl=new TH1F("h1_Ecal_Reso_pipl","",100,-1.,1.);
	TH1F *h1_Ecal_Reso_pimi=new TH1F("h1_Ecal_Reso_pimi","",100,-1.,1.);
	TH2F *h2_Ecal_Etrue_pipl=new TH2F("h2_Ecal_Etrue_pipl","",100,0,10.,100,0,10.);
	TH2F *h2_Ecal_Etrue_pimi=new TH2F("h2_Ecal_Etrue_pimi","",100,0,10.,100,0,10.);
	TH2F *h2_Etrue_Ecal_pipl=new TH2F("h2_Etrue_Ecal_pipl","",100,0,10.,100,0,10.);
	TH2F *h2_Etrue_Ecal_pimi=new TH2F("h2_Etrue_Ecal_pimi","",100,0,10.,100,0,10.);
	TH2F *h2_EresEcalratio_Eres_pipl=new TH2F("h2_EresEcalratio_Eres_pipl","",600,0,5,300,0,2);
	TH2F *h2_EresEcalratio_Eres_pimi=new TH2F("h2_EresEcalratio_Eres_pimi","",600,0,5,300,0,2);
	TH2F *h2_EresEcaldiff_Eres_pipl=new TH2F("h2_EresEcaldiff_Eres_pipl","",600,0,5,300,-3,3);
	TH2F *h2_EresEcaldiff_Eres_pimi=new TH2F("h2_EresEcaldiff_Eres_pimi","",600,0,5,300,-3,3);
	TH2F *h2_N_prot_pi=new TH2F("h2_N_prot_pi","",10,0,5,10,0,5);
	TH2F *h2_N_prot_pi_phot=new TH2F("h2_N_prot_pi_phot","",10,0,5,10,0,5);
	TH2F *h2_N_prot_pi_phot_nonrad=new TH2F("h2_N_prot_pi_phot_nonrad","",10,0,5,10,0,5);
//	TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,0,180);
	TH2F *h2_el_theta_phi = new TH2F("h2_el_theta_phi","",200,0,360,200,10,60);
	TH2F *h2_el_mom_diff = new TH2F("h2_el_mom_diff","",500,0.,1.,500,-0.1,0.1);

	int NBinsNu = 300, NBinsQ2 = 300;
	double MinNu = 0., MaxNu = 4.; double MinQ2 = 0., MaxQ2 = 6.;
	TH2F *h2_Q2_nu = new TH2F("h2_Q2_nu","",NBinsNu,MinNu,MaxNu,NBinsQ2,MinQ2,MaxQ2);
	TH2F *h2_Q2_nu_weight_pipl = new TH2F("h2_Q2_nu_weight_pipl","",NBinsNu,MinNu,MaxNu,NBinsQ2,MinQ2,MaxQ2);
	TH2F *h2_Q2_nu_weight_pimi = new TH2F("h2_Q2_nu_weight_pimi","",NBinsNu,MinNu,MaxNu,NBinsQ2,MinQ2,MaxQ2);
	TH2F *h2_Q2_nu_weight_FirstSector_pipl = new TH2F("h2_Q2_nu_weight_FirstSector_pipl","",0.7*NBinsNu,MinNu,MaxNu,0.7*NBinsQ2,MinQ2,MaxQ2);
	TH2F *h2_Q2_nu_weight_FirstSector_pimi = new TH2F("h2_Q2_nu_weight_FirstSector_pimi","",0.7*NBinsNu,MinNu,MaxNu,0.7*NBinsQ2,MinQ2,MaxQ2);

	TH2F *h2_Q2_xbjk_weight = new TH2F("h2_Q2_xbjk_weight","",200,0,3,200,0,5);
	TH2F *h2_Q2_W=new TH2F("h2_Q2_W","",200,0,3,200,0,5);
	TH2F *h2_xB_W=new TH2F("h2_xB_W","",200,0,3,200,0,3);
	TH2F *h2_Q2_W_weight=new TH2F("h2_Q2_W_weight","",200,0,3,200,0,5);
	TH2F *h2_el_pcorr_puncorr = new TH2F("h2_el_pcorr_puncorr","",100,0,1,100,0,3);
	TH2F *h2_Erec_pperp_pipl = new TH2F("h2_Erec_pperp_pipl","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_pimi = new TH2F("h2_Erec_pperp_pimi","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_newcut2 = new TH2F("h2_Erec_pperp_newcut2","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_cut3 = new TH2F("h2_Erec_pperp_cut3","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_2p = new TH2F("h2_Erec_pperp_2p","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_321p = new TH2F("h2_Erec_pperp_321p","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_31p = new TH2F("h2_Erec_pperp_31p","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_4321p = new TH2F("h2_Erec_pperp_4321p","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_431p = new TH2F("h2_Erec_pperp_431p","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_421p = new TH2F("h2_Erec_pperp_421p","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_41p = new TH2F("h2_Erec_pperp_41p","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_1p1pi = new TH2F("h2_Erec_pperp_1p1pi","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_1p2pi_1p1pi_pipl = new TH2F("h2_Erec_pperp_1p2pi_1p1pi_pipl","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_1p2pi_1p1pi_pimi = new TH2F("h2_Erec_pperp_1p2pi_1p1pi_pimi","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_2p1pi_1p1pi_pipl = new TH2F("h2_Erec_pperp_2p1pi_1p1pi_pipl","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_2p1pi_1p1pi_pimi = new TH2F("h2_Erec_pperp_2p1pi_1p1pi_pimi","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_1p3pi_pipl = new TH2F("h2_Erec_pperp_1p3pi_pipl","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_1p3pi_pimi = new TH2F("h2_Erec_pperp_1p3pi_pimi","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_2p2pi_pipl = new TH2F("h2_Erec_pperp_2p2pi_pipl","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_2p2pi_pimi = new TH2F("h2_Erec_pperp_2p2pi_pimi","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_3p1pi_pipl = new TH2F("h2_Erec_pperp_3p1pi_pipl","",400,0,1,400,0,6.);
	TH2F *h2_Erec_pperp_3p1pi_pimi = new TH2F("h2_Erec_pperp_3p1pi_pimi","",400,0,1,400,0,6.);
	TH2F *h2_pperp_W_pipl=new TH2F("h2_pperp_W_pipl","",200,0,3,200,0,2);
	TH2F *h2_pperp_W_pimi=new TH2F("h2_pperp_W_pimi","",200,0,3,200,0,2);
	TH2F *h2_Etot_pperp_pipl = new TH2F("h2_Etot_pperp_pipl","",400,0,1,400,0,6.);
	TH2F *h2_Etot_pperp_pimi = new TH2F("h2_Etot_pperp_pimi","",400,0,1,400,0,6.);
	TH2F *h2_phot_e_angle_Erec= new TH2F("h2_phot_e_angle_Erec","",400,0,4.7,300,0,180);

	TH2F* h2_QVector_theta_phi = new TH2F("h2_QVector_theta_phi","",200,0,360,200,0,80);

	//Binning for energy reconstruction histograms
	//	int n_bins;
	//	double *x_values;
	//
	//	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
	//		n_bins=38;
	//		x_values=new double[n_bins+1];
	//		for (int i=0;i<=17;i++) x_values[i]=0.4+i*0.04;
	//		for (int i=0;i<=20;i++) x_values[i+18]=1.08+(i+1)*0.02;
	//	}
	//
	//	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
	//		n_bins=54;
	//		x_values=new double[n_bins+1];
	//		for (int i=0;i<=23;i++) x_values[i]=i*0.09;
	//		for (int i=0;i<=30;i++) x_values[i+24]=2.07+(i+1)*0.03;
	//	}
	//
	//	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
	//		n_bins=38;
	//		x_values=new double[n_bins+1];
	//		for (int i=0;i<=21;i++)	x_values[i]=i*0.2;
	//		for (int i=0;i<=16;i++)	x_values[i+22]=4.2+(i+1)*0.05;
	//	}

	//Definitions of further Histograms
	TH1F *h1_E_rec_2p_det = new TH1F("h1_E_rec_2p_det","",n_bins,x_values);
	TH1F *h1_E_tot_2p_det = new TH1F("h1_E_tot_2p_det","",n_bins,x_values);
	TH1F *h1_E_tot_p_bkgd = new TH1F("h1_E_tot_p_bkgd","",n_bins,x_values);
	TH1F *h1_E_rec_p_bkgd = new TH1F("h1_E_rec_p_bkgd","",n_bins,x_values);
	TH1F *h1_E_tot_3pto1p = new TH1F("h1_E_tot_3pto1p","",n_bins,x_values);
	TH1F *h1_E_rec_3pto1p = new TH1F("h1_E_rec_3pto1p","",n_bins,x_values);
	TH1F *h1_E_tot_43pto1p = new TH1F("h1_E_tot_43pto1p","",n_bins,x_values);
	TH1F *h1_E_rec_43pto1p = new TH1F("h1_E_rec_43pto1p","",n_bins,x_values);
	TH1F *h1_E_tot_3pto2p =new TH1F("h1_E_tot_3pto2p","",n_bins,x_values);
	TH1F *h1_E_rec_3pto2p = new TH1F("h1_E_rec_3pto2p","",n_bins,x_values);
	TH1F *h1_E_tot_4pto1p =new TH1F("h1_E_tot_4pto1p","",n_bins,x_values);
	TH1F *h1_E_rec_4pto1p = new TH1F("h1_E_rec_4pto1p","",n_bins,x_values);
	TH1F *h1_E_tot_4pto3p =new TH1F("h1_E_tot_4pto3p","",n_bins,x_values);
	TH1F *h1_E_rec_4pto3p = new TH1F("h1_E_rec_4pto3p","",n_bins,x_values);
	TH1F *h1_E_tot_4pto2p =new TH1F("h1_E_tot_4pto2p","",n_bins,x_values);
	TH1F *h1_E_rec_4pto2p = new TH1F("h1_E_rec_4pto2p","",n_bins,x_values);
	TH1F *h1_E_rec_pipl = new TH1F("h1_E_rec_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_pimi = new TH1F("h1_E_rec_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_1pi = new TH1F("h1_E_rec_1pi","",n_bins,x_values);
	TH1F *h1_E_rec_1pi_weight = new TH1F("h1_E_rec_1pi_weight","",n_bins,x_values);
	TH1F *h1_E_rec_2pi_weight = new TH1F("h1_E_rec_2pi_weight","",n_bins,x_values);
	TH1F *h1_E_rec_3pi_weight = new TH1F("h1_E_rec_3pi_weight","",n_bins,x_values);
	TH1F *h1_E_rec_4pi_weight = new TH1F("h1_E_rec_4pi_weight","",n_bins,x_values);
	TH1F *h1_E_rec_21pi = new TH1F("h1_E_rec_21pi","",n_bins,x_values);
	TH1F *h1_E_rec_1prot_pipl	= new TH1F("h1_E_rec_1prot_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_1prot_pimi	= new TH1F("h1_E_rec_1prot_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_1prot_pipl	= new TH1F("h1_E_tot_1prot_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_1prot_pimi	= new TH1F("h1_E_tot_1prot_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_cutpi1_pipl = new TH1F("h1_E_rec_cutpi1_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_cutpi1_pimi = new TH1F("h1_E_rec_cutpi1_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_cutpi1_pipl = new TH1F("h1_E_tot_cutpi1_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_cutpi1_pimi = new TH1F("h1_E_tot_cutpi1_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_pipl = new TH1F("h1_E_tot_pipl","",n_bins,x_values); //Total events
	TH1F *h1_E_tot_pimi = new TH1F("h1_E_tot_pimi","",n_bins,x_values);
//	TH1F *h1_E_tot_pimi_const_bin = new TH1F("h1_E_tot_pimi_const_bin","",54,0,8.);
//	TH1F *h1_E_tot_pipl_const_bin = new TH1F("h1_E_tot_pipl_const_bin","", 54,0,8.);
	TH1F *h1_E_rec_cut2_new = new TH1F("h1_E_rec_cut2_new","",n_bins,x_values);
	TH1F *h1_E_tot_cut2 = new TH1F("h1_E_tot_cut2","",n_bins,x_values);
	TH1F *h1_E_rec_cut005_newcut3 = new TH1F("h1_E_rec_cut005_newcut3","",n_bins,x_values);
	TH1F *h1_E_rec_undetfactor	= new TH1F("h1_E_rec_undetfactor","",n_bins,x_values);
	TH1F *h1_E_tot_undetfactor	= new TH1F("h1_E_tot_undetfactor","",n_bins,x_values);
	TH1F *h1_E_tot_1p2pi_pipl	= new TH1F("h1_E_tot_1p2pi_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_1p2pi_pimi	= new TH1F("h1_E_tot_1p2pi_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_1p2pi_pipl	= new TH1F("h1_E_rec_1p2pi_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_1p2pi_pimi	= new TH1F("h1_E_rec_1p2pi_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_1p3pi_pipl	= new TH1F("h1_E_tot_1p3pi_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_1p3pi_pimi	= new TH1F("h1_E_tot_1p3pi_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_1p3pi_pipl	= new TH1F("h1_E_rec_1p3pi_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_1p3pi_pimi	= new TH1F("h1_E_rec_1p3pi_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_2p2pi_pipl	= new TH1F("h1_E_tot_2p2pi_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_2p2pi_pimi	= new TH1F("h1_E_tot_2p2pi_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_2p2pi_pipl	= new TH1F("h1_E_rec_2p2pi_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_2p2pi_pimi	= new TH1F("h1_E_rec_2p2pi_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_3p1pi_pipl	= new TH1F("h1_E_tot_3p1pi_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_3p1pi_pimi	= new TH1F("h1_E_tot_3p1pi_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_3p1pi_pipl	= new TH1F("h1_E_rec_3p1pi_pipl","",n_bins,x_values);
  	TH1F *h1_E_rec_3p1pi_pimi	= new TH1F("h1_E_rec_3p1pi_pimi","",n_bins,x_values);
	TH1F *h1_E_tot_2p1pi_1p1pi_pipl	= new TH1F("h1_E_tot_2p1pi_1p1pi_pipl","",n_bins,x_values);
	TH1F *h1_E_tot_2p1pi_1p1pi_pimi	= new TH1F("h1_E_tot_2p1pi_1p1pi_pimi","",n_bins,x_values);
	TH1F *h1_E_rec_2p1pi_1p1pi_pipl	= new TH1F("h1_E_rec_2p1pi_1p1pi_pipl","",n_bins,x_values);
	TH1F *h1_E_rec_2p1pi_1p1pi_pimi= new TH1F("h1_E_rec_2p1pi_1p1pi_pimi","",n_bins,x_values);
	//TH1F *h1_E_tot_1p2pi_pipl_const_bin	= new TH1F("h1_E_tot_1p2pi_pipl_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_1p2pi_pimi_const_bin	= new TH1F("h1_E_tot_1p2pi_pimi_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_2p2pi_pipl_const_bin	= new TH1F("h1_E_tot_2p2pi_pipl_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_2p2pi_pimi_const_bin	= new TH1F("h1_E_tot_2p2pi_pimi_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_3p1pi_pipl_const_bin	= new TH1F("h1_E_tot_3p1pi_pipl_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_3p1pi_pimi_const_bin	= new TH1F("h1_E_tot_3p1pi_pimi_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_2p1pi_1p1pi_pipl_const_bin	= new TH1F("h1_E_tot_2p1pi_1p1pi_pipl_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_2p1pi_1p1pi_pimi_const_bin	= new TH1F("h1_E_tot_2p1pi_1p1pi_pimi_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_1p3pi_pipl_const_bin	= new TH1F("h1_E_tot_3p1pi_pipl_const_bin","",54,0,8.);
	//TH1F *h1_E_tot_1p3pi_pimi_const_bin	= new TH1F("h1_E_tot_3p1pi_pimi_const_bin","",54,0,8.);
	//TH1F *h1_E_rec_2p2pi_pipl	= new TH1F("h1_E_rec_2p2pi_pipl","",n_bins,x_values);
	//TH1F *h1_E_rec_2p2pi_pimi	= new TH1F("h1_E_rec_2p2pi_pimi","",n_bins,x_values);
	//TH1F *h1_E_rec_3p1pi_pipl	= new TH1F("h1_E_rec_3p1pi_pipl","",n_bins,x_values);
    //TH1F *h1_E_rec_3p1pi_pimi	= new TH1F("h1_E_rec_3p1pi_pimi","",n_bins,x_values);
	//TH1F *h1_E_rec_2p1pi_1p1pi_pipl	= new TH1F("h1_E_rec_2p1pi_1p1pi_pipl","",n_bins,x_values);
	//TH1F *h1_E_rec_2p1pi_1p1pi_pimi= new TH1F("h1_E_rec_2p1pi_1p1pi_pimi","",n_bins,x_values);
	
	//Plots for testing
	TH1F *h1_E_tot_1p3pi_pipl_test	= new TH1F("h1_E_tot_1p3pi_pipl_test","",n_bins,x_values);
	TH1F *h1_E_tot_1p3pi_pimi_test	= new TH1F("h1_E_tot_1p3pi_pimi_test","",n_bins,x_values);
	//Plots for E_Cal
	//11.3.21 EDIT: All events with at least 1p1pi
	TH1F *h1_E_cal_1p1pi_pimi_tot = new TH1F("h1_E_cal_1p1pi_pimi_tot","",n_bins,x_values);
	TH1F *h1_E_cal_1p1pi_pipl_tot = new TH1F("h1_E_cal_1p1pi_pipl_tot","",n_bins,x_values);
	TH1F *h1_E_cal_1p1pi_pimi_ONLY = new TH1F("h1_E_cal_1p1pi_pimi_ONLY","",n_bins,x_values);
	TH1F *h1_E_cal_1p1pi_pipl_ONLY = new TH1F("h1_E_cal_1p1pi_pipl_ONLY","",n_bins,x_values);

	//genie_truth analysis histograms and true undetected histograms
	TH1F* h1_E_cal_1p1pi_pimi_TRUE = new TH1F("h1_E_cal_1p1pi_pimi_TRUE","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_TRUE = new TH1F("h1_E_cal_1p1pi_pipl_TRUE","",n_bins, x_values);
	//True reaction Mechanisms
	TH1F* h1_E_cal_1p1pi_pimi_TRUE_qel = new TH1F("h1_E_cal_1p1pi_pimi_TRUE_qel","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_qel = new TH1F("h1_E_cal_1p1pi_pipl_TRUE_qel","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_TRUE_mec = new TH1F("h1_E_cal_1p1pi_pimi_TRUE_mec","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_mec = new TH1F("h1_E_cal_1p1pi_pipl_TRUE_mec","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_TRUE_res = new TH1F("h1_E_cal_1p1pi_pimi_TRUE_res","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_res = new TH1F("h1_E_cal_1p1pi_pipl_TRUE_res","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_TRUE_dis = new TH1F("h1_E_cal_1p1pi_pimi_TRUE_dis","",n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_TRUE_dis = new TH1F("h1_E_cal_1p1pi_pipl_TRUE_dis","",n_bins, x_values);
	//true undetected histograms
	TH1F* h1_E_cal_1p1pi_pimi_2pNotDetected 		= new TH1F("h1_E_cal_1p1pi_pimi_2pNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_2piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pimi_2piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_3piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pimi_3piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pimi_2p2piNotDetected 	= new TH1F("h1_E_cal_1p1pi_pimi_2p2piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_2pNotDetected 		= new TH1F("h1_E_cal_1p1pi_pipl_2pNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_2piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pipl_2piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_3piNotDetected 		= new TH1F("h1_E_cal_1p1pi_pipl_3piNotDetected", "", n_bins, x_values);
	TH1F* h1_E_cal_1p1pi_pipl_2p2piNotDetected 	= new TH1F("h1_E_cal_1p1pi_pipl_2p2piNotDetected", "", n_bins, x_values);
	
	//Added detected histograms
	TH1F *h1_E_cal_2p1pi_pipl_det = new TH1F("h1_E_cal_2p1pi_pipl_det","",n_bins,x_values);
	TH1F *h1_E_cal_2p1pi_pimi_det = new TH1F("h1_E_cal_2p1pi_pimi_det","",n_bins,x_values);

	// Unweighted plots for the number of events

	TH1F *h1_MissMomentum_NoWeight = new TH1F("MissMomentum_NoWeight","",100,0.,1.);

	TH1F *h1_ECal_Slice0_NoWeight = new TH1F("epRecoEnergy_slice_0_NoWeight","",n_bins,x_values);
	TH1F *h1_ECal_Slice1_NoWeight = new TH1F("epRecoEnergy_slice_1_NoWeight","",n_bins,x_values);
	TH1F *h1_ECal_Slice2_NoWeight = new TH1F("epRecoEnergy_slice_2_NoWeight","",n_bins,x_values);
	TH1F *h1_ECal_Slice3_NoWeight = new TH1F("epRecoEnergy_slice_3_NoWeight","",n_bins,x_values);

	TH1F *h1_EQE_Slice0_NoWeight = new TH1F("eRecoEnergy_slice_0_NoWeight","",n_bins,x_values);
	TH1F *h1_EQE_Slice1_NoWeight = new TH1F("eRecoEnergy_slice_1_NoWeight","",n_bins,x_values);
	TH1F *h1_EQE_Slice2_NoWeight = new TH1F("eRecoEnergy_slice_2_NoWeight","",n_bins,x_values);
	TH1F *h1_EQE_Slice3_NoWeight = new TH1F("eRecoEnergy_slice_3_NoWeight","",n_bins,x_values);

	//Histograms for testing 
	TH1F *h1_test_E_cal_2p1pi_pimi = new TH1F("h1_test_E_cal_2p1pi_pimi", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_2p1pi_pipl = new TH1F("h1_test_E_cal_2p1pi_pipl", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_1p2pi_pimi = new TH1F("h1_test_E_cal_1p2pi_pimi", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_1p2pi_pipl = new TH1F("h1_test_E_cal_1p2pi_pipl", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_2p2pi_pimi = new TH1F("h1_test_E_cal_2p2pi_pimi", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_2p2pi_pipl = new TH1F("h1_test_E_cal_2p2pi_pipl", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_3p1pi_pimi = new TH1F("h1_test_E_cal_3p1pi_pimi", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_3p1pi_pipl = new TH1F("h1_test_E_cal_3p1pi_pipl", "", n_bins, x_values);
	TH1F *h1_test_E_cal_1p3pi_pimi = new TH1F("h1_test_E_cal_1p3pi_pimi", "", n_bins, x_values); 
	TH1F *h1_test_E_cal_1p3pi_pipl = new TH1F("h1_test_E_cal_1p3pi_pipl", "", n_bins, x_values); 
	//Histograms for truth analysis 
	TH1F *h1_truth_E_cal_2p1pi_pimi = new TH1F("h1_truth_E_cal_2p1pi_pimi", "", n_bins, x_values);
	TH1F *h1_truth_E_cal_2p1pi_pipl = new TH1F("h1_truth_E_cal_2p1pi_pipl", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_1p2pi_pimi = new TH1F("h1_truth_E_cal_1p2pi_pimi", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_1p2pi_pipl = new TH1F("h1_truth_E_cal_1p2pi_pipl", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_2p2pi_pimi = new TH1F("h1_truth_E_cal_2p2pi_pimi", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_2p2pi_pipl = new TH1F("h1_truth_E_cal_2p2pi_pipl", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_3p1pi_pimi = new TH1F("h1_truth_E_cal_3p1pi_pimi", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_3p1pi_pipl = new TH1F("h1_truth_E_cal_3p1pi_pipl", "", n_bins, x_values);
	TH1F *h1_truth_E_cal_1p3pi_pimi = new TH1F("h1_truth_E_cal_1p3pi_pimi", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_1p3pi_pipl = new TH1F("h1_truth_E_cal_1p3pi_pipl", "", n_bins, x_values);
	//Histogram for slice analysis
	TH1F *h1_test_E_cal_1p2pi_pipl_piplSame = new TH1F("h1_test_E_cal_1p2pi_pipl_piplSame", "", n_bins, x_values);
	TH1F *h1_test_E_cal_1p2pi_pipl_pimiSame = new TH1F("h1_test_E_cal_1p2pi_pipl_pimiSame", "", n_bins, x_values);
	TH1F *h1_test_E_cal_1p2pi_pipl_diff 	= new TH1F("h1_test_E_cal_1p2pi_pipl_diff", "", n_bins, x_values);
	TH1F *h1_test_E_cal_1p2pi_pimi_piplSame = new TH1F("h1_test_E_cal_1p2pi_pimi_piplSame", "", n_bins, x_values);
	TH1F *h1_test_E_cal_1p2pi_pimi_pimiSame = new TH1F("h1_test_E_cal_1p2pi_pimi_pimiSame", "", n_bins, x_values);
	TH1F *h1_test_E_cal_1p2pi_pimi_diff 	= new TH1F("h1_test_E_cal_1p2pi_pimi_diff", "", n_bins, x_values);
	//Histogram for truth slice analysis 
	TH1F *h1_truth_E_cal_1p2pi_pipl_piplSame = new TH1F("h1_truth_E_cal_1p2pi_pipl_piplSame", "", n_bins, x_values);
	TH1F *h1_truth_E_cal_1p2pi_pipl_pimiSame = new TH1F("h1_truth_E_cal_1p2pi_pipl_pimiSame", "", n_bins, x_values);
	TH1F *h1_truth_E_cal_1p2pi_pipl_diff 	 = new TH1F("h1_truth_E_cal_1p2pi_pipl_diff", "", n_bins, x_values); 
	TH1F *h1_truth_E_cal_1p2pi_pimi_piplSame = new TH1F("h1_truth_E_cal_1p2pi_pimi_piplSame", "", n_bins, x_values);
	TH1F *h1_truth_E_cal_1p2pi_pimi_pimiSame = new TH1F("h1_truth_E_cal_1p2pi_pimi_pimiSame", "", n_bins, x_values);
	TH1F *h1_truth_E_cal_1p2pi_pimi_diff 	 = new TH1F("h1_truth_E_cal_1p2pi_pimi_diff", "", n_bins, x_values); 



	//Defintions of Histogram for each slice
	for(int h = 0; h < n_slice; h++)
	{
		h1_Erec_p_bkgd_slice[h]= new TH1F(Form("h1_Erec_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice[h]= new TH1F(Form("h1_Etot_p_bkgd_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_3pto1p_slice[h]= new TH1F(Form("h1_Erec_3pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_3pto1p_slice[h]= new TH1F(Form("h1_Etot_3pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_3pto2p_slice[h]= new TH1F(Form("h1_Erec_3pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_3pto2p_slice[h]= new TH1F(Form("h1_Etot_3pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_3p1pi_slice_pipl[h]= new TH1F(Form("h1_Erec_3p1pi_slice_pipl_%d",h+1),"",n_bins,x_values);
		h1_Erec_3p1pi_slice_pimi[h]= new TH1F(Form("h1_Erec_3p1pi_slice_pimi_%d",h+1),"",n_bins,x_values);
		h1_Etot_3p1pi_slice_pipl[h]= new TH1F(Form("h1_Etot_3p1pi_slice_pipl_%d",h+1),"",n_bins,x_values);
		h1_Etot_3p1pi_slice_pimi[h]= new TH1F(Form("h1_Etot_3p1pi_slice_pimi_%d",h+1),"",n_bins,x_values);
		h1_Etot_43pto1p_slice[h]= new TH1F(Form("h1_Etot_43pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_43pto1p_slice[h]= new TH1F(Form("h1_Erec_43pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_4pto3p_slice[h]= new TH1F(Form("h1_Erec_4pto3p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_4pto3p_slice[h]= new TH1F(Form("h1_Etot_4pto3p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_4pto2p_slice[h]= new TH1F(Form("h1_Erec_4pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_4pto2p_slice[h]= new TH1F(Form("h1_Etot_4pto2p_slice_%d",h+1),"",n_bins,x_values);
		h1_Erec_4pto1p_slice[h]= new TH1F(Form("h1_Erec_4pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_4pto1p_slice[h]= new TH1F(Form("h1_Etot_4pto1p_slice_%d",h+1),"",n_bins,x_values);
		h1_Etot_Npi0[h] = new TH1F(Form("h1_Etot_Npi0_%d",h+1),"",n_bins,x_values);
		h1_Erec_Npi0[h] = new TH1F(Form("h1_Erec_Npi0_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_pipl_pimi_fact[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_pipl_pimi_fact_pipl[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pipl_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_pipl_pimi_fact_pimi[h]= new TH1F(Form("h1_Etot_bkgd_pipl_pimi_fact_pimi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_pipl_pimi_new_fact[h]= new TH1F(Form("h1_Erec_bkgd_pipl_pimi_new_fact_%d",h+1),"",n_bins,x_values);
		h1_Etot_Npi1[h] = new TH1F(Form("h1_Etot_Npi1_%d",h+1),"",n_bins,x_values);
		h1_Erec_Npi1[h] = new TH1F(Form("h1_Erec_Npi1_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_1p2pi_pipl[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_pipl_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_1p2pi_pimi[h] = new TH1F(Form("h1_Etot_bkgd_1p2pi_pimi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_1p2pi_pipl[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_pipl_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_1p2pi_pimi[h] = new TH1F(Form("h1_Erec_bkgd_1p2pi_pimi_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pipl[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pipl_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pimi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pimi_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pipl[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pipl_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pimi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pimi_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p2pi_pipl[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p2pi_pipl_%d",h+1),"",n_bins,x_values);
		h1_Erec_p_bkgd_slice_2p2pi_pimi[h] = new TH1F(Form("h1_Erec_p_bkgd_slice_2p2pi_pimi_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p2pi_pipl[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p2pi_pipl_%d",h+1),"",n_bins,x_values);
		h1_Etot_p_bkgd_slice_2p2pi_pimi[h] = new TH1F(Form("h1_Etot_p_bkgd_slice_2p2pi_pimi_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_1p3pi_pipl[h] = new TH1F(Form("h1_Etot_bkgd_1p3pi_pipl_%d",h+1),"",n_bins,x_values);
		h1_Etot_bkgd_1p3pi_pimi[h] = new TH1F(Form("h1_Etot_bkgd_1p3pi_pimi_%d",h+1),"",n_bins,x_values);
		//h1_Erec_bkgd_1p3pi[h] = new TH1F(Form("h1_Erec_bkgd_1p3pi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_1p3pi_pimi[h] = new TH1F(Form("h1_Erec_bkgd_1p3pi_pimi_%d",h+1),"",n_bins,x_values);
		h1_Erec_bkgd_1p3pi_pipl[h] = new TH1F(Form("h1_Erec_bkgd_1p3pi__pipl_%d",h+1),"",n_bins,x_values);
	}

	for(int h=0;h<20;h++)
	{
		h2_N_pi_phot[h]=new TH2F(Form("h2_N_pi_phot_%d",h),"",10,0,5,10,0,5);
	}

	// Plots for interaction break down for GENIE samples
	const int NInt = 6; // All Interactions = 0, QE = 1, MEC = 2, RES = 3, DIS = 4, Other = 6
	TH1D* h1_ECal_BreakDown_pipl[NInt];
	TH1D* h1_ECal_BreakDown_pimi[NInt];
	TH1D* h1_Eres_BreakDown_pipl[NInt];
	TH1D* h1_Eres_BreakDown_pimi[NInt];
	TH1D* h1_InclusiveEQE_BreakDown[NInt];
	TH1D* h1_Pmiss_BreakDown_pipl[NInt];
	TH1D* h1_Pmiss_BreakDown_pimi[NInt];
	TH1D* h1_Q2_BreakDown_pipl[NInt];
	TH1D* h1_Q2_BreakDown_pimi[NInt];
	TH1D* h1_Nu_BreakDown_pipl[NInt];
	TH1D* h1_Nu_BreakDown_pimi[NInt];
	TH1D* h1_Pe_BreakDown_pipl[NInt];
	TH1D* h1_Pe_BreakDown_pimi[NInt];

	for (int WhichInt = 0; WhichInt < NInt; WhichInt++)
	{
		h1_ECal_BreakDown_pipl[WhichInt] = new TH1D(Form("ECal_Int_%d",WhichInt),";E^{Cal} (GeV)",n_bins,x_values);
		h1_ECal_BreakDown_pimi[WhichInt] = new TH1D(Form("ECal_Int_%d",WhichInt),";E^{Cal} (GeV)",n_bins,x_values);
		h1_Eres_BreakDown_pipl[WhichInt] = new TH1D(Form("EQE_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);
		h1_Eres_BreakDown_pimi[WhichInt] = new TH1D(Form("EQE_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);
		h1_InclusiveEQE_BreakDown[WhichInt] = new TH1D(Form("InclusiveEQE_Int_%d",WhichInt),";E^{QE} (GeV)",n_bins,x_values);
		h1_Pmiss_BreakDown_pipl[WhichInt] = new TH1D(Form("Pmiss_Int_%d",WhichInt),";P_{miss}^{#perp} [GeV/c]",80,0.,1.);
		h1_Pmiss_BreakDown_pimi[WhichInt] = new TH1D(Form("Pmiss_Int_%d",WhichInt),";P_{miss}^{#perp} [GeV/c]",80,0.,1.);
		h1_Q2_BreakDown_pipl[WhichInt] = new TH1D(Form("Q2_Int_%d",WhichInt),";Q^{2} [GeV^{2}/c^{2}]",400,0,6);
		h1_Q2_BreakDown_pimi[WhichInt] = new TH1D(Form("Q2_Int_%d",WhichInt),";Q^{2} [GeV^{2}/c^{2}]",400,0,6);
		h1_Nu_BreakDown_pipl[WhichInt] = new TH1D(Form("Nu_Int_%d",WhichInt),";Energy Transfer [GeV]",400,0,4);
		h1_Nu_BreakDown_pimi[WhichInt] = new TH1D(Form("Nu_Int_%d",WhichInt),";Energy Transfer [GeV]",400,0,4);
		h1_Pe_BreakDown_pipl[WhichInt] = new TH1D(Form("Pe_Int_%d",WhichInt),";P_{e} [GeV/c]",100,0.,5.);
		h1_Pe_BreakDown_pimi[WhichInt] = new TH1D(Form("Pe_Int_%d",WhichInt),";P_{e} [GeV/c]",100,0.,5.);
	}

	// Vector containing kinematic variables using Ecal
	vector<double> CalKineVars{};
	// Weight to fill the plots mentioned above
	double LocalWeight;

	// Signal Event Counter -> 1e1p0pi events (everything lese is bkg)
	int SignalEvents = 0;
	int QESignalEvents = 0;
	int MECSignalEvents = 0;
	int RESSignalEvents = 0;
	int DISSignalEvents = 0;
	int OtherSignalEvents = 0;

	//----- Counts for events to match up -----//
	//----- 1p1pi -----//
	int C1p1piALL = 0;
	int C1p1piPIMI = 0;
	int C1p1piPIPL = 0;
	//----- 2p1pi -----//
	int C2p1piALL = 0;
	int C2p1piPIMI = 0;
	int C2p1piPIPL = 0;
	//----- 3p1pi -----//
	int C3p1piALL = 0;
	int C3p1piPIMI = 0;
	int C3p1piPIPL = 0;
	//----- 1p2pi -----//
	int C1p2piALL = 0;
	int C1p2piPIMI = 0;
	int C1p2piPIPL = 0;
	//----- 2p2pi -----//
	int C2p2piALL = 0;
	int C2p2piPIMI = 0;
	int C2p2piPIPL = 0;
	//----- 1p3pi -----//
	int C1p3piALL = 0;
	int C1p3piPIMI = 0;
	int C1p3piPIPL = 0;


	int EQESignalEventsWithin5Perc = 0, EQESignalEventsWithin5Perc_FirstSlice = 0, EQESignalEventsWithin5Perc_SecondSlice = 0, EQESignalEventsWithin5Perc_ThirdSlice = 0;
	int ECalSignalEventsWithin5Perc = 0, ECalSignalEventsWithin5Perc_FirstSlice = 0, ECalSignalEventsWithin5Perc_SecondSlice = 0, ECalSignalEventsWithin5Perc_ThirdSlice = 0;
	int PMiss_FirstSlice = 0, PMiss_SecondSlice = 0, PMiss_ThirdSlice = 0;

	//Implementation of 1p1pi true stuff//
	//TString FileName = "";
	//if(choice == 1)
	//{
	//	FileName = Form("genie_e2a_ep_%s_%s_neutrino6_united4_radphot_test.root",ftarget.c_str(),fbeam_en.c_str());
	//}
    //	TFile*FilteredFile = new TFile(FileName, "Recreate");
    //	TTree *maketreenew = fChain->GetTree()->CloneTree(0);

	// ---------------------------------------------------------------------------------------------------------------

	// Get the number of events to run overall

    //	int Nentries = TMath::Min(Ntfileentries,NtweightsEntries);

	// ---------------------------------------------------------------------------------------------------------------
	Long64_t firstEntry = 0;
	Long64_t ientry = LoadTree(firstEntry);
	if(firstEntry == 0)
	{ //first entry to initialize TorusCurrent, Fiducials and Subtraction classes
		//The TorusField has to be set before the Fiducialcut parameters are initialized
		if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2. ) //1.1 GeV, we are not using the 1.1 GeV data with 1500 current field
		{
				fTorusCurrent = 750;
		}
		else if( (en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.) || (en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.) ) //2.2 GeV	or 4.4 GeV
		{
				fTorusCurrent = 2250;
		}
		else { std::cout << "genie_analysis::Loop(): fTorusCurrent could not be assigned" << std::endl;}

		fiducialcut->SetConstants(fTorusCurrent, target_name, en_beam);
		fiducialcut->SetFiducialCutParameters(fbeam_en);
		std::cout << " EventLoop: Finished setting up fiducial cut class " << std::endl;
		rotation->InitSubtraction(fbeam_en, target_name, bind_en, N_tot, fiducialcut);
		std::cout << " EventLoop: Finished setting up rotation initialize " << std::endl;
	}


	/** Beginning of Event Loop **/
	//for (Long64_t jentry=0; jentry<nentries;jentry++)
	for(Long64_t jentry = 0; jentry < runEntries; jentry++)
	{
       //	for (Long64_t jentry=0; jentry<Nentries;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		//Read Entry
		int nb = GetEntry(jentry);
		if (nb == 0) { std::cout <<"Event loop: 0 byte read for entry " << jentry << ". Indicate failure in reading the file" <<	std::endl;}

		if (jentry%10000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/fChain->GetEntries()*100. << " %"<< std::endl;}

		if( jentry%200000 == 0 )
		{
			gDirectory->Write("hist_Files", TObject::kOverwrite);
			//cout<<jentry<<endl;
		}

		// -----------------------------------------------------------------------------------------------------------------------------------------------------------

		double SmearedPe;
		double SmearedEe;
		double e_acc_ratio = 1.;	//will be 1 for CLAS data

		// Outgoing e',	Uncorr and corrected are the same read from root file.
		//V4_el and V3_el will be changed by smearing for GENIE simulation data
		TLorentzVector V4_el(pxl,pyl,pzl,El);
		TLorentzVector V4_el_uncorr(pxl,pyl,pzl,El);
		//Comment above line Ali look here 
		TLorentzVector V4_el_unsmeared(pxl,pyl,pzl,El);
		TVector3 V3_el(pxl,pyl,pzl);

		double el_momentum = V3_el.Mag();
		double el_theta = V3_el.Theta();
		double el_theta_deg;

		if (choice == 1)
		{ //smearing, fiducials and acceptance ratio for GENIE simulation data

			//Smearing of Electron Vector from Simulation
			SmearedPe = gRandom->Gaus(pl,reso_e*pl);
			SmearedEe = sqrt( SmearedPe*SmearedPe + e_mass * e_mass );
			//Ali Look here -> Maybe just set V4 and then make V3 from there?
			V3_el.SetXYZ(SmearedPe/pl * pxl,SmearedPe/pl * pyl,SmearedPe/pl * pzl);
			double phi_ElectronOut = V3_el.Phi(); //in Radians

			V3_el.SetPhi(phi_ElectronOut + TMath::Pi() ); // Vec.Phi() is between (-180,180), GENIE coordinate system flipped with respect to CLAS
			V4_el.SetPxPyPzE(V3_el.X(),V3_el.Y(),V3_el.Z(), SmearedEe);
			V4_el.SetPhi(V3_el.Phi());
			//Fiducial Cuts with the smeared values
			if (!EFiducialCut(fbeam_en,V3_el) ) continue; // Electron theta & phi fiducial cuts

			phi_ElectronOut += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS
			el_momentum = V3_el.Mag(); //Momentum after smearing
			el_theta = V3_el.Theta(); //Angle after smearing

			//acceptance_c takes phi in radians and here unmodified by 30 degree.
			//11 - might be particle type (electron)
			e_acc_ratio = acceptance_c(el_momentum, cos(el_theta), phi_ElectronOut, 11,file_acceptance);
			if ( fabs(e_acc_ratio) != e_acc_ratio ) { continue; }

			// --------------------------------------------------------------------------------------------------

			// GENIE Systematic Uncertainties

			//			tweights->GetEntry(jentry);
			//			float* ArrayWeights = weights->GetArray();

			////			double TuningWeight = ArrayWeights[0]; // - 1 sigma variation
			////			double TuningWeight = ArrayWeights[1]; // 0 sigma variation
			//			double TuningWeight = ArrayWeights[2]; // + 1 sigma variation

			//			e_acc_ratio = e_acc_ratio * TuningWeight;

			// --------------------------------------------------------------------------------------------------

		}

		// Explicit cuts on electron momentum
		if (fbeam_en=="1161" && el_momentum < 0.4) { continue; }
		if (fbeam_en=="2261" && el_momentum < 0.55) { continue; }
		if (fbeam_en=="4461" && el_momentum < 1.1) { continue; }

		//Definition as for data. It is also correct for GENIE simulation data since V3_el is rotated above by 180 degree in phi
		double el_phi_mod = V3_el.Phi()*TMath::RadToDeg()  + 30; //Add 30 degree for plotting and photon phi cut
		if(el_phi_mod<0)  el_phi_mod  = el_phi_mod+360; //Add 360 so that electron phi is between 0 and 360 degree


		//Calculated Mott Cross Section and Weights for Inclusive Histograms
		//Wght and e_acc_ratio is 1 for CLAS data
		//double Mott_cross_sec = ( pow(fine_struc_const,2.)*(cos(el_theta)+1))/(2*pow(El,2.)*pow((1-cos(el_theta)),2.));

		double reco_Q2 = -(V4_el-V4_beam).Mag2();
		double Q4 = reco_Q2 * reco_Q2;
		double Mott_cross_sec = (1./Q4);  //not really the Mott_cross_sec


		// ---------------------------------------------------------------------------------------------------------------------

		double WeightIncl = wght*e_acc_ratio / Mott_cross_sec;//wght -> genie weight

		// Securing ourselves against infinities
		if ( fabs(WeightIncl) != WeightIncl ) { continue; }

		//Calculation of Reconstructed Energy from ELectron only
		//using the same value of single nucleon separation E Ecal and Eqe
		double m_delta = 1.232;
		double E_rec = (m_delta*m_delta-(m_prot-bind_en[ftarget])*(m_prot-bind_en[ftarget])+2*(m_prot-bind_en[ftarget])*V4_el.E())/(2*(m_prot-bind_en[ftarget]-V4_el.E()+V4_el.Rho()*cos(el_theta)));

		//Calculation of kinematic quantities (nu, Q2, x bjorken, q and W)
		double nu = -(V4_el-V4_beam).E();
		double x_bjk = reco_Q2/(2*m_prot*nu);

		// QE selection
		//if ( fabs(x_bjk - 1.) > 0.2) { continue; }

		// ---------------------------------------------------------------------------------------------------------------------

		TVector3 V3_q = (V4_beam-V4_el).Vect();
		double V3_q_theta_deg = V3_q.Theta() * 180. / TMath::Pi();
		double V3_q_phi_deg = V3_q.Phi() * 180. / TMath::Pi() + 30.;
		if (V3_q_phi_deg > 360) { V3_q_phi_deg = V3_q_phi_deg - 360.; }
		if (V3_q_phi_deg < 0) { V3_q_phi_deg = V3_q_phi_deg + 360.; }
		double W_var = TMath::Sqrt((m_prot+nu)*(m_prot+nu)-V3_q*V3_q);

		//converting theta to degrees
		el_theta_deg = el_theta*TMath::RadToDeg();

		//Cuts on Q2 and W, only keep events with Q2 > Q2cut and W < Wcut
		if ( reco_Q2 < Q2cut || W_var > Wcut) continue;

		//Set q vector for the following rotations for the subtraction procedure
		rotation->SetQVector(V3_q);
		//rotation->PrintQVector();

		h1_el_mom->Fill(V4_el_uncorr.Rho(),WeightIncl);
		h1_el_mom_ratio->Fill(V4_el.Rho()/V4_el_uncorr.Rho(),WeightIncl);
		h2_el_pcorr_puncorr->Fill(V4_el.Rho(),V4_el.Rho()/V4_el_uncorr.Rho(),WeightIncl);
		h2_el_mom_diff->Fill(V4_el.Rho(),V4_el.Rho()-V4_el_uncorr.Rho(),WeightIncl);

		//Filling Histogram for electron kinematics
		h1_xbjk->Fill(x_bjk);
		h1_Q2->Fill(reco_Q2);
		h1_Wvar->Fill(W_var);
		h1_el_Mott_crosssec->Fill(Mott_cross_sec);
		h2_el_theta_phi->Fill(el_phi_mod,el_theta_deg,WeightIncl);
		h1_el_theta->Fill(el_theta_deg);
		h2_Q2_nu->Fill(nu,reco_Q2);
		h2_Q2_xbjk_weight->Fill(x_bjk,reco_Q2,WeightIncl);
		h2_Q2_W->Fill(W_var,reco_Q2);
		h2_xB_W->Fill(W_var,x_bjk);
		h2_Q2_W_weight->Fill(W_var,reco_Q2,WeightIncl);

		//Now we are done with the selection of electrons. Next step is looking for other hadrons in the events

		//Index variables for hadrons (p and pions)
		int index_p_det[20]; //index for each proton
		int index_pi_det[20]; //index for each pion
		int ind_pi_phot[20]; //index for pions and photons
		int index_pipl[20]; //index for each pi plus
		int index_pimi[20]; //index for each pi minus

		int charge_pi_det[20]; //Charge for the pions and photons
		int charge_pi_true[20]; //Charge for pions/photons (true)
		//Smeared Momentum and Energy values for GENIE (simulation) data
		double Smeared_Pp[20]; //smeared momentum values for protons
		double Smeared_Ep[20]; //smeared energy values for protons
		double Smeared_Ppi[20]; //smeared momentum values for pions
		double Smeared_Epi[20]; //smeared energy values for pions

		// - - - Number of hadrons - - - //
		int num_p_det = 0;
		int num_pi_det = 0;
		int num_pi_phot = 0; //couting all pions and photons
		int num_pimi_det = 0;
		int num_pipl_det = 0;
		int num_pi_phot_nonrad = 0; //counting all pions and non-radiation photons
		int num_phot_rad = 0; //counting radiation photons
		//counters for genie_truth analysis
		vector<int> index_pi_true;
		vector<int> index_p_true;
		vector<int> index_pipl_true;
		vector<int> index_pimi_true;
		int num_p_true = 0;
		int num_pi_true = 0;
		int num_pipl_true = 0;
		int num_pimi_true = 0;
		int num_gamma_true = 0;
		//Index and number variables for neutral particles
		int ec_index_n[20];
		int ec_num_n = 0;
		bool ec_radstat_n[20];

		//Array initialize to -1 or false
		for (int i = 0; i < 20; i++)
		{
			index_p_det[i] = -1;   index_pi_det[i] = -1;   index_pipl[i] = -1;   index_pimi[i] = -1;   ind_pi_phot[i] = -1;
			ec_index_n[i] = -1;   ec_radstat_n[i] = false;
			charge_pi_det[i] = -2; //default number should be not a possible real charge
			charge_pi_true[i] = -2;
			Smeared_Pp[i]  = 0; Smeared_Ep[i]  = 0;  //default 0 momentum and energy after smearing
			Smeared_Ppi[i] = 0; Smeared_Epi[i] = 0;  //default 0 momentum and energy after smearing
		}

		const double phot_rad_cut = 40;
		const double phot_e_phidiffcut = 30; //electron - photon phi difference cut

		// Creating vectors to store id of particles in the array
		vector <int> ProtonID; vector <int> PiPlusID; vector <int> PiMinusID; vector <int> PhotonID;
		ProtonID.clear(); PiPlusID.clear(); PiMinusID.clear();  PhotonID.clear();

		//Loop for Hadrons
		for (int i = 0; i < nf; i++)
		{
			//Start of proton selection
			if (pdgf[i] == 2212  && pf[i] > 0.3)
			{
				if ( choice == 1)
				{ //GENIE data

					//genie_truth analysis
					num_p_true++;
					index_p_true.push_back(i);

					//Smearing of proton
					double temp_smear_P = gRandom->Gaus(pf[i],reso_p*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_prot * m_prot );

					TVector3 V3_prot_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_prot = V3_prot_corr.Phi();
					V3_prot_corr.SetPhi(phi_prot + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS
					if (!PFiducialCut(fbeam_en, V3_prot_corr) ) { continue; } // Proton theta & phi fiducial cuts

					num_p_det = num_p_det + 1;
					index_p_det[num_p_det - 1] = i;
					ProtonID.push_back(i);
					Smeared_Pp[num_p_det - 1] = temp_smear_P;
					Smeared_Ep[num_p_det - 1] = temp_smear_E;
				}
				else if(choice == 0)
				{ //CLAS data does not need Fiducial Cut again
					num_p_det = num_p_det + 1;
					index_p_det[num_p_det - 1] = i;
					ProtonID.push_back(i);
					//temp stuff bc the true will be true since its data 
					index_p_true.push_back(i); 
					num_p_true++;
				}
			}

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if (pdgf[i] == -211  && pf[i] > 0.15)
			{ //PI minus
				if ( choice == 1)
				{ //GENIE data
					//genie_truth analysis
					num_pimi_true++;
					num_pi_true++;
					index_pimi_true.push_back(i);
					index_pi_true.push_back(i);
					charge_pi_true[i] = -1;

					//Smearing of pi minus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					// Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
					if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, -1) ){continue;}
					num_pimi_det = num_pimi_det + 1;
					num_pi_det = num_pi_det + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pimi[num_pi_det - 1] = i;
					index_pi_det[num_pi_det - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiMinusID.push_back(i);
					charge_pi_det[num_pi_det - 1] = -1;
					Smeared_Ppi[num_pi_det - 1] = temp_smear_P;
					Smeared_Epi[num_pi_det - 1] = temp_smear_E;
				}
				else
				{ //CLAS data does not need Fiducial Cut again
					num_pimi_det = num_pimi_det + 1;
					num_pi_det = num_pi_det + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pimi[num_pi_det - 1] = i;
					index_pi_det[num_pi_det - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiMinusID.push_back(i);
					charge_pi_det[num_pi_det - 1] = -1;

					//temp stuff bc data true will be true 
					num_pimi_true++;
					num_pi_true++;
					index_pimi_true.push_back(i);
					index_pi_true.push_back(i);
					charge_pi_true[i] = -1;
				}
			}//end of pimi

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if ( pdgf[i] == 211  && pf[i] > 0.15)
			{
				if (choice == 1)
				{ //GENIE data
					//genie_truth analysis
					num_pipl_true++;
					num_pi_true++;
					index_pipl_true.push_back(i);
					index_pi_true.push_back(i);
					charge_pi_true[i] = 1;

					//Smearing of pi plus
					double temp_smear_P = gRandom->Gaus(pf[i],reso_pi*pf[i]);
					double temp_smear_E = sqrt( temp_smear_P*temp_smear_P + m_pion * m_pion );

					TVector3 V3_pi_corr(temp_smear_P/pf[i] * pxf[i],temp_smear_P/pf[i] * pyf[i],temp_smear_P/pf[i] * pzf[i]);
					double phi_pion = V3_pi_corr.Phi();
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					// Pi_phot_fid_united with +1 is for Piplus and Pi_phot_fid_united with -1 is for Piminus
					if ( !Pi_phot_fid_united(fbeam_en, V3_pi_corr, 1) )     {  continue; }

					num_pipl_det = num_pipl_det + 1;
					num_pi_det  = num_pi_det + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pipl[num_pi_det - 1] = i;
					index_pi_det[num_pi_det - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiPlusID.push_back(i);
					charge_pi_det[num_pi_det - 1] = 1;
					Smeared_Ppi[num_pi_det - 1] = temp_smear_P;
					Smeared_Epi[num_pi_det - 1] = temp_smear_E;
				}
				else
				{ //CLAS data does not need Fiducial Cut again
					num_pipl_det = num_pipl_det + 1;
					num_pi_det  = num_pi_det + 1;
					num_pi_phot = num_pi_phot + 1;
					num_pi_phot_nonrad = num_pi_phot_nonrad + 1;
					index_pipl[num_pi_det - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					ind_pi_phot[num_pi_phot - 1] = i;
					PiPlusID.push_back(i);
					charge_pi_det[num_pi_det - 1] = 1;

					//temp stuff bc data true will be true 
					num_pipl_true++;
					num_pi_true++;
					index_pipl_true.push_back(i);
					index_pi_true.push_back(i);
					charge_pi_true[i] = 1;
					index_pi_det[num_pi_det - 1] = i;

				}
			} //end of pipl

			// -----------------------------------------------------------------------------------------------------------------------------------------------

			if (pdgf[i] == 22  && pf[i] > 0.3)
			{
				continue; //6.29.22 No photons!
			}//end of gamma
		} //end of hadron loop

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------

		//Skip event if there is at least one radiation photon
		if (num_phot_rad > 0) {
		  continue;
		}

		//Skip event if there are more than 3 protons or more than 3 pions 
		//Also skip if true exceeds this number bc we are only trying to fix the lower order mults rn 
		/*if(num_p_det > 3 || num_pi_det > 3 || num_p_true > 3 || num_pi_true > 3 || (num_p_true == 3 && num_pi_true == 2) || (num_p_true == 2 && num_pi_true == 3)) 
		{
			continue;
		} Omit for nwo bc thesis */

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------
		// For GENIE samples, identify the interaction type

		int Interaction = -1;
		if (choice == 1) {
			if (qel) { Interaction = 1; }
			if (mec) { Interaction = 2; }
			if (res) { Interaction = 3; }
			if (dis) { Interaction = 4; }
		}

		// -------------------------------------------------------------------------------------------------------------------------------------------------------------
		
		//Filling Histograms with multiplicities
		h1_Npi->Fill(num_pi_det);
		h1_Nprot->Fill(num_p_det);

		/*
		Commenting out code to speed up testing 
		//genie_truth analysis:
		if(num_p_det == 1 && num_pi_det == 1) //detected 1p and detected 1pi
		{ //define E_cal up here
			double E_cal = Ef[index_p_det[0]] + El + Ef[index_pi_det[0]] - m_prot;
			if(num_p_true == 2 && num_pi_true == 1) //2nd proton undetected
			{
				if(num_pipl_true == 1 && num_pimi_true == 0)
				{
					h1_E_cal_1p1pi_pipl_2pNotDetected->Fill(E_cal, histoweight);
				}
				else if(num_pipl_true == 0 && num_pimi_true == 1)
				{
					h1_E_cal_1p1pi_pimi_2pNotDetected->Fill(E_cal, histoweight);
				}
				else
					cout<<"This should not happen! 2p1pi undetect with extra pion! ln 1252" << endl;
			}
			else if(num_p_true == 1 && num_pi_true == 2) //2nd pi undetected
			{
				//check if the first pion is positively charged
				if(charge_pi_det[0] > 0) //detected pion is pipl
				{
					if(num_pipl_true == 2 || num_pimi_true == 1) //extra charged pion is pipl or pimi
					{
						h1_E_cal_1p1pi_pipl_2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Line 1270" << endl;
					}
				}
				else if(charge_pi_det[0] < 0) //detected pion is pimi
				{
					if(num_pipl_true == 1 || num_pimi_true == 2)//extra charged pion is pipl or PiMinusID
					{
						h1_E_cal_1p1pi_pimi_2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Line 1282" << endl;
					}
				}
			}
			else if(num_p_true == 1 && num_pi_true == 3)
			{
				if(charge_pi_det[0] == 1)
				{
					h1_E_cal_1p1pi_pipl_3piNotDetected->Fill(E_cal, histoweight);
				}
				else if(charge_pi_det[0] == -1)
				{
					h1_E_cal_1p1pi_pimi_3piNotDetected->Fill(E_cal, histoweight);
				}
			}
			else if(num_p_true == 2 && num_pi_true == 2)
			{
				if(charge_pi_det[0] > 0 ) //detected pion is pipl
				{
					if(num_pipl_true == 2 || num_pimi_true == 1) //extra pi is either pipl or pimi
					{
						h1_E_cal_1p1pi_pipl_2p2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Ln 1302. " << endl;
					}
				}
				else if(charge_pi_det[0] < 0) //detected pi is pimi
				{
					if(num_pipl_true == 1 || num_pimi_true == 2)
					{
						h1_E_cal_1p1pi_pimi_2p2piNotDetected->Fill(E_cal, histoweight);
					}
					else
					{
						cout << "This should not happen. Ln 1313." << endl;
					}
				}
				else
				{
					cout << "Somehow the charge of the first pion is 0, this is not possible! Line 1318" << endl;
					cout << "The charge of the first pion is: " << charge_pi_det[0] << endl;
				}
			}
		}

		if(choice == 1 && num_pipl_true == 1 && num_p_true == 1)
		{
			double E_cal = Ef[index_p_true[0]] + El + Ef[index_pipl_true[0]]-m_prot;
			if(num_p_true == 1 && num_pi_true == 1 && num_pipl_true == 1)
			{
				h1_E_cal_1p1pi_pipl_TRUE->Fill(E_cal, histoweight);
			}
		}

		//genie_truth analysis
		if(choice == 1 && num_p_true == 1 && num_pimi_true == 1)
		{
			double E_cal = Ef[index_p_true[0]] + El + Ef[index_pimi_true[0]]-m_prot;
			if(num_p_true == 1 && num_pi_true == 1 && num_pimi_true == 1)
			{
				h1_E_cal_1p1pi_pimi_TRUE->Fill(E_cal, histoweight);
			}
		}

		*/

		//genie_truth test suite 
		if(choice == 1 && num_p_det == 1 && num_pi_det == 1)
		{
			//Checking for 1p1pi events that are secretly npnpi events instead 
			TLorentzVector V4_prot_uncorr(pxf[index_p_det[0]],pyf[index_p_det[0]],pzf[index_p_det[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p_det[0]]*pf[index_p_det[0]]));
			TLorentzVector V4_pi_corr(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[0]]*pf[ind_pi_phot[0]]));
			double E_cal = V4_el.E() + V4_prot_uncorr.E() - m_prot + V4_pi_corr.E();
			//2p1pi
			if(num_p_true == 2 && num_pi_true == 1)
			{
				if(charge_pi_true[index_pi_det[0]] == 1)
					h1_truth_E_cal_2p1pi_pipl->Fill(E_cal, histoweight);
				else if(charge_pi_true[index_pi_det[0] == -1])
					h1_truth_E_cal_2p1pi_pimi->Fill(E_cal, histoweight);
			}
			else if(num_p_true == 1 && num_pi_true == 2)
			{
				//ensure that we are calling the detected and undetected pions separately
				int detPionIndex = index_pi_det[0];
				int undetPionIndex = -2; 
				if(index_pi_det[0] == index_pi_true[0]) 
					undetPionIndex = index_pi_true[1]; 
				else
					undetPionIndex = index_pi_true[0];
				//std::cout << "Undetected Index, Charge: " << undetPionIndex <<  " " << charge_pi_true[undetPionIndex] << " Detected Index, Charge: " << detPionIndex  << " " << charge_pi_true[detPionIndex] << " \n";
				if(charge_pi_true[detPionIndex] == 1)
				{
					h1_truth_E_cal_1p2pi_pipl->Fill(E_cal, histoweight);

					if(charge_pi_true[undetPionIndex] == 1)
					{
						h1_truth_E_cal_1p2pi_pipl_piplSame->Fill(E_cal, histoweight);
					}
					else if(charge_pi_true[undetPionIndex] == -1)
						h1_truth_E_cal_1p2pi_pipl_diff->Fill(E_cal, histoweight);
				}
				else if(charge_pi_true[detPionIndex] == -1)
				{
					h1_truth_E_cal_1p2pi_pimi->Fill(E_cal, histoweight);
					if(charge_pi_true[undetPionIndex] == 1)
						h1_truth_E_cal_1p2pi_pimi_diff->Fill(E_cal, histoweight); 
					else if(charge_pi_true[undetPionIndex] == -1)
						h1_truth_E_cal_1p2pi_pimi_pimiSame->Fill(E_cal, histoweight); 
				}
			} //Carry through the undet/det pion index fix (ALI LOOK HERE)
			else if(num_p_true == 2 && num_pi_true == 2)
			{
				if(charge_pi_det[index_pi_det[0]] == 1)
					h1_truth_E_cal_2p2pi_pipl->Fill(E_cal, histoweight);
				else if(charge_pi_det[index_pi_det[0]] == -1)
					h1_truth_E_cal_2p2pi_pimi->Fill(E_cal, histoweight);
			}
			else if(num_p_true == 3 && num_pi_true == 1)
			{
				if(charge_pi_det[index_pi_det[0]] == 1)
					h1_truth_E_cal_3p1pi_pipl->Fill(E_cal, histoweight);
				else if(charge_pi_det[index_pi_det[0]] == -1)
					h1_truth_E_cal_3p1pi_pimi->Fill(E_cal, histoweight);
			}
			else if(num_p_true == 1 && num_pi_true == 3)
			{
				if(charge_pi_det[index_pi_det[0]] == 1)
					h1_truth_E_cal_1p3pi_pipl->Fill(E_cal, histoweight);
				else if(charge_pi_det[index_pi_det[0]] == -1)
					h1_truth_E_cal_1p3pi_pimi->Fill(E_cal, histoweight);
			}
			else if(num_p_true == 1 && num_pi_true == 1)
			{
				//do nothing
			}
			//else
				//std::cout << "Weird number of p/pi " << num_p_true << " " << num_pi_true << std::endl;
		}
		if (num_p_det > 0)
		{
			h1_Nprot_NonZeroProt->Fill(num_p_det );
			h1_Npi_NonZeroProt->Fill(num_pi_det);
			h2_QVector_theta_phi->Fill(V3_q_phi_deg,V3_q_theta_deg,WeightIncl);
		}

		h1_Nphot->Fill(ec_num_n);
		h1_Npipl->Fill(num_pipl_det);
		h1_Npimi->Fill(num_pimi_det);
		h1_Npiphot->Fill(num_pi_phot);
		h1_Npiphot_norad->Fill(num_pi_phot_nonrad);
		h2_N_prot_pi->Fill(num_pi_det, num_p_det);
		h2_N_prot_pi_phot->Fill(num_pi_det + ec_num_n,num_p_det);
		h2_N_prot_pi_phot_nonrad->Fill(num_pi_phot_nonrad,num_p_det );
		h2_N_pi_phot[num_p_det]->Fill(ec_num_n, num_pi_det);

		bool detTruthFlag = false; 

		//Events with exactly 2 protons
		if(num_p_det == 2)
		{
			//Implement unsmeared
			TLorentzVector V4_2prot_unsmeared[2]; 
			TLorentzVector V4_2prot_smeared[2];
			for(int i = 0; i < 2; i++)
			{
				V4_2prot_unsmeared[i].SetPxPyPzE(pxf[index_p_det[i]], pyf[index_p_det[i]], pzf[index_p_det[i]], TMath::Sqrt(m_prot*m_prot+pf[index_p_det[i]]*pf[index_p_det[i]]));
			}


			//LorentzVectors for protons without momentum smearing or corrections
			TLorentzVector V4_prot_uncorr1(pxf[index_p_det[0]],pyf[index_p_det[0]],pzf[index_p_det[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p_det[0]]*pf[index_p_det[0]]));
			TLorentzVector V4_prot_uncorr2(pxf[index_p_det[1]],pyf[index_p_det[1]],pzf[index_p_det[1]],TMath::Sqrt(m_prot*m_prot+pf[index_p_det[1]]*pf[index_p_det[1]]));
			//LorentzVectors for protons with momentum smearing or corrections
			TVector3 V3_prot_corr1;
			TVector3 V3_prot_corr2;

			double p_acc_ratio1 = 1; //will be 1 for CLAS data
			double p_acc_ratio2 = 1; //will be 1 for CLAS data

			if (choice == 0) { //CLAS data
				V3_prot_corr1.SetXYZ(pxf[index_p_det[0]+60],pyf[index_p_det[0]+60],pzf[index_p_det[0]+60]);
				V3_prot_corr2.SetXYZ(pxf[index_p_det[1]+60],pyf[index_p_det[1]+60],pzf[index_p_det[1]+60]);
			}

			if (choice == 1)
			{ //GENIE data, fiducials are done in hadron loop
				V3_prot_corr1.SetXYZ(Smeared_Pp[0]/pf[index_p_det[0]] * pxf[index_p_det[0]],Smeared_Pp[0]/pf[index_p_det[0]] * pyf[index_p_det[0]],Smeared_Pp[0]/pf[index_p_det[0]] * pzf[index_p_det[0]]);
				double phi_prot1 = V3_prot_corr1.Phi();
				V3_prot_corr1.SetPhi(phi_prot1 + TMath::Pi()); // Vec.Phi() is between (-180,180), // GENIE coordinate system flipped with respect to CLAS
				phi_prot1 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				double p_theta1 = V3_prot_corr1.Theta();
				double prot_mom_corr1 = V3_prot_corr1.Mag();
				//Proton 1 weight
				p_acc_ratio1 = acceptance_c(prot_mom_corr1, cos(p_theta1), phi_prot1, 2212,file_acceptance_p);
				if ( fabs(p_acc_ratio1) != p_acc_ratio1 ) { continue; }

				V3_prot_corr2.SetXYZ(Smeared_Pp[1]/pf[index_p_det[1]] * pxf[index_p_det[1]],Smeared_Pp[1]/pf[index_p_det[1]] * pyf[index_p_det[1]],Smeared_Pp[1]/pf[index_p_det[1]] * pzf[index_p_det[1]]);
				double phi_prot2 = V3_prot_corr2.Phi();
 				V3_prot_corr2.SetPhi(phi_prot2 + TMath::Pi()); // Vec.Phi() is between (-180,180) // GENIE coordinate system flipped with respect to CLAS
				phi_prot2 += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				double p_theta2 = V3_prot_corr2.Theta();
				double prot_mom_corr2 = V3_prot_corr2.Mag();
				//Proton 2 weight
				p_acc_ratio2 = acceptance_c(prot_mom_corr2, cos(p_theta2), phi_prot2, 2212,file_acceptance_p);
				if ( fabs(p_acc_ratio2) != p_acc_ratio2 ) { continue; }

			}

			//Total proton weight
			double weight_protons = p_acc_ratio1 * p_acc_ratio2;

			TVector3 V3_2prot_uncorr[2];
			V3_2prot_uncorr[0] = V4_prot_uncorr1.Vect();
			V3_2prot_uncorr[1] = V4_prot_uncorr2.Vect();
			//Look here -> Make an array to begin with (bad coding)
			TVector3 V3_2prot_corr[2];
			V3_2prot_corr[0] = V3_prot_corr1;
			V3_2prot_corr[1] = V3_prot_corr2;

			TLorentzVector V4_2prot_corr[2];
			V4_2prot_corr[0] = V4_prot_uncorr1;
			V4_2prot_corr[1] = V4_prot_uncorr2;
			//Added Corrected Px,Py,Pz values for V4_2prot since it is used for prot2_pi1 calculations 7.28.22
			V4_2prot_corr[0].SetPxPyPzE(V3_2prot_corr[0].X(), V3_2prot_corr[0].Y(), V3_2prot_corr[0].Z(), Smeared_Ep[0]);
			V4_2prot_corr[1].SetPxPyPzE(V3_2prot_corr[1].X(), V3_2prot_corr[1].Y(), V3_2prot_corr[1].Z(), Smeared_Ep[1]);
			//---------------------------------- 2p 1pi   ----------------------------------------------
			//Const int can be placed somewhere up after if for 2 protons F.H. 05.09.19
			const int N_2prot=2;
			//Variable might/could be placed in a more local context F.H. 05.09.19
			double Ecal_2p1pi_to2p0pi[N_2prot]={0};
			double p_miss_perp_2p1pi_to2p0pi[N_2prot]={0};
			
			if (num_pi_det == 1)
			{
				C2p1piALL++;

				TVector3 V3_1pi_corr;
				TLorentzVector V4_1pi_corr;
				//Unsmeared add Ali Look here 
				TLorentzVector V4_1pi_unsmeared(pxf[index_pi_det[0]], pyf[index_pi_det[0]], pzf[index_pi_det[0]], Ef[index_pi_det[0]]);
				double pion_acc_ratio = 1;

				if (choice == 0)
				{ //CLAS data
					V3_1pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);
					V4_1pi_corr.SetPxPyPzE(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[0]]*pf[ind_pi_phot[0]]));
				}

				if (choice == 1)
				{ //GENIE data
					pion_acc_ratio = 0;//reset to 0 just to be save
					V3_1pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);
					V4_1pi_corr.SetPxPyPzE(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]], TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[0]]*pf[ind_pi_phot[0]]));

					double phi_pion = V3_1pi_corr.Phi();
					V3_1pi_corr.SetPhi(phi_pion + TMath::Pi()); // Vec.Phi() is between (-180,180)
					V4_1pi_corr.SetPhi(V3_1pi_corr.Phi()); //Set new Phi value to reflect V3_1pi_corr's phi
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_1pi_corr.Theta();
					double pion_mom_corr = V3_1pi_corr.Mag();

					if (charge_pi_true[index_pi_det[0]] == 1)
					{ //acceptance for pi plus
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi_true[index_pi_det[0]] == -1)
					{    //acceptance for pi minus. using electron acceptance map
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi_true[index_pi_det[0]] == 0)
					{    //acceptance for neutral, setting to 1 for now F.H. 09/24/19
						pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 2proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }
				}

				double P_2p1pito2p0pi[2] = {0};
				double P_2p1pito1p1pi[2] = {0};
				double P_2p1pito1p0pi[2] = {0};
				double Ptot[2] = {0};
				double E_tot_2p[2] = {0};
				double p_perp_tot_2p[2] = {0};
				//rotation->prot2_pi1_rot_func(V3_2prot_corr, V3_2prot_uncorr, V3_1pi_corr, V4_2prot_corr, V4_1pi_corr, charge_pi_det[0], V4_el, E_tot_2p, p_perp_tot_2p, Ptot);
				rotation->prot2_pi1_rot_func(V4_2prot_unsmeared, V4_el_unsmeared, V4_1pi_unsmeared, V3_q, charge_pi_det[0], E_tot_2p,p_perp_tot_2p, P_2p1pito1p1pi);

				//rotation->prot2_pi1_rot_func(V3_2prot_corr, V3_2prot_uncorr, V3_1pi_corr, V4_2prot_corr, V4_1pi_corr, charge_pi_det[0], V4_el, E_tot_2p,p_perp_tot_2p, P_2p1pito1p1pi);
				//double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec;
				//Is this correct in the following loop? F.H. 09/01/19
				detTruthFlag = (((num_pi_true == 1) && (num_p_true == 2)) && ((num_pi_det == 1) && (num_p_det == 2)));
				for(int z=0; z < N_2prot; z++)
				{ //looping over two protons
					if(charge_pi_true[index_pi_true[0]] == 1)
					{
					 	C2p1piPIPL++;
						if(detTruthFlag || choice == 0)
						{
							h1_test_E_cal_2p1pi_pipl->Fill(E_tot_2p[z], P_2p1pito1p1pi[z]*histoweight);
						}
						//11.3.21 EDIT: Added histogram fill for total events with at LEAST 1p1pi
						h1_E_cal_1p1pi_pipl_tot->Fill(E_tot_2p[z], histoweight);
						//---------------------------------- 2p 1pi ->2p 0pi ----------------------------------------------
						h2_Etot_pperp_pipl->Fill(p_miss_perp_2p1pi_to2p0pi[z],Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
						h2_pperp_W_pipl->Fill(W_var,p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
						h1_theta0_pipl->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z]) *TMath::RadToDeg(),P_2p1pito2p0pi[z]*histoweight);
						h2_Ecal_Eres_pipl->Fill(E_rec,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
						h1_Ecal_pipl->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
						h1_Ecal_Reso_pipl->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
						h2_Ecal_Etrue_pipl->Fill(Ecal_2p1pi_to2p0pi[z],Ev,P_2p1pito2p0pi[z]*histoweight);
						h2_Etrue_Ecal_pipl->Fill(Ev,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
						h2_EresEcalratio_Eres_pipl->Fill(E_rec,E_rec/Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
						h2_EresEcaldiff_Eres_pipl->Fill(E_rec,E_rec-Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);

						h1_xbjk_weight_pipl->Fill(x_bjk,P_2p1pito2p0pi[z]*histoweight);
						h1_Q2_weight_pipl->Fill(reco_Q2,P_2p1pito2p0pi[z]*histoweight);
						h1_Wvar_weight_pipl->Fill(W_var,P_2p1pito2p0pi[z]*histoweight);
						h1_nu_weight_pipl->Fill(nu,P_2p1pito2p0pi[z]*histoweight);
						h1_el_mom_corr_pipl->Fill(V4_el.Rho(),P_2p1pito2p0pi[z]*histoweight);
						h1_prot_mom_pipl->Fill(V3_2prot_corr[z].Mag(),P_2p1pito2p0pi[z]*histoweight);
						h1_MissMomentum_pipl->Fill(p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(Ecal_2p1pi_to2p0pi[z],V4_el);
						LocalWeight = P_2p1pito2p0pi[z]*histoweight;

						h1_nuCal_weight_pipl->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pipl->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pipl->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pipl->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pipl->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pipl->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pipl[0]->Fill(Ecal_2p1pi_to2p0pi[z],LocalWeight);
						h1_Eres_BreakDown_pipl[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pipl[0]->Fill(p_miss_perp_2p1pi_to2p0pi[z],LocalWeight);
						h1_Q2_BreakDown_pipl[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pipl[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pipl[0]->Fill(V4_el.Rho(),LocalWeight);

	 					if (choice == 1)
						{
							h1_ECal_BreakDown_pipl[Interaction]->Fill(Ecal_2p1pi_to2p0pi[z],LocalWeight);
							h1_Eres_BreakDown_pipl[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pipl[Interaction]->Fill(p_miss_perp_2p1pi_to2p0pi[z],LocalWeight);
							h1_Q2_BreakDown_pipl[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pipl[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pipl[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						//---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

						h1_E_tot_2p1pi_1p1pi_pipl->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					    //	h1_E_tot_2p1pi_1p1pi_pipl_const_bin->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h1_E_rec_2p1pi_1p1pi_pipl->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
						h2_Erec_pperp_2p1pi_1p1pi_pipl->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p1pi[z]*histoweight);
						h2_Etot_pperp_pipl->Fill(p_perp_tot_2p[z],E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h1_E_tot_2p1pi_1p1pi_fracfeed_pipl->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
						h1_E_rec_2p1pi_1p1pi_fracfeed_pipl->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
						h2_pperp_W_pipl->Fill(W_var,p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h1_theta0_pipl->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_2p1pito1p1pi[z]*histoweight);
						h2_Ecal_Eres_pipl->Fill(E_rec,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h1_Ecal_pipl->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h1_Ecal_Reso_pipl->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
						h2_Ecal_Etrue_pipl->Fill(E_tot_2p[z],Ev,P_2p1pito1p1pi[z]*histoweight);
						h2_Etrue_Ecal_pipl->Fill(Ev,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h2_EresEcalratio_Eres_pipl->Fill(E_rec,E_rec/E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
						h2_EresEcaldiff_Eres_pipl->Fill(E_rec,E_rec-E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

						h1_xbjk_weight_pipl->Fill(x_bjk,P_2p1pito1p1pi[z]*histoweight);
						h1_Q2_weight_pipl->Fill(reco_Q2,P_2p1pito1p1pi[z]*histoweight);
						h1_Wvar_weight_pipl->Fill(W_var,P_2p1pito1p1pi[z]*histoweight);
						h1_nu_weight_pipl->Fill(nu,P_2p1pito1p1pi[z]*histoweight);
						h1_el_mom_corr_pipl->Fill(V4_el.Rho(),P_2p1pito1p1pi[z]*histoweight);
						h1_prot_mom_pipl->Fill(V3_2prot_corr[z].Mag(),P_2p1pito1p1pi[z]*histoweight);
						h1_MissMomentum_pipl->Fill(p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(E_tot_2p[z],V4_el);
						LocalWeight = P_2p1pito1p1pi[z]*histoweight;

						h1_nuCal_weight_pipl->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pipl->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pipl->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pipl->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pipl->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pipl->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pipl[0]->Fill(E_tot_2p[z],LocalWeight);
						h1_Eres_BreakDown_pipl[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pipl[0]->Fill(p_perp_tot_2p[z],LocalWeight);
						h1_Q2_BreakDown_pipl[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pipl[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pipl[0]->Fill(V4_el.Rho(),LocalWeight);

	 					if (choice == 1)
						{
							h1_ECal_BreakDown_pipl[Interaction]->Fill(E_tot_2p[z],LocalWeight);
							h1_Eres_BreakDown_pipl[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pipl[Interaction]->Fill(p_perp_tot_2p[z],LocalWeight);
							h1_Q2_BreakDown_pipl[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pipl[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pipl[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++)
						{
							if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i])
							{
								h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pipl[i]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
								h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pipl[i]->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
							}
						}

						//---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------
						h2_Etot_pperp_pipl->Fill(p_perp_tot_2p[z],E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
						h2_pperp_W_pipl->Fill(W_var,p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
						h1_theta0_pipl->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),-P_2p1pito1p0pi[z]*histoweight);
						h2_Ecal_Eres_pipl->Fill(E_rec,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
						h1_Ecal_pipl->Fill(E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
						h1_Ecal_Reso_pipl->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-P_2p1pito1p0pi[z]*histoweight);
						h2_Ecal_Etrue_pipl->Fill(E_tot_2p[z],Ev,-P_2p1pito1p0pi[z]*histoweight);
						h2_Etrue_Ecal_pipl->Fill(Ev,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
						h2_EresEcalratio_Eres_pipl->Fill(E_rec,E_rec/E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
						h2_EresEcaldiff_Eres_pipl->Fill(E_rec,E_rec-E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

						h1_xbjk_weight_pipl->Fill(x_bjk,-P_2p1pito1p0pi[z]*histoweight);
						h1_Q2_weight_pipl->Fill(reco_Q2,-P_2p1pito1p0pi[z]*histoweight);
						h1_Wvar_weight_pipl->Fill(W_var,-P_2p1pito1p0pi[z]*histoweight);
						h1_nu_weight_pipl->Fill(nu,-P_2p1pito1p0pi[z]*histoweight);
						h1_el_mom_corr_pipl->Fill(V4_el.Rho(),-P_2p1pito1p0pi[z]*histoweight);
						h1_prot_mom_pipl->Fill(V3_2prot_corr[z].Mag(),-P_2p1pito1p0pi[z]*histoweight);
						h1_MissMomentum_pipl->Fill(p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(E_tot_2p[z],V4_el);
						LocalWeight = -P_2p1pito1p0pi[z]*histoweight;

						h1_nuCal_weight_pipl->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pipl->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pipl->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pipl->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pipl->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pipl->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pipl[0]->Fill(E_tot_2p[z],LocalWeight);
						h1_Eres_BreakDown_pipl[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pipl[0]->Fill(p_perp_tot_2p[z],LocalWeight);
						h1_Q2_BreakDown_pipl[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pipl[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pipl[0]->Fill(V4_el.Rho(),LocalWeight);

	 					if (choice == 1)
						{
							h1_ECal_BreakDown_pipl[Interaction]->Fill(E_tot_2p[z],LocalWeight);
							h1_Eres_BreakDown_pipl[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pipl[Interaction]->Fill(p_perp_tot_2p[z],LocalWeight);
							h1_Q2_BreakDown_pipl[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pipl[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pipl[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}
					}//Minor bug fix: charge_pi_det for pimi
				else if(charge_pi_true[index_pi_det[0]] == -1)
				{
				  	C2p1piPIMI++;
					if(detTruthFlag || choice == 0)
					{
						h1_test_E_cal_2p1pi_pimi->Fill(E_tot_2p[z], P_2p1pito1p1pi[z]*histoweight);
					}
					//11.3.21 EDIT: Added 1p1pi total event counter.
					h1_E_cal_1p1pi_pimi_tot->Fill(E_tot_2p[z],histoweight);

					h2_Etot_pperp_pimi->Fill(p_miss_perp_2p1pi_to2p0pi[z],Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h2_pperp_W_pimi->Fill(W_var,p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_theta0_pimi->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z]) *TMath::RadToDeg(),P_2p1pito2p0pi[z]*histoweight);
					h2_Ecal_Eres_pimi->Fill(E_rec,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_Ecal_pimi->Fill(Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h1_Ecal_Reso_pimi->Fill((Ecal_2p1pi_to2p0pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito2p0pi[z]*histoweight);
					h2_Ecal_Etrue_pimi->Fill(Ecal_2p1pi_to2p0pi[z],Ev,P_2p1pito2p0pi[z]*histoweight);
					h2_Etrue_Ecal_pimi->Fill(Ev,Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h2_EresEcalratio_Eres_pimi->Fill(E_rec,E_rec/Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					h2_EresEcaldiff_Eres_pimi->Fill(E_rec,E_rec-Ecal_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);

					h1_xbjk_weight_pimi->Fill(x_bjk,P_2p1pito2p0pi[z]*histoweight);
					h1_Q2_weight_pimi->Fill(reco_Q2,P_2p1pito2p0pi[z]*histoweight);
					h1_Wvar_weight_pimi->Fill(W_var,P_2p1pito2p0pi[z]*histoweight);
					h1_nu_weight_pimi->Fill(nu,P_2p1pito2p0pi[z]*histoweight);
					h1_el_mom_corr_pimi->Fill(V4_el.Rho(),P_2p1pito2p0pi[z]*histoweight);
					h1_prot_mom_pimi->Fill(V3_2prot_corr[z].Mag(),P_2p1pito2p0pi[z]*histoweight);
					h1_MissMomentum_pimi->Fill(p_miss_perp_2p1pi_to2p0pi[z],P_2p1pito2p0pi[z]*histoweight);
					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(Ecal_2p1pi_to2p0pi[z],V4_el);
					LocalWeight = P_2p1pito2p0pi[z]*histoweight;

					h1_nuCal_weight_pimi->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight_pimi->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight_pimi->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight_pimi->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight_pimi->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pimi->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					h1_ECal_BreakDown_pimi[0]->Fill(Ecal_2p1pi_to2p0pi[z],LocalWeight);
					h1_Eres_BreakDown_pimi[0]->Fill(E_rec,LocalWeight);
					h1_Pmiss_BreakDown_pimi[0]->Fill(p_miss_perp_2p1pi_to2p0pi[z],LocalWeight);
					h1_Q2_BreakDown_pimi[0]->Fill(reco_Q2,LocalWeight);
					h1_Nu_BreakDown_pimi[0]->Fill(nu,LocalWeight);
					h1_Pe_BreakDown_pimi[0]->Fill(V4_el.Rho(),LocalWeight);

					if (choice == 1)
					{
						h1_ECal_BreakDown_pimi[Interaction]->Fill(Ecal_2p1pi_to2p0pi[z],LocalWeight);
						h1_Eres_BreakDown_pimi[Interaction]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pimi[Interaction]->Fill(p_miss_perp_2p1pi_to2p0pi[z],LocalWeight);
						h1_Q2_BreakDown_pimi[Interaction]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pimi[Interaction]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pimi[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					//---------------------------------- 2p 1pi ->1p 1pi   ----------------------------------------------

					h1_E_tot_2p1pi_1p1pi_pimi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
				//	h1_E_tot_2p1pi_1p1pi_pimi_const_bin->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_E_rec_2p1pi_1p1pi_pimi->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
					h2_Erec_pperp_2p1pi_1p1pi_pimi->Fill(p_perp_tot_2p[z],E_rec,P_2p1pito1p1pi[z]*histoweight);
					h2_Etot_pperp_pimi->Fill(p_perp_tot_2p[z],E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_E_tot_2p1pi_1p1pi_fracfeed_pimi->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
					h1_E_rec_2p1pi_1p1pi_fracfeed_pimi->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
					h2_pperp_W_pimi->Fill(W_var,p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_theta0_pimi->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_2p1pito1p1pi[z]*histoweight);
					h2_Ecal_Eres_pimi->Fill(E_rec,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_Ecal_pimi->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h1_Ecal_Reso_pimi->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_2p1pito1p1pi[z]*histoweight);
					h2_Ecal_Etrue_pimi->Fill(E_tot_2p[z],Ev,P_2p1pito1p1pi[z]*histoweight);
					h2_Etrue_Ecal_pimi->Fill(Ev,E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h2_EresEcalratio_Eres_pimi->Fill(E_rec,E_rec/E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
					h2_EresEcaldiff_Eres_pimi->Fill(E_rec,E_rec-E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

					h1_xbjk_weight_pimi->Fill(x_bjk,P_2p1pito1p1pi[z]*histoweight);
					h1_Q2_weight_pimi->Fill(reco_Q2,P_2p1pito1p1pi[z]*histoweight);
					h1_Wvar_weight_pimi->Fill(W_var,P_2p1pito1p1pi[z]*histoweight);
					h1_nu_weight_pimi->Fill(nu,P_2p1pito1p1pi[z]*histoweight);
					h1_el_mom_corr_pimi->Fill(V4_el.Rho(),P_2p1pito1p1pi[z]*histoweight);
					h1_prot_mom_pimi->Fill(V3_2prot_corr[z].Mag(),P_2p1pito1p1pi[z]*histoweight);
					h1_MissMomentum_pimi->Fill(p_perp_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot_2p[z],V4_el);
					LocalWeight = P_2p1pito1p1pi[z]*histoweight;

					h1_nuCal_weight_pimi->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight_pimi->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight_pimi->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight_pimi->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight_pimi->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pimi->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					h1_ECal_BreakDown_pimi[0]->Fill(E_tot_2p[z],LocalWeight);
					h1_Eres_BreakDown_pimi[0]->Fill(E_rec,LocalWeight);
					h1_Pmiss_BreakDown_pimi[0]->Fill(p_perp_tot_2p[z],LocalWeight);
					h1_Q2_BreakDown_pimi[0]->Fill(reco_Q2,LocalWeight);
					h1_Nu_BreakDown_pimi[0]->Fill(nu,LocalWeight);
					h1_Pe_BreakDown_pimi[0]->Fill(V4_el.Rho(),LocalWeight);

					if (choice == 1)
					{
						h1_ECal_BreakDown_pimi[Interaction]->Fill(E_tot_2p[z],LocalWeight);
						h1_Eres_BreakDown_pimi[Interaction]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pimi[Interaction]->Fill(p_perp_tot_2p[z],LocalWeight);
						h1_Q2_BreakDown_pimi[Interaction]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pimi[Interaction]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pimi[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0; i < n_slice; i++){
						if (p_perp_tot_2p[z]<pperp_max[i] && p_perp_tot_2p[z]>pperp_min[i]){
							h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pimi[i]->Fill(E_tot_2p[z],P_2p1pito1p1pi[z]*histoweight);
							h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pimi[i]->Fill(E_rec,P_2p1pito1p1pi[z]*histoweight);
						}
					}

					//---------------------------------- 2p 1pi ->1p 0pi   ----------------------------------------------
					h2_Etot_pperp_pimi->Fill(p_perp_tot_2p[z],E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h2_pperp_W_pimi->Fill(W_var,p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h1_theta0_pimi->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),-P_2p1pito1p0pi[z]*histoweight);
					h2_Ecal_Eres_pimi->Fill(E_rec,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h1_Ecal_pimi->Fill(E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h1_Ecal_Reso_pimi->Fill((E_tot_2p[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],-P_2p1pito1p0pi[z]*histoweight);
					h2_Ecal_Etrue_pimi->Fill(E_tot_2p[z],Ev,-P_2p1pito1p0pi[z]*histoweight);
					h2_Etrue_Ecal_pimi->Fill(Ev,E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h2_EresEcalratio_Eres_pimi->Fill(E_rec,E_rec/E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);
					h2_EresEcaldiff_Eres_pimi->Fill(E_rec,E_rec-E_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

					h1_xbjk_weight_pimi->Fill(x_bjk,-P_2p1pito1p0pi[z]*histoweight);
					h1_Q2_weight_pimi->Fill(reco_Q2,-P_2p1pito1p0pi[z]*histoweight);
					h1_Wvar_weight_pimi->Fill(W_var,-P_2p1pito1p0pi[z]*histoweight);
					h1_nu_weight_pimi->Fill(nu,-P_2p1pito1p0pi[z]*histoweight);
					h1_el_mom_corr_pimi->Fill(V4_el.Rho(),-P_2p1pito1p0pi[z]*histoweight);
					h1_prot_mom_pimi->Fill(V3_2prot_corr[z].Mag(),-P_2p1pito1p0pi[z]*histoweight);
					h1_MissMomentum_pimi->Fill(p_perp_tot_2p[z],-P_2p1pito1p0pi[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(E_tot_2p[z],V4_el);
					LocalWeight = -P_2p1pito1p0pi[z]*histoweight;

					h1_nuCal_weight_pimi->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight_pimi->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight_pimi->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight_pimi->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight_pimi->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pimi->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					h1_ECal_BreakDown_pimi[0]->Fill(E_tot_2p[z],LocalWeight);
					h1_Eres_BreakDown_pimi[0]->Fill(E_rec,LocalWeight);
					h1_Pmiss_BreakDown_pimi[0]->Fill(p_perp_tot_2p[z],LocalWeight);
					h1_Q2_BreakDown_pimi[0]->Fill(reco_Q2,LocalWeight);
					h1_Nu_BreakDown_pimi[0]->Fill(nu,LocalWeight);
					h1_Pe_BreakDown_pimi[0]->Fill(V4_el.Rho(),LocalWeight);

					if (choice == 1)
					{
						h1_ECal_BreakDown_pimi[Interaction]->Fill(E_tot_2p[z],LocalWeight);
						h1_Eres_BreakDown_pimi[Interaction]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pimi[Interaction]->Fill(p_perp_tot_2p[z],LocalWeight);
						h1_Q2_BreakDown_pimi[Interaction]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pimi[Interaction]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pimi[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}
				}


				}//filling the histograms for 2protons

			}//1pi requirement

			//---------------------------------- 2p 2pi ----------------------------------------------

			const int N_2pi=2;
			double Ecal_2p2pi[N_2prot];
			double p_miss_perp_2p2pi[N_2prot];
			double Ptot_2p[2][2]={0};
			TLorentzVector V4_2pi_corr[2];
			if (num_pi_det == 2) {
			  C2p2piALL++;

				TVector3 V3_2pi_corr[N_2pi];
				double pion_acc_ratio[N_2pi] = {1};
				for (int i = 0; i < num_pi_det; i++)
				{
					if (choice == 0)
					{ //CLAS data
						V3_2pi_corr[i].SetXYZ(pxf[ind_pi_phot[i]],pyf[ind_pi_phot[i]],pzf[ind_pi_phot[i]]);
						V4_2pi_corr[i].SetPxPyPzE(pxf[ind_pi_phot[i]],pyf[ind_pi_phot[i]],pzf[ind_pi_phot[i]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[i]]*pf[ind_pi_phot[i]]));
						pion_acc_ratio[i] = 1; //Acceptance 1 for CLAS data
					}

					if (choice == 1) { //GENIE data
						pion_acc_ratio[i] = 0; //reset to 0 just to be same
						V3_2pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
								Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						V4_2pi_corr[i].SetPxPyPzE(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
								Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]], TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[i]]*pf[ind_pi_phot[i]]));
						double phi_pion = V3_2pi_corr[i].Phi();
						V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_2pi_corr[i].Theta();
						double pion_mom_corr = V3_2pi_corr[i].Mag();

						if (charge_pi_true[index_pi_det[i]] == 1) { //acceptance for pi plus
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi_true[index_pi_det[i]] == -1) {		//acceptance for pi minus. using electron acceptance map
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi_true[index_pi_det[i]] == 0) {		//acceptance for photon set to 1 for now F.H. 09/24/19
							pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 2proton and 2 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue;
						}
					}
				}
				double Ecal2p2pi[2][2];
				double p_miss_perp2p2pi[2][2];
				double P_tot_2p[2][2];
				//	double P_tot_2p_1step[2][2]; //for 1step subtraction of 2p2pi
				//	double P_tot_2p_2step[2][2]; //for 2step subtraction of 2p2pi

				//Allows me to select how many steps of subtraction to do
				int selection[3] = {0, 1, 2};
				//rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,V4_2prot_corr, V4_2pi_corr, charge_pi_det,V4_el, Ecal2p2pi, p_miss_perp2p2pi, P_tot_2p, selection[0]);
				rotation->prot2_pi2_rot_func(V4_2prot_unsmeared, V4_el, V4_2pi_corr, V3_q, charge_pi_det, Ecal2p2pi, p_miss_perp2p2pi, P_tot_2p);
				//	rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,V4_2prot_corr, V4_2pi_corr, charge_pi_det,V4_el, Ecal2p2pi, p_miss_perp2p2pi, P_tot_2p_1step, selection[1]);
				//	rotation->prot2_pi2_rot_func(V3_2prot_corr,V3_2prot_uncorr,V3_2pi_corr,V4_2prot_corr, V4_2pi_corr, charge_pi_det,V4_el, Ecal2p2pi, p_miss_perp2p2pi, P_tot_2p_2step, selection[2]);
				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
				//double histoweight = weight_pions * weight_protons * e_acc_ratio * wght/Mott_cross_sec;
				//Is this correct in the following loop? F.H. 09/01/19

				detTruthFlag = (((num_pi_true == 2) && (num_p_true == 2)) && ((num_pi_det == 2) && (num_p_det == 2)));

				for(int z = 0; z < N_2prot; z++){ //looping over two protons
					for(int j=0; j<2;j++) {
						if(charge_pi_true[index_pi_det[j]] == 1) //MOVE STUFF DOWN HERE
						{
						  	C2p2piPIPL++;
							if(detTruthFlag || choice == 0)
							{
								h1_test_E_cal_2p2pi_pipl->Fill(Ecal2p2pi[z][j], P_tot_2p[z][j]*histoweight);
							}
							//11.3.21 EDIT: Added 1p1pi tot h1 fill
							h1_E_cal_1p1pi_pipl_tot->Fill(Ecal2p2pi[z][j],histoweight);
							//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------
							h1_E_tot_2p2pi_pipl->Fill(Ecal2p2pi[z][j], P_tot_2p[z][j]*histoweight);
							//	h1_E_tot_2p2pi_pipl_1step->Fill(Ecal2p2pi[z][j], P_tot_2p_1step[z][j]*histoweight);
							//	h1_E_tot_2p2pi_pipl_2step->Fill(Ecal2p2pi[z][j], P_tot_2p_2step[z][j]*histoweight);
							//	h1_E_tot_2p2pi_pipl_const_bin->Fill(Ecal2p2pi[z][j], P_tot_2p[z][j]*histoweight);
							h1_E_rec_2p2pi_pipl->Fill(E_rec,P_tot_2p[z][j]*histoweight);
							h2_Erec_pperp_2p2pi_pipl->Fill(p_miss_perp2p2pi[z][j],E_rec,P_tot_2p[z][j]*histoweight);
							h2_Etot_pperp_pipl->Fill(p_miss_perp2p2pi[z][j],Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_E_tot_2p2pi_fracfeed_pipl->Fill((Ecal2p2pi[z][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_tot_2p[z][j]*histoweight);
							h1_E_rec_2p2pi_fracfeed_pipl->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_2p[z][j]*histoweight);
							h2_pperp_W_pipl->Fill(W_var,p_miss_perp2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_theta0_pipl->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_tot_2p[z][j]*histoweight);
							h2_Ecal_Eres_pipl->Fill(E_rec,Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_Ecal_pipl->Fill(Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_Ecal_Reso_pipl->Fill((Ecal2p2pi[z][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_tot_2p[z][j]*histoweight);
							h2_Ecal_Etrue_pipl->Fill(Ecal2p2pi[z][j],Ev,P_tot_2p[z][j]*histoweight);
							h2_Etrue_Ecal_pipl->Fill(Ev,Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h2_EresEcalratio_Eres_pipl->Fill(E_rec,E_rec/Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h2_EresEcaldiff_Eres_pipl->Fill(E_rec,E_rec-Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);

							h1_xbjk_weight_pipl->Fill(x_bjk,P_tot_2p[z][j]*histoweight);
							h1_Q2_weight_pipl->Fill(reco_Q2,P_tot_2p[z][j]*histoweight);
							h1_Wvar_weight_pipl->Fill(W_var,P_tot_2p[z][j]*histoweight);
							h1_nu_weight_pipl->Fill(nu,P_tot_2p[z][j]*histoweight);
							h1_el_mom_corr_pipl->Fill(V4_el.Rho(),P_tot_2p[z][j]*histoweight);
							h1_prot_mom_pipl->Fill(V3_2prot_corr[z].Mag(),P_tot_2p[z][j]*histoweight);
							h1_MissMomentum_pipl->Fill(p_miss_perp2p2pi[z][j],P_tot_2p[z][j]*histoweight);

							// -----------------------------------------------------------------------------------------------
							// Reconstruct xB, W, Q2 using Ecal instead of Etrue

							CalKineVars = CalculateCalKineVars(Ecal2p2pi[z][j],V4_el);
							LocalWeight = P_tot_2p[z][j]*histoweight;

							h1_nuCal_weight_pipl->Fill(CalKineVars.at(0),LocalWeight);
							h1_Q2Cal_weight_pipl->Fill(CalKineVars.at(1),LocalWeight);
							h1_xbjkCal_weight_pipl->Fill(CalKineVars.at(2),LocalWeight);
							h1_WvarCal_weight_pipl->Fill(CalKineVars.at(3),LocalWeight);

							h2_Q2_nu_weight_pipl->Fill(nu,reco_Q2,LocalWeight);
							if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pipl->Fill(nu,reco_Q2,LocalWeight); }

							// Fill plots based on underlying interactions

							h1_ECal_BreakDown_pipl[0]->Fill(Ecal2p2pi[z][j],LocalWeight);
							h1_Eres_BreakDown_pipl[0]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pipl[0]->Fill(p_miss_perp2p2pi[z][j],LocalWeight);
							h1_Q2_BreakDown_pipl[0]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pipl[0]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pipl[0]->Fill(V4_el.Rho(),LocalWeight);

		 					if (choice == 1)
							{
								h1_ECal_BreakDown_pipl[Interaction]->Fill(Ecal2p2pi[z][j],LocalWeight);
								h1_Eres_BreakDown_pipl[Interaction]->Fill(E_rec,LocalWeight);
								h1_Pmiss_BreakDown_pipl[Interaction]->Fill(p_miss_perp2p2pi[z][j],LocalWeight);
								h1_Q2_BreakDown_pipl[Interaction]->Fill(reco_Q2,LocalWeight);
								h1_Nu_BreakDown_pipl[Interaction]->Fill(nu,LocalWeight);
								h1_Pe_BreakDown_pipl[Interaction]->Fill(V4_el.Rho(),LocalWeight);
							}

							// -----------------------------------------------------------------------------------------------

							for(int i = 0; i < n_slice; i++)
							{
								if (p_miss_perp2p2pi[z][j]<pperp_max[i] && p_miss_perp2p2pi[z][j]>pperp_min[i]){
									h1_Etot_p_bkgd_slice_2p2pi_pipl[i]->Fill(Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
									h1_Erec_p_bkgd_slice_2p2pi_pipl[i]->Fill(E_rec,P_tot_2p[z][j]*histoweight);
								}
							}
						}//Minor bug fix: added charge conditional
						else if(charge_pi_true[index_pi_det[j]] == -1)
						{
							C2p2piPIMI++;
							if(detTruthFlag)
							{
								h1_test_E_cal_2p2pi_pimi->Fill(Ecal2p2pi[z][j], P_tot_2p[z][j]*histoweight);
							}
							h1_E_cal_1p1pi_pimi_tot->Fill(Ecal2p2pi[z][j], histoweight);
							//---------------------------------- 2p 2pi ->1p 0pi   ----------------------------------------------
							h1_E_tot_2p2pi_pimi->Fill(Ecal2p2pi[z][j], P_tot_2p[z][j]*histoweight);
							//h1_E_tot_2p2pi_pimi_1step->Fill(Ecal2p2pi[z][j], P_tot_2p_1step[z][j]*histoweight);
							//h1_E_tot_2p2pi_pimi_2step->Fill(Ecal_2p2pi[z][j], P_tot_2p_2step[z][j]*histoweight);
							//h1_E_tot_2p2pi_pimi_const_bin->Fill(Ecal2p2pi[z][j], P_tot_2p[z][j]*histoweight);
							h1_E_rec_2p2pi_pimi->Fill(E_rec,P_tot_2p[z][j]*histoweight);
							h2_Erec_pperp_2p2pi_pimi->Fill(p_miss_perp2p2pi[z][j],E_rec,P_tot_2p[z][j]*histoweight);
							h2_Etot_pperp_pimi->Fill(p_miss_perp2p2pi[z][j],Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_E_tot_2p2pi_fracfeed_pimi->Fill((Ecal2p2pi[z][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_tot_2p[z][j]*histoweight);
							h1_E_rec_2p2pi_fracfeed_pimi->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_2p[z][j]*histoweight);
							h2_pperp_W_pimi->Fill(W_var,p_miss_perp2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_theta0_pimi->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_2prot_uncorr[z])*TMath::RadToDeg(),P_tot_2p[z][j]*histoweight);
							h2_Ecal_Eres_pimi->Fill(E_rec,Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_Ecal_pimi->Fill(Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h1_Ecal_Reso_pimi->Fill((Ecal2p2pi[z][j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_tot_2p[z][j]*histoweight);
							h2_Ecal_Etrue_pimi->Fill(Ecal2p2pi[z][j],Ev,P_tot_2p[z][j]*histoweight);
							h2_Etrue_Ecal_pimi->Fill(Ev,Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h2_EresEcalratio_Eres_pimi->Fill(E_rec,E_rec/Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
							h2_EresEcaldiff_Eres_pimi->Fill(E_rec,E_rec-Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);

							h1_xbjk_weight_pimi->Fill(x_bjk,P_tot_2p[z][j]*histoweight);
							h1_Q2_weight_pimi->Fill(reco_Q2,P_tot_2p[z][j]*histoweight);
							h1_Wvar_weight_pimi->Fill(W_var,P_tot_2p[z][j]*histoweight);
							h1_nu_weight_pimi->Fill(nu,P_tot_2p[z][j]*histoweight);
							h1_el_mom_corr_pimi->Fill(V4_el.Rho(),P_tot_2p[z][j]*histoweight);
							h1_prot_mom_pimi->Fill(V3_2prot_corr[z].Mag(),P_tot_2p[z][j]*histoweight);
							h1_MissMomentum_pimi->Fill(p_miss_perp2p2pi[z][j],P_tot_2p[z][j]*histoweight);

							// -----------------------------------------------------------------------------------------------
							// Reconstruct xB, W, Q2 using Ecal instead of Etrue

							CalKineVars = CalculateCalKineVars(Ecal2p2pi[z][j],V4_el);
							LocalWeight = P_tot_2p[z][j]*histoweight;

							h1_nuCal_weight_pimi->Fill(CalKineVars.at(0),LocalWeight);
							h1_Q2Cal_weight_pimi->Fill(CalKineVars.at(1),LocalWeight);
							h1_xbjkCal_weight_pimi->Fill(CalKineVars.at(2),LocalWeight);
							h1_WvarCal_weight_pimi->Fill(CalKineVars.at(3),LocalWeight);

							h2_Q2_nu_weight_pimi->Fill(nu,reco_Q2,LocalWeight);
							if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pimi->Fill(nu,reco_Q2,LocalWeight); }

							// Fill plots based on underlying interactions

							h1_ECal_BreakDown_pimi[0]->Fill(Ecal2p2pi[z][j],LocalWeight);
							h1_Eres_BreakDown_pimi[0]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pimi[0]->Fill(p_miss_perp2p2pi[z][j],LocalWeight);
							h1_Q2_BreakDown_pimi[0]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pimi[0]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pimi[0]->Fill(V4_el.Rho(),LocalWeight);

							if (choice == 1)
							{
								h1_ECal_BreakDown_pimi[Interaction]->Fill(Ecal2p2pi[z][j],LocalWeight);
								h1_Eres_BreakDown_pimi[Interaction]->Fill(E_rec,LocalWeight);
								h1_Pmiss_BreakDown_pimi[Interaction]->Fill(p_miss_perp2p2pi[z][j],LocalWeight);
								h1_Q2_BreakDown_pimi[Interaction]->Fill(reco_Q2,LocalWeight);
								h1_Nu_BreakDown_pimi[Interaction]->Fill(nu,LocalWeight);
								h1_Pe_BreakDown_pimi[Interaction]->Fill(V4_el.Rho(),LocalWeight);
							}

							// -----------------------------------------------------------------------------------------------

							for(int i = 0; i < n_slice; i++) {

								if (p_miss_perp2p2pi[z][j]<pperp_max[i] && p_miss_perp2p2pi[z][j]>pperp_min[i]){
									h1_Etot_p_bkgd_slice_2p2pi_pimi[i]->Fill(Ecal2p2pi[z][j],P_tot_2p[z][j]*histoweight);
									h1_Erec_p_bkgd_slice_2p2pi_pimi[i]->Fill(E_rec,P_tot_2p[z][j]*histoweight);
								}
							}
						} //pimi
					} //Filling the histogram for two pions
				}//end loop over 2 protons
			} //2pi requirement
		}//2prot requirement */

		// -------------------------------------------------------------------------------------------------------------------------------------

		//Events with exactly 3 protons

		
		if(num_p_det == 3)
		{
			const int N_3p = 3;
			TLorentzVector V4_p_uncorr[N_3p], V4_p_corr[N_3p],V4_prot_el[N_3p];
			TVector3 V3_prot_uncorr[N_3p],V3_prot_corr[N_3p],V3_3p_rot[N_3p];
			double E_cal[N_3p],p_miss_perp[N_3p],P_3pto1p[N_3p];
			double N_p1[N_3p]={0};
			double N_p_three=0;
			double E_cal_3pto1p[3]={0};
			double p_miss_perp_3pto1p[3]={0};
			int N_comb = 3;
			const int N_2p = 2;
			double E_cal_3pto2p[3][N_2p]={0};
			double p_miss_perp_3pto2p[3][N_2p]={0};
			double P_3pto2p[3][N_2p]={0};
			TVector3 V3_2p_rot[N_2p], V3_prot_el[N_3p][N_2p];
			double p_acc_ratio[N_3p] = {1};

			for(int i = 0; i < N_3p; i++) {

				N_p1[i] = 0;

				V4_p_uncorr[i].SetPxPyPzE(pxf[index_p_det[i]],pyf[index_p_det[i]],pzf[index_p_det[i]],TMath::Sqrt(m_prot*m_prot+pf[index_p_det[i]]*pf[index_p_det[i]]));
				V3_prot_uncorr[i] = V4_p_uncorr[i].Vect();

				if (choice == 0) { //CLAS data

					V3_prot_corr[i].SetXYZ(pxf[index_p_det[i]+60], pyf[index_p_det[i]+60], pzf[index_p_det[i]+60]);
					V4_p_corr[i].SetPxPyPzE(pxf[index_p_det[i]+60], pyf[index_p_det[i]+60], pzf[index_p_det[i]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p_det[i]+60]*pf[index_p_det[i]+60]));
					p_acc_ratio[i] = 1; //Acceptance is 1 for CLAS datafile
				}

				if (choice == 1) { //GENIE data

					p_acc_ratio[i] = 0; //Reset just to be sure
					V3_prot_corr[i].SetXYZ(Smeared_Pp[i]/pf[index_p_det[i]] * pxf[index_p_det[i]],Smeared_Pp[i]/pf[index_p_det[i]] * pyf[index_p_det[i]],
						Smeared_Pp[i]/pf[index_p_det[i]] * pzf[index_p_det[i]]);
					V4_p_corr[i].SetPxPyPzE(Smeared_Pp[i]/pf[index_p_det[i]] * pxf[index_p_det[i]],Smeared_Pp[i]/pf[index_p_det[i]] * pyf[index_p_det[i]],
						Smeared_Pp[i]/pf[index_p_det[i]] * pzf[index_p_det[i]],Smeared_Ep[i]);

					double phi_prot = V3_el.Phi(); //in Radians
					V3_prot_corr[i].SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double p_theta = V3_prot_corr[i].Theta();
					double prot_mom_corr = V3_prot_corr[i].Mag();
					//Proton acceptance weight
					p_acc_ratio[i] = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);
					if ( fabs(p_acc_ratio[i]) != p_acc_ratio[i] ) { continue; }
				}

				V4_prot_el[i] = V4_p_corr[i] + V4_el;
				E_cal[i] = V4_el.E()+ V4_p_corr[i].E() - m_prot + bind_en[ftarget];
				p_miss_perp[i] = TMath::Sqrt(V4_prot_el[i].Px()*V4_prot_el[i].Px() + V4_prot_el[i].Py()*V4_prot_el[i].Py());

			} //end loop over N_3p

			V3_prot_el[0][0]=V4_el.Vect()+V3_prot_uncorr[0];
			V3_prot_el[0][1]=V4_el.Vect()+V3_prot_uncorr[1];
			V3_prot_el[1][0]=V4_el.Vect()+V3_prot_uncorr[0];
			V3_prot_el[1][1]=V4_el.Vect()+V3_prot_uncorr[2];
			V3_prot_el[2][0]=V4_el.Vect()+V3_prot_uncorr[1];
			V3_prot_el[2][1]=V4_el.Vect()+V3_prot_uncorr[2];

			//acceptance weight for all three protons ( = 1 for CLAS data)
			double weight_protons =	p_acc_ratio[0] * p_acc_ratio[1] * p_acc_ratio[2];
			//----------------------------------3p 1pi ----------------------------------------------------------

			if (num_pi_det == 1)
			{
			  	C3p1piALL++;
				double P_tot_3p[N_3p]={0};
				double Ecal_3p1pi[N_3p]={0};
				double p_miss_perp_3p1pi[N_3p]={0};
				TVector3 V3_pi_corr;
				TLorentzVector V4_pi_corr;
				double pion_acc_ratio = 1;

				if (choice == 0) { //CLAS data
					V3_pi_corr.SetXYZ(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]]);
					V4_pi_corr.SetPxPyPzE(pxf[ind_pi_phot[0]],pyf[ind_pi_phot[0]],pzf[ind_pi_phot[0]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[0]]*pf[ind_pi_phot[0]]));
					pion_acc_ratio = 1; //Acceptance is 1 for CLAS datafile
				}

				if (choice == 1){ //GENIE data

					pion_acc_ratio = 0; //Reset to 0 just to be sure
					V3_pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);
					V4_pi_corr.SetPxPyPzE(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[0]]*pf[ind_pi_phot[0]]));
					double phi_pion = V3_pi_corr.Phi(); //in Radians
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_pi_corr.Theta();
					double pion_mom_corr = V3_pi_corr.Mag();

					if (charge_pi_det[0] == 1) { //acceptance for pi plus
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi_det[0] == -1) {		//acceptance for pi minus. using electron acceptance map
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi_det[0] == 0) {		//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						 pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 3proton and 1 Pion loop. pion_acc_ratio is still 0. Continue with next event " << std::endl;	continue; }

				}
				int selection[3] = {0, 1, 2};
				//	Double P_tot_3p_1step[3] = {0.0};
				//	Double P_tot_3p_2step[3] = {0.0};
				rotation->prot3_pi1_rot_func(V4_p_uncorr, V4_el, V4_pi_corr, V3_q, charge_pi_det[0], Ecal_3p1pi, p_miss_perp_3p1pi, P_tot_3p, V3_prot_corr, V3_pi_corr, V4_p_corr);
				//rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot_uncorr, V3_pi_corr, V4_p_corr, V4_pi_corr, charge_pi_det[0], V4_el, Ecal_3p1pi, p_miss_perp_3p1pi, P_tot_3p, selection[0]);
				//	rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot_uncorr, V3_pi_corr, V4_p_corr, V4_pi_corr, charge_pi_det[0], V4_el, Ecal_3p1pi, p_miss_perp_3p1pi, P_tot_3p_1step, selection);
				//	rotation->prot3_pi1_rot_func(V3_prot_corr,V3_prot_uncorr, V3_pi_corr, V4_p_corr, V4_pi_corr, charge_pi_det[0], V4_el, Ecal_3p1pi, p_miss_perp_3p1pi, P_tot_3p_2step, selection);
				//for CLAS data is histoweight = 1/Mott_cross_sec
				//double histoweight = pion_acc_ratio * weight_protons * e_acc_ratio * wght/Mott_cross_sec;
				//Weight for 3protons, 1 pion, 1 electron, GENIE weight and Mott cross section
				detTruthFlag = (((num_p_det == 3) && (num_p_true == 3)) && ((num_pi_det == 1) && (num_pi_true == 1)));
				for(int j = 0; j < N_3p; j++) { //loop over 3 protons
					if(charge_pi_det[0] == 1)
					{
						if(detTruthFlag)
						{
							h1_test_E_cal_3p1pi_pipl->Fill(E_cal[j], P_tot_3p[j]*histoweight);
						}
						//11.3.21 EDIT: Added 1p1pi_tot h1 fill
						h1_E_cal_1p1pi_pipl_tot->Fill(E_cal[j],histoweight);
					    C3p1piPIPL++;
						h1_E_tot_3p1pi_pipl->Fill(E_cal[j], P_tot_3p[j]*histoweight);
						//	h1_E_tot_3p1pi_pipl_1step->Fill(E_cal[j], P_tot_3p_1step[j]*histoweight);
						//	h1_E_tot_3p1pi_pipl_2step->Fill(E_cal[j], P_tot_3p_2step[j]*histoweight);
						//	h1_E_tot_3p1pi_pipl_const_bin->Fill(E_cal[j], P_tot_3p[j]*histoweight);
						h1_E_rec_3p1pi_pipl->Fill(E_rec,P_tot_3p[j]*histoweight);
						h2_Erec_pperp_3p1pi_pipl->Fill(p_miss_perp[j],E_rec,P_tot_3p[j]*histoweight);
						h2_Etot_pperp_pipl->Fill(p_miss_perp[j],E_cal[j],P_tot_3p[j]*histoweight);
						h1_E_tot_3p1pi_fracfeed_pipl->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_tot_3p[j]*histoweight);
					    h1_E_rec_3p1pi_fracfeed_pipl->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_3p[j]*histoweight);
						h2_pperp_W_pipl->Fill(W_var,p_miss_perp[j],P_tot_3p[j]*histoweight);
						h1_theta0_pipl->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_tot_3p[j]*histoweight);
						h2_Ecal_Eres_pipl->Fill(E_rec,E_cal[j],P_tot_3p[j]*histoweight);
						h1_Ecal_pipl->Fill(E_cal[j],P_tot_3p[j]*histoweight);
						h1_Ecal_Reso_pipl->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_tot_3p[j]*histoweight);
						h2_Ecal_Etrue_pipl->Fill(E_cal[j],Ev,P_tot_3p[j]*histoweight);
						h2_Etrue_Ecal_pipl->Fill(Ev,E_cal[j],P_tot_3p[j]*histoweight);
						h2_EresEcalratio_Eres_pipl->Fill(E_rec,E_rec/E_cal[j],P_tot_3p[j]*histoweight);
						h2_EresEcaldiff_Eres_pipl->Fill(E_rec,E_rec-E_cal[j],P_tot_3p[j]*histoweight);

						h1_xbjk_weight_pipl->Fill(x_bjk,P_tot_3p[j]*histoweight);
						h1_Q2_weight_pipl->Fill(reco_Q2,P_tot_3p[j]*histoweight);
						h1_Wvar_weight_pipl->Fill(W_var,P_tot_3p[j]*histoweight);
						h1_nu_weight_pipl->Fill(nu,P_tot_3p[j]*histoweight);
						h1_el_mom_corr_pipl->Fill(V4_el.Rho(),P_tot_3p[j]*histoweight);
						h1_prot_mom_pipl->Fill(V3_prot_corr[j].Mag(),P_tot_3p[j]*histoweight);
						h1_MissMomentum_pipl->Fill(p_miss_perp[j],P_tot_3p[j]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(E_cal[j],V4_el);
						LocalWeight = P_tot_3p[j]*histoweight;

						h1_nuCal_weight_pipl->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pipl->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pipl->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pipl->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pipl->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pipl->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pipl[0]->Fill(E_cal[j],LocalWeight);
						h1_Eres_BreakDown_pipl[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pipl[0]->Fill(p_miss_perp[j],LocalWeight);
						h1_Q2_BreakDown_pipl[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pipl[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pipl[0]->Fill(V4_el.Rho(),LocalWeight);

						if (choice == 1)
						{
							h1_ECal_BreakDown_pipl[Interaction]->Fill(E_cal[j],LocalWeight);
							h1_Eres_BreakDown_pipl[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pipl[Interaction]->Fill(p_miss_perp[j],LocalWeight);
							h1_Q2_BreakDown_pipl[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pipl[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pipl[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++)
						{
							if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
								h1_Etot_3p1pi_slice_pipl[i]->Fill(E_cal[j],P_tot_3p[j]*histoweight);
								h1_Erec_3p1pi_slice_pipl[i]->Fill(E_rec,P_tot_3p[j]*histoweight);
							}

						}

					}//11.3.21 Minor bug fix: Added pimi charge fix charge_pi_det<0 conditional
					else if(charge_pi_det[0] == -1)
					{
					  	C3p1piPIMI++;
						if(detTruthFlag)
						{
							h1_test_E_cal_3p1pi_pimi->Fill(E_cal[j], P_tot_3p[j]*histoweight);
						}
						h1_E_cal_1p1pi_pimi_tot->Fill(E_cal[j], histoweight);
						h1_E_tot_3p1pi_pimi->Fill(E_cal[j], P_tot_3p[j]*histoweight);
						//		h1_E_tot_3p1pi_pimi_1step->Fill(E_cal[j], P_tot_3p_1step[j]*histoweight);
						//		h1_E_tot_3p1pi_pimi_2step->Fill(E_cal[j], P_tot_3p_2step[j]*histoweight);
						//h1_E_tot_3p1pi_pimi_const_bin->Fill(E_cal[j], P_tot_3p[j]*histoweight);
						h1_E_rec_3p1pi_pimi->Fill(E_rec,P_tot_3p[j]*histoweight);
						h2_Erec_pperp_3p1pi_pimi->Fill(p_miss_perp[j],E_rec,P_tot_3p[j]*histoweight);
						h2_Etot_pperp_pimi->Fill(p_miss_perp[j],E_cal[j],P_tot_3p[j]*histoweight);
						h1_E_tot_3p1pi_fracfeed_pimi->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en], P_tot_3p[j]*histoweight);
						h1_E_rec_3p1pi_fracfeed_pimi->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_tot_3p[j]*histoweight);
						h2_pperp_W_pimi->Fill(W_var,p_miss_perp[j],P_tot_3p[j]*histoweight);
						h1_theta0_pimi->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr[j])*TMath::RadToDeg(),P_tot_3p[j]*histoweight);
						h2_Ecal_Eres_pimi->Fill(E_rec,E_cal[j],P_tot_3p[j]*histoweight);
						h1_Ecal_pimi->Fill(E_cal[j],P_tot_3p[j]*histoweight);
						h1_Ecal_Reso_pimi->Fill((E_cal[j]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_tot_3p[j]*histoweight);
						h2_Ecal_Etrue_pimi->Fill(E_cal[j],Ev,P_tot_3p[j]*histoweight);
						h2_Etrue_Ecal_pimi->Fill(Ev,E_cal[j],P_tot_3p[j]*histoweight);
						h2_EresEcalratio_Eres_pimi->Fill(E_rec,E_rec/E_cal[j],P_tot_3p[j]*histoweight);
						h2_EresEcaldiff_Eres_pimi->Fill(E_rec,E_rec-E_cal[j],P_tot_3p[j]*histoweight);

						h1_xbjk_weight_pimi->Fill(x_bjk,P_tot_3p[j]*histoweight);
						h1_Q2_weight_pimi->Fill(reco_Q2,P_tot_3p[j]*histoweight);
						h1_Wvar_weight_pimi->Fill(W_var,P_tot_3p[j]*histoweight);
						h1_nu_weight_pimi->Fill(nu,P_tot_3p[j]*histoweight);
						h1_el_mom_corr_pimi->Fill(V4_el.Rho(),P_tot_3p[j]*histoweight);
						h1_prot_mom_pimi->Fill(V3_prot_corr[j].Mag(),P_tot_3p[j]*histoweight);
						h1_MissMomentum_pimi->Fill(p_miss_perp[j],P_tot_3p[j]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(E_cal[j],V4_el);
						LocalWeight = P_tot_3p[j]*histoweight;

						h1_nuCal_weight_pimi->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pimi->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pimi->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pimi->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pimi->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pimi->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pimi[0]->Fill(E_cal[j],LocalWeight);
						h1_Eres_BreakDown_pimi[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pimi[0]->Fill(p_miss_perp[j],LocalWeight);
						h1_Q2_BreakDown_pimi[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pimi[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pimi[0]->Fill(V4_el.Rho(),LocalWeight);

		 				if (choice == 1)
						{
							h1_ECal_BreakDown_pimi[Interaction]->Fill(E_cal[j],LocalWeight);
							h1_Eres_BreakDown_pimi[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pimi[Interaction]->Fill(p_miss_perp[j],LocalWeight);
							h1_Q2_BreakDown_pimi[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pimi[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pimi[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++)
						{
							if (p_miss_perp[j]<pperp_max[i] && p_miss_perp[j]>pperp_min[i]){
								h1_Etot_3p1pi_slice_pimi[i]->Fill(E_cal[j],P_tot_3p[j]*histoweight);
								h1_Erec_3p1pi_slice_pimi[i]->Fill(E_rec,P_tot_3p[j]*histoweight);
							}

						}

					}

				} //end loop over N_3p

			} // 1 pi requirement ends

		} //end if num_p_det == 3  3proton requirement 

		//Events with exactly one proton

		if(num_p_det == 1)
		{

			//Vector for proton without momentum smearing
			TLorentzVector V4_prot_uncorr(pxf[index_p_det[0]],pyf[index_p_det[0]],pzf[index_p_det[0]],TMath::Sqrt(m_prot*m_prot+pf[index_p_det[0]]*pf[index_p_det[0]]));
			TVector3 V3_prot_uncorr = V4_prot_uncorr.Vect();

			//Vector for proton with momentum smearing or correction (energy loss)
			TVector3 V3_prot_corr;
			TLorentzVector V4_prot_corr;

			double p_acc_ratio = 1; //acceptance is 1 for CLAS data

			if (choice == 0) { //CLAS data
				V3_prot_corr.SetXYZ(pxf[index_p_det[0]+60], pyf[index_p_det[0]+60], pzf[index_p_det[0]+60]);
				V4_prot_corr.SetPxPyPzE(pxf[index_p_det[0]+60], pyf[index_p_det[0]+60], pzf[index_p_det[0]+60],TMath::Sqrt(m_prot*m_prot+pf[index_p_det[0]+60]*pf[index_p_det[0]+60]));
			}
			if (choice == 1) { //GENIE data
				p_acc_ratio = 0; //Reset just to be sure
				//Fiducial cuts are done in the hadron loop
				//Vector for proton with momentum smearing
				V3_prot_corr.SetXYZ(Smeared_Pp[0]/pf[index_p_det[0]] * pxf[index_p_det[0]],Smeared_Pp[0]/pf[index_p_det[0]] * pyf[index_p_det[0]],Smeared_Pp[0]/pf[index_p_det[0]] * pzf[index_p_det[0]]);
				V4_prot_corr.SetPxPyPzE(Smeared_Pp[0]/pf[index_p_det[0]] * pxf[index_p_det[0]],Smeared_Pp[0]/pf[index_p_det[0]] * pyf[index_p_det[0]],
							Smeared_Pp[0]/pf[index_p_det[0]] * pzf[index_p_det[0]],Smeared_Ep[0]);

				double phi_prot = V3_prot_corr.Phi(); //in Radians
				V3_prot_corr.SetPhi(phi_prot + TMath::Pi() ); // Vec.Phi() is between (-180,180)
				phi_prot += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

				//Proton kinematic variables
				double p_theta = V3_prot_corr.Theta();
				double prot_mom_corr = V3_prot_corr.Mag();
				//Proton weight
				p_acc_ratio = acceptance_c(prot_mom_corr, cos(p_theta), phi_prot, 2212,file_acceptance_p);
				if ( fabs(p_acc_ratio) != p_acc_ratio ) { continue; }

			}

			//---------------------------------- 1p 1pi   ----------------------------------------------
			//The following code is commented out to speed up execution for testing 
			if(num_pi_det == 1){C1p1piALL++;}
			if(num_pi_det == 1)
			{
			  C1p1piALL++;

				double N_piphot_det;
				double N_piphot_undet;
				TVector3 V3_pi_corr;
				double pion_acc_ratio = 1;
				TLorentzVector V4_pi_corr;
				if (choice == 0) { //CLAS data

					pion_acc_ratio = 1; //Acceptance is 1 for CLAS data
					V3_pi_corr.SetXYZ(pxf[index_pi_det[0]], pyf[index_pi_det[0]], pzf[index_pi_det[0]]);
					V4_pi_corr.SetPxPyPzE(pxf[index_pi_det[0]], pyf[index_pi_det[0]], pzf[index_pi_det[0]],TMath::Sqrt(m_pion*m_pion+pf[index_pi_det[0]]*pf[index_pi_det[0]]));
				}

				if (choice == 1) { //GENIE data

					pion_acc_ratio = 1; //Reset to 0 just to be sure
					V3_pi_corr.SetXYZ(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]]);
				  V4_pi_corr.SetPxPyPzE(Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pxf[ind_pi_phot[0]],Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pyf[ind_pi_phot[0]],
							Smeared_Ppi[0]/pf[ind_pi_phot[0]] * pzf[ind_pi_phot[0]] ,TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[0]]*pf[ind_pi_phot[0]]));
					double phi_pion = V3_pi_corr.Phi(); //in Radians
					V3_pi_corr.SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
					phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

					double pion_theta = V3_pi_corr.Theta();
					double pion_mom_corr = V3_pi_corr.Mag();

					if (charge_pi_det[0] == 1) { //acceptance for pi plus

						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi_det[0] == -1) { //acceptance for pi minus. using electron acceptance map
						pion_acc_ratio = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
						if ( fabs(pion_acc_ratio) != pion_acc_ratio ) { continue; }
					}
					else if (charge_pi_det[0] == 0) { //acceptance for photon/pi0 is 1 for now F.H. 09/24/19
						pion_acc_ratio = 1;
					}
					else { std::cout << "WARNING: 1 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }

				}
				double Ecal = 0;
				double p_perp_tot = 0;
				rotation->prot1_pi1_en_calc(V4_prot_uncorr, V4_pi_corr, charge_pi_det[0], V4_el, &Ecal, &p_perp_tot);

				//histoweight is 1/Mott_cross_sec for CLAS data
				//double histoweight = pion_acc_ratio * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;
				//1proton, 1 Pion, 1 electron acceptance, GENIE weight and Mott

				if(charge_pi_true[index_pi_det[0]] == 1)
				{
				  C1p1piPIPL++;
					// 11.3.21 EDIT: Added histogram filling for 1p1pi_pipl ONLY and tot events
					h1_E_cal_1p1pi_pipl_ONLY->Fill(Ecal,histoweight);
					h1_E_cal_1p1pi_pipl_tot->Fill(Ecal, histoweight);
					h1_E_tot_pipl->Fill(Ecal,histoweight);
				//	h1_E_tot_pipl_const_bin->Fill(Ecal,histoweight);
					h1_Ecal_pipl->Fill(Ecal,histoweight);
					h1_E_rec_pipl->Fill(E_rec,histoweight);
					h1_E_rec_1prot_pipl->Fill(E_rec,histoweight);
					h1_E_tot_1prot_pipl->Fill(Ecal,histoweight);
					h2_Erec_pperp_pipl->Fill(p_perp_tot,E_rec);
					h2_Etot_pperp_pipl->Fill(p_perp_tot,Ecal,histoweight);
					h1_E_rec_cutpi1_pipl->Fill(E_rec,histoweight);
					h1_E_tot_cutpi1_pipl->Fill(Ecal,histoweight);
				}//11.3.21 Minor Bug fix: added conditional for charge < 0
				else if(charge_pi_true[index_pi_det[0]] == -1)
				{
				  C1p1piPIMI++;
					// 11.3.21 EDIT: Added histogram filling for 1p1pi ONLY  and tot events
					h1_E_cal_1p1pi_pimi_ONLY->Fill(Ecal,histoweight);
					h1_E_cal_1p1pi_pimi_tot->Fill(Ecal, histoweight);
					h1_E_tot_pimi->Fill(Ecal,histoweight);
				//	h1_E_tot_pimi_const_bin->Fill(Ecal,histoweight);
					h1_Ecal_pimi->Fill(Ecal,histoweight);
					h1_E_rec_pimi->Fill(E_rec,histoweight);
					h1_E_rec_1prot_pimi->Fill(E_rec,histoweight);
					h1_E_tot_1prot_pimi->Fill(Ecal,histoweight);
					h2_Erec_pperp_pimi->Fill(p_perp_tot,E_rec);
					h2_Etot_pperp_pimi->Fill(p_perp_tot,Ecal,histoweight);

					h1_E_rec_cutpi1_pimi->Fill(E_rec,histoweight);
					h1_E_tot_cutpi1_pimi->Fill(Ecal,histoweight);
				}
			 }//end of 1p 1pi requirement

			//---------------------------------- 1p 2pi   ----------------------------------------------

			if(num_pi_det == 2)
			{
			  	C1p2piALL++;

				const int N_2pi=2;
				TVector3 V3_2pi_corr[N_2pi],V3_2pi_rot[N_2pi],V3_p_rot;
				TLorentzVector V4_2pi_corr[2];
				double P_1p0pi=0;
				//double P_1p1pi[N_2pi]={0};
				double P_1p1pi[3] = {0};
				double pion_acc_ratio[N_2pi] = {0};

				for (int i = 0; i < num_pi_det; i++)
				{
					if (choice == 0) { //CLAS data
						V3_2pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
						V4_2pi_corr[i].SetPxPyPzE(pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]], TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[i]]*pf[ind_pi_phot[i]]));
						pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
					}

					if (choice == 1) { //GENIE data

						V3_2pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						V4_2pi_corr[i].SetPxPyPzE(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[i]]*pf[ind_pi_phot[i]]));
						double phi_pion = V3_2pi_corr[i].Phi(); //in Radians
						V3_2pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_2pi_corr[i].Theta();
						double pion_mom_corr = V3_2pi_corr[i].Mag();

						if (charge_pi_true[index_pi_det[i]] == 1) { //acceptance for pi plus
							 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
							  if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi_true[index_pi_det[i]] == -1) {	//acceptance for pi minus. using electron acceptance map
							 pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
							  if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi_true[index_pi_det[i]] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
							 pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 1 Proton 2 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }
					}

				} //end loop over num_pi_det
				double Ecal[3] = {0};
				double p_miss_perp[3] = {0};

				rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi_corr, V4_prot_corr, V4_2pi_corr, charge_pi_det, V4_el, Ecal, p_miss_perp, P_1p1pi, 0);

				//weight_pions is 1 for CLAS data
				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1];
				//histoweight is 1/Mott_cross_sec for CLAS data
				//	double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;
				//1proton, 2 Pion, 1 electron acceptance, GENIE weight and Mott
				//test case to see if calc works
				/*double P_1p1pi_pipl[3] = {0};
				double Ecal_pipl[3] = {0};
				double Ecal_pimi[3] = {0};
				double P_1p1pi_pimi[3] = {0};
				double p_miss_perp_pipl[3] = {0};
				double p_miss_perp_pimi[3] = {0};*/
				double P_1p1pi_pipl[4] = {0};
				double P_1p1pi_pimi[4] = {0};
				double Ecal_pipl[4] = {0};
				double Ecal_pimi[4] = {0};
				//double Ecal_pimi[3] = {0};
				//double P_1p1pi_pimi[3] = {0};
				double p_miss_perp_pipl[3] = {0};
				double p_miss_perp_pimi[3] = {0};
				rotation->prot1_pi2_rot_funcNew(V3_prot_uncorr,V3_2pi_corr, V4_prot_corr, V4_2pi_corr, charge_pi_det, V4_el, Ecal_pipl, p_miss_perp_pipl, P_1p1pi_pipl, 1);
				rotation->prot1_pi2_rot_funcNew(V3_prot_uncorr,V3_2pi_corr, V4_prot_corr, V4_2pi_corr, charge_pi_det, V4_el, Ecal_pimi, p_miss_perp_pimi, P_1p1pi_pimi, -1);
				//rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi_corr, V4_prot_corr, V4_2pi_corr, charge_pi_det, V4_el, Ecal_pipl, p_miss_perp_pipl, P_1p1pi_pipl, 1);
				//rotation->prot1_pi2_rot_func(V3_prot_uncorr,V3_2pi_corr, V4_prot_corr, V4_2pi_corr, charge_pi_det, V4_el, Ecal_pimi, p_miss_perp_pimi, P_1p1pi_pimi, -1);

				for(int k = 0; k < 2; k++)
				{
					h1_E_tot_1p2pi_pipl->Fill(Ecal_pipl[k],P_1p1pi_pipl[k]*histoweight);
					h1_E_tot_1p2pi_pimi->Fill(Ecal_pimi[k],P_1p1pi_pimi[k]*histoweight);
				}
				//---------------------------------- 1p 2pi->1p1pi   ----------------------------------------------

				detTruthFlag = (((num_pi_true == 2) && (num_p_true == 1)) && ((num_pi_det == 2) && (num_p_det == 1)));
				if(detTruthFlag)
				{
					for(int z = 0; z < 4; z++)
					{
						if(charge_pi_true[index_pi_det[z%2]] == 1)
						{
							h1_test_E_cal_1p2pi_pipl->Fill(Ecal_pipl[z], P_1p1pi_pipl[z]);
							if(z <= 1) {h1_test_E_cal_1p2pi_pipl_piplSame->Fill(Ecal_pipl[z], P_1p1pi_pipl[z]*histoweight);}
							if(z > 1)
							{
								h1_test_E_cal_1p2pi_pipl_diff->Fill(Ecal_pipl[z], P_1p1pi_pipl[z]*histoweight);
							}
							// ^ I am a dumbass, Ecal[2] is literally 0 >:(
						}
						else if(charge_pi_true[index_pi_det[z%2]] == -1)
						{
							h1_test_E_cal_1p2pi_pimi->Fill(Ecal_pimi[z], P_1p1pi_pimi[z]);
							if(z <= 1)
							{
								h1_test_E_cal_1p2pi_pimi_pimiSame->Fill(Ecal_pimi[z], P_1p1pi_pimi[z]*histoweight);
							}
							else if(z > 1)
							{
								h1_test_E_cal_1p2pi_pimi_diff->Fill(Ecal_pimi[z], P_1p1pi_pimi[z]*histoweight);
							}
						}
					}
				}
				for(int z = 0; z < N_2pi; z++)
				{  //to consider 2 diff. 1pi states
					if(charge_pi_true[index_pi_det[z]] == 1)
					{
						C1p2piPIPL++;
						//11.3.21 EDIT: Added 1p1pi tot h1
						h1_E_cal_1p1pi_pipl_tot->Fill(Ecal[z],histoweight);

						//	h1_E_tot_1p2pi_pipl->Fill(Ecal[z],P_1p1pi[z]*histoweight);
						h1_E_rec_1p2pi_pipl->Fill(E_rec,P_1p1pi[z]*histoweight);
						h2_Erec_pperp_1p2pi_1p1pi_pipl->Fill(p_miss_perp[z],E_rec,P_1p1pi[z]*histoweight);
						h2_Etot_pperp_pipl->Fill(p_miss_perp[z],Ecal[z],P_1p1pi[z]*histoweight);
						h1_E_tot_1p2pi_fracfeed_pipl->Fill((Ecal[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
						h1_E_rec_1p2pi_fracfeed_pipl->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p1pi[z]*histoweight);
						h2_pperp_W_pipl->Fill(W_var,p_miss_perp[z],P_1p1pi[z]*histoweight);
						h1_theta0_pipl->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p1pi[z]*histoweight);
						h2_Ecal_Eres_pipl->Fill(E_rec,Ecal[z],P_1p1pi[z]*histoweight);
						h1_Ecal_pipl->Fill(Ecal[z],P_1p1pi[z]*histoweight);
						h1_Ecal_Reso_pipl->Fill((Ecal[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
						h2_Ecal_Etrue_pipl->Fill(Ecal[z],Ev,P_1p1pi[z]*histoweight);
						h2_Etrue_Ecal_pipl->Fill(Ev,Ecal[z],P_1p1pi[z]*histoweight);
						h2_EresEcalratio_Eres_pipl->Fill(E_rec,E_rec/Ecal[z],P_1p1pi[z]*histoweight);
						h2_EresEcaldiff_Eres_pipl->Fill(E_rec,E_rec-Ecal[z],P_1p1pi[z]*histoweight);

						h1_xbjk_weight_pipl->Fill(x_bjk,P_1p1pi[z]*histoweight);
						h1_Q2_weight_pipl->Fill(reco_Q2,P_1p1pi[z]*histoweight);
						h1_Wvar_weight_pipl->Fill(W_var,P_1p1pi[z]*histoweight);
						h1_nu_weight_pipl->Fill(nu,P_1p1pi[z]*histoweight);
						h1_el_mom_corr_pipl->Fill(V4_el.Rho(),P_1p1pi[z]*histoweight);
						h1_prot_mom_pipl->Fill(V3_prot_corr.Mag(),P_1p1pi[z]*histoweight);
						h1_MissMomentum_pipl->Fill(p_miss_perp[z],P_1p1pi[z]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(Ecal[z],V4_el);
						LocalWeight = P_1p1pi[z]*histoweight;

						h1_nuCal_weight_pipl->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pipl->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pipl->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pipl->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pipl->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pipl->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pipl[0]->Fill(Ecal[z],LocalWeight);
						h1_Eres_BreakDown_pipl[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pipl[0]->Fill(p_miss_perp[z],LocalWeight);
						h1_Q2_BreakDown_pipl[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pipl[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pipl[0]->Fill(V4_el.Rho(),LocalWeight);

						if (choice == 1)
						{
							h1_ECal_BreakDown_pipl[Interaction]->Fill(Ecal[z],LocalWeight);
							h1_Eres_BreakDown_pipl[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pipl[Interaction]->Fill(p_miss_perp[z],LocalWeight);
							h1_Q2_BreakDown_pipl[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pipl[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pipl[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++){

							if (p_miss_perp[z]<pperp_max[i] && p_miss_perp[z]>pperp_min[i]){

								h1_Etot_bkgd_1p2pi_pipl[i]->Fill(Ecal[z],P_1p1pi[z]*histoweight);
								h1_Erec_bkgd_1p2pi_pipl[i]->Fill(E_rec,P_1p1pi[z]*histoweight);
							}
						}

					}
					else if(charge_pi_true[index_pi_det[z]] == -1)
					{
					  	C1p2piPIMI++;
						//11.3.21 EDIT: Added 1p1pi tot h1
						h1_E_cal_1p1pi_pimi_tot->Fill(Ecal[z], histoweight);

						//h1_E_tot_1p2pi_pimi->Fill(Ecal[z],P_1p1pi[z]*histoweight);
						//	h1_E_tot_1p2pi_pimi_const_bin->Fill(Ecal[z],P_1p1pi[z]*histoweight);
						h1_E_rec_1p2pi_pimi->Fill(E_rec,P_1p1pi[z]*histoweight);
						h2_Erec_pperp_1p2pi_1p1pi_pimi->Fill(p_miss_perp[z],E_rec,P_1p1pi[z]*histoweight);
						h2_Etot_pperp_pimi->Fill(p_miss_perp[z],Ecal[z],P_1p1pi[z]*histoweight);
						h1_E_tot_1p2pi_fracfeed_pimi->Fill((Ecal[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
						h1_E_rec_1p2pi_fracfeed_pimi->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p1pi[z]*histoweight);
						h2_pperp_W_pimi->Fill(W_var,p_miss_perp[z],P_1p1pi[z]*histoweight);
						h1_theta0_pimi->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p1pi[z]*histoweight);
						h2_Ecal_Eres_pimi->Fill(E_rec,Ecal[z],P_1p1pi[z]*histoweight);
						h1_Ecal_pimi->Fill(Ecal[z],P_1p1pi[z]*histoweight);
						h1_Ecal_Reso_pimi->Fill((Ecal[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p1pi[z]*histoweight);
						h2_Ecal_Etrue_pimi->Fill(Ecal[z],Ev,P_1p1pi[z]*histoweight);
						h2_Etrue_Ecal_pimi->Fill(Ev,Ecal[z],P_1p1pi[z]*histoweight);
						h2_EresEcalratio_Eres_pimi->Fill(E_rec,E_rec/Ecal[z],P_1p1pi[z]*histoweight);
						h2_EresEcaldiff_Eres_pimi->Fill(E_rec,E_rec-Ecal[z],P_1p1pi[z]*histoweight);

						h1_xbjk_weight_pimi->Fill(x_bjk,P_1p1pi[z]*histoweight);
						h1_Q2_weight_pimi->Fill(reco_Q2,P_1p1pi[z]*histoweight);
						h1_Wvar_weight_pimi->Fill(W_var,P_1p1pi[z]*histoweight);
						h1_nu_weight_pimi->Fill(nu,P_1p1pi[z]*histoweight);
						h1_el_mom_corr_pimi->Fill(V4_el.Rho(),P_1p1pi[z]*histoweight);
						h1_prot_mom_pimi->Fill(V3_prot_corr.Mag(),P_1p1pi[z]*histoweight);
						h1_MissMomentum_pimi->Fill(p_miss_perp[z],P_1p1pi[z]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(Ecal[z],V4_el);
						LocalWeight = P_1p1pi[z]*histoweight;

						h1_nuCal_weight_pimi->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pimi->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pimi->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pimi->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pimi->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pimi->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pimi[0]->Fill(Ecal[z],LocalWeight);
						h1_Eres_BreakDown_pimi[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pimi[0]->Fill(p_miss_perp[z],LocalWeight);
						h1_Q2_BreakDown_pimi[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pimi[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pimi[0]->Fill(V4_el.Rho(),LocalWeight);

						if (choice == 1)
						{
							h1_ECal_BreakDown_pimi[Interaction]->Fill(Ecal[z],LocalWeight);
							h1_Eres_BreakDown_pimi[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pimi[Interaction]->Fill(p_miss_perp[z],LocalWeight);
							h1_Q2_BreakDown_pimi[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pimi[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pimi[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++){

							if (p_miss_perp[z]<pperp_max[i] && p_miss_perp[z]>pperp_min[i]){

								h1_Etot_bkgd_1p2pi_pimi[i]->Fill(Ecal[z],P_1p1pi[z]*histoweight);
								h1_Erec_bkgd_1p2pi_pimi[i]->Fill(E_rec,P_1p1pi[z]*histoweight);
							}
						}

					} //end if 
				} //end loop over N_2pi

			}//1p 2pi statetment ends

			//---------------------------------- 1p 3pi   ----------------------------------------------

			if(num_pi_det == 3)
			{
			 	C1p3piALL++;
				const int N_3pi = 3;
				TVector3 V3_3pi_corr[N_3pi],V3_3pi_rot[N_3pi],V3_p_rot;
				TLorentzVector V4_3pi_corr[3];
				double P_1p3pi[3] = {0};
				double pion_acc_ratio[N_3pi] = {1};

				for (int i = 0; i < num_pi_det; i++)
				{

					if (choice == 0) { //CLAS data
						V3_3pi_corr[i].SetXYZ( pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]]);
						V4_3pi_corr[i].SetPxPyPzE(pxf[ind_pi_phot[i]], pyf[ind_pi_phot[i]], pzf[ind_pi_phot[i]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[i]]*pf[ind_pi_phot[i]]));
						pion_acc_ratio[i] = 1; //Acceptance is 1 for CLAS data
					}

					if (choice == 1) { //GENIE data
						pion_acc_ratio[i] = 0; //Reset to 0 just to be sure
						V3_3pi_corr[i].SetXYZ(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Pp[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Pp[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]]);
						V4_3pi_corr[i].SetPxPyPzE(Smeared_Ppi[i]/pf[ind_pi_phot[i]] * pxf[ind_pi_phot[i]],Smeared_Pp[i]/pf[ind_pi_phot[i]] * pyf[ind_pi_phot[i]],
									Smeared_Pp[i]/pf[ind_pi_phot[i]] * pzf[ind_pi_phot[i]],TMath::Sqrt(m_pion*m_pion+pf[ind_pi_phot[i]]*pf[ind_pi_phot[i]]));
						double phi_pion = V3_3pi_corr[i].Phi(); //in Radians
						V3_3pi_corr[i].SetPhi(phi_pion + TMath::Pi() ); // Vec.Phi() is between (-180,180)
						phi_pion += TMath::Pi(); // GENIE coordinate system flipped with respect to CLAS

						double pion_theta = V3_3pi_corr[i].Theta();
						double pion_mom_corr = V3_3pi_corr[i].Mag();

						if (charge_pi_true[index_pi_det[i]] == 1) { //acceptance for pi plus
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, 211, file_acceptance_pip);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi_true[index_pi_det[i]] == -1) {	//acceptance for pi minus. using electron acceptance map
							pion_acc_ratio[i] = acceptance_c(pion_mom_corr, cos(pion_theta), phi_pion, -211, file_acceptance);
							if ( fabs(pion_acc_ratio[i]) != pion_acc_ratio[i] ) { continue; }
						}
						else if (charge_pi_true[index_pi_det[i]] == 0) {	//acceptance for photon/pi0 is 1 for now F.H. 09/24/19
							pion_acc_ratio[i] = 1;
						}
						else { std::cout << "WARNING: 3 Pion Events. pion_acc_ratio is still 0. Continue with next event " << std::endl;  continue; }
					}

				} //end loop over num_p_deti
				double Ecal1p3pi[3] = {0};
				double p_perp1p3pi[3] = {0};
				rotation->prot1_pi3_rot_func(V3_prot_uncorr, V3_3pi_corr, V4_prot_corr, V4_3pi_corr, charge_pi_det, V4_el, Ecal1p3pi, p_perp1p3pi, P_1p3pi, 0);
		 		//weight_pions is 1 for CLAS data
				double weight_pions = pion_acc_ratio[0] * pion_acc_ratio[1] * pion_acc_ratio[2];
				//histoweight is 1/Mott_cross_sec for CLAS data
				//double histoweight = weight_pions * p_acc_ratio * e_acc_ratio * wght/Mott_cross_sec;
				//1proton, 3 Pion, 1 electron acceptance, GENIE weight and Mott
				//Testing new subtraction function 
				double Ecal1p3pi_pimi[3] = {0};
				double p_perp1p3pi_pimi[3] = {0};
				double P_1p3pi_pimi[3] = {0};
				double Ecal1p3pi_pipl[3] = {0};
				double p_perp1p3pi_pipl[3] = {0};
				double P_1p3pi_pipl[3] = {0};
				rotation->prot1_pi3_rot_func(V3_prot_uncorr, V3_3pi_corr, V4_prot_corr, V4_3pi_corr, charge_pi_det, V4_el, Ecal1p3pi_pimi, p_perp1p3pi_pimi, P_1p3pi_pimi, -1);
				rotation->prot1_pi3_rot_func(V3_prot_uncorr, V3_3pi_corr, V4_prot_corr, V4_3pi_corr, charge_pi_det, V4_el, Ecal1p3pi_pipl, p_perp1p3pi_pipl, P_1p3pi_pipl, 1);

				//---------------------------------- 1p 3pi->1p 0pi  total ?? F.H. 08/13/19 check logic here compared to 1p 2pi case ----------------------------
				detTruthFlag = (((num_pi_true == 3) && (num_p_true == 1)) && ((num_pi_det == 3) && (num_p_det == 1)));

				for(int z=0;z<3;z++)
				{
					if(charge_pi_true[index_pi_det[z]] == 1)
					{
						C1p3piPIPL++;
						if(detTruthFlag)
						{
							h1_test_E_cal_1p3pi_pipl->Fill(Ecal1p3pi_pipl[z], P_1p3pi_pipl[z]*histoweight);
						}
						//11.3.21 EDIT: Added 1p1pi tot h1
						h1_E_cal_1p1pi_pipl_tot->Fill(Ecal1p3pi[z], histoweight);
						h1_E_tot_1p3pi_pipl->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
						//h1_E_tot_1p3pi_pipl_const_bin->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
						h1_E_rec_1p3pi_pipl->Fill(E_rec,P_1p3pi[z]*histoweight);
						h2_Erec_pperp_1p3pi_pipl->Fill(p_perp1p3pi[z],E_rec,P_1p3pi[z]*histoweight);
						h2_Etot_pperp_pipl->Fill(p_perp1p3pi[z],Ecal1p3pi[z],P_1p3pi[z]*histoweight);
						h1_E_tot_1p3pi_fracfeed_pipl->Fill((Ecal1p3pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi[z]*histoweight);
						h1_E_rec_1p3pi_fracfeed_pipl->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p3pi[z]*histoweight);
						h2_pperp_W_pipl->Fill(W_var,p_perp1p3pi[z],P_1p3pi[z]*histoweight);
						h1_theta0_pipl->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p3pi[z]*histoweight);
						h2_Ecal_Eres_pipl->Fill(E_rec,Ecal1p3pi[z],P_1p3pi[z]*histoweight);
						h1_Ecal_pipl->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
						h1_Ecal_Reso_pipl->Fill((Ecal1p3pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi[z]*histoweight);
						h2_Ecal_Etrue_pipl->Fill(Ecal1p3pi[z],Ev,P_1p3pi[z]*histoweight);
						h2_Etrue_Ecal_pipl->Fill(Ev,Ecal1p3pi[z],P_1p3pi[z]*histoweight);
						h2_EresEcalratio_Eres_pipl->Fill(E_rec,E_rec/Ecal1p3pi[z],P_1p3pi[z]*histoweight);
						h2_EresEcaldiff_Eres_pipl->Fill(E_rec,E_rec-Ecal1p3pi[z],P_1p3pi[z]*histoweight);

						h1_xbjk_weight_pipl->Fill(x_bjk,P_1p3pi[z]*histoweight);
						h1_Q2_weight_pipl->Fill(reco_Q2,P_1p3pi[z]*histoweight);
						h1_Wvar_weight_pipl->Fill(W_var,P_1p3pi[z]*histoweight);
						h1_nu_weight_pipl->Fill(nu,P_1p3pi[z]*histoweight);
						h1_el_mom_corr_pipl->Fill(V4_el.Rho(),P_1p3pi[z]*histoweight);
						h1_prot_mom_pipl->Fill(V3_prot_corr.Mag(),P_1p3pi[z]*histoweight);
						h1_MissMomentum_pipl->Fill(p_perp1p3pi[z],P_1p3pi[z]*histoweight);

						// -----------------------------------------------------------------------------------------------
						// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

						CalKineVars = CalculateCalKineVars(Ecal1p3pi[z],V4_el);
						LocalWeight = P_1p3pi[z]*histoweight;

						h1_nuCal_weight_pipl->Fill(CalKineVars.at(0),LocalWeight);
						h1_Q2Cal_weight_pipl->Fill(CalKineVars.at(1),LocalWeight);
						h1_xbjkCal_weight_pipl->Fill(CalKineVars.at(2),LocalWeight);
						h1_WvarCal_weight_pipl->Fill(CalKineVars.at(3),LocalWeight);

						h2_Q2_nu_weight_pipl->Fill(nu,reco_Q2,LocalWeight);
						if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pipl->Fill(nu,reco_Q2,LocalWeight); }

						// Fill plots based on underlying interactions

						h1_ECal_BreakDown_pipl[0]->Fill(Ecal1p3pi[z],LocalWeight);
						h1_Eres_BreakDown_pipl[0]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pipl[0]->Fill(p_perp1p3pi[z],LocalWeight);
						h1_Q2_BreakDown_pipl[0]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pipl[0]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pipl[0]->Fill(V4_el.Rho(),LocalWeight);

						if (choice == 1)
						{
							h1_ECal_BreakDown_pipl[Interaction]->Fill(Ecal1p3pi[z],LocalWeight);
							h1_Eres_BreakDown_pipl[Interaction]->Fill(E_rec,LocalWeight);
							h1_Pmiss_BreakDown_pipl[Interaction]->Fill(p_perp1p3pi[z],LocalWeight);
							h1_Q2_BreakDown_pipl[Interaction]->Fill(reco_Q2,LocalWeight);
							h1_Nu_BreakDown_pipl[Interaction]->Fill(nu,LocalWeight);
							h1_Pe_BreakDown_pipl[Interaction]->Fill(V4_el.Rho(),LocalWeight);
						}

						// -----------------------------------------------------------------------------------------------

						for(int i = 0; i < n_slice; i++){
							if (p_perp1p3pi[z]<pperp_max[i] && p_perp1p3pi[z]>pperp_min[i]){
								h1_Etot_bkgd_1p3pi_pipl[i]->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
								h1_Erec_bkgd_1p3pi_pipl[i]->Fill(E_rec,P_1p3pi[z]*histoweight);
							}
						}
					}//11.3.21 Minor Bug fix: did the charge stuff :)
				else if(charge_pi_true[index_pi_det[z]] == -1)
				{
				  	C1p3piPIMI++;
					if(detTruthFlag)
					{
						//to-do implement this
					}
					//11.3.21 EDIT: Added 1p1pi tot h1
					h1_E_cal_1p1pi_pipl_tot->Fill(Ecal1p3pi[z], histoweight);

					h1_E_tot_1p3pi_pimi->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
					//h1_E_tot_1p3pi_pimi_const_bin->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
					h1_E_rec_1p3pi_pimi->Fill(E_rec,P_1p3pi[z]*histoweight);
					h2_Erec_pperp_1p3pi_pimi->Fill(p_perp1p3pi[z],E_rec,P_1p3pi[z]*histoweight);
					h2_Etot_pperp_pimi->Fill(p_perp1p3pi[z],Ecal1p3pi[z],P_1p3pi[z]*histoweight);
					h1_E_tot_1p3pi_fracfeed_pimi->Fill((Ecal1p3pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi[z]*histoweight);
					h1_E_rec_1p3pi_fracfeed_pimi->Fill((E_rec-en_beam_Eqe[fbeam_en])/en_beam_Eqe[fbeam_en],P_1p3pi[z]*histoweight);
					h2_pperp_W_pimi->Fill(W_var,p_perp1p3pi[z],P_1p3pi[z]*histoweight);
					h1_theta0_pimi->Fill((V4_beam.Vect()).Angle(V4_el.Vect()+V3_prot_uncorr)*TMath::RadToDeg(),P_1p3pi[z]*histoweight);
					h2_Ecal_Eres_pimi->Fill(E_rec,Ecal1p3pi[z],P_1p3pi[z]*histoweight);
					h1_Ecal_pimi->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
					h1_Ecal_Reso_pimi->Fill((Ecal1p3pi[z]-en_beam_Ecal[fbeam_en])/en_beam_Ecal[fbeam_en],P_1p3pi[z]*histoweight);
					h2_Ecal_Etrue_pimi->Fill(Ecal1p3pi[z],Ev,P_1p3pi[z]*histoweight);
					h2_Etrue_Ecal_pimi->Fill(Ev,Ecal1p3pi[z],P_1p3pi[z]*histoweight);
					h2_EresEcalratio_Eres_pimi->Fill(E_rec,E_rec/Ecal1p3pi[z],P_1p3pi[z]*histoweight);
					h2_EresEcaldiff_Eres_pimi->Fill(E_rec,E_rec-Ecal1p3pi[z],P_1p3pi[z]*histoweight);

					h1_xbjk_weight_pimi->Fill(x_bjk,P_1p3pi[z]*histoweight);
					h1_Q2_weight_pimi->Fill(reco_Q2,P_1p3pi[z]*histoweight);
					h1_Wvar_weight_pimi->Fill(W_var,P_1p3pi[z]*histoweight);
					h1_nu_weight_pimi->Fill(nu,P_1p3pi[z]*histoweight);
					h1_el_mom_corr_pimi->Fill(V4_el.Rho(),P_1p3pi[z]*histoweight);
					h1_prot_mom_pimi->Fill(V3_prot_corr.Mag(),P_1p3pi[z]*histoweight);
					h1_MissMomentum_pimi->Fill(p_perp1p3pi[z],P_1p3pi[z]*histoweight);

					// -----------------------------------------------------------------------------------------------
					// apapadop: Reconstruct xB, W, Q2 using Ecal instead of Etrue

					CalKineVars = CalculateCalKineVars(Ecal1p3pi[z],V4_el);
					LocalWeight = P_1p3pi[z]*histoweight;

					h1_nuCal_weight_pimi->Fill(CalKineVars.at(0),LocalWeight);
					h1_Q2Cal_weight_pimi->Fill(CalKineVars.at(1),LocalWeight);
					h1_xbjkCal_weight_pimi->Fill(CalKineVars.at(2),LocalWeight);
					h1_WvarCal_weight_pimi->Fill(CalKineVars.at(3),LocalWeight);

					h2_Q2_nu_weight_pimi->Fill(nu,reco_Q2,LocalWeight);
					if (el_phi_mod > 0 && el_phi_mod< 60) {h2_Q2_nu_weight_FirstSector_pimi->Fill(nu,reco_Q2,LocalWeight); }

					// Fill plots based on underlying interactions

					h1_ECal_BreakDown_pimi[0]->Fill(Ecal1p3pi[z],LocalWeight);
					h1_Eres_BreakDown_pimi[0]->Fill(E_rec,LocalWeight);
					h1_Pmiss_BreakDown_pimi[0]->Fill(p_perp1p3pi[z],LocalWeight);
					h1_Q2_BreakDown_pimi[0]->Fill(reco_Q2,LocalWeight);
					h1_Nu_BreakDown_pimi[0]->Fill(nu,LocalWeight);
					h1_Pe_BreakDown_pimi[0]->Fill(V4_el.Rho(),LocalWeight);

					if (choice == 1)
					{
						h1_ECal_BreakDown_pimi[Interaction]->Fill(Ecal1p3pi[z],LocalWeight);
						h1_Eres_BreakDown_pimi[Interaction]->Fill(E_rec,LocalWeight);
						h1_Pmiss_BreakDown_pimi[Interaction]->Fill(p_perp1p3pi[z],LocalWeight);
						h1_Q2_BreakDown_pimi[Interaction]->Fill(reco_Q2,LocalWeight);
						h1_Nu_BreakDown_pimi[Interaction]->Fill(nu,LocalWeight);
						h1_Pe_BreakDown_pimi[Interaction]->Fill(V4_el.Rho(),LocalWeight);
					}

					// -----------------------------------------------------------------------------------------------

					for(int i = 0; i < n_slice; i++){
						if (p_perp1p3pi[z]<pperp_max[i] && p_perp1p3pi[z]>pperp_min[i]){
							h1_Etot_bkgd_1p3pi_pimi[i]->Fill(Ecal1p3pi[z],P_1p3pi[z]*histoweight);
							h1_Erec_bkgd_1p3pi_pimi[i]->Fill(E_rec,P_1p3pi[z]*histoweight);
						}
					}
				}
				}

			}//end of 1p 3pi requirement

		} // 1proton ends

	} //end of event loop (jentry)

	gStyle->SetOptFit(1);
	TH1F* h1_E_cal_pimi_sub = (TH1F*) h1_E_tot_pimi->Clone("h1_E_cal_pimi_sub");
	//Creates new histogram filled with floats and makes it a clone of h1_Ecal_pimi
	//Quotation marks are title of histogram when it is formed -> Should be same as variable name
	h1_E_cal_pimi_sub->Add(h1_E_tot_2p1pi_1p1pi_pimi, 1);
	h1_E_cal_pimi_sub->Write("process1");
	h1_E_cal_pimi_sub->Add(h1_E_tot_1p2pi_pimi, 1);
	h1_E_cal_pimi_sub->Write("process2");
	h1_E_cal_pimi_sub->Add(h1_E_tot_2p2pi_pimi, 1);
	h1_E_cal_pimi_sub->Write("process3");
	h1_E_cal_pimi_sub->Add(h1_E_tot_3p1pi_pimi, 1);
	h1_E_cal_pimi_sub->Write("process4");
	h1_E_cal_pimi_sub->Add(h1_E_tot_1p3pi_pimi, 1);
	h1_E_cal_pimi_sub->Write("process1");

	TH1F* h1_E_cal_pimi_sub_thesis = (TH1F*) h1_E_tot_pimi->Clone("h1_E_cal_pimi_sub_thesis");
	h1_E_cal_pimi_sub_thesis->Add(h1_test_E_cal_1p2pi_pimi, 1); 
	h1_E_cal_pimi_sub_thesis->Add(h1_test_E_cal_2p1pi_pimi, 1);


	TH1F* h1_E_cal_pipl_sub = (TH1F*) h1_E_tot_pipl->Clone("h1_E_cal_pipl_sub");

	//2.24.22 -> If the weights are NEGATIVE -> These must be ADDED to reflect proper SUBTRACTION
	h1_E_cal_pipl_sub->Add(h1_E_tot_3p1pi_pipl, 1);
	h1_E_cal_pipl_sub->Add(h1_E_tot_2p1pi_1p1pi_pipl, 1);
	h1_E_cal_pipl_sub->Add(h1_E_tot_2p2pi_pipl, 1);
	h1_E_cal_pipl_sub->Add(h1_E_tot_1p2pi_pipl, 1);
	h1_E_cal_pipl_sub->Add(h1_E_tot_1p3pi_pipl, 1);
	//Weights of +1 to add positively

	TH1F* h1_E_cal_pipl_sub_thesis = (TH1F*) h1_E_tot_pipl->Clone("h1_E_cal_pipl_sub_thesis");
	h1_E_cal_pipl_sub_thesis->Add(h1_test_E_cal_1p2pi_pipl, 1); 
	h1_E_cal_pipl_sub_thesis->Add(h1_test_E_cal_2p1pi_pipl, 1);

/*	for(int i = 0; i <= n_slice-1; i++) {

		//------------------------------------using the ratio of the pi- to pi+  --------------------------------------

		h1_Etot_piplpimi_subtruct_fact[i] = (TH1F*)  h1_Etot_Npi0[i]->Clone(Form("h1_Etot_piplpimi_subtruct_fact_%d",i+1));
		h1_Etot_piplpimi_subtruct_fact[i]->Add(h1_Etot_bkgd_pipl_pimi_fact[i],-1);
		h1_Erec_piplpimi_subtruct_fact[i]=(TH1F*)  h1_Erec_Npi0[i]->Clone(Form("h1_Erec_piplpimi_subtruct_fact_%d",i+1));
		h1_Erec_piplpimi_subtruct_fact[i]->Add(h1_Erec_bkgd_pipl_pimi_new_fact[i],-1);

		//------------------------------------subtracting 2p contribution from 1p events  --------------------------------------

		h1_Etot_p_bkgd_slice_sub[i]=(TH1F*) h1_Etot_piplpimi_subtruct_fact[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub_%d",i+1));
		h1_Etot_p_bkgd_slice_sub[i]->Add(h1_Etot_p_bkgd_slice[i],-1);
		h1_Erec_p_bkgd_slice_sub[i]=(TH1F*) h1_Erec_piplpimi_subtruct_fact[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub_%d",i+1));
		h1_Erec_p_bkgd_slice_sub[i]->Add(h1_Erec_p_bkgd_slice[i],-1);

		//------------------------------------undetected 3 to 2 proton subtraction --------------------------------------

		h1_Etot_p_bkgd_slice_sub32[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub32_%d",i+1));
		h1_Etot_p_bkgd_slice_sub32[i]->Add(h1_Etot_3pto2p_slice[i]);
		h1_Erec_p_bkgd_slice_sub32[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub32_%d",i+1));
		h1_Erec_p_bkgd_slice_sub32[i]->Add(h1_Erec_3pto2p_slice[i]);

		//------------------------------------undetected 3 to 1 proton addition --------------------------------------

		h1_Etot_p_bkgd_slice_sub31[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub32[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub31_%d",i+1));
		h1_Etot_p_bkgd_slice_sub31[i]->Add(h1_Etot_3pto1p_slice[i],-1);
		h1_Erec_p_bkgd_slice_sub31[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub32[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub31_%d",i+1));
		h1_Erec_p_bkgd_slice_sub31[i]->Add(h1_Erec_3pto1p_slice[i],-1);

		//------------------------------------undetected 4 to 3->2->1 proton addition --------------------------------------

		h1_Etot_p_bkgd_slice_sub43[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub31[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub43_%d",i+1));
		h1_Etot_p_bkgd_slice_sub43[i]->Add(h1_Etot_4pto3p_slice[i],-1);
		h1_Erec_p_bkgd_slice_sub43[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub31[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub43_%d",i+1));
		h1_Erec_p_bkgd_slice_sub43[i]->Add(h1_Erec_4pto3p_slice[i],-1);

		//------------------------------------undetected 4 to 3->1 proton addition --------------------------------------

		h1_Etot_p_bkgd_slice_sub431[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub43[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub431_%d",i+1));
		h1_Etot_p_bkgd_slice_sub431[i]->Add(h1_Etot_43pto1p_slice[i]);
		h1_Erec_p_bkgd_slice_sub431[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub43[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub431_%d",i+1));
		h1_Erec_p_bkgd_slice_sub431[i]->Add(h1_Erec_43pto1p_slice[i]);

		//------------------------------------undetected 4 to 2 proton addition --------------------------------------

		h1_Etot_p_bkgd_slice_sub42[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub431[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub42_%d",i+1));
		h1_Etot_p_bkgd_slice_sub42[i]->Add(h1_Etot_4pto2p_slice[i]);
		h1_Erec_p_bkgd_slice_sub42[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub431[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub42_%d",i+1));
		h1_Erec_p_bkgd_slice_sub42[i]->Add(h1_Erec_4pto2p_slice[i]);

		//------------------------------------undetected 4 to 1 proton addition --------------------------------------

		h1_Etot_p_bkgd_slice_sub41[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub42[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub41_%d",i+1));
		h1_Etot_p_bkgd_slice_sub41[i]->Add(h1_Etot_4pto1p_slice[i],-1);
		h1_Erec_p_bkgd_slice_sub41[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub42[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub41_%d",i+1));
		h1_Erec_p_bkgd_slice_sub41[i]->Add(h1_Erec_4pto1p_slice[i],-1);

		//------------------------------------undetected 1p 2pi  ------ --------------------------------------

		h1_Etot_p_bkgd_slice_sub1p2pi[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub41[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub1p2pi_%d",i+1));
		h1_Etot_p_bkgd_slice_sub1p2pi[i]->Add(h1_Etot_bkgd_1p2pi_pipl[i]);
		h1_Etot_p_bkgd_slice_sub1p2pi[i]->Add(h1_Etot_bkgd_1p2pi_pimi[i]);
		h1_Erec_p_bkgd_slice_sub1p2pi[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub41[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub1p2pi_%d",i+1));
		h1_Erec_p_bkgd_slice_sub1p2pi[i]->Add(h1_Erec_bkgd_1p2pi_pimi[i]);
		h1_Erec_p_bkgd_slice_sub1p2pi[i]->Add(h1_Erec_bkgd_1p2pi_pipl[i]);

		//------------------------------------undetected 2p 1pi ->1p 1pi  ------ --------------------------------------

		h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]=(TH1F*) h1_Etot_p_bkgd_slice_sub2p1pi_2p[i]->Clone(Form("h1_Etot_p_bkgd_slice_sub2p1pi_1p_%d",i+1));
		h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pipl[i]);
		h1_Etot_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Etot_p_bkgd_slice_2p1pi_to1p1pi_pimi[i]);
		h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]=(TH1F*) h1_Erec_p_bkgd_slice_sub2p1pi_2p[i]->Clone(Form("h1_Erec_p_bkgd_slice_sub2p1pi_1p_%d",i+1));
		h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pipl[i]);
		h1_Erec_p_bkgd_slice_sub2p1pi_1p[i]->Add(h1_Erec_p_bkgd_slice_2p1pi_to1p1pi_pimi[i]);

	}

	//------------------------------------fractional energy reconstruction plots --------------------------------------

	//------------------------------------using the ratio of the pi- to pi+  ---------------------------------------

//	TH1F *h_Erec_subtruct_piplpimi_factor =(TH1F*)  h1_E_rec_cut2_new->Clone("eRecoEnergy_slice_0");
	TH1F *h_Erec_subtruct_piplpimi_factor =(TH1F*)	h1_E_rec_cut2_new->Clone("h_Erec_subtruct_piplpimi_factor");
	h_Erec_subtruct_piplpimi_factor->Add(h1_E_rec_undetfactor,-1);

//	TH1F *h_Etot_subtruct_piplpimi_factor =(TH1F*)  h1_E_tot_cut2->Clone("epRecoEnergy_slice_0");
	TH1F *h_Etot_subtruct_piplpimi_factor=(TH1F*)  h1_E_tot_cut2->Clone("h_Etot_subtruct_piplpimi_factor");
	h_Etot_subtruct_piplpimi_factor->Add(h1_E_tot_undetfactor,-1);

	TH2F *h2_Erec_pperp_1p1pisub=(TH2F*) h2_Erec_pperp_newcut2->Clone("h2_Erec_pperp_1p1pisub");
	h2_Erec_pperp_1p1pisub->Add(h2_Erec_pperp_1p1pi,-1);

	TH1F *h_Erec_subtruct_piplpimi_factor_fracfeed =(TH1F*)  h1_E_rec_cut2_new_fracfeed->Clone("h_Erec_subtruct_piplpimi_factor_fracfeed");
	h_Erec_subtruct_piplpimi_factor_fracfeed->Add(h1_E_rec_undetfactor_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_factor_fracfeed=(TH1F*)  h1_E_tot_cut2_fracfeed->Clone("h_Etot_subtruct_piplpimi_factor_fracfeed");
	h_Etot_subtruct_piplpimi_factor_fracfeed->Add(h1_E_tot_undetfactor_fracfeed,-1);

	//-----------------------------------undetected 2 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_prot=(TH1F*)	h_Erec_subtruct_piplpimi_factor->Clone("h_Erec_subtruct_piplpimi_prot");
	h_Erec_subtruct_piplpimi_prot->Add(h1_E_rec_p_bkgd,-1);

	TH1F *h_Etot_subtruct_piplpimi_prot=(TH1F*)	h_Etot_subtruct_piplpimi_factor->Clone("h_Etot_subtruct_piplpimi_prot");
	h_Etot_subtruct_piplpimi_prot->Add(h1_E_tot_p_bkgd,-1);

	TH2F *h2_Erec_pperp_2psub=(TH2F*) h2_Erec_pperp_1p1pisub->Clone("h2_Erec_pperp_2psub");
	h2_Erec_pperp_2psub->Add(h2_Erec_pperp_2p,-1);

	TH1F *h_Erec_subtruct_piplpimi_prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_factor_fracfeed->Clone("h_Erec_subtruct_piplpimi_prot_fracfeed");
	h_Erec_subtruct_piplpimi_prot_fracfeed->Add(h1_E_rec_p_bkgd_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_factor_fracfeed->Clone("h_Etot_subtruct_piplpimi_prot_fracfeed");
	h_Etot_subtruct_piplpimi_prot_fracfeed->Add(h1_E_tot_p_bkgd_fracfeed,-1);

	 //-----------------------------------undetected 3 to 2 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_32prot=(TH1F*)	h_Erec_subtruct_piplpimi_prot->Clone("h_Erec_subtruct_piplpimi_32prot");
	h_Erec_subtruct_piplpimi_32prot->Add(h1_E_rec_3pto2p);

	TH1F *h_Etot_subtruct_piplpimi_32prot=(TH1F*)	h_Etot_subtruct_piplpimi_prot->Clone("h_Etot_subtruct_piplpimi_32prot");
	h_Etot_subtruct_piplpimi_32prot->Add(h1_E_tot_3pto2p);

	TH2F *h2_Erec_pperp_32psub=(TH2F*) h2_Erec_pperp_2psub->Clone("h2_Erec_pperp_32psub");
	h2_Erec_pperp_32psub->Add(h2_Erec_pperp_321p);

	TH1F *h_Erec_subtruct_piplpimi_32prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_32prot_fracfeed");
	h_Erec_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_rec_3pto2p_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_32prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_32prot_fracfeed");
	h_Etot_subtruct_piplpimi_32prot_fracfeed->Add(h1_E_tot_3pto2p_fracfeed);

	 //-----------------------------------undetected 3 to 1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_31prot=(TH1F*)	h_Erec_subtruct_piplpimi_32prot->Clone("h_Erec_subtruct_piplpimi_31prot");
	h_Erec_subtruct_piplpimi_31prot->Add(h1_E_rec_3pto1p,-1);

	TH1F *h_Etot_subtruct_piplpimi_31prot=(TH1F*)	h_Etot_subtruct_piplpimi_32prot->Clone("h_Etot_subtruct_piplpimi_31prot");
	h_Etot_subtruct_piplpimi_31prot->Add(h1_E_tot_3pto1p,-1);

	TH2F *h2_Erec_pperp_31psub=(TH2F*) h2_Erec_pperp_32psub->Clone("h2_Erec_pperp_31psub");
	h2_Erec_pperp_31psub->Add(h2_Erec_pperp_31p,-1);

	TH1F *h_Erec_subtruct_piplpimi_31prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_32prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_31prot_fracfeed");
	h_Erec_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_rec_3pto1p_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_31prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_32prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_31prot_fracfeed");
	h_Etot_subtruct_piplpimi_31prot_fracfeed->Add(h1_E_tot_3pto1p_fracfeed,-1);

	 //-----------------------------------undetected 4 to 3->2->1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_43prot=(TH1F*)	h_Erec_subtruct_piplpimi_31prot->Clone("h_Erec_subtruct_piplpimi_43prot");
	h_Erec_subtruct_piplpimi_43prot->Add(h1_E_rec_4pto3p,-1);

	TH1F *h_Etot_subtruct_piplpimi_43prot=(TH1F*)	h_Etot_subtruct_piplpimi_31prot->Clone("h_Etot_subtruct_piplpimi_43prot");
	h_Etot_subtruct_piplpimi_43prot->Add(h1_E_tot_4pto3p,-1);

	TH2F *h2_Erec_pperp_43psub=(TH2F*) h2_Erec_pperp_31psub->Clone("h2_Erec_pperp_43psub");
	h2_Erec_pperp_43psub->Add(h2_Erec_pperp_4321p,-1);

	TH1F *h_Erec_subtruct_piplpimi_43prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_31prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_43prot_fracfeed");
	h_Erec_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_rec_4pto3p_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_43prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_31prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_43prot_fracfeed");
	h_Etot_subtruct_piplpimi_43prot_fracfeed->Add(h1_E_tot_4pto3p_fracfeed,-1);

	 //-----------------------------------undetected 4 to 3->1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_431prot=(TH1F*)	h_Erec_subtruct_piplpimi_43prot->Clone("h_Erec_subtruct_piplpimi_431prot");
	h_Erec_subtruct_piplpimi_431prot->Add(h1_E_rec_43pto1p);

	TH1F *h_Etot_subtruct_piplpimi_431prot=(TH1F*)	h_Etot_subtruct_piplpimi_43prot->Clone("h_Etot_subtruct_piplpimi_431prot");
	h_Etot_subtruct_piplpimi_431prot->Add(h1_E_tot_43pto1p);

	TH2F *h2_Erec_pperp_431psub=(TH2F*) h2_Erec_pperp_43psub->Clone("h2_Erec_pperp_431psub");
	h2_Erec_pperp_431psub->Add(h2_Erec_pperp_431p);

	TH1F *h_Erec_subtruct_piplpimi_431prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_43prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_431prot_fracfeed");
	h_Erec_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_rec_43pto1p_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_431prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_43prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_431prot_fracfeed");
	h_Etot_subtruct_piplpimi_431prot_fracfeed->Add(h1_E_tot_43pto1p_fracfeed);

	//-----------------------------------undetected 4 to 2 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_42prot=(TH1F*)	h_Erec_subtruct_piplpimi_431prot->Clone("h_Erec_subtruct_piplpimi_42prot");
	h_Erec_subtruct_piplpimi_42prot->Add(h1_E_rec_4pto2p);

	TH1F *h_Etot_subtruct_piplpimi_42prot=(TH1F*) h_Etot_subtruct_piplpimi_431prot->Clone("h_Etot_subtruct_piplpimi_42prot");
	h_Etot_subtruct_piplpimi_42prot->Add(h1_E_tot_4pto2p);

	TH2F *h2_Erec_pperp_42psub=(TH2F*) h2_Erec_pperp_431psub->Clone("h2_Erec_pperp_42psub");
	h2_Erec_pperp_42psub->Add(h2_Erec_pperp_421p);

	TH1F *h_Erec_subtruct_piplpimi_42prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_431prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_42prot_fracfeed");
	h_Erec_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_rec_4pto2p_fracfeed);

	TH1F *h_Etot_subtruct_piplpimi_42prot_fracfeed=(TH1F*) h_Etot_subtruct_piplpimi_431prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_42prot_fracfeed");
	h_Etot_subtruct_piplpimi_42prot_fracfeed->Add(h1_E_tot_4pto2p_fracfeed);

	 //-----------------------------------undetected 4 to 1 proton subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_41prot=(TH1F*)	h_Erec_subtruct_piplpimi_42prot->Clone("h_Erec_subtruct_piplpimi_41prot");
	h_Erec_subtruct_piplpimi_41prot->Add(h1_E_rec_4pto1p,-1);

	TH1F *h_Etot_subtruct_piplpimi_41prot=(TH1F*)	h_Etot_subtruct_piplpimi_42prot->Clone("h_Etot_subtruct_piplpimi_41prot");
	h_Etot_subtruct_piplpimi_41prot->Add(h1_E_tot_4pto1p,-1);

	TH2F *h2_Erec_pperp_41psub=(TH2F*) h2_Erec_pperp_42psub->Clone("h2_Erec_pperp_41psub");
	h2_Erec_pperp_41psub->Add(h2_Erec_pperp_41p,-1);

	TH1F *h_Erec_subtruct_piplpimi_41prot_fracfeed=(TH1F*)	h_Erec_subtruct_piplpimi_42prot_fracfeed->Clone("h_Erec_subtruct_piplpimi_41prot_fracfeed");
	h_Erec_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_rec_4pto1p_fracfeed,-1);

	TH1F *h_Etot_subtruct_piplpimi_41prot_fracfeed=(TH1F*)	h_Etot_subtruct_piplpimi_42prot_fracfeed->Clone("h_Etot_subtruct_piplpimi_41prot_fracfeed");
	h_Etot_subtruct_piplpimi_41prot_fracfeed->Add(h1_E_tot_4pto1p_fracfeed,-1);

	//------------------------------------undetected 1p 2pi ->1 p1pi ------ --------------------------------------

	/*TH1F *h_Erec_subtruct_pipl_1p2pi=(TH1F*)	h_Erec_subtruct_pipl_41prot->Clone("h_Erec_subtruct_piplpimi_1p2pi");
	h_Erec_subtruct_pipl_1p2pi->Add(h1_E_rec_1p2pi_pipl);

	TH1F *h_Etot_subtruct_pipl_1p2pi=(TH1F*)	h_Etot_subtruct_pipl_41prot->Clone("h_Etot_subtruct_piplpimi_1p2pi");
	h_Etot_subtruct_pipl_1p2pi->Add(h1_E_tot_1p2pi_pipl);

	TH2F *h2_Erec_pperp_sub_1p2pi_1p1pi_pipl=(TH2F*) h2_Erec_pperp_41psub->Clone("h2_Erec_pperp_sub_1p2pi_1p1pi_pipl");
	h2_Erec_pperp_sub_1p2pi_1p1pi_pipl->Add(h2_Erec_pperp_1p2pi_1p1pi_pipl);

	TH1F *h_Erec_subtruct_pipl_1p2pi_fracfeed=(TH1F*)	h_Erec_subtruct_pipl_41prot_fracfeed->Clone("h_Erec_subtruct_pipl_1p2pi_fracfeed");
	h_Erec_subtruct_pipl_1p2pi_fracfeed->Add(h1_E_rec_1p2pi_fracfeed_pipl);

	TH1F *h_Etot_subtruct_pipl_1p2pi_fracfeed=(TH1F*)	h_Etot_subtruct_pipl_41prot_fracfeed->Clone("h_Etot_subtruct_pipl_1p2pi_fracfeed");
	h_Etot_subtruct_pipl_1p2pi_fracfeed->Add(h1_E_tot_1p2pi_fracfeed_pipl);

	TH1F *h_Erec_subtruct_pimi_1p2pi=(TH1F*)	h_Erec_subtruct_pimi_41prot->Clone("h_Erec_subtruct_pimi_1p2pi");
	h_Erec_subtruct_pimi_1p2pi->Add(h1_E_rec_1p2pi_pimi);

	TH1F *h_Etot_subtruct_pimi_1p2pi=(TH1F*)	h_Etot_subtruct_pimi_41prot->Clone("h_Etot_subtruct_pimi_1p2pi");
	h_Etot_subtruct_pimi_1p2pi->Add(h1_E_tot_1p2pi_pimi);

	TH2F *h2_Erec_pperp_sub_1p2pi_1p1pi_pimi=(TH2F*) h2_Erec_pperp_41psub->Clone("h2_Erec_pperp_sub_1p2pi_1p1pi_pimi");
	h2_Erec_pperp_sub_1p2pi_1p1pi_pimi->Add(h2_Erec_pperp_1p2pi_1p1pi_pimi);

	TH1F *h_Erec_subtruct_pimi_1p2pi_fracfeed=(TH1F*)	h_Erec_subtruct_pimi_41prot_fracfeed->Clone("h_Erec_subtruct_pimi_1p2pi_fracfeed");
	h_Erec_subtruct_pimi_1p2pi_fracfeed->Add(h1_E_rec_1p2pi_fracfeed_pimi);

	TH1F *h_Etot_subtruct_pimi_1p2pi_fracfeed=(TH1F*)	h_Etot_subtruct_pimi_41prot_fracfeed->Clone("h_Etot_subtruct_pimi_1p2pi_fracfeed");
	h_Etot_subtruct_pimi_1p2pi_fracfeed->Add(h1_E_tot_1p2pi_fracfeed_pimi);*/

	 //-----------------------------------looking only at e-, 2pi undetected pion subtraction  ---------------------------------------
/* Commented out bc dependent on 0pi Histograms
	TH1F *h_Erec_subtruct_piplpimi_noprot_2pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot->Clone("h_Erec_subtruct_piplpimi_noprot_2pi");
	h_Erec_subtruct_piplpimi_noprot_2pi->Add(h1_E_rec_2pi_weight);

	TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed2pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_frac_feed->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed2pi");
	h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Add(h1_E_rec_2pi_weight_frac_feed);

	 //-----------------------------------looking only at e-, 3pi, undetected pion subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_noprot_3pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_2pi->Clone("h_Erec_subtruct_piplpimi_noprot_3pi");
	h_Erec_subtruct_piplpimi_noprot_3pi->Add(h1_E_rec_3pi_weight);

	TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed3pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_frac_feed2pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed3pi");
	h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Add(h1_E_rec_3pi_weight_frac_feed);

	 //-----------------------------------looking only at e-, 4pi, undetected pion subtraction  ---------------------------------------

	TH1F *h_Erec_subtruct_piplpimi_noprot_4pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_3pi->Clone("h_Erec_subtruct_piplpimi_noprot_4pi");
	h_Erec_subtruct_piplpimi_noprot_4pi->Add(h1_E_rec_4pi_weight);

	TH1F *h_Erec_subtruct_piplpimi_noprot_frac_feed4pi = (TH1F*)	h_Erec_subtruct_piplpimi_noprot_frac_feed3pi->Clone("h_Erec_subtruct_piplpimi_noprot_frac_feed4pi");
	h_Erec_subtruct_piplpimi_noprot_frac_feed4pi->Add(h1_E_rec_4pi_weight_frac_feed);*/

	gDirectory->Write("hist_Files", TObject::kOverwrite);
	// skim_tree->AutoSave();

	// --------------------------------------------------------------------------------------------------------

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Total Number of Events: " << nentries << std::endl;
	std::cout << std::endl << "One Pion Events = " << std::endl;
	std::cout << std::endl << "1p1pi ALL: " << C1p1piALL << std::endl;
	std::cout << std::endl << "1p1pi PIPL: " << C1p1piPIPL << std::endl;
	std::cout << std::endl << "1p1pi PIMI: " << C1p1piPIMI << std::endl;
	std::cout << std::endl << "2p1pi ALL: " <<C2p1piALL << std::endl;
	std::cout << std::endl << "2p1pi PIPL: " <<C2p1piPIPL << std::endl;
	std::cout << std::endl << "2p1pi PIMI: " <<C2p1piPIMI << std::endl;
	std::cout << std::endl << "3p1pi ALL: " <<C3p1piALL << std::endl;
	std::cout << std::endl << "3p1pi PIPL: " <<C3p1piPIPL << std::endl;
	std::cout << std::endl << "3p1pi PIMI: " <<C3p1piPIMI << std::endl;
	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Two Pion Events = " << std::endl;
	std::cout << std::endl << "1p2pi ALL: " <<C1p2piALL << std::endl;
	std::cout << std::endl << "1p2pi PIPL: " <<C1p2piPIPL << std::endl;
	std::cout << std::endl << "1p2pi PIMI: " <<C1p2piPIMI << std::endl;
	std::cout << std::endl << "2p2pi ALL: " << C2p2piALL << std::endl;
	std::cout << std::endl << "2p2pi PIPL: " << C2p2piPIPL << std::endl;
	std::cout << std::endl << "2p2pi PIMI: " <<C2p2piPIMI << std::endl;
	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Three Pion Events = " << std::endl;
	std::cout << std::endl << "1p3pi ALL: " <<C1p3piALL << std::endl;
	std::cout << std::endl << "1p3pi PIPL: " << C1p3piPIPL << std::endl;
	std::cout << std::endl << "1p3pi PIMI: " << C1p3piPIMI << std::endl;
	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "Initial # Events = " << fChain->GetEntries() << std::endl;
	std::cout << std::endl << "1e1p0pi Signal # Events = " << SignalEvents << std::endl;
	std::cout << std::endl << "Passing Rate = " << int(double(SignalEvents) / double(fChain->GetEntries())*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "PMiss Fraction 1st Slice = " << int(double(PMiss_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "PMiss Fraction 2nd Slice = " << int(double(PMiss_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "PMiss Fraction 3rd Slice = " << int(double(PMiss_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Events With ECal Within 5\% of ETrue = " << ECalSignalEventsWithin5Perc << std::endl;
	std::cout << std::endl << "ECal 5% Fraction = " << int(double(ECalSignalEventsWithin5Perc) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 1st Slice = " << int(double(ECalSignalEventsWithin5Perc_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 2nd Slice = " << int(double(ECalSignalEventsWithin5Perc_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "ECal 5% Fraction 3rd Slice = " << int(double(ECalSignalEventsWithin5Perc_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << std::endl << "# Events With EQE Within 5\% of ETrue = " << EQESignalEventsWithin5Perc << std::endl;
	std::cout << std::endl << "EQE 5% Fraction = " << int(double(EQESignalEventsWithin5Perc) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 1st Slice = " << int(double(EQESignalEventsWithin5Perc_FirstSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 2nd Slice = " << int(double(EQESignalEventsWithin5Perc_SecondSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;
	std::cout << std::endl << "EQE 5% Fraction 3rd Slice = " << int(double(EQESignalEventsWithin5Perc_ThirdSlice) / double(SignalEvents)*100.) << " \%"<< std::endl << std::endl;

	if (choice == 1) {

		std::cout << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "MEC Fractional Contribution = " << int(double(MECSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "RES Fractional Contribution = " << int(double(RESSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "DIS Fractional Contribution = " << int(double(DISSignalEvents) / double(SignalEvents)*100.) << " \%" << std::endl;
		std::cout << std::endl << "-----------------------------------------------------------------------------------------------------" << std::endl;

	}
	
}

//End Loop function


double genie_analysis::acceptance_c(double p, double cost, double phi, int particle_id,TFile* file_acceptance)
{

	//Redefinition of the phi angle
	// because the acceptance maps are defined between (-30,330)

	// Check that phi is between (0,360)

	//int redef = -30;
	int redef = 0;

	TH3D * acc;
	TH3D * gen;

	acc = (TH3D*)file_acceptance->Get("Accepted Particles");
	gen = (TH3D*)file_acceptance->Get("Generated Particles");

	//map 330 till 360 to [-30:0] for the acceptance map histogram
	if(phi > (2*TMath::Pi() - TMath::Pi()/6.) ) { phi -= 2*TMath::Pi(); }
	//Find number of generated events

	double pbin_gen = gen->GetXaxis()->FindBin(p);
	double tbin_gen = gen->GetYaxis()->FindBin(cost);
	double phibin_gen = gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	double num_gen = gen->GetBinContent(pbin_gen, tbin_gen, phibin_gen);

	//Find number of accepted events

	double pbin_acc = acc->GetXaxis()->FindBin(p);
	double tbin_acc = acc->GetYaxis()->FindBin(cost);
	double phibin_acc = acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	double num_acc = acc->GetBinContent(pbin_acc, tbin_acc, phibin_acc);

	double acc_ratio = (double)num_acc / (double)num_gen;
	double acc_err = (double)sqrt(acc_ratio*(1-acc_ratio)) / (double)num_gen;


	return acc_ratio;

}
