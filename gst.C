#define gst_cxx
#include "gst.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TVector3.h>
#include <TF1.h>

#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

#include "Constants.h"

using namespace std;

void gst::Loop() {

	std::cout<<"Beginning Loop..."<<std::endl;

	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;

	// - - - - - - Getting number of bins for the histograms

	//Getting number of bins for stuff
	int n_bins;
	double *x_values;
	double *x_qe;

	std::map<std::string,double> en_beam;
	//std::map<std::string,double>bind_en;
	std::map<std::string, double> Ecal_offset;
/*
	Ecal_offset["3He"]=0.004;
	Ecal_offset["4He"]=0.005;
	Ecal_offset["C12"]=0.005;
	Ecal_offset["56Fe"]=0.011;

	bind_en["3He"] = He3_bind_en-D2_bind_en + Ecal_offset["3He"]; //the offset is used to shift the peak to be at 0
	bind_en["4He"] = He4_bind_en-H3_bind_en + Ecal_offset["4He"];
	bind_en["C12"] = C12_bind_en-B_bind_en	+ Ecal_offset["C12"];
	bind_en["56Fe"]= Fe_bind_en-Mn_bind_en	+ Ecal_offset["56Fe"];
	bind_en["CH2"] = C12_bind_en-B_bind_en;*/

	en_beam["1161"]=1.161;
	en_beam["2261"]=2.261;
	en_beam["4461"]=4.461;

	if(en_beam[fbeam_en]>1. && en_beam[fbeam_en]<2.){
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=17;i++) { x_values[i]=0.4+i*0.04; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=20;i++) { x_values[i+18]=1.08+(i+1)*0.02; x_qe[i+18] = (x_values[i+18] - en_beam[fbeam_en]) / en_beam[fbeam_en]; }
	}

	if(en_beam[fbeam_en]>2. && en_beam[fbeam_en]<3.){
		n_bins=54;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=23;i++) { x_values[i]=i*0.09; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=30;i++) { x_values[i+24]=2.07+(i+1)*0.03; x_qe[i+24] = (x_values[i+24] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}

	if(en_beam[fbeam_en]>4. && en_beam[fbeam_en]<5.){
		n_bins=38;
		x_values=new double[n_bins+1]; x_qe=new double[n_bins+1];
		for (int i=0;i<=21;i++)	{ x_values[i]=i*0.2; x_qe[i] = (x_values[i] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
		for (int i=0;i<=16;i++)	{ x_values[i+22]=4.2+(i+1)*0.05; x_qe[i+22] = (x_values[i+22] - en_beam[fbeam_en]) / en_beam[fbeam_en];}
	}
	//

// - - - - - Number of bins obtained



	// --------------------------------------------------------------

	TString FileName = Form("FilteredGenie_%s_%sGeV_trueEvents.root",ftarget.c_str(), fbeam_en.c_str());
	//TString FileName = Form("/w/hallb-scshelf2102/clas/claseg2/apapadop/eresmaid_%s_%s_hA2018_LFG_FSI_NoRadCorr_3M.root", ftarget.c_str(), fbeam_en.c_str());
	TFile* FilteredFile = new TFile(FileName,"recreate");
	TTree* maketreenew = fChain->GetTree()->CloneTree(0);
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
	// --------------------------------------------------------------

	int SelectedEvents = 0;

	// --------------------------------------------------------------

	TF1 *myElectronFit = new TF1("myElectronFit","[0]+[1]/x",0.,5.);
	myElectronFit->SetParameters(13.5,15);

	TF1* myPiMinusFit = new TF1("myPiMinusFit","(x<0.35)*(25.+7./TMath::Power(x,1.)) + (x>0.35)*(16.+10/TMath::Power(x,1.))",0,5.);

	// --------------------------------------------------------------

//	for (Long64_t jentry=0; jentry<nentries;jentry++) {
	for (Long64_t jentry=0; jentry<10000000;jentry++) {

		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

		int Interaction = -1;
		if (qel) { Interaction = 1; }
		if (mec) { Interaction = 2; }
		if (res) { Interaction = 3; }
		if (dis) { Interaction = 4; }

		//- - - - - Creation of Vectors - - - - -//
		 vector<int> PionIndexCounter;
		 vector<int> ProtonIndexCounter;
		 vector<int> piplIndexCounter;
		 vector<int> pimiIndexCounter;
		// --------------------------------------------------------------

		if (jentry%1000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/fChain->GetEntries()*100. << " %"<< std::endl;}

		// --------------------------------------------------------------

		// Electron threshods for 4.4 GeV

		if (El < 1.1) { continue; }
		if (Q2 < 0.8) { continue; }
		if (W > 2) { continue; }

		TVector3 ElectronV3(pxl,pyl,pzl);
		double E_el = El;
		double el_momentum = ElectronV3.Mag();
		double el_theta = ElectronV3.Theta();

		double theta_min = myElectronFit->Eval(el_momentum);
		if (el_theta*180./TMath::Pi() < theta_min) { continue; }

		// --------------------------------------------------------------

		int ProtonCounter = 0;
		int ChargedPionCounter = 0;
		int PiPlusCounter = 0;
		int PiMinusCounter = 0;
		int GammaCounter = 0;
		//Loop for Hadrons

		for (int i = 0; i < nf; i++) {

			double theta = TMath::ACos(cthf[i]) * 180./TMath::Pi();

			if (pdgf[i] == 2212  && pf[i] > 0.3 && theta > MinThetaProton)
			{
				ProtonCounter++;
				ProtonIndexCounter.push_back(i);
			}
			if (pdgf[i] == 211  && pf[i] > 0.15 && theta > MinThetaPiPlus)
			{
				PiPlusCounter++;
				ChargedPionCounter++;
				piplIndexCounter.push_back(i);
				PionIndexCounter.push_back(i);
			}
			if (pdgf[i] == -211  && pf[i] > 0.15) {

				TVector3 PiMinusV3(pxf[i],pyf[i],pzf[i]);
				double piminus_momentum = PiMinusV3.Mag();
				double piminus_theta = PiMinusV3.Theta();

				double piminus_theta_min = myPiMinusFit->Eval(piminus_momentum);
				if (piminus_theta*180./TMath::Pi() > piminus_theta_min)
				{

					PiMinusCounter++;
					ChargedPionCounter++;
					pimiIndexCounter.push_back(i);
					PionIndexCounter.push_back(i);

				}
			}

			if (pdgf[i] == 22  && pf[i] > 0.3 && theta > MinThetaGamma) { GammaCounter++; }


		}

		// --------------------------------------------------------------

                // Chosen topology

                if (
                    !(
											(ProtonCounter == 1 & ChargedPionCounter == 1)
/*
                        (ProtonCounter == 1 && ChargedPionCounter == 0)
                     || (ProtonCounter == 1 && ChargedPionCounter == 1)
                     || (ProtonCounter == 1 && ChargedPionCounter == 2)
                     || (ProtonCounter == 1 && ChargedPionCounter == 3)
                     || (ProtonCounter == 2 && ChargedPionCounter == 0)
                     || (ProtonCounter == 2 && ChargedPionCounter == 1)
                     || (ProtonCounter == 1 && ChargedPionCounter == 2)
                     || (ProtonCounter == 3 && ChargedPionCounter == 0)
                     || (ProtonCounter == 3 && ChargedPionCounter == 1)*/

                    )

                ) { continue; }

					//Energy of prot -> Ef[ProtonIndexCounter[0]]
					//Energy of electron -> V4_el.E();
					//Energy of pion -> Ef[PionIndexCounter[0]]
					double Ecal = 0.;
					if(PiMinusCounter == 1)
					{
						if(pimiIndexCounter.size() != 1 || ProtonIndexCounter.size() != 1)
						{
							std::cout<<"This should not happen. Sanity check failed."<<std::endl;
						}
						Ecal = Ef[ProtonIndexCounter[0]] + E_el + Ef[pimiIndexCounter[0]] - m_prot;
						h1_E_cal_1p1pi_pimi_TRUE->Fill(Ecal);
						switch (Interaction) {
							case 1:
								h1_E_cal_1p1pi_pimi_TRUE_qel->Fill(Ecal);
								break;
							case 2:
								h1_E_cal_1p1pi_pimi_TRUE_mec->Fill(Ecal);
								break;
							case 3:
								h1_E_cal_1p1pi_pimi_TRUE_res->Fill(Ecal);
								break;
							case 4:
								h1_E_cal_1p1pi_pimi_TRUE_dis->Fill(Ecal);
								break;
							default:
								break;
						}
					}
					else if(PiPlusCounter == 1)
					{
						if(piplIndexCounter.size() != 1)
						{
							std::cout<<"This should not happen. Sanity check failed."<<std::endl;
						}
						Ecal = Ef[ProtonIndexCounter[0]] + E_el + Ef[piplIndexCounter[0]] - m_prot;
						h1_E_cal_1p1pi_pipl_TRUE->Fill(Ecal);
						switch (Interaction) {
							case 1:
								h1_E_cal_1p1pi_pipl_TRUE_qel->Fill(Ecal);
								break;
							case 2:
								h1_E_cal_1p1pi_pipl_TRUE_mec->Fill(Ecal);
								break;
							case 3:
								h1_E_cal_1p1pi_pipl_TRUE_res->Fill(Ecal);
								break;
							case 4:
								h1_E_cal_1p1pi_pipl_TRUE_dis->Fill(Ecal);
								break;
							default:
								break;
						}
					}
		// --------------------------------------------------------------

		// 1p0pi + 2p0pi
//		if ( !( (ProtonCounter == 1 && ChargedPionCounter == 0) || (ProtonCounter == 2 && ChargedPionCounter == 0) ) ) { continue; }
//		if (GammaCounter != 0) { continue; }

		// --------------------------------------------------------------

		SelectedEvents++;
		maketreenew->Fill();

		//Make Histogram for 1p1piTrue stuff
		// --------------------------------------------------------------


   } // End of the loop over the events

   std::cout << "Selected events = " << SelectedEvents << std::endl;
   std::cout << "Created file " << FileName << std::endl;

   maketreenew->Write();
	 FilteredFile->Write("hist_files",TObject::kOverwrite);
   FilteredFile->Close();

} // End of the Loop()
