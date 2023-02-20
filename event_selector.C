#define event_selector_cxx
#include "event_selector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TString.h>
#include <TVector3.h>
#include <TF1.h>
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

#include "Constants.h"

using namespace std; 

void event_selector::Loop() {
    //start
    std::cout << "Beginning Loop..." << std::endl;

    Long64_t nbytes = 0, nb = 0; 

    std::map<std::string, double> en_beam; 
    std::map<std::string, double> Ecal_offset; 
    std::map<std::string, double> target_mass; 
    std::map<std::string, double> residual_target_mass; 
    
    
    en_beam["1161"] = 1.161; 
    en_beam["2261"] = 2.261; 
    en_beam["4461"] = 4.461;

    int numSelectedEvents = 0; 

    if (fChain == 0) return; 

    /**
     * @brief Data file setup and headers
     */
    std::ofstream dataFile; 
    dataFile.open("data1p2pi.csv");
    dataFile << "Event Type, Detected Particle charge, Same/Diff Charges (true), E0, Erec, Ee', pex, pey, pez, Ep, ppx, ppy, ppz \n"; 

    Long64_t nentries = fChain->GetEntries();
    std::cout << "Number of Entries: " << nentries << std::endl;


    /**
     * @brief Event loop, loop over the entries in the file
     */
    for(Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		//Read Entry
		int nb = GetEntry(jentry);
		if (nb == 0) { std::cout <<"Event loop: 0 byte read for entry " << jentry << ". Indicate failure in reading the file" <<	std::endl;}

		if (jentry%10000 == 0) {std::cout << jentry/1000 << " k " << std::setprecision(3) << double(jentry)/fChain->GetEntries()*100. << " %"<< std::endl;}

        TLorentzVector V4_el(pxl,pyl,pzl, El);
        TVector3 V3_el(pxl,pyl,pzl);
        double el_momentum = V3_el.Mag();
        double el_theta = V3_el.Theta(); 

        //Electron cuts
        if((el_theta * TMath::RadToDeg()) < 15.0)
        {
            continue;
        }
        //Cuts on electron momentum
        if (fbeam_en == "1161" && el_momentum < 0.4) {continue;}
        if (fbeam_en=="2261" && el_momentum < 0.55) { continue; }
		if (fbeam_en=="4461" && el_momentum < 1.1) { continue; }

        int true_ProtonCounter = 0; 
        int true_PionCounter = 0;
        int true_PiplCounter = 0;
        int true_PimiCounter = 0;
        int det_ProtonCounter = 0; 
        int det_PionCounter = 0;
        int det_PiplCounter = 0; 
        int det_PimiCounter = 0; 

        vector<int> IndexProton; 
        vector<int> IndexPion;
        /**
         * @brief Looping for hadrons
         */
        for(int i = 0; i < nf; i++) 
        {
            //Proton Selection
            if(pdgf[i] == 2212 && pf[i] > 0.3)
            {
                true_ProtonCounter++; 

                TVector3 V3_prot(pxf[i], pyf[i], pzf[i]);
                if((V3_prot.Theta() * TMath::RadToDeg()) < 10) {continue;}
                if(V3_prot.Mag() < 0.3) {continue;}
                det_ProtonCounter++;
                IndexProton.push_back(i);
            }
            if(pdgf[i] == -211 && pf[i] > 0.15)
            {
                //Pi minus subtraction 
                true_PionCounter++;
                true_PimiCounter++;
                TVector3 V3_pion(pxf[i],pyf[i], pzf[i]);
                if((V3_pion.Theta() *TMath::RadToDeg()) < 15) {continue;}
                if(V3_pion.Mag() < 0.15) {continue;}
                det_PionCounter++;
                det_PimiCounter++;
                IndexPion.push_back(i);
            }
            if(pdgf[i] == 211 && pf[i] > 0.15)
            {
                //Pi plus detection
                true_PionCounter++;
                true_PiplCounter++;
                TVector3 V3_pion(pxf[i], pyf[i], pzf[i]);
                if((V3_pion.Theta() * TMath::RadToDeg()) < 15.0) {continue;}
                if(V3_pion.Mag() < 0.15) {continue;}
                det_PionCounter++;
                det_PiplCounter++;
                IndexPion.push_back(i);
            }

        }

        if(true_ProtonCounter == 1 && true_PionCounter == 2 && det_ProtonCounter == 1 && det_PionCounter == 1)
        {
            double E_cal = Ef[IndexProton[0]] + El - m_prot;
            int charge = 0; 
            string chargeType = "";
            if(det_PiplCounter == 1)
                charge = 1; 
            else if(det_PimiCounter == 1)
                charge = -1;
            if(true_PimiCounter == true_PiplCounter == 1)
                chargeType = "Different Charges";
            else if((true_PimiCounter == 2 && true_PiplCounter == 0) || (true_PiplCounter == 2 && true_PimiCounter == 0))
                chargeType = "Same Charges";
            dataFile << true_ProtonCounter << "p" << true_PiplCounter << "pipl" << true_PimiCounter << "pimi," << charge << "," << chargeType << "," << Ev << "," <<  E_cal << "," << El << "," << pxl << "," << pyl << "," << pzl << "," << Ef[IndexProton[0]] << "," << pxf[IndexProton[0]] << "," << pyf[IndexProton[0]] << "," << pzf[IndexProton[0]] << "\n";
        }


    }

    dataFile.close();
}