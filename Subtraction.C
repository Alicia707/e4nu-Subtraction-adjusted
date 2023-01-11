#ifndef SUBTRACTION_CXX
#define SUBTRACTION_CXX

#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVectorT.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TGraph.h>
#include <vector>
#include "Subtraction.h"

void Subtraction::prot1_pi2_rot_func(TVector3  V3prot, TVector3 V3pi[2], TLorentzVector V4prot, TLorentzVector V4pi[2], int q_pi[2], TLorentzVector V4_el, double Ecal[3], double p_miss_perp[3], double P_1p1pi[], int targetCharge)
{
  const int N_pi = 2;
  double rotation_ang = 0.0;
  TVector3 V3_rot_pi[2], V3_p_rot;
  bool status_pi[2] = {true};
  bool status_prot = true;

  double N_all = 0;
  double N_1p1pi[2] = {0};
  double N_1p1pi_diff[2] = {0}; //for when there is 1 pimi and 1 pipl

  std::vector<int> piplIndexCounter;
  std::vector<int> pimiIndexCounter;
  int PiPlusCounter = 0;
  int PiMinusCounter = 0;

  for(int i = 0; i < N_pi; i++)
  {
    if(q_pi[i] > 0)
    {
      piplIndexCounter.push_back(i);
      PiPlusCounter++;
    }
    else if(q_pi[i] < 0)
    {
      pimiIndexCounter.push_back(i);
      PiMinusCounter++;
    }
  }

  //Check if putting for inside if statement for combo check speeds up the code at all (use linux "time" command)
  for(int g = 0; g < N_tot; g++) // Get number of charges, set up vectors, etc.
  {
    rotation_ang = gRandom->Uniform(0,2*TMath::Pi());
    V3_p_rot = V3prot;

    V3_p_rot.Rotate(rotation_ang, V3q);

    status_prot = PFiducialCut(fbeam_en, V3_p_rot);

    for(int i = 0; i < N_pi; i++)// Get number of charges, set up vectors, etc.
    {
      V3_rot_pi[i]=V3pi[i];
      V3_rot_pi[i].Rotate(rotation_ang, V3q);
      status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
    }

    if(status_prot  && status_pi[0]  && status_pi[1] ) N_all++; //both pion det

    if(targetCharge > 0 && PiPlusCounter == 2 && PiMinusCounter == 0)
    {
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  && !status_pi[piplIndexCounter[1]] ) N_1p1pi[0]++; //1st pion det
      if(status_prot  && !status_pi[piplIndexCounter[0]]  &&  status_pi[piplIndexCounter[1]] ) N_1p1pi[1]++; //2nd pion det
    }
    else if(targetCharge < 0 && PiPlusCounter == 0 && PiMinusCounter == 2)
    {
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  && !status_pi[pimiIndexCounter[1]] ) N_1p1pi[0]++; //1st pion det
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  &&  status_pi[pimiIndexCounter[1]] ) N_1p1pi[1]++; //2nd pion det
    }
    else if(PiPlusCounter == 1 && PiMinusCounter == 1)
    {
      if(status_prot &&  status_pi[piplIndexCounter[0]] && !status_pi[pimiIndexCounter[0]]) N_1p1pi_diff[piplIndexCounter[0]]++; //the pipl is det
      if(status_prot && !status_pi[piplIndexCounter[0]] &&  status_pi[pimiIndexCounter[0]]) N_1p1pi_diff[pimiIndexCounter[0]]++; //the pimi is det
    }
    else if(targetCharge == 0) //Previous logic to avoid seg fault errors 6.29.22 -> Previous logic did not account for detection of diff pion charges
    {
      if(status_prot  &&  status_pi[0]  && !status_pi[1] ) N_1p1pi[0]++;
      if(status_prot  && !status_pi[0]  &&  status_pi[1] ) N_1p1pi[1]++;
    }
  }
  //Check if this is to contribute to pipl or pimi histos
  if(targetCharge == 1) // contribute to pipl
  {
    if(PiPlusCounter == 2 && PiMinusCounter == 0)
    {
      for(int h = 0; h<N_pi; h++) //Go over two pions that are the same (both pipl)
      {
        prot1_pi1_en_calc(V4prot, V4pi[h], q_pi[h], V4_el, &Ecal[h], &p_miss_perp[h]);
        if(N_all!=0)
        {
          P_1p1pi[h] = -(N_1p1pi[h]/N_all);
        }
        else
        P_1p1pi[h]=0; //N_all!=0 statement
      }
    }
    else if(PiPlusCounter == 1 && PiMinusCounter == 1)
    {
      prot1_pi1_en_calc(V4prot, V4pi[piplIndexCounter[0]], q_pi[piplIndexCounter[0]], V4_el, &Ecal[2], &p_miss_perp[2]);
      if(N_all!=0)
      {
        P_1p1pi[2] = -(N_1p1pi_diff[piplIndexCounter[0]]/N_all);
      }
      else
        P_1p1pi[2] = 0;
    }
  }

  else if(targetCharge == -1)
  {
    if(PiPlusCounter == 0 && PiMinusCounter == 2)
    {
      for(int h = 0; h<N_pi; h++)
      {
        prot1_pi1_en_calc(V4prot, V4pi[h], q_pi[h], V4_el, &Ecal[h], &p_miss_perp[h]);
        if(N_all!=0)
        P_1p1pi[h] = -(N_1p1pi[h]/N_all);
        else
        P_1p1pi[h]=0; //N_all!=0 statement
      }
    }
    else if(PiPlusCounter == 1 && PiMinusCounter == 1)
    {
      prot1_pi1_en_calc(V4prot, V4pi[pimiIndexCounter[0]], q_pi[pimiIndexCounter[0]], V4_el, &Ecal[2], &p_miss_perp[2]);
      if(N_all!=0)
        P_1p1pi[2] = -(N_1p1pi_diff[pimiIndexCounter[0]]/N_all);
      else
        P_1p1pi[2] = 0;
    }
  }
  else if(targetCharge == 0)
  {
    //----------------------1p2pi->1p1pi
    for(int h=0;h<N_pi;h++)
    {
      prot1_pi1_en_calc(V4prot, V4pi[h], q_pi[h], V4_el, &Ecal[h], &p_miss_perp[h]);
      if(N_all!=0)
      P_1p1pi[h] = -(N_1p1pi[h]/N_all);
      else
      P_1p1pi[h]=0; //N_all!=0 statement
    }
  }
  else
  {
    std::cout << "This should not happen!" << std::endl;
  }
}


void Subtraction::prot1_pi3_rot_func(TVector3  V3prot, TVector3 V3pi[3], TLorentzVector V4prot, TLorentzVector V4pi[3], int q_pi[3], TLorentzVector V4_el, double Ecal[3], double p_miss_perp[3], double P_tot[3], int targetCharge)
{
  const int N_pi = 3;
  double P_1p3pito1p1pi[3] = {0};
  double rotation_ang;
  TVector3 V3_rot_pi[N_pi], V3_p_rot;
  bool status_pi[N_pi]={true};
  bool status_prot = true;
  double N_all = 0;
  double N_1p1pi[3]={0},N_1p2pi[3]={0}, N_1p1pi_diff[3] = {0};
  double N_1p2pi_diff[2][2] = {0.0};


  //Imeplmenting charge selection
  std::vector<int> piplIndexCounter; 
  std::vector<int> pimiIndexCounter;
  int PiPlusCounter = 0; 
  int PiMinusCounter = 0;
  //Implementing charge selection
  for(int i = 0; i < N_pi; i++)
  {
    if(q_pi[i] == 1)
    {
      piplIndexCounter.push_back(i); 
      PiPlusCounter++;
    }
    else if(q_pi[i] == -1)
    {
      pimiIndexCounter.push_back(i); 
      PiMinusCounter++;
    }
  }

  for(int g=0; g < N_tot; g++)
  {
    rotation_ang=gRandom->Uniform(0,2*TMath::Pi());
    V3_p_rot = V3prot;

    V3_p_rot.Rotate(rotation_ang,V3q);
    
    status_prot = PFiducialCut(fbeam_en, V3_p_rot);

    for(int i=0;i<N_pi;i++)
    {
      V3_rot_pi[i]=V3pi[i];
      V3_rot_pi[i].Rotate(rotation_ang,V3q);
      status_pi[i]=Pi_phot_fid_united(fbeam_en, V3_rot_pi[i],q_pi[i]);
    }
    //All pions are detected
    if(status_prot && status_pi[0] && status_pi[1] && status_pi[2]) N_all++;
    //Case that there are 3 pions that are pipl
    if(targetCharge == 1 && PiPlusCounter == 3)
    {
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  && !status_pi[piplIndexCounter[1]] && !status_pi[piplIndexCounter[2]]) N_1p1pi[0]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  &&  status_pi[piplIndexCounter[1]] && !status_pi[piplIndexCounter[2]]) N_1p1pi[1]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  && !status_pi[piplIndexCounter[1]] &&  status_pi[piplIndexCounter[2]]) N_1p1pi[2]++;
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  &&  status_pi[piplIndexCounter[1]] && !status_pi[piplIndexCounter[2]]) N_1p2pi[0]++;
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  && !status_pi[piplIndexCounter[1]] &&  status_pi[piplIndexCounter[2]]) N_1p2pi[1]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  &&  status_pi[piplIndexCounter[1]] &&  status_pi[piplIndexCounter[2]]) N_1p2pi[2]++;
    }
    else if(targetCharge == 1 && PiPlusCounter == 1 && PiMinusCounter == 2)
    {
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  && !status_pi[pimiIndexCounter[0]] && !status_pi[pimiIndexCounter[1]]) N_1p1pi_diff[piplIndexCounter[0]]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  &&  status_pi[pimiIndexCounter[0]] && !status_pi[pimiIndexCounter[1]]) N_1p1pi_diff[pimiIndexCounter[0]]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  && !status_pi[pimiIndexCounter[0]] &&  status_pi[pimiIndexCounter[1]]) N_1p1pi_diff[pimiIndexCounter[1]]++;
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  &&  status_pi[pimiIndexCounter[0]] && !status_pi[pimiIndexCounter[1]]) N_1p2pi_diff[0][1]++;
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  && !status_pi[pimiIndexCounter[0]] &&  status_pi[pimiIndexCounter[1]]) N_1p2pi_diff[0][1]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  &&  status_pi[pimiIndexCounter[0]] &&  status_pi[pimiIndexCounter[1]]) N_1p2pi_diff[1][1]++;
    }
    else if(targetCharge == 1 && PiPlusCounter == 2 && PiMinusCounter == 1)
    {
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  && !status_pi[piplIndexCounter[1]] && !status_pi[pimiIndexCounter[0]]) N_1p1pi_diff[piplIndexCounter[0]]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  &&  status_pi[piplIndexCounter[1]] && !status_pi[pimiIndexCounter[0]]) N_1p1pi_diff[piplIndexCounter[1]]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  && !status_pi[piplIndexCounter[1]] &&  status_pi[pimiIndexCounter[0]]) N_1p1pi_diff[pimiIndexCounter[0]]++;
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  &&  status_pi[piplIndexCounter[1]] && !status_pi[pimiIndexCounter[0]]) N_1p2pi_diff[0][0]++;
      if(status_prot  &&  status_pi[piplIndexCounter[0]]  && !status_pi[piplIndexCounter[1]] &&  status_pi[pimiIndexCounter[0]]) N_1p2pi_diff[0][1]++;
      if(status_prot  && !status_pi[piplIndexCounter[0]]  &&  status_pi[piplIndexCounter[1]] &&  status_pi[pimiIndexCounter[0]]) N_1p2pi_diff[0][1]++;
    }
    else if(targetCharge == -1 && PiMinusCounter == 3)
    {
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  && !status_pi[pimiIndexCounter[1]] && !status_pi[pimiIndexCounter[2]]) N_1p1pi[0]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  &&  status_pi[pimiIndexCounter[1]] && !status_pi[pimiIndexCounter[2]]) N_1p1pi[1]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  && !status_pi[pimiIndexCounter[1]] &&  status_pi[pimiIndexCounter[2]]) N_1p1pi[2]++;
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  &&  status_pi[pimiIndexCounter[1]] && !status_pi[pimiIndexCounter[2]]) N_1p2pi[0]++;
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  && !status_pi[pimiIndexCounter[1]] &&  status_pi[pimiIndexCounter[2]]) N_1p2pi[1]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  &&  status_pi[pimiIndexCounter[1]] &&  status_pi[pimiIndexCounter[2]]) N_1p2pi[2]++;
    }
    else if(targetCharge == -1 && PiPlusCounter == 2 && PiMinusCounter == 1)
    {
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  && !status_pi[piplIndexCounter[0]] && !status_pi[piplIndexCounter[1]]) N_1p1pi_diff[pimiIndexCounter[0]]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  &&  status_pi[piplIndexCounter[0]] && !status_pi[piplIndexCounter[1]]) N_1p1pi_diff[piplIndexCounter[0]]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  && !status_pi[piplIndexCounter[0]] &&  status_pi[piplIndexCounter[1]]) N_1p1pi_diff[piplIndexCounter[1]]++;
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  &&  status_pi[piplIndexCounter[0]] && !status_pi[piplIndexCounter[1]]) N_1p2pi_diff[1][0]++;
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  && !status_pi[piplIndexCounter[0]] &&  status_pi[piplIndexCounter[1]]) N_1p2pi_diff[1][0]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  &&  status_pi[piplIndexCounter[0]] &&  status_pi[piplIndexCounter[1]]) N_1p2pi_diff[0][0]++;
    }
    else if(targetCharge == 1 && PiPlusCounter == 1 && PiMinusCounter == 2)
    {
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  && !status_pi[pimiIndexCounter[1]] && !status_pi[piplIndexCounter[0]]) N_1p1pi_diff[pimiIndexCounter[0]]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  &&  status_pi[pimiIndexCounter[1]] && !status_pi[piplIndexCounter[0]]) N_1p1pi_diff[pimiIndexCounter[1]]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  && !status_pi[pimiIndexCounter[1]] &&  status_pi[piplIndexCounter[0]]) N_1p1pi_diff[piplIndexCounter[0]]++;
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  &&  status_pi[pimiIndexCounter[1]] && !status_pi[piplIndexCounter[0]]) N_1p2pi_diff[1][1]++;
      if(status_prot  &&  status_pi[pimiIndexCounter[0]]  && !status_pi[pimiIndexCounter[1]] &&  status_pi[piplIndexCounter[0]]) N_1p2pi_diff[1][0]++;
      if(status_prot  && !status_pi[pimiIndexCounter[0]]  &&  status_pi[pimiIndexCounter[1]] &&  status_pi[piplIndexCounter[0]]) N_1p2pi_diff[1][0]++;
    }
    else if(targetCharge == 0) // Prevents seg fault errors
    {
      if(status_prot  &&  status_pi[0]  && !status_pi[1] && !status_pi[2]) N_1p1pi[0]++;
      if(status_prot  && !status_pi[0]  &&  status_pi[1] && !status_pi[2]) N_1p1pi[1]++;
      if(status_prot  && !status_pi[0]  && !status_pi[1] &&  status_pi[2]) N_1p1pi[2]++;
      if(status_prot  &&  status_pi[0]  &&  status_pi[1] && !status_pi[2]) N_1p2pi[0]++;
      if(status_prot  &&  status_pi[0]  && !status_pi[1] &&  status_pi[2]) N_1p2pi[1]++;
      if(status_prot  && !status_pi[0]  &&  status_pi[1] &&  status_pi[2]) N_1p2pi[2]++;
    }

  }

  if(N_all != 0)
  {
    //Check contribution
    if(targetCharge == 1)
    {
      if(PiPlusCounter == 3)
      {
        //Loop over three pions
        for (int i = 0; i < N_pi; i++)
        {
          if(N_all != 0)
          {
            P_1p3pito1p1pi[i] = -(N_1p1pi[i]/N_all);
          }
          else 
          {
            P_1p3pito1p1pi[i] = 0.0;
          }
        }     
      }
      else
      {
        TVector3 V3pi2[2];
        TLorentzVector V4pi2[2];
        int q_pi2[2] = {0};
        double P_1p1pi[2] = {0};
        double Ecal2[2] = {0};
        double p_miss_perp2[2] = {0};

        std::string strArray[4] = {"++", "+-", "-+", "--"};
        if(N_all != 0)
        {
          for(int i = 0; i < N_pi; i++)
          {
            for(int j = 0; j < N_pi; j++)
            {
              P_1p1pi[0] = P_1p1pi[1] = 0;
              V3pi2[0] = V3pi[i];
              V3pi2[1] = V3pi[j];
              V4pi2[0] = V4pi[i];
              V4pi2[1] = V4pi[j];
              q_pi2[0] = q_pi[i];
              q_pi2[1] = q_pi[j];
              prot1_pi2_rot_func(V3prot, V3pi2, V4prot, V4pi2, q_pi2, V4_el, Ecal2, p_miss_perp2, P_1p1pi, targetCharge);
              Ecal[i] = Ecal2[0];
              Ecal[j] = Ecal2[1];
              p_miss_perp[i] = p_miss_perp2[0];
              p_miss_perp[j] = p_miss_perp2[1];
              P_1p3pito1p1pi[i] += -P_1p1pi[0]*(N_1p2pi_diff[0][0] + N_1p2pi_diff[0][1])/N_all;
              P_1p3pito1p1pi[j] += -P_1p1pi[1]*(N_1p2pi_diff[1][0] + N_1p2pi_diff[1][1])/N_all;
            }
          }
        }
      }

    }
    else if(targetCharge == -1)
    {
      if(PiMinusCounter == 3)
      {
        //Loop over three pions
        for (int i = 0; i < N_pi; i++)
        {
          if(N_all != 0)
          {
            P_1p3pito1p1pi[i] = -(N_1p1pi[i]/N_all);
          }
          else 
          {
            P_1p3pito1p1pi[i] = 0.0;
          }
        }     
      }
      else
      {
        int count = 0;
        TVector3 V3pi2[2];
        TLorentzVector V4pi2[2];
        int q_pi2[2] = {0};
        double P_1p1pi[2] = {0};
        double Ecal2[2] = {0};
        double p_miss_perp2[2] = {0};

        std::string strArray[4] = {"++", "+-", "-+", "--"};
        if(N_all != 0)
        {
          for(int i = 0; i < N_pi; i++)
          {
            for(int j = 0; j < N_pi; j++)
            {
              P_1p1pi[0] = P_1p1pi[1] = 0;
              V3pi2[0] = V3pi[i];
              V3pi2[1] = V3pi[j];
              V4pi2[0] = V4pi[i];
              V4pi2[1] = V4pi[j];
              q_pi2[0] = q_pi[i];
              q_pi2[1] = q_pi[j];
              prot1_pi2_rot_func(V3prot, V3pi2, V4prot, V4pi2, q_pi2, V4_el, Ecal2, p_miss_perp2, P_1p1pi, targetCharge);
              Ecal[i] = Ecal2[0];
              Ecal[j] = Ecal2[1];
              p_miss_perp[i] = p_miss_perp2[0];
              p_miss_perp[j] = p_miss_perp2[1];
              P_1p3pito1p1pi[i] += -P_1p1pi[0]*(N_1p2pi_diff[0][0] + N_1p2pi_diff[0][1])/N_all;
              P_1p3pito1p1pi[j] += -P_1p1pi[1]*(N_1p2pi_diff[1][0] + N_1p2pi_diff[1][1])/N_all;
              count = count+1;
            }
          }
        }
      }

    }
    //Double-check that this is in-line with charge dist. 
    P_tot[0] = P_1p3pito1p1pi[0] + P_1p3pito1p1pi[0] + P_1p3pito1p1pi[1];
    P_tot[1] = P_1p3pito1p1pi[1] + P_1p3pito1p1pi[0] + P_1p3pito1p1pi[2];
    P_tot[2] = P_1p3pito1p1pi[2] + P_1p3pito1p1pi[1] + P_1p3pito1p1pi[2];
  }

  else if(targetCharge == 0)
  {
    for(int i=0; i< N_pi; i++)
    {
      if(N_all!=0)
      {
        P_1p3pito1p1pi[i] = -(N_1p1pi[i]/N_all);
      }
      else
      {
        P_1p3pito1p1pi[i] = 0;
      }
    }
    //--------------1p2pi->1p1pi------
    int count = 0;
    TVector3 V3pi2[2];
    TLorentzVector V4pi2[2];
    int q_pi2[2] = {0};
    double P_1p1pi[2] = {0};
    double Ecal2[2] = {0};
    double p_miss_perp2[2] = {0};
    if(N_all!=0)
    {
      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          if(i<j)
          {
            P_1p1pi[0] = P_1p1pi[1] = 0;
            V3pi2[0] = V3pi[i];
            V3pi2[1] = V3pi[j];
            V4pi2[0] = V4pi[i];
            V4pi2[1] = V4pi[j];
            q_pi2[0] = q_pi[i];
            q_pi2[1] = q_pi[j];
            prot1_pi2_rot_func(V3prot, V3pi2, V4prot, V4pi2, q_pi2, V4_el, Ecal2, p_miss_perp2, P_1p1pi, 0);
            Ecal[i] = Ecal2[0];
            Ecal[j] = Ecal2[1];
            p_miss_perp[i] = p_miss_perp2[0];
            p_miss_perp[j] = p_miss_perp2[1];
            P_1p3pito1p1pi[i] += -P_1p1pi[0]*N_1p2pi[count]/N_all;
            P_1p3pito1p1pi[j] += -P_1p1pi[1]*N_1p2pi[count]/N_all;
            count = count+1;
          }
        }
      }
      //Ask about this
      P_tot[0] = P_1p3pito1p1pi[0] + P_1p3pito1p1pi[0] + P_1p3pito1p1pi[1];
      P_tot[1] = P_1p3pito1p1pi[1] + P_1p3pito1p1pi[0] + P_1p3pito1p1pi[2];
      P_tot[2] = P_1p3pito1p1pi[2] + P_1p3pito1p1pi[1] + P_1p3pito1p1pi[2];
    }
  }
}


void Subtraction::prot2_pi1_rot_func(TLorentzVector V4_2prot_uncorr[2], TLorentzVector V4_el, TLorentzVector V4_1pi, TVector3 V3_q, int q_pi, double Ecal[2], double p_miss_perp[2], double P_tot[2])
{
  //8.30.22 Removing arguments that don't even get used (V3_2prot_corr, )
  /*
    8.30.22 Changed corrected to uncorrected 
     - Removed V3 dependence -> We only want 4-vectors 
     - Implemented good programming practices aka consistent naming schema
  */
 //Check charge of pion matches what we want to fill :)
 const int num_prot = 2;
 TVector3 V3_2p_rotated[2], V3_1pi_rotated;
 bool status_pi = true;
 bool status_prot[2] = {true};
 double N_all = 0, N_1p_1pi[2] = {0};
 double P_2pto1p[2] = {0},N_2p_det = 0;
 double rot_angle = 0;

 for(int g = 0; g < N_tot; g++)
 {
     rot_angle = gRandom->Uniform(0,2*TMath::Pi());

     for(int i = 0; i < num_prot; i++)
     {//get rid of corrections (FOR NOW)
       V3_2p_rotated[i] = V4_2prot_uncorr[i].Vect();
       V3_2p_rotated[i].Rotate(rot_angle,V3_q);
       status_prot[i] = PFiducialCut(fbeam_en, V3_2p_rotated[i]);
     }

     V3_1pi_rotated = V4_1pi.Vect();
     //We rotate the particle around the momentum transfer vector (q)
     V3_1pi_rotated.Rotate(rot_angle,V3_q);
     status_pi = Pi_phot_fid_united(fbeam_en, V3_1pi_rotated, q_pi);

     if( status_prot[0] && !status_prot[1] && status_pi) N_1p_1pi[0]++;
     if(!status_prot[0] && status_prot[1] && status_pi) N_1p_1pi[1]++;
     if( status_prot[0] && status_prot[1] && status_pi) N_all++;
  }

  for(int h = 0; h < 2; h++)
  {
     prot1_pi1_en_calc(V4_2prot_uncorr[h], V4_1pi, q_pi, V4_el, &Ecal[h], &p_miss_perp[h]);
     if(N_all > 0)
        P_tot[h] = -N_1p_1pi[h]/N_all;
      else
        P_tot[h] = 0;
  }
}


void Subtraction::prot2_pi2_rot_func(TLorentzVector V4_2prot_uncorr[2], TLorentzVector V4_el, TLorentzVector V4_2pi[2], TVector3 V3_q, int q_pi[2], double Ecal[2][2], double p_miss_perp[2][2], double P_tot_2p[2][2])
{
  const int num_prot=2,num_pi=2;
  TVector3 V3_2p_rotated[num_prot],V3_2pirot[num_pi];
  bool pi2_stat[num_pi]={true};
  double rot_angle;
  double N_2p_1pi[num_pi]={0},N_1p_2pi[num_prot]={0},N_all=0,N_1p_1pi[num_prot][num_pi]={0};
  double N_pidet=0,N_piundet=0;
  double P_2pto1p[num_prot]={0},N_2p_det=0;
  double P_1p1pi[num_pi]={0};
  double P_2p1pito1p1pi[2]={0},Ptot=0;
  double P_2p2pito1p1pi[num_prot]={0},P_2p2pito1p2pi[num_prot]={0},P_2p2pito2p1pi[num_prot]={0};
  P_tot_2p[0][0] = P_tot_2p[0][1] = P_tot_2p[1][0] = P_tot_2p[1][1] = 0;

  for(int g=0; g<N_tot; g++)
  {
    rot_angle=gRandom->Uniform(0,2*TMath::Pi());

    for(int k=0; k<num_pi; k++)
    {
      V3_2p_rotated[k]=V4_2prot_uncorr[k].Vect();
      V3_2p_rotated[k].Rotate(rot_angle,V3q);

      V3_2pirot[k]=V4_2pi[k].Vect();
      V3_2pirot[k].Rotate(rot_angle,V3q);
      pi2_stat[k]=Pi_phot_fid_united(fbeam_en, V3_2pirot[k], q_pi[k]);
    }

    if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_2p_1pi[0]++;
    if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_2p_1pi[1]++;
    if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[0]++;
    if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_1p_2pi[1]++;
    if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[0][0]++;
    if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && !PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[0][1]++;
    if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && !pi2_stat[1])  N_1p_1pi[1][0]++;
    if(!PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && !pi2_stat[0]  && pi2_stat[1])  N_1p_1pi[1][1]++;
    if(PFiducialCut(fbeam_en, V3_2p_rotated[0]) && PFiducialCut(fbeam_en, V3_2p_rotated[1]) && pi2_stat[0]  && pi2_stat[1])  N_all++;
  }
  if(N_all!=0)
  {    
    //---------------------------------------------------2p2pi->1p2pi-------------------------------------------------------
    double prob2p2pito1p1pi[2][2] = {0};
    double P_tot[2] = {0};
    double Ecal2[2] = {0};
    double p_miss_perp2[2] = {0};
    //Energy calculation is const between all subtraction types
    TVector3 V3_2pi[2]; 
    for(int i = 0; i < 2; i++){
      V3_2pi[i] = V4_2pi[i].Vect();
    }
    for(int i=0; i < 2; i++) //Helps avoid segmentation fault error :)
    {
      prot1_pi2_rot_func(V4_2prot_uncorr[i].Vect(), V3_2pi, V4_2prot_uncorr[i], V4_2pi, q_pi, V4_el, Ecal2, p_miss_perp2, P_tot, 0);
      Ecal[i][0] = Ecal2[0];
      Ecal[i][1] = Ecal2[1];
      p_miss_perp[i][0] = p_miss_perp2[0];
      p_miss_perp[i][1] = p_miss_perp2[1];
    }

    for(int i=0;i<2;i++)
    {
      //Changed below to reflect the rotation function being performed. From N_2p_1pi -> N_1p_2pi
      prob2p2pito1p1pi[i][0]= -(N_1p_2pi[i]/N_all)*P_tot[0]; //-> Already made negative when P_Tot is calc
      prob2p2pito1p1pi[i][1]= -(N_1p_2pi[i]/N_all)*P_tot[1];
    }
    //---------------------------------------------------2p2pi->2p1pi->1p1pi-------------------------------------------------------
    P_tot[0] = P_tot[1] = 0;
    //Loop twice because we can get this two ways
    for(int i=0;i<2;i++)
    {
      prot2_pi1_rot_func(V4_2prot_uncorr, V4_el, V4_2pi[i], V3_q, q_pi[i], Ecal[i], p_miss_perp2, P_tot);
      p_miss_perp[0][i] = p_miss_perp2[0];
      p_miss_perp[1][i] = p_miss_perp2[1];
      prob2p2pito1p1pi[0][i] = -(N_2p_1pi[i]/N_all)*P_tot[0]; //Already made negative in 2p_1pi calc
      prob2p2pito1p1pi[1][i] = -(N_2p_1pi[i]/N_all)*P_tot[1];
      //Positive bc it is doubly subtracted ^
    }
    //----------------------------------------------------2p2pi->1p1pi---------------------------------------------------------
    for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
      {
        prob2p2pito1p1pi[i][j] += -N_1p_1pi[i][j]/N_all; //probability = 2stepWeight - 1stepWeight (prob = prob - N1...)
        P_tot_2p[i][j] = prob2p2pito1p1pi[i][j];
      }
    }
  }
  else if(N_all == 0)
  {
    P_tot_2p[0][0] = P_tot_2p[0][1] = P_tot_2p[1][0] = P_tot_2p[1][1] = 0;
  }
  else
  {
    //This should not happen 
    std::cout << "Error - N_all  == 0 && != 0" << std::endl;
  }

}

void Subtraction::prot3_pi1_rot_func(TLorentzVector V4_3prot_uncorr[3], TLorentzVector V4_el, TLorentzVector V4_1pi, TVector3 V3_q, int q_pi, double Ecal[3], double p_miss_perp[3], double P_tot_3p[3], TVector3 V3_3prot_uncorr[3],TVector3 V3_pi, TLorentzVector V4_3prot_corr[3])
{
    const int N_3prot=3;
    TVector3 V3_3p_rotated[N_3prot],V3_pirot;
    bool pi_stat=true;
    double rot_angle;

    double N_all=0,N_1p1pi[N_3prot]={0},N_2p1pi[N_3prot]={0};
    double  P_3p1pito1p1pi[N_3prot]={0};
    double   N_pidet=0,N_piundet=0;
    TVector3 V3_2p_corr[N_3prot],V3_2p_uncorr[N_3prot];
    double P_2pto1p[2]={0},N_2p_det=0;
    int count=0;
    double N_p1[N_3prot]={0},N_p_three=0;
    double P_3pto2p[N_3prot][2]={0};
    double P_3p1pito2p1pi[N_3prot]={0};
    double P_2p1pito1p1pi[2]={0},Ptot=0;

    for(int g=0; g<N_tot; g++)
    {
       rot_angle=gRandom->Uniform(0,2*TMath::Pi());
       for(int k=0; k<N_3prot; k++)
       {
         V3_3p_rotated[k]=V3_3prot_uncorr[k];
         V3_3p_rotated[k].Rotate(rot_angle,V3q);
       }

       V3_pirot=V3_pi;
       V3_pirot.Rotate(rot_angle,V3q);
       pi_stat=Pi_phot_fid_united(fbeam_en, V3_pirot, q_pi);

       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[0]=N_1p1pi[0]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[1]=N_1p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_1p1pi[2]=N_1p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && !PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[0]=N_2p1pi[0]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && !PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[1]=N_2p1pi[1]+1;
       if(!PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_2p1pi[2]=N_2p1pi[2]+1;
       if(PFiducialCut(fbeam_en, V3_3p_rotated[0]) && PFiducialCut(fbeam_en, V3_3p_rotated[1]) && PFiducialCut(fbeam_en, V3_3p_rotated[2])  && pi_stat)  N_all=N_all+1;
    }
    if(N_all!=0)
    {
      for(int z=0;z<N_3prot;z++)
      {
        //---------------------------------- 3p 1pi ->1p 1pi   ----------------------------------------------
        P_3p1pito1p1pi[z] = -(N_1p1pi[z]/N_all);
       //---------------------------------- 3p 1pi ->2p 1pi   ----------------------------------------------
       TVector3 V3_prot_uncorr[2];
       TLorentzVector V4_prot_corr[2];
       TLorentzVector V4_prot_uncorr[2];
       double Ecal2[2] = {0};
       double p_miss_perp2[2] = {0};
       for(int i=0;i<N_3prot;i++)
       {
         //looping through 2p combinations  out of 3p
         if(z!=i && z<i)
         {
              // 3 pairs of 2proton combinations with z, i indexes(z<i)
              P_2p1pito1p1pi[0]=P_2p1pito1p1pi[1]=0;
              Ptot=0;
              V3_prot_uncorr[0] = V4_3prot_uncorr[z].Vect();
              V3_prot_uncorr[1] = V4_3prot_uncorr[i].Vect();
              V4_prot_uncorr[0] = V4_3prot_uncorr[z];
              V4_prot_uncorr[1] = V4_3prot_uncorr[i];
              V4_prot_corr[0] = V4_3prot_corr[z];
              V4_prot_corr[1] = V4_3prot_corr[z];

              
              //Ali look here -> Implement prot2_pi1 correctly
              prot2_pi1_rot_func(V4_prot_uncorr, V4_el, V4_1pi, V3_q, q_pi, Ecal2, p_miss_perp2, P_2p1pito1p1pi);
              //prot2_pi1_rot_func(V3_prot_corr,V3_prot_uncorr,V3_pi,V4_prot_corr, V4_pi, q_pi, V4_el, Ecal2, p_miss_perp2, P_2p1pito1p1pi);
              Ecal[z] = Ecal2[0];
              Ecal[i] = Ecal2[1];
              p_miss_perp[z] = p_miss_perp2[0];
              p_miss_perp[i] = p_miss_perp2[1];

              P_3p1pito2p1pi[z] += -(N_2p1pi[count]/N_all)*(P_2p1pito1p1pi[0]); //Probability will automatically be + in case that it is 2p2pi->1p1pi
              P_3p1pito2p1pi[i] += -(N_2p1pi[count]/N_all)*(P_2p1pito1p1pi[1]);

      	      count=count+1;
      		 }
      		}
      }//looping through 3p
      P_tot_3p[0]=P_3p1pito2p1pi[0]+P_3p1pito1p1pi[0];
      P_tot_3p[1]=P_3p1pito2p1pi[1]+P_3p1pito1p1pi[1];
      P_tot_3p[2]=P_3p1pito2p1pi[2]+P_3p1pito1p1pi[2];

    }
    else
    {
      P_tot_3p[0]= P_tot_3p[1]=P_tot_3p[2]=0;
    }

  }

void Subtraction::prot1_pi1_en_calc(TLorentzVector V4prot, TLorentzVector V4pi, int q_pi, TLorentzVector V4_el, double *Ecal, double *p_miss_perp)
{
  double m_prot=0.9382720813;
  TVector3 V3_total = V4prot.Vect() + V4pi.Vect() + V4_el.Vect();
  *Ecal = V4_el.E() + V4prot.E() - m_prot + V4pi.E();
  *p_miss_perp = TMath::Sqrt(V3_total.Px()*V3_total.Px()+V3_total.Py()*V3_total.Py());
}
#endif
