#include <TChain.h>


void Make_class()
{

  TChain *ch = new TChain("ch", "energy_reconstruction");   //here you write the name of the code you want to create "e2a_ep_neutrino6_errors_united2"
 ch->Add("*_skim.root/h10");  //here we have the path to data files "/work/clas/clase2/Mariana/data/e2a_56Fe_2261_v1/"
 ch->MakeClass("energy_reconstruction");  //again the name of the code to be created "e2a_ep_neutrino6_errors_united2"

}
