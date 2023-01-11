#ifndef runevent_selector
#define runevent_selector
#include "event_selector.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{


  if( argc < 4 ){
    std::cout<<"Please specify the target (3He, 56Fe, C12, 4He), the beam energy (2261 or 4461), the filepath"<<std::endl;
    std::cout<<"================= Usage ==============="<<std::endl;
    std::cout<<"./genie_analysis target beam_energy filepath"<<std::endl;
    exit(1);
  }


  std::string target  = argv[1];
  std::string beam_en = argv[2];
  std::string file_name = argv[3];

  event_selector  t(target,beam_en, file_name);
  t.Loop();


  return 0;
}
#endif
