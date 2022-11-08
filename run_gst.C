#ifndef rungst
#define rungst
#include "gst.h"

#include <iostream>

using namespace std;

int main(int argc, char **argv)
{


  if( argc < 3 ){
    std::cout<<"Please specify the target (3He, 56Fe, C12, 4He), the beam energy (2261 or 4461)"<<std::endl;
    std::cout<<"================= Usage ==============="<<std::endl;
    std::cout<<"./genie_analysis target beam_energy"<<std::endl;
    exit(1);
  }


  std::string target  = argv[1];
  std::string beam_en = argv[2];

  gst  t(target,beam_en);
  t.Loop();


  return 0;
}
#endif
