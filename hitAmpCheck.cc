/*************************************************************
 *************************************************************/


//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Hit.h"

int main(int argv, char** argc) {

  if(argv!=3){
    std::cout << "ERROR: Usage is hitCheckAmp <input_file_list> <output_file>" << std::endl;
    return -1;
  }
  
  TFile f_output(argc[2],"RECREATE");

  TH2F *h_U_ph_ch = new TH2F("h_U_ph_ch","U Plane; Channel; Hit Pulse Amplitude",2400,-0.5,2399.5,200,0,200);
  TH2F *h_V_ph_ch = new TH2F("h_V_ph_ch","V Plane; Channel; Hit Pulse Amplitude",2400,2399.5,4799.5,200,0,200);
  TH2F *h_Y_ph_ch = new TH2F("h_Y_ph_ch","Y Plane; Channel; Hit Pulse Amplitude",3456,4799.5,8255.5,200,0,200);
  
  std::vector<std::string> filenames;

  std::string file_name;
  std::ifstream input_file(argc[1]);
  while(getline(input_file,file_name))
    filenames.push_back(file_name);


  art::InputTag hit_tag { "gaushit" };

  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    
    auto const& hit_handle = ev.getValidHandle< std::vector<recob::Hit> >(hit_tag);
    auto const& hitVec = *hit_handle;

    for(auto const& hit : hitVec){
      if(hit.WireID().Plane==0)
	h_U_ph_ch->Fill(hit.Channel(),hit.PeakAmplitude());
      else if(hit.WireID().Plane==1)
	h_V_ph_ch->Fill(hit.Channel(),hit.PeakAmplitude());
      else if(hit.WireID().Plane==2)
	h_Y=h_ph_ch->Fill(hit.Channel(),hit.PeakAmplitude());
    }
    
  } //end loop over events!


  //and ... write to file!
  f_output.Write();
  f_output.Close();

}
