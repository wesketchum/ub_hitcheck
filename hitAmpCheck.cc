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
#include <cmath>

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
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "Regions.hh"

void FillShortedRegionPlot(TH2F *h)
{
    for(int i_y=1; i_y <= h->GetNbinsY(); ++i_y)
      for(int i_z=1; i_z <= h->GetNbinsX(); ++i_z)	
	h->SetBinContent(i_z,i_y,
			 GetShortedRegionType(h->GetYaxis()->GetBinCenter(i_y),
					      h->GetXaxis()->GetBinCenter(i_z)));
}

int main(int argv, char** argc) {

  if(argv!=4){
    std::cout << "ERROR: Usage is hitCheckAmp <config_file> <input_file_list> <output_file>" << std::endl;
    return -1;
  }
  
  fhicl::ParameterSet pset;
  putenv(const_cast<char*>("FHICL_FILE_PATH=.:$FHICL_FILE_PATH"));
  cet::filepath_lookup policy("FHICL_FILE_PATH");
  fhicl::make_ParameterSet(argc[1],policy,pset);

  art::InputTag digit_tag = pset.get<art::InputTag>("DigitModuleLabel");
  art::InputTag hit_tag = pset.get<art::InputTag>("HitModuleLabel");
  bool makeShortedRegionPlot = pset.get<bool>("MakeShortedRegionPlot",false);
  bool make2DHistos = pset.get<bool>("Make2DHistos",false);

  int hitMultiplicityCut = pset.get<int>("HitMultiplicityCut");
  int hitRMSCut = pset.get<int>("HitRMSCut");
  
  TFile f_output(argc[3],"RECREATE");
  
  
  if(makeShortedRegionPlot){
    TH2F *h_shortedregions = new TH2F("h_shortedregions","Shorted Regions;Z (cm); Y (cm)",
				      1100,0,1100.,
				      250,-125.,125.);
    FillShortedRegionPlot(h_shortedregions);    
  }

  TH2F *h_U_ph_ch;
  TH2F *h_V_ph_ch;
  TH2F *h_Y_ph_ch;
  if(make2DHistos){
    h_U_ph_ch = new TH2F("h_U_ph_ch","U Plane; Channel; Hit Pulse Amplitude",2400,-0.5,2399.5,200,0,200);
    h_V_ph_ch = new TH2F("h_V_ph_ch","V Plane; Channel; Hit Pulse Amplitude",2400,2399.5,4799.5,200,0,200);
    h_Y_ph_ch = new TH2F("h_Y_ph_ch","Y Plane; Channel; Hit Pulse Amplitude",3456,4799.5,8255.5,200,0,200);
  }
  
  TNtuple *nt_hits = new TNtuple("nt_hits","Hits Ntuple","run:subrun:ev:timestamp:ch:hit_amp:hit_integral:hit_sumadc:hit_rms:hit_mult");
  TNtuple *nt_ch = new TNtuple("nt_ch","Channel Ntuple","run:subrun:ev:timestamp:ch:pedestal:rms:n_hits");


  
  std::vector<std::string> filenames;
  std::string file_name;
  std::ifstream input_file(argc[2]);
  while(getline(input_file,file_name))
    filenames.push_back(file_name);

  std::unordered_map<raw::ChannelID_t,size_t> nhits_map;  

  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    
    auto const& hit_handle = ev.getValidHandle< std::vector<recob::Hit> >(hit_tag);
    auto const& hitVec = *hit_handle;

    for(auto const& hit : hitVec){

      if(hit.Multiplicity() > hitMultiplicityCut &&
	 hit.RMS() < hitRMSCut ){
	nt_hits->Fill(ev.eventAuxiliary().run(),
		      ev.eventAuxiliary().subRun(),
		      ev.eventAuxiliary().event(),
		      ev.eventAuxiliary().time().timeHigh(),
		      hit.Channel(),
		      hit.PeakAmplitude(),
		      hit.Integral(),
		      hit.SummedADC(),
		      hit.RMS(),
		      hit.Multiplicity());
	nhits_map[hit.Channel()] += 1;
      }
      if(make2DHistos){

	if(hit.WireID().Plane==0)
	  h_U_ph_ch->Fill(hit.Channel(),hit.PeakAmplitude());
	else if(hit.WireID().Plane==1)
	  h_V_ph_ch->Fill(hit.Channel(),hit.PeakAmplitude());
	else if(hit.WireID().Plane==2)
	  h_Y_ph_ch->Fill(hit.Channel(),hit.PeakAmplitude());
      }
      
    }//end loop over hits

    auto const& digit_handle = ev.getValidHandle< std::vector<raw::RawDigit> >(digit_tag);
    auto const& digitVec = *digit_handle;
    for(auto const& digit : digitVec){

      float ped = (float)std::accumulate(digit.ADCs().begin(),digit.ADCs().end(),0) / (float)(digit.ADCs().size());
      float rms = 0.0;
      std::for_each (digit.ADCs().begin(), digit.ADCs().end(), [&](const float d) {
	  rms += (d - ped) * (d - ped);
	});
      rms = std::sqrt(rms / digit.ADCs().size());

      nt_ch->Fill(ev.eventAuxiliary().run(),
		  ev.eventAuxiliary().subRun(),
		  ev.eventAuxiliary().event(),
		  ev.eventAuxiliary().time().timeHigh(),
		  digit.Channel(),
		  ped,
		  rms,
		  nhits_map[digit.Channel()]);
    }
    
  } //end loop over events!


  //and ... write to file!
  f_output.Write();
  f_output.Close();

}
