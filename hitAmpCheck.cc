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
#include "TVirtualFFT.h"

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

  int inEventsFFTPerRun = pset.get<int>("EventsFFTPerRun");
  size_t nEventsFFTPerRun = (inEventsFFTPerRun>0) ? inEventsFFTPerRun : 999999;
    
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

  TVirtualFFT *fftr2c;
  std::vector<double> doubleVecInput,doubleVecOutput;
  // = TVirtualFFT::FFT(1, &NDIM, "R2C EX K");

  TNtuple *nt_hits = new TNtuple("nt_hits","Hits Ntuple","run:subrun:ev:timestamp:ch:hit_amp:hit_integral:hit_sumadc:hit_rms:hit_mult");

  typedef struct chTreeObj{
    unsigned int run;
    unsigned int subrun;
    unsigned int ev;
    unsigned int timestamp;
    unsigned int ch;
    float pedestal;
    float rms;
    float median;
    float rms_trunc;
    int n_hits;
    float fft_m_1q;
    float fft_m_2q;
    float fft_m_3q;
    float fft_m_4q;
    float fft_m_1q_trunc;
    float fft_m_2q_trunc;
    float fft_m_3q_trunc;
    float fft_m_4q_trunc; 
  } chTreeObj_t;

  chTreeObj_t chTreeData;

  TTree *nt_ch = new TTree("nt_ch","Channel Ntuple");
  nt_ch->Branch("",&chTreeData,"run/i:subrun/i:ev/i:timestamp/i:ch/i:pedestal/F:rms/F:median/F:rms_trunc/F:n_hits/I:fft_m_1q/F:fft_m_2q/F:fft_m_3q/F:fft_m_4q/F:fft_m_1q_trunc/F:fft_m_2q_trunc/F:fft_m_3q_trunc/F:fft_m_4q_trunc/F");


  
  std::vector<std::string> filenames;
  std::string file_name;
  std::ifstream input_file(argc[2]);
  while(getline(input_file,file_name))
    filenames.push_back(file_name);

  std::unordered_map<raw::ChannelID_t,size_t> nhits_map;  

  size_t n_evts=0;
  
  bool fft_Initialize = true;
  unsigned int prev_run = 0;
  size_t n_evts_thisrun = 0;
  
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {

    std::cout << "Processing event " << n_evts << std::endl;
    nhits_map.clear();

    if(ev.eventAuxiliary().run()!=prev_run){
      n_evts_thisrun=0;
      prev_run = ev.eventAuxiliary().run();
    }
    
    auto const& hit_handle = ev.getValidHandle< std::vector<recob::Hit> >(hit_tag);
    auto const& hitVec = *hit_handle;

    for(auto const& hit : hitVec){

      if(hit.Multiplicity() <= hitMultiplicityCut &&
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

    chTreeData.run = ev.eventAuxiliary().run();
    chTreeData.subrun = ev.eventAuxiliary().subRun();
    chTreeData.ev = ev.eventAuxiliary().event();
    chTreeData.timestamp = ev.eventAuxiliary().time().timeHigh();
    chTreeData.fft_m_1q=-1;
    chTreeData.fft_m_2q=-1;
    chTreeData.fft_m_3q=-1;
    chTreeData.fft_m_4q=-1;
    chTreeData.fft_m_1q_trunc=-1;
    chTreeData.fft_m_2q_trunc=-1;
    chTreeData.fft_m_3q_trunc=-1;
    chTreeData.fft_m_4q_trunc=-1;

    
    auto const& digit_handle = ev.getValidHandle< std::vector<raw::RawDigit> >(digit_tag);
    auto const& digitVec = *digit_handle;
    for(auto const& digit : digitVec){

      chTreeData.ch = digit.Channel();
      doubleVecInput.assign(digit.ADCs().begin(),digit.ADCs().end());

      //sort to get median
      std::nth_element(doubleVecInput.begin(),doubleVecInput.begin()+(doubleVecInput.size()/2),doubleVecInput.end());
      chTreeData.median = *(doubleVecInput.begin()+(doubleVecInput.size()/2));

      //sort to get first and third quartile
      std::nth_element(doubleVecInput.begin(),doubleVecInput.begin()+(doubleVecInput.size()/4),doubleVecInput.begin()+(doubleVecInput.size()/2));
      std::nth_element(doubleVecInput.begin()+(doubleVecInput.size()/2),doubleVecInput.begin()+(3*doubleVecInput.size()/4),doubleVecInput.end());

      float q1 = *(doubleVecInput.begin()+(doubleVecInput.size()/4));
      float q3 = *(doubleVecInput.begin()+(3*doubleVecInput.size()/4));

      chTreeData.rms_trunc = std::sqrt(0.5*((q1-chTreeData.median)*(q1-chTreeData.median) +
					    (q3-chTreeData.median)*(q3-chTreeData.median))) / 0.6745;
      
      chTreeData.pedestal = (float)std::accumulate(digit.ADCs().begin(),digit.ADCs().end(),0) / (float)(digit.ADCs().size());
      chTreeData.rms = 0.0;
      std::for_each (digit.ADCs().begin(), digit.ADCs().end(), [&](const float d) {
	  chTreeData.rms += (d - chTreeData.pedestal) * (d - chTreeData.pedestal);
	});
      chTreeData.rms = std::sqrt(chTreeData.rms / digit.ADCs().size());

      
      
      if(n_evts_thisrun < nEventsFFTPerRun){
	
	if(fft_Initialize){
	  int ndim = digit.ADCs().size();
	  fftr2c = TVirtualFFT::FFT(1, &ndim, "R2C EX K");
	  doubleVecOutput.resize(ndim);
	  fft_Initialize = false;
	}

	doubleVecInput.assign(digit.ADCs().begin(),digit.ADCs().end());
	fftr2c->SetPoints(doubleVecInput.data());
	fftr2c->Transform();

	TH1 *hfft_m = 0;
	hfft_m = TH1::TransformHisto(fftr2c,hfft_m,"MAG");

	chTreeData.fft_m_1q = hfft_m->Integral(2,hfft_m->GetNbinsX()/8);
	chTreeData.fft_m_2q = hfft_m->Integral(hfft_m->GetNbinsX()/8 + 1,2*hfft_m->GetNbinsX()/8);
	chTreeData.fft_m_3q = hfft_m->Integral(2*hfft_m->GetNbinsX()/8 + 1,3*hfft_m->GetNbinsX()/8);
	chTreeData.fft_m_4q = hfft_m->Integral(3*hfft_m->GetNbinsX()/8 + 1,4*hfft_m->GetNbinsX()/8);

	delete hfft_m;

	doubleVecInput.assign(digit.ADCs().begin(),digit.ADCs().end());
	for(auto & val : doubleVecInput)
	  if( std::abs(val-chTreeData.median)>(3.5*chTreeData.rms_trunc) ) val = chTreeData.median;
	
	fftr2c->SetPoints(doubleVecInput.data());
	fftr2c->Transform();
	//fftr2c->GetPoints(doubleVecOutput.data());

	TH1 *hfft_m_trunc = 0;
	hfft_m_trunc = TH1::TransformHisto(fftr2c,hfft_m_trunc,"MAG");

	chTreeData.fft_m_1q_trunc = hfft_m_trunc->Integral(2,hfft_m_trunc->GetNbinsX()/8);
	chTreeData.fft_m_2q_trunc = hfft_m_trunc->Integral(hfft_m_trunc->GetNbinsX()/8 + 1,2*hfft_m_trunc->GetNbinsX()/8);
	chTreeData.fft_m_3q_trunc = hfft_m_trunc->Integral(2*hfft_m_trunc->GetNbinsX()/8 + 1,3*hfft_m_trunc->GetNbinsX()/8);
	chTreeData.fft_m_4q_trunc = hfft_m_trunc->Integral(3*hfft_m_trunc->GetNbinsX()/8 + 1,4*hfft_m_trunc->GetNbinsX()/8);

	delete hfft_m_trunc;

      }
      
      nt_ch->Fill();
    }

    ++n_evts;
    ++n_evts_thisrun;

  } //end loop over events!

  std::cout << "Processed " << n_evts << " events." << std::endl;
  

  //and ... write to file!
  f_output.Write();
  f_output.Close();

}
