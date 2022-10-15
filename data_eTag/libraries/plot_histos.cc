#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
using std::cout;
using std::endl;

#include "infos.h"

// ROOT header files
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"

// settings
#define nbins 50
#define max_time 30. //[ns]

void print_histos(char * fileName, int openCloseFile){

  // open a new file
  TFile* myFile = TFile::Open(fileName);

  /******* READERS & HISTOS INITIALIZATION *******/

  TTreeReader* reader = new TTreeReader("eTagDataTuples", myFile);
  std::vector< TTreeReaderValue<Double_t> > times;

  std::vector <TH1F> histograms;

  for (int channel = 0; channel < numberOfChannels; ++channel){

    std::string ChannelName = std::to_string(channel+1);
    
    std::string dirNameS = "PhotonsTime[" + ChannelName + "]";
    const char* dirName = dirNameS.c_str();
    times.push_back(TTreeReaderValue<Double_t>(*reader,dirName));

    // histograms initialization
    std::string histoNameS = "histogram[" + ChannelName + "]";
    const char* histoName = histoNameS.c_str();
    std::string histoTitleS = "noOfPhotons[" + ChannelName + "]";
    const char* histoTitle = histoTitleS.c_str();
    histograms.push_back(TH1F(histoName, histoTitle, nbins, 0, max_time));

  }

  /************ REMOVE ZEROS AND FILL HISTOS ************/

  double_t time = 0;

  while(reader->Next()){

    for (int channel = 0; channel < numberOfChannels; channel++){

      time = *times[channel];      
      if(time != 0) histograms[channel].Fill(time);
      
    }
  }

  double_t max_count = 0, val = 0;

  for (int channel = 0; channel < numberOfChannels; channel++){

    val = histograms[channel].GetMaximum();
    if (val > max_count) max_count = val;
      
  }
  
  /************ SUPERIMPOSED PLOT ************/
  
  TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3200, 3600);
  canva->SetGrid();
  const int color[8] = {kRed+2, kRed+2, kRed+2, kRed+2, kTeal+3, kTeal+3, kTeal+3, kTeal+3};
  const int colFill[8] = {kRed-3, kRed-3, kRed-3, kRed-3, kTeal+5, kTeal+5, kTeal+5, kTeal+5};
  
  histograms[0].GetXaxis()->SetTitle("Time [ns]");
  histograms[0].GetYaxis()->SetTitle("Number of events");

  TLegend* legend = new TLegend(0.7, 0.6, 1.1, 0.9);

  for (int channel = 0; channel < numberOfChannels; channel++){

    histograms[channel].GetYaxis()->SetRangeUser(0.,max_count*1.05);
    histograms[channel].Draw("C SAME");
    histograms[channel].SetLineColor(color[channel]);
    
    std::string histoTitleS = "noOfPhotons[" + std::to_string(channel+1) + "]";
    const char* histoTitle = histoTitleS.c_str();
    legend->AddEntry(&histograms[channel], histoTitle, "lf");
    canva->Update();

  }

  legend->Draw();
  if(openCloseFile % 2 == 0)  canva->Print("images/signal_waveforms.pdf(","pdf");
  else  canva->Print("images/signal_waveforms.pdf","pdf");
  
  canva->Clear();

  /************ DIVIDED PLOT ************/

  //1
  canva->Divide(1,4);

  for (int channel = 0; channel < 4; channel++)
  {
    canva->cd(channel+1);
    histograms[channel].Draw("SAMES E1 C");
    histograms[channel].GetXaxis()->SetTitle("Time [ns]");
    histograms[channel].GetYaxis()->SetTitle("Number of events");
    histograms[channel].SetFillColor(colFill[channel]);
    histograms[channel].SetLineColor(color[channel]);
    canva->Update();
  }

  canva->Print("images/signal_waveforms.pdf","pdf");
  canva->Clear();

  //2
  canva->Divide(1,4);

  for (int channel = 0; channel < 4; channel++)
  {
    canva->cd(channel+1);
    histograms[channel+4].Draw("SAMES E1 C");
    histograms[channel+4].GetXaxis()->SetTitle("Time [ns]");
    histograms[channel+4].GetYaxis()->SetTitle("Number of events");
    histograms[channel+4].SetLineColor(color[channel+4]);
    canva->Update();
  }

  if(openCloseFile == 1 || openCloseFile == 2)  canva->Print("images/signal_waveforms.pdf)","pdf");
  else  canva->Print("images/signal_waveforms.pdf","pdf");

  /************ CLOSE FILES ************/

  myFile->Close();

  delete legend;
  delete canva;
  delete reader;
  delete myFile;

}