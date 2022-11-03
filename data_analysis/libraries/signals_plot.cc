#include "infos.h"

// settings
#define nbins 50
#define max_Time 50. //[ns]

void print_signals(char * fileName, int openCloseFile){

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
    histograms.push_back(TH1F(histoName, histoTitle, nbins, 0, max_Time));

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
  
  /************ SUPERIMPOSED PLOTS ************/
  
  TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3200, 3600);
  canva->SetGrid();
  const int color[8] = {kRed+2, kRed+2, kRed+2, kRed+2, kCyan+3, kCyan+3, kCyan+3, kCyan+3};
  const int colorShades[8] = {kRed+3, kRed+2, kRed+1, kRed, kCyan+1, kCyan+2, kCyan+3, kCyan+4};
  const int colFill[8] = {kRed-3, kRed-3, kRed-3, kRed-3, kCyan-6, kCyan-6, kCyan-6, kCyan-6};
  
  histograms[0].GetXaxis()->SetTitle("Time [ns]");
  histograms[0].GetYaxis()->SetTitle("Number of events");

  TLegend* legend = new TLegend(0.7, 0.6, 1.1, 0.9);

  for (int channel = 0; channel < numberOfChannels; channel++){

    histograms[channel].GetYaxis()->SetRangeUser(0.,max_count*1.05);
    histograms[channel].Draw("C SAME");
    histograms[channel].SetLineColor(colorShades[channel]);
    
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
    histograms[channel].GetXaxis()->SetTitle("Time [ns]");
    histograms[channel].GetYaxis()->SetTitle("Number of events");
    histograms[channel].SetLineColor(kBlack);
    gStyle->SetEndErrorSize(8);
    histograms[channel].Draw("E1");
    histograms[channel].SetFillColorAlpha(colFill[channel], 0.3);
    histograms[channel].SetFillStyle(3002);
    histograms[channel].Draw("SAME");
    canva->Update();
  }

  canva->Print("images/signal_waveforms.pdf","pdf");
  canva->Clear();

  //2
  canva->Divide(1,4);

  for (int channel = 0; channel < 4; channel++)
  {
    canva->cd(channel+1);
    histograms[channel+4].GetXaxis()->SetTitle("Time [ns]");
    histograms[channel+4].GetYaxis()->SetTitle("Number of events");
    histograms[channel+4].SetLineColor(kBlack);
    gStyle->SetEndErrorSize(8);
    histograms[channel+4].Draw("E1");
    histograms[channel+4].SetFillColorAlpha(colFill[channel+4], 0.3);
    histograms[channel+4].SetFillStyle(3002);
    histograms[channel+4].Draw("SAME");
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