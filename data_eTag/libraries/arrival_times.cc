#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
using std::cout;
using std::endl;

// ROOT header files
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TF1.h"

#include "infos.h"

// settings
#define nbins 30
#define max_time 3.5 //[ns]

void print_arrivalTimes(char * fileName, int openCloseFile, std::vector <std::vector <Double_t> >* means, std::vector <std::vector <Double_t> >* stdDevs){

  cout << "Inspecting file: " << fileName << endl;

  /********** READ FILE AND INITIALIZE READERS **********/

  // reading objects
  TFile* myFile = TFile::Open(fileName);
  TTreeReader reader = TTreeReader();
  reader.SetTree("firstTimes", myFile);
  std::vector< TTreeReaderValue<Double_t> > times;

  // 8 objects: one for each SiPM
  std::vector <TH1F> histograms;
  std::vector <std::vector <Double_t> > min_times = {{}, {}, {}, {}, {}, {}, {}, {}};

  for (int channel = 0; channel < numberOfChannels; ++channel)
  {
    // readers
    std::string ChannelName = std::to_string(channel+1);
    std::string dirNameS = "1stTimeSiPM[" + ChannelName + "]";
    const char* dirName = dirNameS.c_str();
    times.push_back(TTreeReaderValue<Double_t>(reader,dirName));

    // histograms
    std::string histoNameS = "histogram[" + ChannelName + "]";
    const char* histoName = histoNameS.c_str();
    std::string histoTitleS = "noOfPhotons[" + ChannelName + "]";
    const char* histoTitle = histoTitleS.c_str();
    histograms.push_back(TH1F(histoName, histoTitle, nbins, 0, max_time));
  }

  /********** READ FILE AND INITIALIZE READERS **********/

  double_t time = 0;

  while(reader.Next())
  {
    for (int channel = 0; channel < numberOfChannels; channel++)
    {

      time = *times[channel];
      if(time != 0){
        min_times[channel].push_back(time);
        histograms[channel].Fill(time);
      }

    }
  }

  /********** PLOT HISTOGRAM OF EACH SiPM **********/

  TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3200, 3800);
  TF1* gauss_fit = new TF1();
  const int color[8] = {kBlue+3, kBlue+3, kBlue+3, kBlue+3, kOrange+9, kOrange+9, kOrange+9, kOrange+9};

  canva->Divide(2,2);

  for (int channel = 0; channel < numberOfChannels; channel++){

    // mean and std deviation
    (*means)[channel].push_back( histograms[channel].GetMean() );
    (*stdDevs)[channel].push_back( histograms[channel].GetStdDev() );

    // histograms
    canva->cd((channel)%4 + 1);
    histograms[channel].Draw("");
    histograms[channel].Draw("SAME E1");
    histograms[channel].GetXaxis()->SetTitle("Time [ns]");
    histograms[channel].GetYaxis()->SetTitle("Number of events");
    histograms[channel].SetLineColor(kBlack);

    gStyle->SetEndErrorSize(8); //4 is the number of pixels

    gauss_fit = new TF1("fitting a line", "gaus", -HALF_LEN_X, HALF_LEN_X);
    histograms[channel].Fit(gauss_fit, "0", "0");
    gauss_fit->SetLineColor(color[channel]);
    gauss_fit->SetLineWidth(1);
    gauss_fit->SetFillStyle(3002);
    gauss_fit->SetFillColorAlpha(color[channel],0.5);
    gauss_fit->Draw("SAME C");


    if(channel == (int)(numberOfChannels/2 -1))
    {
      if (openCloseFile % 2 == 0) canva->Print("images/InitialTimes.pdf(","pdf");
      else canva->Print("images/InitialTimes.pdf","pdf");
      canva->Clear();
      canva->Divide(2,2);
    }
    else if(channel == numberOfChannels - 1)
    {
      if (openCloseFile == 1 || openCloseFile == 2) canva->Print("images/InitialTimes.pdf)","pdf");
      else canva->Print("images/InitialTimes.pdf","pdf");
      canva->Clear();
    }

  }

  myFile->Close();

  delete gauss_fit;
  delete canva;
  delete myFile;

};