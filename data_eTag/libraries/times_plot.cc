#include "infos.h"

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

// settings
#define nbins 30
#define max_time 3.5 //[ns]

void plotHisto_arrivalTimes(char * fileName, int openCloseFile, std::vector <std::vector <Double_t> >* means, std::vector <std::vector <Double_t> >* stdDevs){

  std::cout << "Inspecting file: " << fileName << std::endl;

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
    histograms[channel].SetLineWidth((Width_t)1.5);

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

}

void plotHisto_arrivalTimes2(char * fileName, int openCloseFile, std::vector <std::vector <Double_t> >* means2, std::vector <std::vector <Double_t> >* stdDevs2){

    /********** READ FILE AND INITIALIZE READERS **********/

    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("firstTimes", myFile);
    std::vector < TTreeReaderValue <Double_t> > times;

    // to manipulate and save data
    const int DX = 0, SX = 1;
    std::vector <std::vector <Double_t> > min_times = {{}, {}};
    std::vector < TH1F > histograms;
    TH1F histogram_differences = TH1F("histogram[DX -SX]", "difference dx - sx", nbins*5, -max_time, max_time);

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "1stTimeSiPM[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      times.push_back(TTreeReaderValue<Double_t>(reader,dirName));
    }

    histograms.push_back(TH1F("histogram[SiPM_DX]", "SiPM_DX", nbins, 0, max_time));
    histograms.push_back(TH1F("histogram[SiPM_SX]", "SiPM_SX", nbins, 0, max_time));

    int counter = 0, entry = 0;
    double_t time = 0, min_timeSX = 0, min_timeDX = 0;
    bool checkSX, checkDX;

    while(reader.Next())
    {
        entry = (int)counter/8;

        for (int channel = 0; channel < numberOfChannels; channel++)
        {
            time = *times[channel];
            if(time != 0)
            {
                checkDX = channel < 4;
                checkSX = !checkDX;

                if(checkDX)
                {
                    if(min_times[DX].size() == entry+1)
                    {
                        if(min_times[DX][entry]>time) min_times[DX][entry] = time;
                    }
                    else min_times[DX].push_back(time);
                }
                else if(checkSX)
                {

                    if(min_times[SX].size() == entry+1)
                    {
                        if(min_times[SX][entry]>time) min_times[SX][entry] = time;
                    }
                    else min_times[SX].push_back(time);
                }
            }
        }
        
        ++counter;
        
        if(counter % 8 == 0)
        {
            histograms[DX].Fill(min_times[DX].back());
            histograms[SX].Fill(min_times[SX].back());
            histogram_differences.Fill(min_times[DX].back() - min_times[SX].back());
        }
    }

    /********** PLOT HISTOGRAM OF EACH SIDE OF SiPMs **********/

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3800, 3800);
    TF1* gauss_fit = new TF1();
    const int color[8] = {kBlue+3, kOrange+9};

    canva->Divide(2,1);

    for(int side = 0; side < 2; ++side)
    {
      // mean and std deviation
      (*means2)[side].push_back( histograms[side].GetMean() );
      (*stdDevs2)[side].push_back( histograms[side].GetStdDev() );

      // histograms
      canva->cd(2-side);
      histograms[side].Draw("");
      histograms[side].Draw("SAME E1");
      histograms[side].GetXaxis()->SetTitle("Time [ns]");
      histograms[side].GetYaxis()->SetTitle("Number of events");
      histograms[side].SetLineColor(kBlack);
      histograms[side].SetLineWidth((Width_t)1.5);

      gStyle->SetEndErrorSize(8); //4 is the number of pixels

      gauss_fit = new TF1("fitting a gaussian", "gaus", -HALF_LEN_X, HALF_LEN_X);
      histograms[side].Fit(gauss_fit, "0", "0");
      gauss_fit->SetLineColor(color[side]);
      gauss_fit->SetLineWidth(1);
      gauss_fit->SetFillStyle(3002);
      gauss_fit->SetFillColorAlpha(color[side],0.5);
      gauss_fit->Draw("SAME C");
        
    }

  if(openCloseFile == 0) canva->Print("images/InitialTimes2.pdf(","pdf");
  if(openCloseFile == 1) canva->Print("images/InitialTimes2.pdf)","pdf");
  if(openCloseFile == 3 || openCloseFile == 2) canva->Print("images/InitialTimes2.pdf","pdf");

  canva->Clear();
  gauss_fit->Clear();
  
  /********** PLOT HISTOGRAM DIFFERENCES BETWEEN SIDES OF SiPMs **********/

  //canva = new TCanvas("canva2", "canvas for plotting", 3800, 3800);
  const int mycolor = kRed+2;

  // histograms
  histogram_differences.Draw("");
  histogram_differences.Draw("SAME E1");
  histogram_differences.GetXaxis()->SetTitle("Time [ns]");
  histogram_differences.GetYaxis()->SetTitle("Number of events");
  histogram_differences.SetLineColor(kBlack);
  histogram_differences.SetLineWidth((Width_t)1.5);

  gStyle->SetEndErrorSize(8); //4 is the number of pixels

  gauss_fit = new TF1("fitting a gaussian", "gaus", -HALF_LEN_X, HALF_LEN_X);
  histogram_differences.Fit(gauss_fit, "0", "0");
  gauss_fit->SetLineColor(mycolor);
  gauss_fit->SetLineWidth(1);
  gauss_fit->SetFillStyle(3002);
  gauss_fit->SetFillColorAlpha(mycolor,0.5);
  gauss_fit->Draw("SAME C");

  if(openCloseFile == 0) canva->Print("images/differences_InitialTimes.pdf(","pdf");
  if(openCloseFile == 1) canva->Print("images/differences_InitialTimes.pdf)","pdf");
  if(openCloseFile == 3 || openCloseFile == 2) canva->Print("images/differences_InitialTimes.pdf","pdf");
  
  /********** DELETE STUFF **********/
  
  delete gauss_fit;
  delete canva;

  myFile->Close();
  delete myFile;

  return;

}