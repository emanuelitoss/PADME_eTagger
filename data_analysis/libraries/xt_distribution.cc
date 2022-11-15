#include "infos.h"

#include "TH1F.h"
#include "TH1.h"
#include "TString.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TStyle.h"

void HistoFillDeltaXperFile(TH1F* histogram, TF1* function, char * fileName,  double true_x){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("arrival_times", myFile);
    std::vector < TTreeReaderValue <Double_t> > times;

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "arrival_times_SiPM[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      times.push_back(TTreeReaderValue<Double_t>(reader,dirName));
    }

    double_t time = 0, min_timeSX = 0, min_timeDX = 0;
    double_t reconstructed_x;

    while(reader.Next())
    {

        min_timeSX = OFFSET;
        min_timeDX = OFFSET;

        for (int channel = 0; channel < numberOfChannels; channel++)
        {
            time = *times[channel];
            if(time != 0)
            {
                if(channel < numberOfChannels/2)
                {
                    if(time < min_timeDX) min_timeDX = time;
                }
                else
                {
                    if(time < min_timeSX) min_timeSX = time;
                }            
            }
        }

        reconstructed_x = function->Eval(min_timeDX-min_timeSX);
        histogram->Fill(reconstructed_x - true_x);

    }

    return;
}

void HistoFillDeltaXRandomFile(TH1F* histogram, TF1* function, TString fileName){
        
    // reading objects
    TFile* myFile = TFile::Open(fileName);

    TTreeReader readerT = TTreeReader();
    readerT.SetTree("arrival_times", myFile);

    TTreeReader readerX = TTreeReader();
    readerX.SetTree("positions", myFile);
    TTreeReaderValue <Double_t> position = TTreeReaderValue <Double_t>(readerX,"X_position");

    std::vector < TTreeReaderValue <Double_t> > times;
    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "arrival_times_SiPM[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      times.push_back(TTreeReaderValue<Double_t>(readerT,dirName));
    }

    double_t time = 0, min_timeSX = 0, min_timeDX = 0;
    double_t reconstructed_x, true_x;

    while(readerT.Next() && readerX.Next())
    {

        min_timeSX = OFFSET;
        min_timeDX = OFFSET;

        true_x = *position;

        for (int channel = 0; channel < numberOfChannels; channel++)
        {
            time = *times[channel];
            if(time != 0)
            {
                if(channel < numberOfChannels/2)
                {
                    if(time < min_timeDX) min_timeDX = time;
                }
                else
                {
                    if(time < min_timeSX) min_timeSX = time;
                }            
            }
        }

        reconstructed_x = function->Eval(min_timeDX-min_timeSX);
        histogram->Fill(reconstructed_x - true_x);

    }

    myFile->Close();
    delete myFile;
}

void PlotHistogramDeltaX(TH1F* histo_deltaX, TString output_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 4000);
    TF1* gauss_fit = new TF1();

    canva->SetGrid();

    histo_deltaX->Draw("EP");
    histo_deltaX->GetXaxis()->SetTitle("Length [mm]");
    histo_deltaX->GetYaxis()->SetTitle("Number of events");
    histo_deltaX->SetLineColor(kBlack);
    histo_deltaX->SetLineWidth((Width_t)1.5);
    
    gStyle->SetEndErrorSize(8);
    
    gauss_fit = new TF1("fitting a gaussian", "gaus", -HALF_LEN_X, HALF_LEN_X);
    histo_deltaX->Fit(gauss_fit, "0", "0");
    gauss_fit->SetLineColor(kAzure-5);
    gauss_fit->SetLineWidth(1);
    gauss_fit->SetFillStyle(3002);
    gauss_fit->SetFillColorAlpha(kAzure-5,0.5);
    gauss_fit->Draw("C SAME");
    
    canva->Print(TString(output_name),"pdf");
    canva->Clear();
  
    delete gauss_fit;
    delete canva;

}