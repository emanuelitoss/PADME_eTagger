#include "infos.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TStyle.h"

void HistoFillDeltaXperFileTimes(TH1F* histogram, TF1* function, char * fileName,  double true_x, vector <Double_t>* deltas, int my_option){
    
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
        
        if(my_option == OPTION_DELTA_POSITION)
        {
            histogram->Fill(reconstructed_x - true_x);
            (*deltas).push_back(reconstructed_x - true_x);
        }
        else if(my_option == OPTION_POSITION)
        {
            histogram->Fill(reconstructed_x);
            (*deltas).push_back(reconstructed_x);
        }

    }

    return;
}

void HistoFillDeltaXRandomFileTimes(TH1F* histogram, TF1* function, TString fileName, vector <Double_t> * deltas){
        
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
        deltas->push_back(reconstructed_x - true_x);

    }

    myFile->Close();
    delete myFile;
}

void PlotHistogramDeltaXTimes(TH1F* histo_deltaX, TString output_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 4000);
    TF1* gauss_fit = new TF1();

    canva->SetGrid();

    histo_deltaX->Draw("EP");
    histo_deltaX->GetXaxis()->SetTitle("Length [mm]");
    histo_deltaX->GetYaxis()->SetTitle("Number of events");
    histo_deltaX->SetLineColor(kBlack);
    histo_deltaX->SetLineWidth((Width_t)1.5);
    
    //gStyle->Reset();
    gStyle->SetEndErrorSize(8);
    gStyle->SetOptFit(1110);
    gStyle->SetOptStat(2210);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetStatFont(42);
    //gStyle->SetLegendFillColor(0);
    //gStyle->SetTitleFillColor(0);
    //gStyle->SetFillColor(0);
    //gStyle->SetPadBorderSize(0.);
    gStyle->SetLineWidth(1);
    //gStyle->SetDrawBorder(false);
    
    gauss_fit = new TF1("fitting a gaussian", "gaus", -HALF_LEN_X, HALF_LEN_X);
    histo_deltaX->Fit(gauss_fit, "L Q", "0");
    gauss_fit->SetLineColor(kAzure-5);
    gauss_fit->SetLineWidth(1);
    gauss_fit->SetFillStyle(3002);
    gauss_fit->SetFillColorAlpha(kAzure-5,0.5);
    gauss_fit->SetObjectStat(true);
    gauss_fit->Draw("C SAME");
    
    canva->Print(TString(output_name),"pdf");
    canva->Clear();
  
    delete gauss_fit;
    delete canva;

}

void PlotPositionsQT(vector <Double_t >* deltasT,vector <Double_t >* deltasQ, int openclosefile, int option_diff_ratio){

    // check that I have couples of points.
    int len = (*deltasT).size();
    bool check_dimensions = (len==(*deltasQ).size());
    if(!check_dimensions){
        PrintColor("ERROR: Not coupled points", OBOLDRED);
        return;
    }

    TH2F histogram = TH2F("hist","#deltax (time analysis) vs #deltax (charge analysis)",70,-HALF_LEN_X,HALF_LEN_X,70,-HALF_LEN_X,HALF_LEN_X);
    
    if(option_diff_ratio == OPTION_Q_DIFFERENCE){
        histogram.SetName("hist_diff");
        histogram.SetTitle("#deltax (time analysis) vs #deltax (charge differences analysis)");
    }
    else if(option_diff_ratio == OPTION_Q_RATIO){
        histogram.SetName("hist_ratio");
        histogram.SetTitle("#deltax (time analysis) vs #deltax (charge ratios analysis)");
    }

    // storing.
    Double_t deltaxT[len], deltaxQ[len];
    for(int i = 0; i<len; ++i)  histogram.Fill( (*deltasT)[i], (*deltasQ)[i] );

    // plot.
    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 2000, 1500);

    histogram.SetXTitle("x (times analysis) [mm]");
    histogram.SetYTitle("x (charges analysis) [mm]");
    gStyle->SetStatH(0.1);
    gStyle->SetStatX(0.42);
    gStyle->SetStatY(0.89);

    histogram.Draw("CONT4Z");
    if(openclosefile == OPEN_OUTPUT || openclosefile == SINGLE_OUTPUT)  canva->Print("images/5delta_positions_correlation.pdf(","pdf");
    else canva->Print("images/5delta_positions_correlation.pdf","pdf");
    canva->Clear();

    gStyle->SetStatX(0.99);
    gStyle->SetStatY(0.85);

    histogram.Draw("LEGO2");
    histogram.SetXTitle("#deltax (times analysis) [mm]");
    histogram.SetYTitle("#deltax (charges analysis) [mm]");
    if(openclosefile == CLOSE_OUTPUT || openclosefile == SINGLE_OUTPUT) canva->Print("images/5delta_positions_correlation.pdf)","pdf");
    else canva->Print("images/5delta_positions_correlation.pdf","pdf");
    canva->Clear();

    delete canva;


}