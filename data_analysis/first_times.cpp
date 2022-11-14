#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "libraries/times_plot.cc"
#include "libraries/times_fit.cc"
#include "libraries/infos.h"

void HistoFillDeltaXperFile(TH1F* histogram, TF1* function, char * fileName,  double true_x);
void PlotHistogramDeltaX(TH1F* histo_deltaX, TString output_name);

int main(int argc, char** argv){
    
    // this because, for the folloxing fit, I nedd two points or more
    if(argc == 2)
    {
        PrintColor("Error: modify .cpp file or insert more input files", OBOLDRED);
        return EXIT_FAILURE;
    }

    char* fileName;

    // 8 arrays, one for each SiPM, contatinins means and errors on initial times
    std::vector <std::vector <double> >* means
        = new std::vector <std::vector <double> > {{}, {}, {}, {}, {}, {}, {}, {}};
    std::vector <std::vector <double> >* stdDevs
        = new std::vector <std::vector <double> > {{}, {}, {}, {}, {}, {}, {}, {}};
    
    // 2 arrays, one for each side of SiPMs, containing means and errors on initial times
    std::vector <std::vector <double> >* means2
        = new std::vector <std::vector <double> > {{}, {}};
    std::vector <std::vector <double> >* stdDevs2
        = new std::vector <std::vector <double> > {{}, {}};
    
    // positions in x axis for each input file
    std::vector <double> positions_x = {-90, -80., -70., -60., -50., -30., 0., 30., 50., 60., 70., 80., 90.};
    for (auto x = positions_x.begin(); x != positions_x.end(); ++x) *x *= HALF_LEN_X/100.;
    
    // loop over input files
    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {

        fileName = argv[file_counter];
        
        // Plot histos of arrival times and append the mean & error to relative vectors
        if (file_counter == 1)
        {
            plotHisto_arrivalTimes(fileName, OPEN_OUTPUT, means, stdDevs);
            plotHisto_arrivalTimes2(fileName, OPEN_OUTPUT, means2, stdDevs2);
        }
        else if (file_counter == argc-1)
        {
            plotHisto_arrivalTimes(fileName, CLOSE_OUTPUT, means, stdDevs);
            plotHisto_arrivalTimes2(fileName, CLOSE_OUTPUT, means2, stdDevs2);
        }
        else
        {
            plotHisto_arrivalTimes(fileName, ADD_OUTPUT, means, stdDevs);
            plotHisto_arrivalTimes2(fileName, ADD_OUTPUT, means2, stdDevs2);
        }

    }

    // Print means and errors on the shell interface
    PrintFitResults(means, stdDevs, argv);
    PrintFitResults2(means2, stdDevs2, argv);

    // The correlation function between time difference dx-sx and the position of incoming particle
    TH1F* histo_deltaX = new TH1F("histo_dx_time","Histogram of x_{rec} - x_{true}", nbins, -HALF_LEN_X, HALF_LEN_X);
    TF1* correlation_function = PlotFitResults2(means2, stdDevs2, positions_x);
    double_t coefficient, quote;
    coefficient = correlation_function->GetParameter(1);
    quote = correlation_function->GetParameter(0);
    // funct is t(x). I invert to obtain x(x):
    // in principle I should propagate errors on the new coefficients but at the moment I do not need it
    correlation_function->SetParameter(0, -quote/coefficient);
    correlation_function->SetParameter(1, 1./coefficient);

    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {
        fileName = argv[file_counter];
        HistoFillDeltaXperFile(histo_deltaX, correlation_function, fileName, positions_x[file_counter-1]);
    }

    PlotHistogramDeltaX(histo_deltaX, "images/t_vs_x.pdf)");

    delete correlation_function;
    delete histo_deltaX;
    delete stdDevs2;
    delete means2;
    delete stdDevs;
    delete means;

    return EXIT_SUCCESS;
}

void PlotHistogramDeltaX(TH1F* histo_deltaX, TString output_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 4000);
    TF1* gauss_fit = new TF1();

    histo_deltaX->Draw("");
    histo_deltaX->Draw("SAME E1");
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
    
    canva->Print("images/t_vs_x.pdf)","pdf");
    canva->Clear();
  
    delete gauss_fit;
    delete canva;

}

void HistoFillDeltaXperFile(TH1F* histogram, TF1* function, char * fileName,  double true_x){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("firstTimes", myFile);
    std::vector < TTreeReaderValue <Double_t> > times;

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "1stTimeSiPM[" + ChannelName + "]";
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