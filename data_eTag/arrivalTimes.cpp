#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "libraries/arrival_times.cc"

// ROOT header files
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#define OPEN_OUTPUT 0
#define CLOSE_OUTPUT 1
#define SINGLE_OUTPUT 2
#define ADD_OUTPUT 3

#define HALF_LEN_X 300 // [mm]
#define HALF_LEN_Y 22.5 // [mm]

void PrintFitResults(std::vector <std::vector <double> >* means, std::vector <std::vector <double> >* stdDevs, char** files);
void PlotFitResults(std::vector <std::vector <double> >* means, std::vector <std::vector <double> >* stdDevs, std::vector <double> positions_x);

int main(int argc, char** argv){

    char* fileName;

    // 8 arrays: one for each SiPm
    std::vector <std::vector <double> >* means
        = new std::vector <std::vector <Double_t> > {{}, {}, {}, {}, {}, {}, {}, {}};
    std::vector <std::vector <double> >* stdDevs
        = new std::vector <std::vector <Double_t> > {{}, {}, {}, {}, {}, {}, {}, {}};
    
    // positions in x axis for each input file
    std::vector <double> positions_x = {-90, -70., -50., -30., 0., 30., 50., 70., 90.};
    for (auto x = positions_x.begin(); x != positions_x.end(); ++x) *x *= HALF_LEN_X/100.;

    if(argc == 2) print_arrivalTimes(argv[1], SINGLE_OUTPUT, means, stdDevs);
    else
    {

        for(int file_counter = 1; file_counter < argc; ++file_counter)
        {

            fileName = argv[file_counter];

            if (file_counter == 1) print_arrivalTimes(fileName, OPEN_OUTPUT, means, stdDevs);
            else if (file_counter == argc-1) print_arrivalTimes(fileName, CLOSE_OUTPUT, means, stdDevs);
            else print_arrivalTimes(fileName, ADD_OUTPUT, means, stdDevs);

        }
    }

    PrintFitResults(means, stdDevs, argv);
    PlotFitResults(means, stdDevs, positions_x);

    return EXIT_SUCCESS;
}

void PrintFitResults(std::vector <std::vector <Double_t> >* means, std::vector <std::vector <Double_t> >* stdDevs, char** files){

    std::cout << OBOLDYELLOW << "Statistical results of initial times:" << ORESET << std::endl;
    
    // loop over the entries (one per file)
    for(int entry = 0; entry < (*means)[0].size(); ++entry)
    {
        std::cout << "Result of the file " << files[entry+1] << std::endl;
        
        for(int channel = 0; channel < numberOfChannels; ++channel)
        {
            std::cout << "\tSiPM # [" << channel+1 << "]\t"
                << "mean: " << (*means)[channel][entry]
                << "\tstd deviation: " << (*stdDevs)[channel][entry] << std::endl;
        }

        std::cout << std::endl;
    }
}

void PlotFitResults(std::vector <std::vector <double> >* means, std::vector <std::vector <double> >* stdDevs, std::vector <double> positions_x){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3800, 3600);
    canva->Divide(2,2);
    const int color[8] = {kBlue+3, kBlue+3, kBlue+3, kBlue+3, kOrange+9, kOrange+9, kOrange+9, kOrange+9};

    auto graph = new TGraphErrors();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    for(int channel = 0; channel < numberOfChannels; ++channel)
    {
        canva->cd((channel%4)+1);

        for(int entry = 0; entry < noOfPoints; ++entry)
        {
            x[entry] = positions_x[entry];
            y[entry] = (*means)[channel][entry];
            dx[entry] = 0;
            dy[entry] = (*stdDevs)[channel][entry];
        }

        graph = new TGraphErrors(noOfPoints, x, y, dx, dy);
        std::string title = "Beam position (x,0) vs time of first detected #gamma. SiPM #" + std::to_string(channel+1);
        TString Ttitle = title;
        graph->SetTitle(Ttitle);
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlack);
        graph->SetMarkerSize(3.);
        graph->GetXaxis()->CenterTitle();
        graph->GetXaxis()->SetTitle("position x [mm]");
        graph->GetYaxis()->SetTitle("Time [ns]");
        graph->SetLineColor(color[channel]);
        graph->Draw("AL*");

        
        canva->Update();

        if (channel == (int)numberOfChannels/2 - 1)
        {
            canva->Print("images/fit_line.pdf(","pdf");
            canva->Clear();
            canva->Divide(2,2);
        } else if (channel == numberOfChannels-1) canva->Print("images/fit_line.pdf)","pdf");

    }

}