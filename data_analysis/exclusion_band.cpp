#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "libraries/infos.h"
#include "libraries/charge_analysis.cc"

int main(int argc, char** argv){

    if(argc > 2)
    {
        PrintColor("Error: Insert ONE input files", OBOLDRED);
        return EXIT_FAILURE;
    }

    char* fileName = argv[1];
    TFile* myFile = TFile::Open(fileName);

    // charges vs arrival time >> for scatterplot
    // channels 0,1 for left and right sides of SiPMs
    // channel 2 for the difference
    vector <vector <double_t> > charge_points = {{}, {}, {}};
    vector <vector <double_t> > time_points = {{}, {}, {}};

    /****************** READING & STORE DATA ******************/
    TTreeReader readerQ = TTreeReader();
    TTreeReader readerT = TTreeReader();
    readerQ.SetTree("charges", myFile);
    readerT.SetTree("arrival_times", myFile);
    vector< TTreeReaderValue<double_t> > charges;
    vector< TTreeReaderValue<double_t> > times;

    string ChannelName, dirNameStr;

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
        ChannelName = to_string(channel+1);

        dirNameStr = "Charges[" + ChannelName + "]";
        const char* dirNameQ = dirNameStr.c_str();
        charges.push_back(TTreeReaderValue<double_t>(readerQ,dirNameQ));

        dirNameStr = "arrival_times_SiPM[" + ChannelName + "]";
        const char* dirNameT = dirNameStr.c_str();
        times.push_back(TTreeReaderValue<double_t>(readerT,dirNameT));
    }

    double_t charge = 0, time = 0;
    double_t charge_dx, charge_sx, time_dx, time_sx;
    const double_t offset = 8;

    while ( readerQ.Next() && readerT.Next() )
    {
        charge_dx = 0;
        charge_sx = 0;
        time_dx = offset;
        time_sx = offset;

        for(int channel = 0; channel < numberOfChannels; channel++)
        {
            charge = *charges[channel];
            time = *times[channel];

            if(channel < numberOfChannels/2){
                charge_dx += charge;
                if(time_dx > time) time_dx = time;
            }
            else{
                charge_sx += charge;
                if(time_sx > time) time_sx = time;
            }
        }

        charge_points[0].push_back(charge_dx);
        charge_points[1].push_back(charge_sx);
        charge_points[2].push_back(charge_dx - charge_sx);
        time_points[0].push_back(time_dx);
        time_points[1].push_back(time_sx);
        time_points[2].push_back(time_dx - time_sx);
    }

    // check the right behavior of the program & data
    bool channel0Bool = charge_points[0].size() == time_points[0].size();
    bool channel1Bool = charge_points[1].size() == time_points[1].size();
    if(!channel0Bool || !channel1Bool)
    {
        cout << ORED << "Storing data problem! Exit failure." << ORESET << endl;
        return EXIT_FAILURE;
    }

    /****************** SCATTER PLOT ******************/

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);
    TF1* fit = nullptr;

    const int len = time_points.size();
    const int color[] = {kGreen+3, kOrange+9, kAzure-5, kAzure-5};
    const TString titles[] = {"Scatterplot: times vs charges (right side)",
                                "Scatterplot: times vs charges (left side)",
                                "Scatterplot: times vs charges (right side - left side)",
                                "Scatterplot: times vs charges (right side + left side)"};

    for(int ch=0; ch<len; ++ch)
    {
        int noOfPoints = charge_points[ch].size();
        double x[noOfPoints], y[noOfPoints];
        for(int entry = 0; entry < noOfPoints; ++entry)
        {
            x[entry] = time_points[ch][entry];
            y[entry] = charge_points[ch][entry];
        }

        canva->SetGrid();

        TGraph* graph = new TGraph(noOfPoints, x, y);

        gStyle->SetOptFit(1110);
        gStyle->SetOptStat(2210);
        gStyle->SetStatX(0.95);
        gStyle->SetStatY(0.87);

        // markers
        graph->SetMarkerStyle(43);
        graph->SetMarkerColor(color[ch]);
        graph->SetMarkerSize(4.);

        // axis
        graph->SetTitle(titles[ch]);
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->CenterTitle();
        graph->GetXaxis()->SetTitle("time [ns]");
        graph->GetYaxis()->SetTitle("Number of detected #gamma");

        graph->Draw("ap");
        
        if(ch<2) fit = new TF1("exponential function", "[1] + [2]*exp(-[3]*x)", graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
        else fit = new TF1("Poly3", "pol3", graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
        graph->Fit(fit, "Q", "0");
        fit->SetLineColor(kBlack);
        fit->SetLineWidth(1);
        fit->Draw("SAME");

        if(ch == 0) canva->Print("images/4time_vs_charge.pdf(","pdf");
        else if (ch == len-1) canva->Print("images/4time_vs_charge.pdf)","pdf");
        else canva->Print("images/4time_vs_charge.pdf","pdf");

        graph->Clear();
        canva->Clear();

        delete graph;
    }

    /****************** CLOSE & EXIT ******************/

    delete fit;
    delete canva;

    myFile->Close();
    delete myFile;

    return EXIT_SUCCESS;
}
