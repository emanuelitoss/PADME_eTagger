#include "infos.h"

// ROOT header files
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"

void plotScatter_arrivalTimes(int argc, std::vector <double> original_posX, char** argv){

    /********** READ FILE AND INITIALIZE READERS **********/

    int file_counter, numOfFiles = argc-1;
    std::vector <Double_t> time_differences = {};
    std::vector <Double_t> positions_x = {};

    for(file_counter = 1; file_counter<numOfFiles; ++file_counter)
    {
        // reading objects
        TFile* myFile = TFile::Open(argv[file_counter]);
        TTreeReader reader = TTreeReader();
        reader.SetTree("firstTimes", myFile);
        std::vector< TTreeReaderValue<Double_t> > times;

        // to manipulate and save data
        const int DX = 0, SX = 1;
        std::vector <std::vector <Double_t> > min_times = {{}, {}};

        // initializing readers
        for (int channel = 0; channel < numberOfChannels; ++channel)
        {
            // readers
            std::string ChannelName = std::to_string(channel+1);
            std::string dirNameS = "1stTimeSiPM[" + ChannelName + "]";
            const char* dirName = dirNameS.c_str();
            times.push_back(TTreeReaderValue<Double_t>(reader,dirName));
        }

        double_t time = 0;
        int counter = 0, entry = 0;
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

            // This because I know how the numbers are stored.
            // Inspect a .root output to understand
            if(counter % 8 == 0){
                time_differences.push_back(min_times[DX].back() - min_times[SX].back());
                positions_x.push_back(original_posX[file_counter-1]);
            }
        }

        myFile->Close();
        delete myFile;
    }

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3800, 3500);

    int noOfPoints = positions_x.size();
    double t[noOfPoints], x[noOfPoints], dt[noOfPoints], dx[noOfPoints];
    for(int i = 0; i<noOfPoints; ++i){
        t[i] = time_differences[i];
        x[i] = positions_x[i];
        dt[i] = 0;
        dx[i] = 0;
    }
    auto graph = new TGraphErrors(noOfPoints, t, x, dx, dt);

    // axis
    std::string title = "Beam position (x,0) vs Delta_time of first detected #gamma (sides method)";
    TString Ttitle = title;
    graph->SetTitle(Ttitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    graph->GetXaxis()->SetTitle("Time [ns]");
    graph->GetYaxis()->SetTitle("position x [mm]");

    // markers
    graph->SetMarkerStyle(kDot);
    graph->SetMarkerColor(kAzure-6);
    graph->SetMarkerSize(2.);

    graph->Draw("AP");

    canva->SetGrid();
    canva->Print("images/scatter_deltatimes.pdf","pdf");

    delete graph;
    delete canva;

}
#include "TAttMarker.h"