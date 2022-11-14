#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "libraries/infos.h"
#include "libraries/charge_analysis.cc"

void AppendMeanStdDev(char * fileName, int openCloseFile, std::vector <std::vector <Double_t> >* means2, std::vector <std::vector <Double_t> >* stdDevs2);

int main(int argc, char** argv){

    if(argc == 2)
    {
        PrintColor("Error: modify .cpp file or insert more input files",OBOLDRED);
        return EXIT_FAILURE;
    }

    char* fileName;
    
    // 2 arrays, one for each side of SiPMs
    std::vector <std::vector <double> >* means2
        = new std::vector <std::vector <Double_t> > {{}, {}};
    std::vector <std::vector <double> >* stdDevs2
        = new std::vector <std::vector <Double_t> > {{}, {}};
    
   
    // positions in x axis for each input file
    std::vector <double> positions_x = {-90, -80., -70., -60., -50., -30., 0., 30., 50., 60., 70., 80., 90.};
    for (auto x = positions_x.begin(); x != positions_x.end(); ++x) *x *= HALF_LEN_X/100.;

    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {
        fileName = argv[file_counter];
        AppendMeanStdDev(fileName, OPEN_OUTPUT, means2, stdDevs2);

    }
    // 2 arrays: one for each side of SiPM
    vector <vector <double_t> >* charges_means
        = new vector <vector <double_t> > {{}, {}};
    vector <vector <double_t> >* charges_stdDevs
        = new vector <vector <double_t> > {{}, {}};

    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {
        fileName = argv[file_counter];
        AddCharges(fileName, charges_means, charges_stdDevs);
    }

    /************* PLOT ************/

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);
    const int color[2] = {kGreen+3, kOrange+9};

    TGraphErrors* graph = new TGraphErrors();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    for(int ch = 0; ch < 2; ++ch)
    {
        canva->SetGrid();

        for(int entry = 0; entry < noOfPoints; ++entry)
        {
            x[entry] = (*means2)[ch][entry];
            y[entry] = (*charges_means)[ch][entry];
            dx[entry] = (*stdDevs2)[ch][entry];
            dy[entry] = (*charges_means)[ch][entry];
        }

        graph = new TGraphErrors(noOfPoints, x, y, dx, dy);

        // markers
        graph->SetMarkerStyle(kCircle);
        graph->SetMarkerColor(color[ch]);
        graph->SetMarkerSize(2.);
        gStyle->SetEndErrorSize(8);

        // axis
        if(ch == 0) graph->SetTitle("Total charges vs initial time (right side)");
        else graph->SetTitle("Total charges vs initial time (left side)");
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->CenterTitle();
        graph->GetXaxis()->SetTitle("initial time t[ns]");
        graph->GetYaxis()->SetTitle("Number of detected #gamma");

        graph->Draw("ap");

        if(ch == 0) canva->Print("images/scatterTvsQ_x.pdf(","pdf");
        if(ch == 1) canva->Print("images/scatterTvsQ_x.pdf)","pdf");

        graph->Clear();
        canva->Clear();

    }

    delete graph;
    delete canva;

    delete charges_stdDevs;
    delete charges_means;
    delete stdDevs2;
    delete means2;

    return EXIT_SUCCESS;
}


void AppendMeanStdDev(char * fileName, int openCloseFile, std::vector <std::vector <Double_t> >* means2, std::vector <std::vector <Double_t> >* stdDevs2){

    /********** READ FILE AND INITIALIZE READERS **********/

    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("firstTimes", myFile);
    std::vector < TTreeReaderValue <Double_t> > times;

    // to manipulate and save data
    const int DX = 0, SX = 1;
    std::vector <std::vector <Double_t> > min_times = {{}, {}};

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "1stTimeSiPM[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      times.push_back(TTreeReaderValue<Double_t>(reader,dirName));
    }

    double_t time = 0, min_timeSX = 0, min_timeDX = 0;

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
    
      min_times[DX].push_back(min_timeDX);
      min_times[SX].push_back(min_timeSX);
    
    }

    (*means2)[DX].push_back(Mean(min_times[DX]));
    (*means2)[SX].push_back(Mean(min_times[SX]));
    (*stdDevs2)[DX].push_back(StdDeviation(min_times[DX]));
    (*stdDevs2)[SX].push_back(StdDeviation(min_times[SX]));

    myFile->Close();
    delete myFile;

    return;

}