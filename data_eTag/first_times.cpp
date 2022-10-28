#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "libraries/times_plot.cc"
#include "libraries/times_fit.cc"
#include "libraries/infos.h"

// ROOT header files
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TObjString.h"

int main(int argc, char** argv){

    char* fileName;

    // 8 arrays: one for each SiPM
    std::vector <std::vector <double> >* means
        = new std::vector <std::vector <Double_t> > {{}, {}, {}, {}, {}, {}, {}, {}};
    std::vector <std::vector <double> >* stdDevs
        = new std::vector <std::vector <Double_t> > {{}, {}, {}, {}, {}, {}, {}, {}};
    
    // 2 arrays, one for each side of SiPMs
    std::vector <std::vector <double> >* means2
        = new std::vector <std::vector <Double_t> > {{}, {}};
    std::vector <std::vector <double> >* stdDevs2
        = new std::vector <std::vector <Double_t> > {{}, {}};
    
   
    // positions in x axis for each input file
    std::vector <double> positions_x = {-90, -70., -50., -30., 0., 30., 50., 70., 90.};
    for (auto x = positions_x.begin(); x != positions_x.end(); ++x) *x *= HALF_LEN_X/100.;

    if(argc == 2)
    {
        PrintColor("Error: modify .cpp file or insert more input files",OBOLDRED);
        return EXIT_FAILURE;
    }
    
    else
    {

        for(int file_counter = 1; file_counter < argc; ++file_counter)
        {

            fileName = argv[file_counter];

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
    }

    PlotFitResults(means, stdDevs, positions_x);
    PlotFitResults2(means2, stdDevs2, positions_x);
    PrintFitResults(means, stdDevs, argv);
    PrintFitResults2(means2, stdDevs2, argv);

    delete stdDevs2;
    delete means2;
    delete stdDevs;
    delete means;

    return EXIT_SUCCESS;
}