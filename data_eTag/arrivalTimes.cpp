#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "libraries/arrival_times.cc"
#include "libraries/fit_plot_lines.cc"
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

    PlotFitResults(means, stdDevs, positions_x);
    PrintFitResults(means, stdDevs, argv);

    delete stdDevs;
    delete means;

    return EXIT_SUCCESS;
}