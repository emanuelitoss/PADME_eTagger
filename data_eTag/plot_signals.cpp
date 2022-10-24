#include "libraries/signals_plot.cc"
#include "libraries/infos.h"

// ROOT header files
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"

int main(int argc, char** argv){

    char* fileName;

    if(argc == 2) print_signals(argv[1], SINGLE_OUTPUT);
    else
    {

        for(int file_counter = 1; file_counter < argc; ++file_counter)
        {

        fileName = argv[file_counter];

        if (file_counter == 1) print_signals(fileName, OPEN_OUTPUT);
        else if (file_counter == argc-1) print_signals(fileName, CLOSE_OUTPUT);
        else print_signals(fileName, ADD_OUTPUT);

        }
    }

    return EXIT_SUCCESS;
}