#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
using std::cout;
using std::endl;

#include "libraries/plot_histos.cc"

// ROOT header files
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"

#define OPEN_OUTPUT 0
#define CLOSE_OUTPUT 1
#define SINGLE_OUTPUT 2
#define ADD_OUTPUT 3

int main(int argc, char** argv){

    char* fileName;

    if(argc == 2) print_histos(argv[1], SINGLE_OUTPUT);
    else
    {

        for(int file_counter = 1; file_counter < argc; ++file_counter)
        {

        fileName = argv[file_counter];

        if (file_counter == 1) print_histos(fileName, OPEN_OUTPUT);
        else if (file_counter == argc-1) print_histos(fileName, CLOSE_OUTPUT);
        else print_histos(fileName, ADD_OUTPUT);

        }
    }

    return EXIT_SUCCESS;
}