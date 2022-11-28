#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "libraries/infos.h"
#include "libraries/charge_analysis.cc"

#define EPE true // if true, it uses event per event method

int main(int argc, char** argv){

    char* fileName;

    // 2 arrays: one for each side of SiPM
    vector <vector <double_t> >* charges_means
        = new vector <vector <double_t> > {{}, {}};
    vector <vector <double_t> >* charges_stdDevs
        = new vector <vector <double_t> > {{}, {}};

    // arrays for 
    vector <vector <double_t> >* charges_functions_means
        = new vector <vector <double_t> > {{}, {}, {}, {}};
    vector <vector <double_t> >* charges_functions_stdDevs
        = new vector <vector <double_t> > {{}, {}, {}, {}};

    // positions in x axis for each input file
    vector <double_t> positions_x = {-90, -80., -70., -60., -50., -30., 0., 30., 50., 60., 70., 80., 90.};
    for (auto x = positions_x.begin(); x != positions_x.end(); ++x) *x *= HALF_LEN_X/100.;
    
    if(argc <= 2)
    {
        PrintColor("Error: Insert input files (more than one)", OBOLDRED);
        return EXIT_FAILURE;
    }
    
    else
    {
        for(int file_counter = 1; file_counter < argc; ++file_counter)
        {
            fileName = argv[file_counter];
            if(EPE) AddChargesEPE(fileName, charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);
            else AddCharges(fileName, charges_means, charges_stdDevs);
        }
    }

    if(!EPE) AddFunctionsOfCharges(charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);

    PlotCharges(charges_means, charges_stdDevs, positions_x, EPE);
    PlotChargesFunctions(charges_functions_means, charges_functions_stdDevs, positions_x, EPE);

    delete charges_functions_stdDevs;
    delete charges_functions_means;
    delete charges_stdDevs;
    delete charges_means;

    return EXIT_SUCCESS;
}