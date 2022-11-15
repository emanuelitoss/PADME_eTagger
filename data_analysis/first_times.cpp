#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "libraries/times_plot.cc"
#include "libraries/times_fit.cc"
#include "libraries/infos.h"
#include "libraries/xt_distribution.cc"

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
    
    /***************** INSPECTING FILES *****************/

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

    /***************** FITTING CORRELATION LINE *****************/

    // Print means and errors on the shell interface
    PrintFitResults(means, stdDevs, argv);
    PrintFitResults2(means2, stdDevs2, argv);

    /***************** RECONSTRUCTION OF THE DISTRIBUTION OF TRUE_X - RECO_X *****************/

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

    PlotHistogramDeltaX(histo_deltaX, "images/t_vs_x.pdf");

    /***************** TESTING OF THE DISTRIBUTION OF TRUE_X - RECO_X OVER RANDOM (x,y) DATA *****************/

    histo_deltaX->Clear();
    HistoFillDeltaXRandomFile(histo_deltaX, correlation_function, "data_eTagRAND.root");
    PlotHistogramDeltaX(histo_deltaX, "images/t_vs_x.pdf)");

    /***************** EXIT PROGRAM *****************/

    delete correlation_function;
    delete histo_deltaX;
    delete stdDevs2;
    delete means2;
    delete stdDevs;
    delete means;

    return EXIT_SUCCESS;
}
