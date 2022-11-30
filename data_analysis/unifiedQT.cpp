#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "libraries/infos.h"
#include "libraries/charge_analysis.cc"
#include "libraries/xt_distribution.cc"
#include "libraries/times_plot.cc"
#include "libraries/times_fit.cc"

#define EPE true // if true, it uses event per event method

int main(int argc, char** argv){

    char* fileName;

    if(argc <= 2)
    {
        PrintColor("Error: Insert input files (more than one)", OBOLDRED);
        return EXIT_FAILURE;
    }

    // positions in x axis for each input file
    vector <double_t> positions_x = {-90, -80., -70., -60., -50., -40., -30., -20., -10., 0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
    for (auto x = positions_x.begin(); x != positions_x.end(); ++x) *x *= HALF_LEN_X/100.;

    /***************** DECLARE ARRAYS FOR ANALYSIS *****************/

    // means of accumulated charges per side,  at each position x
    // Two entries, one for side of eTag
    vector <vector <double_t> >* charges_means
        = new vector <vector <double_t> > {{}, {}};
    vector <vector <double_t> >* charges_stdDevs
        = new vector <vector <double_t> > {{}, {}};

    // 4 functions of previous charges f(Q_dx,Q,sx) at each position x
    vector <vector <double_t> >* charges_functions_means
        = new vector <vector <double_t> > {{}, {}, {}, {}};
    vector <vector <double_t> >* charges_functions_stdDevs
        = new vector <vector <double_t> > {{}, {}, {}, {}};

    // 2 arrays, one for each side of SiPMs, containing means and errors on initial times at each position x
    std::vector <std::vector <double> >* times_means2
        = new std::vector <std::vector <double> > {{}, {}};
    std::vector <std::vector <double> >* times_stdDevs2
        = new std::vector <std::vector <double> > {{}, {}};

    // 8 arrays, one for each SiPM, contatinins means and errors on initial times at each position x
    std::vector <std::vector <double> >* times_means8
        = new std::vector <std::vector <double> > {{}, {}, {}, {}, {}, {}, {}, {}};
    std::vector <std::vector <double> >* times_stdDevs8
        = new std::vector <std::vector <double> > {{}, {}, {}, {}, {}, {}, {}, {}};

    // fitting functions for (x,chargeDX-chargeSX) and (x,timeDX-timeSX)
    TF1* correlation_functionQ;
    TF1* correlation_functionT;

    /***************** FILLING ARRAYS *****************/

    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {

        fileName = argv[file_counter];

        // filling charges. If EPE is true, the functions of charges are evaluated event per event.
        // If not, they are evaluated at the end through error propagation.
        if(EPE) AddChargesEPE(fileName, charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);
        else AddCharges(fileName, charges_means, charges_stdDevs);

        // Similar for time analysis.
        // OPEN/CLOSE/ADD_OUTPUT variable tell to the function how to treat printing on external pdf and whenever close it or create a new one.
        if (file_counter == 1)
        {
            plotHisto_arrivalTimes(fileName, OPEN_OUTPUT, times_means8, times_stdDevs8);
            plotHisto_arrivalTimes2(fileName, OPEN_OUTPUT, times_means2, times_stdDevs2);
        }
        else if (file_counter == argc-1)
        {
            plotHisto_arrivalTimes(fileName, CLOSE_OUTPUT, times_means8, times_stdDevs8);
            plotHisto_arrivalTimes2(fileName, CLOSE_OUTPUT, times_means2, times_stdDevs2);
        }
        else
        {
            plotHisto_arrivalTimes(fileName, ADD_OUTPUT, times_means8, times_stdDevs8);
            plotHisto_arrivalTimes2(fileName, ADD_OUTPUT, times_means2, times_stdDevs2);
        }

    }

    /***************** PRINTING PROGRESS AND FIRST PLOTS *****************/

    // if I'm not acting through "event per event" method, I evaluate functions of charges at the end
    if(!EPE) AddFunctionsOfCharges(charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);
    PlotCharges(charges_means, charges_stdDevs, positions_x, EPE);
    correlation_functionQ = PlotChargesFunctions(charges_functions_means, charges_functions_stdDevs, positions_x, EPE);

    // Print means and errors for time analysis on the shell interface
    PrintFitResults(times_means8, times_stdDevs8, argv);
    PrintFitResults2(times_means2, times_stdDevs2, argv);
    correlation_functionT = PlotFitResults2(times_means2, times_stdDevs2, positions_x);

    /***************** RECONSTRUCTION OF THE DISTRIBUTION OF TRUE_X - RECO_X ("delta_x") *****************/

    // I prepare a vector of doublet delta_x(charge analysis), delta_x(time analysis)
    vector <Double_t> * deltasT = new vector <Double_t> {};
    vector <Double_t> * deltasQ = new vector <Double_t> {};

    TH1F* histo_deltaXT = new TH1F("histo_delta_x(T)","Histogram of x_{rec} - x_{true} from time analysis", nbins, -HALF_LEN_X, HALF_LEN_X);
    TH1F* histo_deltaXQ = new TH1F("histo_delta_x(Q)","Histogram of x_{rec} - x_{true} from charge analysis", nbins, -HALF_LEN_X, HALF_LEN_X);

    // inverting a line function >> Applied only to time-position function.
    // for the charge-position function I use an other method
    double_t coefficient, quote;
    coefficient = correlation_functionT->GetParameter(1);
    quote = correlation_functionT->GetParameter(0);
    // funct is t(x). I invert to obtain x(t):
    correlation_functionT->SetParameter(0, -quote/coefficient);
    correlation_functionT->SetParameter(1, 1./coefficient);

    // fill histograms
    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {
        fileName = argv[file_counter];
        HistoFillDeltaXperFileCharges(histo_deltaXQ, correlation_functionQ, fileName, positions_x[file_counter-1], deltasQ);
        HistoFillDeltaXperFileTimes(histo_deltaXT, correlation_functionT, fileName, positions_x[file_counter-1], deltasT);
    }

    // Plot histograms
    if(EPE) PlotHistogramDeltaXCharges(histo_deltaXQ, "images/chargesEpE.pdf");
    else PlotHistogramDeltaXCharges(histo_deltaXQ, "images/charges.pdf");
    PlotHistogramDeltaXTimes(histo_deltaXT, "images/t_vs_x.pdf");

    /***************** DELTAX_CHARGES VS DELTAX_TIMES *****************/

    // Plot delta_x(times) vs delta_x(charges) 2D histogram
    PlotDeltaPositionsQT(deltasT, deltasQ, OPEN_OUTPUT);

    /***************** TESTING OF THE DISTRIBUTION OF TRUE_X - RECO_X OVER RANDOM (x,y) DATA *****************/

    deltasQ->clear();
    histo_deltaXQ->Clear();
    histo_deltaXQ = new TH1F("histo_dx_charge","Histogram of x_{rec} - x_{true}", nbins, -HALF_LEN_X, HALF_LEN_X);
    HistoFillDeltaXRandomFileCharges(histo_deltaXQ, correlation_functionQ, "data_eTagRAND.root", deltasQ);
    if(EPE) PlotHistogramDeltaXCharges(histo_deltaXQ, "images/chargesEpE.pdf)");
    else PlotHistogramDeltaXCharges(histo_deltaXQ, "images/charges.pdf)");

    deltasT->clear();
    histo_deltaXT->Clear();
    histo_deltaXT = new TH1F("histo_dx_time","Histogram of x_{rec} - x_{true}", nbins, -HALF_LEN_X, HALF_LEN_X);
    HistoFillDeltaXRandomFileTimes(histo_deltaXT, correlation_functionT, "data_eTagRAND.root", deltasT);
    PlotHistogramDeltaXTimes(histo_deltaXT, "images/t_vs_x.pdf)");
    
    /***************** DELTAX_CHARGES VS DELTAX_TIMES OVER RANDOM (x,y) DATA *****************/

    // Plot delta_x(times) vs delta_x(charges) 2D histogram
    PlotDeltaPositionsQT(deltasT, deltasQ, CLOSE_OUTPUT);

    /***************** DELETE STUFF & EXIT *****************/

    delete histo_deltaXQ;
    delete histo_deltaXT;
    delete deltasQ;
    delete deltasT;
    delete correlation_functionT;
    delete correlation_functionQ;
    delete times_stdDevs8;
    delete times_means8;
    delete times_stdDevs2;
    delete times_means2;
    delete charges_functions_stdDevs;
    delete charges_functions_means;
    delete charges_stdDevs;
    delete charges_means;

    return EXIT_SUCCESS;
}