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
    TF1* correlation_functionQ_difference;
    TF1* correlation_functionQ_ratio;
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

    /***************** PRINTING PROGRESS AND FIRST PLOTS, SAVING FITTING FUNCTIONS *****************/

    /***** CHARGES *****/

    // if I'm not acting through "event per event" method, I evaluate functions of charges at the end
    if(!EPE) AddFunctionsOfCharges(charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);
    PlotCharges(charges_means, charges_stdDevs, positions_x, EPE);
    correlation_functionQ_difference = PlotChargesFunctions(charges_functions_means, charges_functions_stdDevs, positions_x, EPE, OPTION_Q_DIFFERENCE);
    correlation_functionQ_ratio = PlotChargesFunctions(charges_functions_means, charges_functions_stdDevs, positions_x, EPE, OPTION_Q_RATIO);

    /***** TIMES *****/

    // Print means and errors for time analysis on the shell interface
    PrintFitResults(times_means8, times_stdDevs8, argv);
    PrintFitResults2(times_means2, times_stdDevs2, argv);
    correlation_functionT = PlotFitResults2(times_means2, times_stdDevs2, positions_x);

    /***************** DEFINITIONS FOR X_T VS X_Q ANALYSIS *****************/

    // I prepare a vector of doublet delta_x(charge analysis), delta_x(time analysis)
    vector <Double_t> * deltasT = new vector <Double_t> {};
    vector <Double_t> * deltasQ_diff = new vector <Double_t> {};
    vector <Double_t> * deltasQ_ratio = new vector <Double_t> {};

    TH1F* histo_deltaXT = new TH1F("histo_delta_x(T)","Histogram of x_{rec} - x_{true} from time analysis", nbins, -HALF_LEN_X, HALF_LEN_X);
    TH1F* histo_deltaXQ_diff = new TH1F("histo_delta_x(Q)_differences","Histogram of x_{rec} - x_{true} from charge differences", nbins, -HALF_LEN_X, HALF_LEN_X);
    TH1F* histo_deltaXQ_ratio = new TH1F("histo_delta_x(Q)_ratios","Histogram of x_{rec} - x_{true} from charge ratios", nbins, -HALF_LEN_X, HALF_LEN_X);

    // inverting a line function >> Applied only to time-position function.
    // for the charge-position function I use an other method
    double_t coefficient, quote;
    coefficient = correlation_functionT->GetParameter(1);
    quote = correlation_functionT->GetParameter(0);
    // funct is t(x). I invert to obtain x(t):
    correlation_functionT->SetParameter(0, -quote/coefficient);
    correlation_functionT->SetParameter(1, 1./coefficient);

    /***************** DISTRIBUTION OF TRUE_X - RECO_X WITH PREVIOUS DATA *****************/

    // evaluate delta_x and fill the histograms
    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {
        fileName = argv[file_counter];
        HistoFillDeltaXperFileCharges(histo_deltaXQ_diff, correlation_functionQ_difference, fileName, positions_x[file_counter-1], deltasQ_diff, OPTION_Q_DIFFERENCE, OPTION_DELTA_POSITION);
        HistoFillDeltaXperFileCharges(histo_deltaXQ_ratio, correlation_functionQ_ratio, fileName, positions_x[file_counter-1], deltasQ_ratio, OPTION_Q_RATIO, OPTION_DELTA_POSITION);
        HistoFillDeltaXperFileTimes(histo_deltaXT, correlation_functionT, fileName, positions_x[file_counter-1], deltasT, OPTION_DELTA_POSITION);
    }

    // plot the histograms
    PlotHistogramDeltaXTimes(histo_deltaXT, "images/2t_vs_x.pdf");
    if(EPE){
        PlotHistogramDeltaXCharges(histo_deltaXQ_diff, "images/3chargesEpE.pdf");
        PlotHistogramDeltaXCharges(histo_deltaXQ_ratio, "images/3chargesEpE.pdf");
    }
    else{
        PlotHistogramDeltaXCharges(histo_deltaXQ_diff, "images/charges.pdf");
        PlotHistogramDeltaXCharges(histo_deltaXQ_ratio, "images/charges.pdf");
    }

    /***************** X_CHARGES VS X_TIMES SCATTERPLOT *****************/

    deltasQ_diff->clear();
    deltasQ_ratio->clear();
    deltasT->clear();

    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {
        fileName = argv[file_counter];
        HistoFillDeltaXperFileCharges(histo_deltaXQ_diff, correlation_functionQ_difference, fileName, positions_x[file_counter-1], deltasQ_diff, OPTION_Q_DIFFERENCE, OPTION_POSITION);
        HistoFillDeltaXperFileCharges(histo_deltaXQ_ratio, correlation_functionQ_ratio, fileName, positions_x[file_counter-1], deltasQ_ratio, OPTION_Q_RATIO, OPTION_POSITION);
        HistoFillDeltaXperFileTimes(histo_deltaXT, correlation_functionT, fileName, positions_x[file_counter-1], deltasT, OPTION_POSITION);
    }

    // Plot delta_x(times) vs delta_x(charges) 2D histogram
    PlotPositionsQT(deltasT, deltasQ_diff, SINGLE_OUTPUT, OPTION_Q_DIFFERENCE);
    // PlotPositionsQT(deltasT, deltasQ_ratio, ADD_OUTPUT, OPTION_Q_RATIO);

    /***************** TESTING OF THE CORRELATION X_T vs X_Q OVER RANDOM (x,y) DATA *****************/

    // Here I append to previous histograms, other tests.
    // I use randomic (x,y) data to obtain a simulated realistic reconstruction

    deltasQ_diff->clear();
    histo_deltaXQ_diff->Clear();
    histo_deltaXQ_diff = new TH1F("histo_dx_charge_differences","Histogram of x_{rec} - x_{true} (charges analysis - differences)", nbins, -HALF_LEN_X, HALF_LEN_X);
    HistoFillDeltaXRandomFileCharges(histo_deltaXQ_diff, correlation_functionQ_difference, "data_eTagRAND.root", deltasQ_diff);
    if(EPE) PlotHistogramDeltaXCharges(histo_deltaXQ_diff, "images/3chargesEpE.pdf");
    else PlotHistogramDeltaXCharges(histo_deltaXQ_diff, "images/charges.pdf");

    deltasQ_ratio->clear();
    histo_deltaXQ_ratio->Clear();
    histo_deltaXQ_ratio = new TH1F("histo_dx_charge_ratios","Histogram of x_{rec} - x_{true} (charges analysis - ratios)", nbins, -HALF_LEN_X, HALF_LEN_X);
    HistoFillDeltaXRandomFileCharges(histo_deltaXQ_ratio, correlation_functionQ_ratio, "data_eTagRAND.root", deltasQ_ratio);
    if(EPE) PlotHistogramDeltaXCharges(histo_deltaXQ_ratio, "images/3chargesEpE.pdf)");
    else PlotHistogramDeltaXCharges(histo_deltaXQ_ratio, "images/charges.pdf)");

    deltasT->clear();
    histo_deltaXT->Clear();
    histo_deltaXT = new TH1F("histo_dx_time","Histogram of x_{rec} - x_{true} from times analysis", nbins, -HALF_LEN_X, HALF_LEN_X);
    HistoFillDeltaXRandomFileTimes(histo_deltaXT, correlation_functionT, "data_eTagRAND.root", deltasT);
    PlotHistogramDeltaXTimes(histo_deltaXT, "images/2t_vs_x.pdf)");
    
    /***************** DELTAX_CHARGES VS DELTAX_TIMES OVER RANDOM (x,y) DATA *****************/

    // COMMENTO MIO: Aggiustare ed inserire le posizioni ricostruite, non le delta!
    // deltasT->clear();
    // deltasQ_diff->clear();
    // deltasQ_ratio->clear();
    // plot delta_x(times) vs delta_x(charges) 2D histogram
    // PlotPositionsQT(deltasT, deltasQ_diff, ADD_OUTPUT);
    // PlotPositionsQT(deltasT, deltasQ_ratio, CLOSE_OUTPUT);

    /***************** DELETE STUFF & EXIT *****************/

    delete histo_deltaXQ_ratio;
    delete histo_deltaXQ_diff;
    delete histo_deltaXT;
    delete deltasQ_ratio;
    delete deltasQ_diff;
    delete deltasT;
    delete correlation_functionT;
    delete correlation_functionQ_ratio;
    delete correlation_functionQ_difference;
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