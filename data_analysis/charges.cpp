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
#define nbins 70

void HistoFillDeltaXperFileCharges(TH1F* histogram, TF1* correlation_function, char* fileName, double true_x);
void PlotHistogramDeltaX(TH1F* histo_deltaX, TString output_name);
void HistoFillDeltaXRandomFileCharges(TH1F* histogram, TF1* function, TString fileName);

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
    vector <double_t> positions_x = {-90, -80., -70., -60., -50., -40., -30., -20., -10., 0., 10., 20., 30., 40., 50., 60., 70., 80., 90.};
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

    /***************** RECONSTRUCTION OF THE DISTRIBUTION OF TRUE_X - RECO_X *****************/

    // The correlation function between time difference dx-sx and the position of incoming particle
    TH1F* histo_deltaX = new TH1F("histo_dx_time","Histogram of x_{rec} - x_{true}", nbins, -HALF_LEN_X, HALF_LEN_X);
    TF1* correlation_function = PlotChargesFunctions(charges_functions_means, charges_functions_stdDevs, positions_x, EPE);
    for(int file_counter = 1; file_counter < argc; ++file_counter)
    {
        fileName = argv[file_counter];
        HistoFillDeltaXperFileCharges(histo_deltaX, correlation_function, fileName, positions_x[file_counter-1]);
    }

    if(EPE) PlotHistogramDeltaX(histo_deltaX, "images/chargesEpE.pdf");
    else PlotHistogramDeltaX(histo_deltaX, "images/charges.pdf");

    /***************** TESTING OF THE DISTRIBUTION OF TRUE_X - RECO_X OVER RANDOM (x,y) DATA *****************/
    histo_deltaX->Clear();
    histo_deltaX = new TH1F("histo_dx_charge","Histogram of x_{rec} - x_{true}", nbins, -HALF_LEN_X, HALF_LEN_X);
    HistoFillDeltaXRandomFileCharges(histo_deltaX, correlation_function, "data_eTagRAND.root");
    if(EPE) PlotHistogramDeltaX(histo_deltaX, "images/chargesEpE.pdf)");
    else PlotHistogramDeltaX(histo_deltaX, "images/charges.pdf)");

    delete correlation_function;
    delete histo_deltaX;
    delete charges_functions_stdDevs;
    delete charges_functions_means;
    delete charges_stdDevs;
    delete charges_means;

    return EXIT_SUCCESS;
}

void HistoFillDeltaXRandomFileCharges(TH1F* histogram, TF1* function, TString fileName){
        
    // reading objects
    TFile* myFile = TFile::Open(fileName);

    TTreeReader readerQ = TTreeReader();
    readerQ.SetTree("charges", myFile);

    TTreeReader readerX = TTreeReader();
    readerX.SetTree("positions", myFile);
    TTreeReaderValue <Double_t> position = TTreeReaderValue <Double_t>(readerX,"X_position");

    std::vector < TTreeReaderValue <Double_t> > charges;
    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "Charges[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      charges.push_back(TTreeReaderValue<Double_t>(readerQ,dirName));
    }

    double_t q = 0, chargeSX = 0, chargeDX = 0;
    double_t reconstructed_x, true_x;

    while(readerQ.Next() && readerX.Next())
    {

        chargeSX = 0;
        chargeDX = 0;

        true_x = *position;

        for (int channel = 0; channel < numberOfChannels; channel++)
        {
            q = *charges[channel];
            
            if(channel < numberOfChannels/2) chargeDX+=q;
            else chargeSX+=q;
        }

        reconstructed_x = function->GetX(chargeDX-chargeSX);
        histogram->Fill(reconstructed_x - true_x);

    }

    myFile->Close();
    delete myFile;
}

void PlotHistogramDeltaX(TH1F* histogram, TString output_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 4000);
    TF1* gauss_fit = new TF1();

    canva->SetGrid();

    histogram->Draw("EP");
    histogram->GetXaxis()->SetTitle("Length [mm]");
    histogram->GetYaxis()->SetTitle("Number of events");
    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth((Width_t)1.5);

    gStyle->Reset();
    gStyle->SetEndErrorSize(8);
    gStyle->SetOptFit(1110);
    gStyle->SetOptStat(2210);
    gStyle->SetStatX(0.89);
    gStyle->SetStatY(0.89);
    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.7);
    
    gauss_fit = new TF1("fitting a gaussian", "gaus", -HALF_LEN_X, HALF_LEN_X);
    histogram->Fit(gauss_fit, "L Q", "0");
    gauss_fit->SetLineColor(kAzure-5);
    gauss_fit->SetLineWidth(1);
    gauss_fit->SetFillStyle(3002);
    gauss_fit->SetFillColorAlpha(kAzure-5,0.5);
    gauss_fit->Draw("C SAME");
    
    canva->Print(TString(output_name),"pdf");
    canva->Clear();
  
    delete gauss_fit;
    delete canva;

}

void HistoFillDeltaXperFileCharges(TH1F* histogram, TF1* correlation_function, char* fileName, double true_x){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("charges", myFile);
    std::vector < TTreeReaderValue <Double_t> > charges;

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "Charges[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      charges.push_back(TTreeReaderValue<Double_t>(reader,dirName));
    }

    double_t q = 0, chargeSX = 0, chargeDX = 0;
    double_t reconstructed_x;

    while(reader.Next())
    {

        chargeSX = OFFSET;
        chargeDX = OFFSET;

        for (int channel = 0; channel < numberOfChannels; channel++)
        {
            q = *charges[channel];

            if(channel < numberOfChannels/2) chargeDX+=q;
            else chargeSX +=q;
            
        }

        reconstructed_x = correlation_function->GetX(chargeDX-chargeSX);
        histogram->Fill(reconstructed_x - true_x);

    }

    return;
}