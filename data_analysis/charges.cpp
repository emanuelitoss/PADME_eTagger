#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TStyle.h"

#include "libraries/infos.h"

#define EVE true // if true, it uses event per event method

void AddCharges(char * filename, vector <vector <double_t> >* means, vector <vector <double_t> >* errs);
void AddChargesEPE(char * filename, vector <vector <double_t> >* means, vector <vector <double_t> >* errs,
     vector <vector <double_t> >* functions, vector <vector <double_t> >* functions_err);
void AddFunctionsOfCharges(vector <vector <double_t> >* means, vector <vector <double_t> >* errs,
    vector <vector <double_t> >* functions, vector <vector <double_t> >* functions_err);
double_t Mean(vector <double_t> vec);
double_t StdDeviation(vector <double_t> vec);
void PlotCharges(vector <vector <double> >* means, vector <vector <double> >* stdDevs, vector <double> positions_x);
void PlotChargesFunctions(vector <vector <double> >* means, vector <vector <double> >* stdDevs, vector <double> positions_x);

int main(int argc, char** argv){

    char* fileName;

    // 2 arrays: one for each side of SiPM
    vector <vector <double> >* charges_means
        = new vector <vector <double_t> > {{}, {}};
    vector <vector <double> >* charges_stdDevs
        = new vector <vector <double_t> > {{}, {}};

    // arrays for 
    vector <vector <double> >* charges_functions_means
        = new vector <vector <double_t> > {{}, {}, {}, {}};
    vector <vector <double> >* charges_functions_stdDevs
        = new vector <vector <double_t> > {{}, {}, {}, {}};

    // positions in x axis for each input file
    vector <double> positions_x = {-90, -80., -70., -60., -50., -30., 0., 30., 50., 60., 70., 80., 90.};
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
            if(EVE) AddChargesEPE(fileName, charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);
            else AddCharges(fileName, charges_means, charges_stdDevs);
        }
    }

    if(!EVE) AddFunctionsOfCharges(charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);

    PlotCharges(charges_means, charges_stdDevs, positions_x);
    PlotChargesFunctions(charges_functions_means, charges_functions_stdDevs, positions_x);

    delete charges_functions_stdDevs;
    delete charges_functions_means;
    delete charges_stdDevs;
    delete charges_means;

    return EXIT_SUCCESS;
}

double_t Mean(vector <double_t> vec){
    
    Double_t mean = 0;
    Int_t num = vec.size();

    for (int i = 0; i<num; ++i) mean += vec[i];

    mean/=num;
    return mean;
}

double_t StdDeviation(vector <double_t> vec){
    
    double_t err = 0, mean = Mean(vec);
    Int_t num = vec.size();

    for (int i = 0; i<num; ++i) err += vec[i]*vec[i];

    err /= num;
    err -= mean*mean;
    err = sqrt(err);

    return err;

}

void AddCharges(char * fileName, vector <vector <double_t> >* charge_means, vector <vector <double_t> >* charge_stdDevs){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("totalCharges", myFile);
    vector< TTreeReaderValue<double_t> > charges;

    cout << "\nInspecting file \"" << fileName << "\":" << endl;

    vector < vector <double_t> > charges_for_side = {{}, {}};

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // initialize readers
      string ChannelName = to_string(channel+1);
      string dirNameStr = "Charges[" + ChannelName + "]";
      const char* dirName = dirNameStr.c_str();
      charges.push_back(TTreeReaderValue<double_t>(reader,dirName));
    }

    double_t charge = 0;
    double_t chargeDX, chargeSX;

    while (reader.Next())
    {
        chargeDX = 0;
        chargeSX = 0;

        for(int channel = 0; channel < numberOfChannels; channel++)
        {
            charge = *charges[channel];
            if(channel < numberOfChannels/2) chargeDX += charge;
            else chargeSX += charge;
        }

        charges_for_side[0].push_back(chargeDX);
        charges_for_side[1].push_back(chargeSX);
    }

    vector <double_t> MeanCharges = {{}, {}};
    vector <double_t> ErrCharges = {{}, {}};

    for (int ch = 0; ch < 2; ++ch)
    {
        MeanCharges[ch] = Mean(charges_for_side[ch]);
        ErrCharges[ch] = StdDeviation(charges_for_side[ch]);

        (*charge_means)[ch].push_back(MeanCharges[ch]);
        (*charge_stdDevs)[ch].push_back(ErrCharges[ch]);

        if (ch == 0) cout << "\n\tRight side:\t";
        else cout << "\tLeft side:\t";
        cout << "Mean = " << MeanCharges[ch] << ",\t StdDev = " << ErrCharges[ch] << endl;
    }

    cout << endl;

}

void AddChargesEPE(char * fileName, vector <vector <double_t> >* charge_means, vector <vector <double_t> >* charge_stdDevs,
    vector <vector <double_t> >* functions, vector <vector <double_t> >* functions_err){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("totalCharges", myFile);
    vector< TTreeReaderValue<double_t> > charges;

    cout << "\nInspecting file \"" << fileName << "\":" << endl;

    vector < vector <double_t> > charges_for_side = {{}, {}};
    vector < vector <double_t> > charges_for_side_err = {{}, {}};
    // sum dx+sx, difference dx-sx, ratio dx/sx, ratio sx/dx
    vector < vector <double_t> > charges_functions = {{}, {}, {}, {}};
    vector < vector <double_t> > charges_functions_err = {{}, {}, {}, {}};

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // initialize readers
      string ChannelName = to_string(channel+1);
      string dirNameStr = "Charges[" + ChannelName + "]";
      const char* dirName = dirNameStr.c_str();
      charges.push_back(TTreeReaderValue<double_t>(reader,dirName));
    }

    double_t charge = 0;
    double_t chargeDX, chargeSX;
    double_t errorSUM, errorRATIO_relative;

    while (reader.Next())
    {
        chargeDX = 0;
        chargeSX = 0;

        for(int channel = 0; channel < numberOfChannels; channel++)
        {
            charge = *charges[channel];
            if(channel < numberOfChannels/2) chargeDX += charge;
            else chargeSX += charge;
        }

        charges_for_side[0].push_back(chargeDX);
        charges_for_side[1].push_back(chargeSX);
        // charges_for_side_err[0].push_back(sqrt(chargeDX));
        // charges_for_side_err[1].push_back(sqrt(chargeSX));

        charges_functions[0].push_back(chargeDX+chargeSX);
        charges_functions[1].push_back(chargeDX-chargeSX);
        charges_functions[2].push_back(chargeDX/chargeSX);
        charges_functions[3].push_back(chargeSX/chargeDX);

        /*
        errorSUM = sqrt(chargeDX + chargeSX);
        errorRATIO_relative = sqrt( 1./chargeDX + 1./chargeSX );

        charges_functions_err[0].push_back(errorSUM);
        charges_functions_err[1].push_back(errorSUM);
        charges_functions_err[2].push_back((chargeDX/chargeSX)*errorRATIO_relative);
        charges_functions_err[3].push_back((chargeSX/chargeDX)*errorRATIO_relative);*/

    }

    vector <double_t> MeanCharges = {{}, {}};
    vector <double_t> ErrCharges = {{}, {}};
    double_t value;
    Int_t size;

    for(int ch = 0; ch < 2; ++ch)
    {
        MeanCharges[ch] = Mean(charges_for_side[ch]);
        ErrCharges[ch] = StdDeviation(charges_for_side[ch]);
    
        (*charge_means)[ch].push_back(MeanCharges[ch]);
        (*charge_stdDevs)[ch].push_back(ErrCharges[ch]);

        if (ch == 0) cout << "\n\tRight side:\t";
        else cout << "\tLeft side:\t";
        cout << "Mean = " << MeanCharges[ch] << ",\t StdDev = " << ErrCharges[ch] << endl;
    }

    for(int func_idx = 0; func_idx < 4; ++func_idx)
    {
        (*functions)[func_idx].push_back(Mean(charges_functions[func_idx]));
        (*functions_err)[func_idx].push_back(StdDeviation(charges_functions[func_idx]));
    }

    cout << endl;

}

void AddFunctionsOfCharges(vector <vector <double_t> >* means, vector <vector <double_t> >* errs, vector <vector <double_t> >* functions, vector <vector <double_t> >* functions_err){

    Int_t numMeasures = (*means)[0].size();
    const int sum_idx = 0, diff_idx = 1, ratio_DS_idx = 2, ratio_SD_idx = 3;
    double_t err = 0;
    double_t val1, val2, err1, err2;

    for(int meas = 0; meas < numMeasures; ++meas)
    {
        val1 = (*means)[0][meas];
        val2 = (*means)[1][meas];
        err1 = (*errs)[0][meas];
        err2 = (*errs)[1][meas];

        // charge sum
        (*functions)[sum_idx].push_back(val1 + val2);
        err = sqrt( err1*err1 + err2*err2 );
        (*functions_err)[sum_idx].push_back(err);

        // charge difference
        (*functions)[diff_idx].push_back(val1 - val2);
        (*functions_err)[diff_idx].push_back(err);

        // charge ratio DX/SX
        (*functions)[ratio_DS_idx].push_back(val1/val2);
        err = (val1/val2)*sqrt( pow(err1/val1,2) + pow(err2/val2,2) );
        (*functions_err)[ratio_DS_idx].push_back(err);

        // charge ratio SX/DX
        (*functions)[ratio_SD_idx].push_back(val2/val1);
        err = (val2/val1)*sqrt( pow(err1/val1,2) + pow(err2/val2,2) );
        (*functions_err)[ratio_SD_idx].push_back(err);

    }
    
}

void PlotCharges(vector <vector <double> >* means, vector <vector <double> >* stdDevs, vector <double> positions_x){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);
    canva->Divide(2,1);
    const int color[2] = {kGreen+3, kOrange+9};

    TGraphErrors* graph = new TGraphErrors();
    TF1* exp_fit = new TF1();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    for(int ch = 0; ch < 2; ++ch)
    {
        canva->cd((ch%2==0)+1);
        canva->SetGrid();

        for(int entry = 0; entry < noOfPoints; ++entry)
        {
            x[entry] = positions_x[entry];
            y[entry] = (*means)[ch][entry];
            dx[entry] = 0;
            dy[entry] = (*stdDevs)[ch][entry];
        }

        graph = new TGraphErrors(noOfPoints, x, y, dx, dy);
        if(ch == 0) graph->SetTitle("DX SiPMs");
        else graph->SetTitle("SX SiPMs");

        // markers
        graph->SetMarkerStyle(kCircle);
        graph->SetMarkerColor(color[ch]);
        graph->SetMarkerSize(2.);
        graph->SetDrawOption("AP");
        gStyle->SetEndErrorSize(8);

        // axis
        string title = "Total charges vs Beam position (x,0)";
        TString Ttitle = title;
        graph->SetTitle(Ttitle);
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->CenterTitle();
        graph->GetXaxis()->SetTitle("position x [mm]");
        graph->GetYaxis()->SetTitle("Number of detected #gamma");

        graph->Draw();

        // exponential fit
        //syntax: TF1(const char* name, const char* formula, Double_t xmin = 0, Double_t xmax = 1)
        exp_fit = new TF1("exponential+offset", "[0] + [1]*exp(x*[2])", -HALF_LEN_X, HALF_LEN_X);
        graph->Fit(exp_fit, "0", "0");
        exp_fit->SetLineColor(color[ch]);
        exp_fit->SetLineWidth(1);
        exp_fit->SetFillStyle(3002);
        exp_fit->SetFillColorAlpha(color[ch],0.5);
        exp_fit->Draw("SAME");

        canva->Update();
        graph->Clear();
    }

    if(EVE) canva->Print("images/chargesEVE.pdf(","pdf");
    else canva->Print("images/charges.pdf(","pdf");

    delete exp_fit;
    delete graph;
    delete canva;

}

void PlotChargesFunctions(vector <vector <double> >* fmeans, vector <vector <double> >* fstdDevs, vector <double> positions_x){
    
    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);

    const int color[4] = {kAzure-5, kAzure-5, kViolet+6, kViolet+6};
    TString name[4] = {"Charges sum", "Charges difference", "Charges ratio", "Charges ratio"};
    TString yAxisName[4] = {"Sum N_{DX} + N_{SX}", "Difference N_{DX} - N_{SX}", "Ratio N_{DX} / N_{SX}", "Ratio N_{SX} / N_{DX}"};

    TGraphErrors* graph = new TGraphErrors();
    TF1* fit = new TF1();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    for(int ch = 0; ch < 4; ++ch)
    {
        canva->SetGrid();

        for(int entry = 0; entry < noOfPoints; ++entry)
        {
            x[entry] = positions_x[entry];
            y[entry] = (*fmeans)[ch][entry];
            dx[entry] = 0;
            dy[entry] = (*fstdDevs)[ch][entry];
        }

        graph = new TGraphErrors(noOfPoints, x, y, dx, dy);
        graph->SetTitle(name[ch]);
        
        // markers
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(color[ch]);
        graph->SetMarkerSize(4.);
        graph->SetDrawOption("AP");
        gStyle->SetEndErrorSize(8);

        // axis
        graph->SetTitle(name[ch]);
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->CenterTitle();
        graph->GetXaxis()->SetTitle("position x [mm]");
        graph->GetYaxis()->SetTitle(yAxisName[ch]);
        
     
        graph->Draw();

        // hyperbolic cosine + offset
        //if(ch == 0) fit = new TF1("f(x)+f(-x)", "2*[0] + [1]*(exp(x*[2])+exp(-x*[2]))", -HALF_LEN_X, HALF_LEN_X);
        fit = new TF1("f(x)+f(-x)", "pol6", -HALF_LEN_X, HALF_LEN_X);
        // hyperbolic cosine + offset
        //else if (ch==1) fit = new TF1("f(x)-f(-x)", "[0]*(exp(x*[1])-exp(-x*[1]))", -HALF_LEN_X, HALF_LEN_X);
        // hyperbolic tangent
        //else fit = new TF1("f(x)/f(-x)", "(1+[0]*exp(x*[1]))/(1+[0]*exp(-x*[1]))", -HALF_LEN_X, HALF_LEN_X);
        graph->Fit(fit, "0", "0");
        fit->SetLineColor(color[ch]);
        fit->SetLineWidth(1);
        fit->SetFillStyle(3002);
        fit->SetFillColorAlpha(color[ch],0.5);
        fit->Draw("SAME");

        canva->Draw();

        if(EVE)
        {
            if (ch == 3) canva->Print("images/chargesEVE.pdf)","pdf");
            else canva->Print("images/chargesEVE.pdf","pdf");
        }
        else{
            if (ch == 3) canva->Print("images/charges.pdf)","pdf");
            else canva->Print("images/charges.pdf","pdf");
        }

        graph->Clear();
        canva->Clear();

    }

    delete fit;
    delete graph;
    delete canva;

}