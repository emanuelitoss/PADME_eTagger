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

void AddCharges(char * filename, vector <vector <double_t> >* means, vector <vector <double_t> >* errs);
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
        = new vector <vector <double_t> > {{}, {}, {}};
    vector <vector <double> >* charges_functions_stdDevs
        = new vector <vector <double_t> > {{}, {}, {}};

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
            AddCharges(fileName, charges_means, charges_stdDevs);
        }
    }

    AddFunctionsOfCharges(charges_means, charges_stdDevs, charges_functions_means, charges_functions_stdDevs);

    PlotCharges(charges_means, charges_stdDevs, positions_x);
    PlotChargesFunctions(charges_functions_means, charges_functions_stdDevs, positions_x);

    delete charges_functions_stdDevs;
    delete charges_functions_means;
    delete charges_stdDevs;
    delete charges_means;

    return EXIT_SUCCESS;
}

double_t Mean(vector <double_t> vec){
    
    double_t mean = 0;

    for (int i = 0; i<vec.size(); ++i){
        mean += vec[i];
    }

    mean/=(vec.size());
    return mean;
}

double_t StdDeviation(vector <double_t> vec){
    
    double_t err = 0;
    double_t mean = Mean(vec);

    for (int i = 0; i<vec.size(); ++i){
        err += (vec[i]-mean)*(vec[i]-mean);
    }

    err = sqrt(err/vec.size());

    return err;

}

void AddCharges(char * fileName, vector <vector <double_t> >* charge_means, vector <vector <double_t> >* charge_stdDevs){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("totalCharges", myFile);
    vector< TTreeReaderValue<double_t> > charges;

    cout << "\nInspectin file \"" << fileName << "\":" << endl;

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
    double_t tot_charge = 0;

    while (reader.Next())
    {
        for(int channel = 0; channel < numberOfChannels; channel++)
        {
            charge = *charges[channel];
            if(charge != 0)
            {
                if(channel < numberOfChannels/2) charges_for_side[0].push_back(charge);
                else{
                    charges_for_side[1].push_back(charge);
                } 
                tot_charge += charge;
            }
        }
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

void AddFunctionsOfCharges(vector <vector <double_t> >* means, vector <vector <double_t> >* errs, vector <vector <double_t> >* functions, vector <vector <double_t> >* functions_err){

    Int_t numMeasures = (*means)[0].size();

    const int sum_idx = 0, diff_idx = 1, ratio_idx = 2;
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

        // charge ratio
        (*functions)[ratio_idx].push_back(val1/val2);
        err = (val1/val2)*sqrt( pow(err1/val1,2) + pow(err2/val2,2) );
        (*functions_err)[ratio_idx].push_back(err);

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
        canva->cd(ch+1);
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

    canva->Print("images/charges.pdf(","pdf");

    delete exp_fit;
    delete graph;
    delete canva;

}

void PlotChargesFunctions(vector <vector <double> >* fmeans, vector <vector <double> >* fstdDevs, vector <double> positions_x){
    
    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);

    const int color[3] = {kAzure-5, kAzure-5, kAzure-5};//kViolet+6};
    TString name[3] = {"Charges sum", "Charges difference", "Charges ratio"};

    TGraphErrors* graph = new TGraphErrors();
    TF1* fit = new TF1();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    for(int ch = 0; ch < 3; ++ch)
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
        graph->SetDrawOption("C AP");
        gStyle->SetEndErrorSize(8);

        // axis
        graph->SetTitle(name[ch]);
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->CenterTitle();
        graph->GetXaxis()->SetTitle("position x [mm]");
        graph->GetYaxis()->SetTitle("Number of detected #gamma");
     
        graph->Draw();

        // hyperbolic cosine + offset
        //if(ch == 0) fit = new TF1("f(x)+f(-x)", "2*[0] + [1]*(exp(x*[2])+exp(-x*[2]))", -HALF_LEN_X, HALF_LEN_X);
        fit = new TF1("f(x)+f(-x)", "pol5", -HALF_LEN_X, HALF_LEN_X);
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

        if (ch == 2) canva->Print("images/charges.pdf)","pdf");
        else canva->Print("images/charges.pdf","pdf");

        graph->Clear();
        canva->Clear();

    }

    delete fit;
    delete graph;
    delete canva;

}