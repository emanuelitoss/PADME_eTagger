#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

#include "libraries/times_plot.cc"
#include "libraries/times_fit.cc"
#include "libraries/infos.h"
#include "libraries/times_scatter.cc"

int main(int argc, char** argv){

    char* fileName;

    // 8 arrays: one for each SiPM
    std::vector <std::vector <double> >* charges_mean
        = new std::vector <std::vector <Double_t> > {{}, {}};
    std::vector <std::vector <double> >* charges_stdDevs
        = new std::vector <std::vector <Double_t> > {{}, {}};
    
    // positions in x axis for each input file
    std::vector <double> positions_x = {-90, -70., -50., -30., 0., 30., 50., 70., 90.};
    for (auto x = positions_x.begin(); x != positions_x.end(); ++x) *x *= HALF_LEN_X/100.;
    
    if(argc == 2)
    {
        PrintColor("Error: modify .cpp file or insert more input files", OBOLDRED);
        return EXIT_FAILURE;
    }
    
    else
    {
        for(int file_counter = 1; file_counter < argc; ++file_counter)
        {
            fileName = argv[file_counter];
            addCharges(fileName, charges_mean, charges_stdDevs);
        }
    }

    PlotCharges(charges_mean, charges_stdDevs, positions_x);

    delete charges_stdDevs;
    delete charges_mean;

    return EXIT_SUCCESS;
}

void addCharges(char * fileName, std::vector <std::vector <Double_t> >* charge_means, std::vector <std::vector <Double_t> >* charge_stdDevs){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("totalCharges", myFile);
    std::vector< TTreeReaderValue<Double_t> > charges;

    std::vector < std::vector <double_t> > charges_for_side = {{}, {}};

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // initialize readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameStr = "totalCharges[" + ChannelName + "]";
      const char* dirName = dirNameStr.c_str();
      charges.push_back(TTreeReaderValue<Double_t>(reader,dirName));
    }

    double_t charge = 0;

    while (reader.Next())
    {
        for(int channel = 0; channel < numberOfChannels; channel++)
        {
            charge = *charges[channel];
            if(charge != 0)
            {
                if(channel < numberOfChannels/2) charges_for_side[0].push_back(charge);
                else charges_for_side[1].push_back(charge);
            }
        }
    }

    std::vector <double_t> MeanCharges = {{}, {}};
    std::vector <double_t> ErrCharges = {{}, {}};

    for (int ch = 0; ch < 2; ++ch)
    {
        MeanCharges[ch] = Mean(charges_for_side[ch]);
        ErrCharges[ch] = StdDeviation(charges_for_side[ch]);
    }
    
    charge_means->push_back(MeanCharges);
    charge_stdDevs->push_back(ErrCharges);

}

double_t Mean(std::vector <double_t> vec){
    
    double_t mean = 0;

    for (int i = 0; i<vec.size(); ++i){
        mean += vec[i];
    }

    mean/=(vec.size());
    return mean;
}

double_t StdDeviation(std::vector <double_t> vec){
    
    double_t err = 0;

    for (int i = 0; i<vec.size(); ++i){
        err += (vec[i])*(vec[i]);
    }

    err -= ( Mean(vec)*Mean(vec) );
    err /= vec.size();
    err = sqrt(err);

    return err;

}

void PlotCharges(std::vector <std::vector <double> >* means, std::vector <std::vector <double> >* stdDevs, std::vector <double> positions_x){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3800, 3500);
    const int color[2] = {kGreen+3, kOrange+9};

    auto graph = new TGraphErrors();
    TF1* line_fit = new TF1();
    TLegend* leg = new TLegend();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    /*********** graph with error bars DX ***********/
    for(int entry = 0; entry < noOfPoints; ++entry)
    {
        y[entry] = positions_x[entry];
        x[entry] = (*means)[0][entry];
        dy[entry] = 0;
        dx[entry] = sqrt((*stdDevs)[0][entry]);
    }
    graph = new TGraphErrors(noOfPoints, x, y, dx, dy);

    // axis
    std::string title = "Total charges vs Beam position (x,0)";
    TString Ttitle = title;
    graph->SetTitle(Ttitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    graph->GetXaxis()->SetTitle("Time [ns]");
    graph->GetYaxis()->SetTitle("position x [mm]");

    // markers
    graph->SetMarkerStyle(kCircle);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerSize(2.);

    graph->Draw("AP");
    gStyle->SetEndErrorSize(8);

    // linear fit
    line_fit = new TF1("fitting a line", "pol1", -HALF_LEN_X, HALF_LEN_X);
    graph->Fit(line_fit, "0", "0");
    line_fit->SetLineColor(color[1]);
    line_fit->SetLineWidth(1);
    line_fit->SetFillStyle(3002);
    line_fit->SetFillColorAlpha(color[0],0.5);

    /*********** graph with error bars SX ***********/
    for(int entry = 0; entry < noOfPoints; ++entry)
    {
        y[entry] = positions_x[entry];
        x[entry] = (*means)[1][entry];
        dy[entry] = 0;
        dx[entry] = sqrt((*stdDevs)[1][entry]);
    }
    graph = new TGraphErrors(noOfPoints, x, y, dx, dy);

    // axis
    std::string title = "Total charges vs Beam position (x,0)";
    TString Ttitle = title;
    graph->SetTitle(Ttitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    graph->GetXaxis()->SetTitle("Time [ns]");
    graph->GetYaxis()->SetTitle("position x [mm]");

    // markers
    graph->SetMarkerStyle(kCircle);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerSize(2.);

    graph->Draw("AP");
    gStyle->SetEndErrorSize(8);

    // linear fit
    line_fit = new TF1("fitting a line", "pol1", -HALF_LEN_X, HALF_LEN_X);
    graph->Fit(line_fit, "0", "0");
    line_fit->SetLineColor(color[1]);
    line_fit->SetLineWidth(1);
    line_fit->SetFillStyle(3002);
    line_fit->SetFillColorAlpha(color[0],0.5);

    // legend
    leg = new TLegend(0.4, 0.75, 0.89, 0.89);
    leg->SetHeader("Fit results:","");
    leg->AddEntry("", Form("coefficient: %.4g +/- %.4g",line_fit->GetParameter(1),line_fit->GetParError(1)),"L");
    leg->AddEntry("", Form("quote: %.4g +/- %.4g",line_fit->GetParameter(0),line_fit->GetParError(0)), "L");
    leg->Draw("SAME");

    line_fit->Draw("SAME");

    canva->Print("images/charges.pdf","pdf");

    delete leg;
    delete line_fit;
    delete graph;
    delete canva;

}