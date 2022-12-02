#include "infos.h"
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

#ifndef charge_analysis_h
#define charge_analysis_h

void AddCharges(char * fileName, vector <vector <double_t> >* charge_means, vector <vector <double_t> >* charge_stdDevs){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("charges", myFile);
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

        charges_for_side[DX].push_back(chargeDX);
        charges_for_side[SX].push_back(chargeSX);
    }

    vector <double_t> MeanCharges = {{}, {}};
    vector <double_t> ErrCharges = {{}, {}};

    for (int ch = 0; ch < 2; ++ch)
    {
        MeanCharges[ch] = Mean(charges_for_side[ch]);
        ErrCharges[ch] = StdDeviation(charges_for_side[ch]);

        (*charge_means)[ch].push_back(MeanCharges[ch]);
        (*charge_stdDevs)[ch].push_back(ErrCharges[ch]);

        if (ch == 0) cout << "\n\tRight side:\t" << "Mean = " << MeanCharges[ch] << ",\t StdDev = " << ErrCharges[ch] << endl;
        else cout << "\tLeft side:\t" << "Mean = " << MeanCharges[ch] << ",\t StdDev = " << ErrCharges[ch] << endl;
    }

    cout << endl;

}

void AddChargesEPE(char * fileName, vector <vector <double_t> >* charge_means, vector <vector <double_t> >* charge_stdDevs,
    vector <vector <double_t> >* functions, vector <vector <double_t> >* functions_err){

    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("charges", myFile);
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

        charges_for_side[DX].push_back(chargeDX);
        charges_for_side[SX].push_back(chargeSX);

        charges_functions[0].push_back(chargeDX+chargeSX);
        charges_functions[1].push_back(chargeDX-chargeSX);
        charges_functions[2].push_back(chargeDX/chargeSX);
        charges_functions[3].push_back(chargeSX/chargeDX);

    }

    vector <double_t> MeanCharges = {0, 0};
    vector <double_t> ErrCharges = {0, 0};
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
    }

    for(int func_idx = 0; func_idx < 4; ++func_idx)
    {
        (*functions)[func_idx].push_back( Mean(charges_functions[func_idx]) );
        (*functions_err)[func_idx].push_back( StdDeviation(charges_functions[func_idx]) );
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
        val1 = (*means)[DX][meas];
        val2 = (*means)[SX][meas];
        err1 = (*errs)[DX][meas];
        err2 = (*errs)[SX][meas];

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

void PlotCharges(vector <vector <double_t> >* means, vector <vector <double_t> >* stdDevs, vector <double_t> positions_x, bool EPE){

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
        gStyle->SetOptFit(1110);
        gStyle->SetOptStat(2210);
        if(ch == 0) gStyle->SetStatX(0.5);
        else gStyle->SetStatX(0.92);

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
        graph->Fit(exp_fit, "Q", "0");
        exp_fit->SetLineColor(color[ch]);
        exp_fit->SetLineWidth(1);
        exp_fit->SetFillStyle(3002);
        exp_fit->SetFillColorAlpha(color[ch],0.5);
        exp_fit->Draw("SAME");

        canva->Update();
        graph->Clear();
    }

    if(EPE) canva->Print("images/3chargesEpE.pdf(","pdf");
    else canva->Print("images/charges.pdf(","pdf");

    delete exp_fit;
    delete graph;
    delete canva;

}

TF1* PlotChargesFunctions(vector <vector <double_t> >* fmeans, vector <vector <double_t> >* fstdDevs, vector <double_t> positions_x, bool EPE){
    
    vector <TF1*> fits = {nullptr, nullptr, nullptr, nullptr};

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);

    TString name[4] = {"Charges sum", "Charges difference", "Charges ratio", "Charges ratio"};
    TString yAxisName[4] = {"Sum N_{DX} + N_{SX}", "Difference N_{DX} - N_{SX}", "Ratio N_{DX} / N_{SX}", "Ratio N_{SX} / N_{DX}"};

    TGraphErrors* graph = new TGraphErrors();

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
        graph->SetMarkerColor(kAzure-5);
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

        gStyle->SetOptFit(1110);
        gStyle->SetOptStat(2210);
        gStyle->SetStatFontSize(2);
        if(ch==0){
            gStyle->SetStatX(0.6);
            gStyle->SetStatY(0.89);

        }
        else if(ch==1 || ch==2){
            gStyle->SetStatX(0.48);
            gStyle->SetStatY(0.89);
        }
        else if(ch==3){
            gStyle->SetStatX(0.89);
        }

        fits[ch] = new TF1("function", "pol5", -1.5*HALF_LEN_X, 1.5*HALF_LEN_X);
        graph->Fit(fits[ch], "Q", "0");
        fits[ch]->SetLineColor(kAzure-5);
        fits[ch]->SetLineWidth(1);
        fits[ch]->SetFillStyle(3002);
        fits[ch]->SetFillColorAlpha(kAzure-5,0.5);
        fits[ch]->Draw("SAME");


        canva->Draw();

        if(EPE) canva->Print("images/3chargesEpE.pdf","pdf");
        else canva->Print("images/charges.pdf","pdf");

        graph->Clear();
        canva->Clear();

    }

    delete graph;
    delete canva;

    return fits[1];

}


void HistoFillDeltaXRandomFileCharges(TH1F* histogram, TF1* function, TString fileName, vector <Double_t> * deltas){
        
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
        deltas->push_back(reconstructed_x - true_x);

    }

    myFile->Close();
    delete myFile;
}

void PlotHistogramDeltaXCharges(TH1F* histogram, TString output_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 4000);
    TF1* gauss_fit = new TF1();

    canva->SetGrid();

    histogram->Draw("E1 P");
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
    histogram->Fit(gauss_fit, "Q", "0");
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

void HistoFillDeltaXperFileCharges(TH1F* histogram, TF1* correlation_function, char* fileName, double true_x, vector <Double_t>* deltas){
    
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
        (*deltas).push_back(reconstructed_x - true_x);

    }

    return;
}

#endif
