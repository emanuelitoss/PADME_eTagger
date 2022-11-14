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

void PlotCharges(vector <vector <double_t> >* means, vector <vector <double_t> >* stdDevs, vector <double_t> positions_x, bool EVE){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);
    canva->Divide(2,1);
    const int color[2] = {kGreen+3, kOrange+9};

    TGraphErrors* graph = new TGraphErrors();
    TF1* exp_fit = new TF1();
    TLegend* leg;

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

        // legend
        if(ch == 0) leg = new TLegend(0.11, 0.78, 0.70, 0.89);
        else leg = new TLegend(0.3, 0.78, 0.89, 0.89);
        leg->SetHeader("Fit results:","");
        leg->AddEntry("-", Form("quote: %.4g +/- %.4g",exp_fit->GetParameter(0), exp_fit->GetParError(0)),"L");
        leg->AddEntry("-", Form("coefficient: %.4g +/- %.4g",exp_fit->GetParameter(1), exp_fit->GetParError(1)), "L");
        leg->AddEntry("-", Form("decayconstant: %.4g +/- %.4g",exp_fit->GetParameter(2), exp_fit->GetParError(2)), "L");
        leg->Draw("SAME");

        canva->Update();
        graph->Clear();
    }

    if(EVE) canva->Print("images/chargesEpE.pdf(","pdf");
    else canva->Print("images/charges.pdf(","pdf");

    delete leg;
    delete exp_fit;
    delete graph;
    delete canva;

}

void PlotChargesFunctions(vector <vector <double_t> >* fmeans, vector <vector <double_t> >* fstdDevs, vector <double_t> positions_x, bool EVE){
    
    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3500);

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

        // hyperbolic cosine + offset
        //if(ch == 0) fit = new TF1("f(x)+f(-x)", "2*[0] + [1]*(exp(x*[2])+exp(-x*[2]))", -HALF_LEN_X, HALF_LEN_X);
        fit = new TF1("f(x)+f(-x)", "pol6", -HALF_LEN_X, HALF_LEN_X);
        // hyperbolic cosine + offset
        //else if (ch==1) fit = new TF1("f(x)-f(-x)", "[0]*(exp(x*[1])-exp(-x*[1]))", -HALF_LEN_X, HALF_LEN_X);
        // hyperbolic tangent
        //else fit = new TF1("f(x)/f(-x)", "(1+[0]*exp(x*[1]))/(1+[0]*exp(-x*[1]))", -HALF_LEN_X, HALF_LEN_X);
        graph->Fit(fit, "0", "0");
        fit->SetLineColor(kAzure-5);
        fit->SetLineWidth(1);
        fit->SetFillStyle(3002);
        fit->SetFillColorAlpha(kAzure-5,0.5);
        fit->Draw("SAME");

        canva->Draw();

        if(EVE)
        {
            if (ch == 3) canva->Print("images/chargesEpE.pdf)","pdf");
            else canva->Print("images/chargesEpE.pdf","pdf");
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

#endif
