#include "infos.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF2.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

void HistoFillDeltaXperFileTimes(TH1F* histogram, TF1* function, TString fileName,  double true_x, vector <Double_t>* deltas, int my_option){
    
    // reading objects
    TFile* myFile = TFile::Open(fileName);
    TTreeReader reader = TTreeReader();
    reader.SetTree("arrival_times", myFile);
    std::vector < TTreeReaderValue <Double_t> > times;

    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "arrival_times_SiPM[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      times.push_back(TTreeReaderValue<Double_t>(reader,dirName));
    }

    double_t time = 0, min_timeSX = 0, min_timeDX = 0;
    double_t reconstructed_x;

    while(reader.Next())
    {

        min_timeSX = OFFSET;
        min_timeDX = OFFSET;

        for (int channel = 0; channel < numberOfChannels; channel++)
        {
            time = *times[channel];
            if(time != 0)
            {
                if(channel < numberOfChannels/2)
                {
                    if(time < min_timeDX) min_timeDX = time;
                }
                else
                {
                    if(time < min_timeSX) min_timeSX = time;
                }            
            }
        }

        reconstructed_x = function->Eval(min_timeDX-min_timeSX);
        
        if(my_option == OPTION_DELTA_POSITION)
        {
            histogram->Fill(reconstructed_x - true_x);
            (*deltas).push_back(reconstructed_x - true_x);
        }
        else if(my_option == OPTION_POSITION)
        {
            histogram->Fill(reconstructed_x);
            (*deltas).push_back(reconstructed_x);
        }

    }

    myFile->Close();
    delete myFile;

    return;
}

void HistoFillDeltaXRandomFileTimes(TH1F* histogram, TF1* function, TString fileName, vector <Double_t> * deltas){
        
    // reading objects
    TFile* myFile = TFile::Open(fileName);

    TTreeReader readerT = TTreeReader();
    readerT.SetTree("arrival_times", myFile);

    TTreeReader readerX = TTreeReader();
    readerX.SetTree("positions", myFile);
    TTreeReaderValue <Double_t> position = TTreeReaderValue <Double_t>(readerX,"X_position");

    std::vector < TTreeReaderValue <Double_t> > times;
    for (int channel = 0; channel < numberOfChannels; ++channel)
    {
      // readers
      std::string ChannelName = std::to_string(channel+1);
      std::string dirNameS = "arrival_times_SiPM[" + ChannelName + "]";
      const char* dirName = dirNameS.c_str();
      times.push_back(TTreeReaderValue<Double_t>(readerT,dirName));
    }

    double_t time = 0, min_timeSX = 0, min_timeDX = 0;
    double_t reconstructed_x, true_x;

    while(readerT.Next() && readerX.Next())
    {

        min_timeSX = OFFSET;
        min_timeDX = OFFSET;

        true_x = *position;

        for (int channel = 0; channel < numberOfChannels; channel++)
        {
            time = *times[channel];
            if(time != 0)
            {
                if(channel < numberOfChannels/2)
                {
                    if(time < min_timeDX) min_timeDX = time;
                }
                else
                {
                    if(time < min_timeSX) min_timeSX = time;
                }            
            }
        }

        reconstructed_x = function->Eval(min_timeDX-min_timeSX);
        histogram->Fill(reconstructed_x - true_x);
        deltas->push_back(reconstructed_x - true_x);

    }

    myFile->Close();
    delete myFile;
}

void PlotHistogramDeltaXTimes(TH1F* histo_deltaX, TString output_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 4000);

    canva->SetGrid();

    //gStyle->Reset();
    gStyle->SetEndErrorSize(8);
    gStyle->SetOptFit(1110);
    gStyle->SetOptStat(2210);
    gStyle->SetStatFont(42);
    gStyle->SetStatFontSize(0.02);
    gStyle->SetLineWidth(1);
    gStyle->SetStatX(0.89);
    gStyle->SetStatY(0.89);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.7);

    histo_deltaX->Draw("E1 P");
    histo_deltaX->GetXaxis()->SetTitle("Length [mm]");
    histo_deltaX->GetYaxis()->SetTitle("Number of events");
    histo_deltaX->SetLineColor(kBlack);
    histo_deltaX->SetLineWidth((Width_t)1.5);
    
    TF1 gauss_fit = TF1("fitting a gaussian", "gaus", -HALF_LEN_X, HALF_LEN_X);
    histo_deltaX->Fit(&gauss_fit, "L Q", "0");
    gauss_fit.SetLineColor(kAzure-5);
    gauss_fit.SetLineWidth(1);
    gauss_fit.SetFillStyle(3002);
    gauss_fit.SetFillColorAlpha(kAzure-5,0.5);
    gauss_fit.SetObjectStat(true);
    gauss_fit.Draw("C SAME");
    
    canva->Print(TString(output_name),"pdf");
    canva->Clear();
  
    delete canva;

}


void PlotHistogramDeltaXTimesSpecial(TH1F* histo_deltaX, TF1* gauss_fit, TString output_name, vector <Double_t> * sigmas, vector <Double_t> * err_sigmas){

    TCanvas canva = TCanvas("canva", "canvas for plotting", 4000, 4000);

    canva.SetGrid();

    histo_deltaX->Draw("E1 P");
    histo_deltaX->GetXaxis()->SetTitle("Length [mm]");
    histo_deltaX->GetYaxis()->SetTitle("Number of events");
    histo_deltaX->SetLineColor(kBlack);
    histo_deltaX->SetLineWidth((Width_t)1.5);
    
    gStyle->SetEndErrorSize(8);
    gStyle->SetOptFit(1110);
    gStyle->SetOptStat(2210);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetStatFont(42);
    gStyle->SetLineWidth(1);
    
    histo_deltaX->Fit(gauss_fit, "L Q", "0");
    gauss_fit->SetLineColor(kAzure-5);
    gauss_fit->SetLineWidth(1);
    gauss_fit->SetFillStyle(3002);
    gauss_fit->SetFillColorAlpha(kAzure-5,0.5);
    gauss_fit->SetObjectStat(true);
    gauss_fit->Draw("C SAME");

    sigmas->push_back(gauss_fit->GetParameter(2));
    err_sigmas->push_back(gauss_fit->GetParError(2));
    
    canva.Print(TString(output_name),"pdf");
    canva.Clear();

}


void PlotPositionsQT(vector <Double_t >* deltasT,vector <Double_t >* deltasQ, int openclosefile, int option_diff_dos){

    // check that I have couples of points.
    int len = (*deltasT).size();
    bool check_dimensions = (len==(*deltasQ).size());
    if(!check_dimensions){
        PrintColor("ERROR: Not coupled points", OBOLDRED);
        return;
    }

    TH2F histogram = TH2F("hist","x_{reco} (time analysis) vs x_{reco} (charge analysis)",70,-HALF_LEN_X,HALF_LEN_X,70,-HALF_LEN_X,HALF_LEN_X);
    
    if(option_diff_dos == OPTION_Q_DIFFERENCE){
        histogram.SetName("hist_diff");
        histogram.SetTitle("x_{reco} from #Deltat vs x_{reco} from #DeltaQ");
    }
    else if(option_diff_dos == OPTION_Q_DOS){
        histogram.SetName("hist_dos");
        histogram.SetTitle("x_{reco} from #Deltat vs x_{reco} from f(Q_{DX},Q_{SX})");
    }

    // storing.
    Double_t deltaxT[len], deltaxQ[len];
    for(int i = 0; i<len; ++i)  histogram.Fill( (*deltasT)[i], (*deltasQ)[i] );

    // plot.
    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 2000, 1500);

    /********** FANCY HISTOGRAM LEVEL CURVES **********/
    histogram.SetXTitle("x_{reco} (times analysis) [mm]");
    histogram.SetYTitle("x_{reco} (charges analysis) [mm]");
    histogram.SetMarkerStyle(kDot);
    histogram.SetMarkerSize(3.);
    histogram.SetMarkerColor(kBlack);
    gStyle->SetStatH(0.22);
    gStyle->SetStatW(0.10);
    gStyle->SetStatX(0.89);
    gStyle->SetStatY(0.53);
    gStyle->SetStatFontSize(0.02);
    histogram.Draw("");

    cout << OBOLDWHITE << "\nCovariance matrix:\n" << ORESET << endl;
    cout << OWHITE << "\tCovariance = " << histogram.GetCovariance() << "\n"
        << "\tVariance x = " << histogram.GetStdDev(1) << "\n"
        << "\tVariance y = " << histogram.GetStdDev(2) << "\n"
        << "\tCorrelation function = " << histogram.GetCorrelationFactor() << ORESET << endl;

    TF1 fitting_function = TF1("Line","pol1",-HALF_LEN_X, HALF_LEN_X);
    histogram.Fit(&fitting_function, "", "");
    fitting_function.SetLineColor(kRed+1);
    fitting_function.SetLineWidth(1);
    fitting_function.Draw("SAME E1");

    if(openclosefile == OPEN_OUTPUT || openclosefile == SINGLE_OUTPUT)  canva->Print("images/5delta_positions_correlation.pdf(","pdf");
    else canva->Print("images/5delta_positions_correlation.pdf","pdf");
    
    canva->Clear();

    /********** FANCY HISTOGRAM LEVEL CURVES **********/
    histogram.SetXTitle("x_{reco} (times analysis) [mm]");
    histogram.SetYTitle("x_{reco} (charges analysis) [mm]");
    gStyle->SetStatH(0.35);
    gStyle->SetStatW(0.1);
    gStyle->SetStatX(0.25);
    gStyle->SetStatY(0.89);
    histogram.Draw("CONT4Z");

    canva->Print("images/5delta_positions_correlation.pdf","pdf");
    canva->Clear();

    /********** FANCY HISTOGRAM 3D **********/
    histogram.SetXTitle("x_{reco} (times analysis) [mm]");
    histogram.SetYTitle("x_{reco} (charges analysis) [mm]");
    gStyle->SetStatX(0.99);
    gStyle->SetStatY(0.85);
    histogram.Draw("LEGO2");
    if(openclosefile == CLOSE_OUTPUT || openclosefile == SINGLE_OUTPUT) canva->Print("images/5delta_positions_correlation.pdf)","pdf");
    else canva->Print("images/5delta_positions_correlation.pdf","pdf");
    canva->Clear();

    delete canva;


}

void PlotSigmaCorrections(vector <Double_t> * sigmaT_vec, vector <Double_t> * sigmaQ_vec, vector <Double_t> * err_sigmaT_vec, vector <Double_t> * err_sigmaQ_vec, vector <double_t> positions_vec, TString filename){

    Int_t size = positions_vec.size();
    double sigma_t[size], sigma_q[size], err_sigma_t[size], err_sigma_q[size], positions[size], zeros[size];

    for(int i=0; i<size; ++i){
        sigma_t[i] = (*sigmaT_vec)[i];
        sigma_q[i] = (*sigmaQ_vec)[i];
        err_sigma_t[i] = (*err_sigmaT_vec)[i];
        err_sigma_q[i] = (*err_sigmaQ_vec)[i];
        positions[i] = positions_vec[i];
        zeros[i] = 0;
    }

    TCanvas canva = new TCanvas("canva", "canvas for plotting", 2000, 1500);
    canva.SetGrid();

    TMultiGraph* mg = new TMultiGraph();
    mg->SetTitle("Error #sigma_{x}(x_{reco}) on reconstruction of x in function of x_{reco}");

    TGraphErrors graphT = TGraphErrors(size, positions, sigma_t, zeros, err_sigma_t);
    TGraphErrors graphQ = TGraphErrors(size, positions, sigma_q, zeros, err_sigma_q);
    
    TF1 fitT = TF1();
    TF1 fitQ = TF1();

    gStyle->SetEndErrorSize(4);
    gStyle->SetLegendTextSize(0.02);

    graphT.SetMarkerStyle(kFullTriangleUp);
    graphT.SetMarkerColor(kBlack);
    graphT.SetMarkerSize(0.6);
    graphT.SetLineColor(kBlue-3);
    graphT.SetTitle("Sigma_x from times analysis");
    mg->Add(&graphT);
    
    graphQ.SetMarkerStyle(kFullSquare);
    graphQ.SetMarkerColor(kBlack);
    graphQ.SetMarkerSize(0.4);
    graphQ.SetLineColor(kRed-3);
    graphQ.SetTitle("Sigma_x from charges analysis");
    mg->Add(&graphQ);

    mg->GetXaxis()->SetTitle("position [mm]");
    mg->GetYaxis()->SetTitle("uncertainty #sigma(x_{reco}) [mm]");
    mg->Draw("ALP");

    canva.BuildLegend(0.55, 0.35, 0.89, 0.50);

    canva.Print(filename,"pdf");

    delete mg;

}