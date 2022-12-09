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

    myFile->Close();
    delete myFile;

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
    vector < vector <double_t> > charges_functions = {{}, {}, {}, {}, {}};
    vector < vector <double_t> > charges_functions_err = {{}, {}, {}, {}, {}};

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
        charges_functions[4].push_back((chargeDX-chargeSX)/(chargeDX+chargeSX));

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

    }

    for(int func_idx = 0; func_idx < 5; ++func_idx)
    {
        (*functions)[func_idx].push_back( Mean(charges_functions[func_idx]) );
        (*functions_err)[func_idx].push_back( StdDeviation(charges_functions[func_idx]) );
    }

    cout << endl;

    myFile->Close();
    delete myFile;

}

void AddFunctionsOfCharges(vector <vector <double_t> >* means, vector <vector <double_t> >* errs, vector <vector <double_t> >* functions, vector <vector <double_t> >* functions_err){

    Int_t numMeasures = (*means)[0].size();
    const int sum_idx = 0, diff_idx = 1, ratio_DS_idx = 2, ratio_SD_idx = 3, diffonsum = 4;
    double_t err = 0., err_app;
    double_t val1, val2, err1, err2, app;

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

        // charge difference on sum
        (*functions)[diffonsum].push_back((val1 - val2)/(val1 + val2));
        app = val1/val2;
        err_app = app*sqrt( pow(err1/val1,2) + pow(err2/val2,2) );
        err = (2/pow((val1/val2 +1),2))*err_app;
        (*functions_err)[diffonsum].push_back(err);

    }

    return;
    
}

void PlotCharges(vector <vector <double_t> >* means, vector <vector <double_t> >* stdDevs, vector <double_t> positions_x, TString output_file_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3200);
    canva->Divide(2,1);
    canva->SetGrid();
    const int color[2] = {kGreen+3, kOrange+9};

    TF1 exp_fitDX = TF1("exponential+offset", "[0] + [1]*exp(x*[2])", -HALF_LEN_X, HALF_LEN_X);
    TF1 exp_fitSX = TF1("exponential+offset", "[0] + [1]*exp(x*[2])", -HALF_LEN_X, HALF_LEN_X);
    vector <TF1> fits = {exp_fitDX, exp_fitSX};

    TGraphErrors graph1 = TGraphErrors();
    TGraphErrors graph2 = TGraphErrors();
    vector <TGraphErrors> graphs = {graph1, graph2};

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

        graphs[ch] = TGraphErrors(noOfPoints, x, y, dx, dy);
        if(ch == 0) graphs[ch].SetTitle("DX SiPMs");
        else graphs[ch].SetTitle("SX SiPMs");

        // markers
        graphs[ch].SetMarkerStyle(kCircle);
        graphs[ch].SetMarkerColor(color[ch]);
        graphs[ch].SetMarkerSize(2.);
        graphs[ch].SetDrawOption("AP");

        gStyle->SetEndErrorSize(8);
        gStyle->SetOptFit(1110);
        gStyle->SetOptStat(2210);
        gStyle->SetStatFontSize(45);
        if(ch == 0) gStyle->SetStatX(0.5);
        else gStyle->SetStatX(0.92);

        // axis
        string title = "Total charges vs Beam position x";
        TString Ttitle = title;
        graphs[ch].SetTitle(Ttitle);
        graphs[ch].GetXaxis()->CenterTitle();
        graphs[ch].GetYaxis()->CenterTitle();
        graphs[ch].GetXaxis()->SetTitle("position x [mm]");
        graphs[ch].GetYaxis()->SetTitle("Number of detected #gamma");

        graphs[ch].Draw();

        // exponential fit
        graphs[ch].Fit(&(fits[ch]), "Q", "0");
        fits[ch].SetLineColor(color[ch]);
        fits[ch].SetLineWidth(1);
        fits[ch].SetFillStyle(3002);
        fits[ch].SetFillColorAlpha(color[ch],0.5);
        fits[ch].Draw("SAME");
        graphs[ch].PaintStats(&fits[ch]);

        canva->Update();
    }

    canva->Print(output_file_name,"pdf");

    delete canva;

}

void PlotChargesFunctions(vector <vector <double_t> >* fmeans, vector <vector <double_t> >* fstdDevs, vector <double_t> positions_x, bool EPE, int ratioQ_or_differenceQ, TF1* returning_function){
    
    vector <TF1*> fits = {};

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 3300);

    TString name[5] = {"Charges sum", "Charges difference", "Charges ratio", "Charges ratio", "Charges sFunction"};
    TString yAxisName[5] = {"Sum N_{DX} + N_{SX} of detected #gamma", "Difference N_{DX} - N_{SX} of detected #gamma", "Ratio N_{DX} / N_{SX} of detected #gamma", "Ratio N_{SX} / N_{DX} of detected #gamma", "function #frac{N_{DX}-N_{SX}}{N_{DX}+N_{SX}}"};

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    for(int ch = 0; ch < 5; ++ch)
    {
        canva->SetGrid();

        for(int entry = 0; entry < noOfPoints; ++entry)
        {
            x[entry] = positions_x[entry];
            y[entry] = (*fmeans)[ch][entry];
            dx[entry] = 0;
            dy[entry] = (*fstdDevs)[ch][entry];
        }

        TGraphErrors graph = TGraphErrors(noOfPoints, x, y, dx, dy);
        graph.SetTitle(name[ch]);
        
        // markers
        graph.SetMarkerStyle(20);
        graph.SetMarkerColor(kAzure-5);
        graph.SetMarkerSize(4.);
        graph.SetDrawOption("AP");
        gStyle->SetEndErrorSize(8);

        // axis
        graph.SetTitle(name[ch]);
        graph.GetXaxis()->CenterTitle();
        graph.GetYaxis()->CenterTitle();
        graph.GetXaxis()->SetTitle("position x [mm]");
        graph.GetYaxis()->SetTitle(yAxisName[ch]);
     
        graph.Draw();

        gStyle->SetOptFit(1110);
        gStyle->SetOptStat(2210);
        gStyle->SetStatFontSize(65);

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
        else if(ch==4){
            gStyle->SetStatX(0.48);
            gStyle->SetStatY(0.89);
        }

        fits.push_back(new TF1("function", "pol5", -2*HALF_LEN_X, 2*HALF_LEN_X));
        //graph.Fit("pol5", "Q", "0", -1.5*HALF_LEN_X, 1.5*HALF_LEN_X);
        graph.Fit(fits[ch], "Q", "0");
        fits[ch]->SetLineColor(kAzure-5);
        fits[ch]->SetLineWidth(1);
        fits[ch]->SetFillStyle(3002);
        fits[ch]->SetFillColorAlpha(kAzure-5,0.5);
        fits[ch]->Draw("SAME");
        
        canva->Draw();

        if(ratioQ_or_differenceQ == OPTION_Q_DIFFERENCE){
            if(EPE) canva->Print("images/3chargesEpE.pdf","pdf");
            else canva->Print("images/charges.pdf","pdf");
        }

        if( (ratioQ_or_differenceQ == OPTION_Q_DIFFERENCE) && ch== 1)    graph.Fit(returning_function, "Q", "0");
        else if( (ratioQ_or_differenceQ == OPTION_Q_DOS) && ch== 4)  graph.Fit(returning_function, "Q", "0");
        
        graph.Clear();
        canva->Clear();

    }

    for(int ch = 0; ch < 5; ++ch) delete fits[ch];
    delete canva;
    
    return;
}


void HistoFillDeltaXRandomFileCharges(TH1F* histogram, TF1* function, TString fileName, vector <Double_t> * deltas, int my_option){
        
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

        if(my_option == OPTION_Q_DIFFERENCE) reconstructed_x = function->GetX(chargeDX-chargeSX);
        else if(my_option == OPTION_Q_DOS) reconstructed_x = function->GetX((chargeDX-chargeSX)/(chargeDX+chargeSX));
        histogram->Fill(reconstructed_x - true_x);
        deltas->push_back(reconstructed_x - true_x);

    }

    myFile->Close();
    delete myFile;
}
 
void PlotHistogramDeltaXCharges(TH1F* histogram, TString output_name){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 4000, 4000);

    canva->SetGrid();

    //gStyle->Reset();
    gStyle->SetEndErrorSize(8);
    gStyle->SetOptFit(1110);
    gStyle->SetOptStat(2210);
    gStyle->SetStatX(0.89);
    gStyle->SetStatY(0.89);
    gStyle->SetStatW(0.1);
    gStyle->SetStatH(0.7);
    gStyle->SetStatFont(42);
    gStyle->SetStatFontSize(0.02);
    gStyle->SetLineWidth(1);
    
    histogram->Draw("E1 P");
    histogram->GetXaxis()->SetTitle("Length [mm]");
    histogram->GetYaxis()->SetTitle("Number of events");
    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth((Width_t)1.5);
    
    TF1 gauss_fit = TF1("fitting a gaussian", "gaus", -HALF_LEN_X, HALF_LEN_X);
    histogram->Fit(&gauss_fit, "Q", "0");
    gauss_fit.SetLineColor(kAzure-5);
    gauss_fit.SetLineWidth(1);
    gauss_fit.SetFillStyle(3002);
    gauss_fit.SetFillColorAlpha(kAzure-5,0.5);
    gauss_fit.Draw("C SAME");
    
    canva->Print(TString(output_name),"pdf");
    canva->Clear();
  
    delete canva;

}

void PlotHistogramDeltaXChargesSpecial(TH1F* histogram, TF1* gauss_fit, TString output_name, vector <Double_t> * sigmas, vector <Double_t> * err_sigmas){

    TCanvas canva = TCanvas("canva", "canvas for plotting", 4000, 4000);

    canva.SetGrid();

    histogram->Draw("E1 P");
    histogram->GetXaxis()->SetTitle("Length [mm]");
    histogram->GetYaxis()->SetTitle("Number of events");
    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth((Width_t)1.5);

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

    histogram->Fit(gauss_fit, "Q", "0");
    // histogram->Fit("gaus", "Q", "0", -HALF_LEN_X, HALF_LEN_X);
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

void HistoFillDeltaXperFileCharges(TH1F* histogram, TF1* correlation_function, char* fileName, double true_x, vector <Double_t>* deltas, int option_Qdiff_or_Qratio, int my_option){
    
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

        if(option_Qdiff_or_Qratio == OPTION_Q_DIFFERENCE) reconstructed_x = correlation_function->GetX(chargeDX-chargeSX);        
        if(option_Qdiff_or_Qratio == OPTION_Q_DOS) reconstructed_x = correlation_function->GetX((chargeDX-chargeSX)/(chargeDX+chargeSX));
        
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
