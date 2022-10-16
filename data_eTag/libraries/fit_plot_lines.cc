#include "infos.h"

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TObjString.h"

void PlotFitResults(std::vector <std::vector <double> >* means, std::vector <std::vector <double> >* stdDevs, std::vector <double> positions_x){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3800, 3600);
    canva->Divide(2,2);
    const int color[8] = {kGreen+3, kGreen+3, kGreen+3, kGreen+3, kOrange+9, kOrange+9, kOrange+9, kOrange+9};

    auto graph = new TGraphErrors();
    TF1* line_fit = new TF1();
    TLegend* leg = new TLegend();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    // loop over channels
    for(int channel = 0; channel < numberOfChannels; ++channel)
    {
        canva->cd((channel%4)+1);

        // graph with error bars
        for(int entry = 0; entry < noOfPoints; ++entry)
        {
            x[entry] = positions_x[entry];
            y[entry] = (*means)[channel][entry];
            dx[entry] = 0;
            dy[entry] = (*stdDevs)[channel][entry];
        }
        graph = new TGraphErrors(noOfPoints, x, y, dx, dy);

        // axis
        std::string title = "Beam position (x,0) vs time of first detected #gamma. SiPM #" + std::to_string(channel+1);
        TString Ttitle = title;
        graph->SetTitle(Ttitle);
        graph->GetXaxis()->CenterTitle();
        graph->GetXaxis()->SetTitle("position x [mm]");
        graph->GetYaxis()->SetTitle("Time [ns]");

        // markers
        graph->SetMarkerStyle(kCircle);
        graph->SetMarkerColor(kBlack);
        graph->SetMarkerSize(2.);

        graph->Draw("AP");
        gStyle->SetEndErrorSize(8);

        // linear fit
        line_fit = new TF1("fitting a line", "pol1", -HALF_LEN_X, HALF_LEN_X);
        graph->Fit(line_fit, "0", "0");
        line_fit->SetLineColor(color[channel]);
        line_fit->SetLineWidth(1);
        line_fit->SetFillStyle(3002);
        line_fit->SetFillColorAlpha(color[channel],0.5);

        // legend
        if(channel < 4) leg = new TLegend(0.4, 0.75, 0.89, 0.89);
        else leg = new TLegend(0.12, 0.75, 0.6, 0.89);
        leg->SetHeader("Fit results:","");
        leg->AddEntry("", Form("coefficient: %.4g +/- %.4g",line_fit->GetParameter(1),line_fit->GetParError(1)),"L");
        leg->AddEntry("", Form("quote: %.4g +/- %.4g",line_fit->GetParameter(0),line_fit->GetParError(0)), "L");
        leg->Draw("SAME");

        line_fit->Draw("SAME");

        if (channel == (int)numberOfChannels/2 - 1)
        {
            canva->Print("images/line_time_position.pdf(","pdf");
            canva->Clear();
            canva->Divide(2,2);
        } else if (channel == numberOfChannels-1) canva->Print("images/line_time_position.pdf)","pdf");

    }

    delete leg;
    delete line_fit;
    delete graph;
    delete canva;

}

void PrintFitResults(std::vector <std::vector <Double_t> >* means, std::vector <std::vector <Double_t> >* stdDevs, char** files){

    std::cout << OBOLDYELLOW << "Statistical results of initial times:" << ORESET << std::endl;
    
    // loop over the entries (one per file)
    for(int entry = 0; entry < (*means)[0].size(); ++entry)
    {
        std::cout << "Result of the file " << files[entry+1] << std::endl;
        
        for(int channel = 0; channel < numberOfChannels; ++channel)
        {
            std::cout << "\tSiPM # [" << channel+1 << "]\t"
                << "mean: " << (*means)[channel][entry]
                << "\tstd deviation: " << (*stdDevs)[channel][entry] << std::endl;
        }

        std::cout << std::endl;
    }
}