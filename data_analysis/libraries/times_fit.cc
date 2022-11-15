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
#include "TStyle.h"

TF1* PlotFitResults2(std::vector <std::vector <double> >* means2, std::vector <std::vector <double> >* stdDevs2, std::vector <double> positions_x){

    TCanvas* canva = new TCanvas("canva", "canvas for plotting", 3800, 3500);
    const int color[2] = {kGreen+3, kOrange+9};

    canva->SetGrid();

    auto graph = new TGraphErrors();
    TF1* line_fit = new TF1();
    TLegend* leg = new TLegend();

    int noOfPoints = positions_x.size();
    double x[noOfPoints], y[noOfPoints], dx[noOfPoints], dy[noOfPoints];

    // graph with error bars
    for(int entry = 0; entry < noOfPoints; ++entry)
    {
        x[entry] = positions_x[entry];
        y[entry] = (*means2)[0][entry] - (*means2)[1][entry];
        dx[entry] = 0;
        dy[entry] = sqrt((*stdDevs2)[0][entry]*(*stdDevs2)[0][entry] + (*stdDevs2)[1][entry]*(*stdDevs2)[1][entry]);
    }
    
    graph = new TGraphErrors(noOfPoints, x, y, dx, dy);

    // axis
    std::string title = "Beam position (x,0) vs D_time of first detected #gamma (of the two sides)";
    TString Ttitle = title;
    graph->SetTitle(Ttitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
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
    line_fit->SetLineColor(color[1]);
    line_fit->SetLineWidth(1);
    line_fit->SetFillStyle(3002);
    line_fit->SetFillColorAlpha(color[0],0.9);
    //line_fit->DrawClone("SAME E3AL");

    // legend
    leg = new TLegend(0.4, 0.75, 0.89, 0.89);
    leg->SetHeader("Fit results:","");
    leg->AddEntry("", Form("coefficient: %.4g +/- %.4g",line_fit->GetParameter(1),line_fit->GetParError(1)),"L");
    leg->AddEntry("", Form("quote: %.4g +/- %.4g",line_fit->GetParameter(0),line_fit->GetParError(0)), "L");
    leg->Draw("SAME");

    line_fit->Draw("SAME");

    canva->Print("images/t_vs_x.pdf(","pdf");

    delete leg;
    delete graph;
    delete canva;

    return line_fit;

}

void PrintFitResults(std::vector <std::vector <Double_t> >* means, std::vector <std::vector <Double_t> >* stdDevs, char** files){

    PrintColor("Statistical results of initial times (\"single SiPM\" method):", OBOLDYELLOW);
    
    // loop over the entries (one per file)
    for(int entry = 0; entry < (*means)[0].size(); ++entry)
    {
        std::cout << "Result for the file " << files[entry+1] << std::endl;
        
        for(int channel = 0; channel < numberOfChannels; ++channel)
        {
            std::cout << "\tSiPM # [" << channel+1 << "]\t"
                << "mean: " << (*means)[channel][entry]
                << "\tstd deviation: " << (*stdDevs)[channel][entry] << std::endl;
        }

        std::cout << std::endl;
    }
}

void PrintFitResults2(std::vector <std::vector <Double_t> >* means2, std::vector <std::vector <Double_t> >* stdDevs2, char** files){

    PrintColor("Statistical results of initial times (\"sides\" method):", OBOLDYELLOW);

    double_t mean = 0;
    double_t stdDev = 0;
    
    // loop over the entries (one per file)
    for(int entry = 0; entry < (*means2)[0].size(); ++entry)
    {
        mean = (*means2)[0][entry] - (*means2)[1][entry];
        stdDev = sqrt((*stdDevs2)[0][entry]*(*stdDevs2)[0][entry] + (*stdDevs2)[1][entry]*(*stdDevs2)[1][entry]);

        std::cout << "Result for the file " << files[entry+1] << std::endl;
        
        std::cout << "\tRight side minus left side\t"
                << "mean: " << mean
                << "\tstd deviation: " << stdDev << "\n" << std::endl;

    }
}