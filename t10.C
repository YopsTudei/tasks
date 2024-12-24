#include <iostream>
#include <cmath>
#include <fstream>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <vector>

void loadData(const char* fileName, std::vector<double>& data) {
    std::ifstream file(fileName);
    double value;
    while (file >> value) {
        data.push_back(value);
    }
}

TH1F* createHistogram(const std::vector<double>& data, const char* name) {
    TH1F* hist = new TH1F(name, "Data Histogram", 100, 500, 600);
    for (double val : data) {
        hist->Fill(val);
    }
    return hist;
}

double func1(double *x, double *par) {
    double c = par[0];
    double A = par[1];
    double mu = par[2];
    double sigma = par[3];
    return c + A * exp(-0.5 * pow((x[0] - mu) / sigma, 2));
}

double func2(double *x, double *par) {
    return par[0];
}

void fitHistograms(TH1F* hist1, TH1F* hist2, TF1*& f1, TF1*& f2) {
    f1 = new TF1("f1", func1, 500, 600, 4);
    f2 = new TF1("f2", func2, 500, 600, 1);
    f1->SetParameters(1.0, 1.0, 550, 10);
    f2->SetParameters(1.0);
    hist1->Fit(f1, "R");
    f2->FixParameter(0, f1->GetParameter(0));
    hist2->Fit(f2, "R");
    std::cout << "Fit parameters for data_1: c = " << f1->GetParameter(0)
              << ", A = " << f1->GetParameter(1)
              << ", mu = " << f1->GetParameter(2)
              << ", sigma = " << f1->GetParameter(3) << std::endl;
    std::cout << "Fit parameters for data_2: c = " << f2->GetParameter(0) << std::endl;
}

void calculateGaussEvents(TF1* f1) {
    double gaussEvents = f1->Integral(500, 600);
    std::cout << "Number of events under the Gaussian: " << gaussEvents << std::endl;
}

void plotAndSaveSeparate(TH1F* hist1, TH1F* hist2, TF1* f1, TF1* f2) {
    TCanvas *canvas1 = new TCanvas("canvas1", "Histogram 1 with Fit", 800, 600);
    hist1->SetLineColor(kBlue);
    hist1->Draw();
    f1->Draw("SAME");
    canvas1->SaveAs("fit_results_data1.pdf");
    TCanvas *canvas2 = new TCanvas("canvas2", "Histogram 2 with Fit", 800, 600);
    hist2->SetLineColor(kRed);
    hist2->Draw();
    f2->Draw("SAME");
    canvas2->SaveAs("fit_results_data2.pdf");
    delete canvas1;
    delete canvas2;
}

void analyzeData() {
    std::vector<double> data1, data2;
    loadData("data_1.dat", data1);
    loadData("data_2.dat", data2);
    TH1F* hist1 = createHistogram(data1, "hist1");
    TH1F* hist2 = createHistogram(data2, "hist2");
    TF1* f1 = nullptr;
    TF1* f2 = nullptr;
    fitHistograms(hist1, hist2, f1, f2);
    calculateGaussEvents(f1);
    plotAndSaveSeparate(hist1, hist2, f1, f2);
    double chi2_1 = f1->GetChisquare();
    double chi2_2 = f2->GetChisquare();
    double total_chi2 = chi2_1 + chi2_2;
    std::cout << "Chi-square for data_1 fit: " << chi2_1 << std::endl;
    std::cout << "Chi-square for data_2 fit: " << chi2_2 << std::endl;
    std::cout << "Total Chi-square: " << total_chi2 << std::endl;
}

int main(int argc, char **argv) {
    TApplication app("app", &argc, argv);
    analyzeData();
    app.Run();
    return 0;
}
