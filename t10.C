#include <iostream>
#include <cmath>
#include <fstream>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TMinuit.h>
#include <vector>

TH1F* gHist1;
TH1F* gHist2;

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

void combinedChi2(int& npar, double* grad, double& f, double* par, int flag) {
    double c = par[0];
    double A = par[1];
    double mu = par[2];
    double sigma = par[3];
    double chi2_1 = 0.0;
    for (int i = 1; i <= gHist1->GetNbinsX(); ++i) {
        double x = gHist1->GetBinCenter(i);
        double observed = gHist1->GetBinContent(i);
        double expected = c + A * std::exp(-0.5 * std::pow((x - mu) / sigma, 2.0));  
        double error = gHist1->GetBinError(i);
        chi2_1 += ((observed - expected) / error) * ((observed - expected) / error);
    }
    double chi2_2 = 0.0;
    for (int i = 1; i <= gHist2->GetNbinsX(); ++i) {
        double observed = gHist2->GetBinContent(i);
        double error = gHist2->GetBinError(i);
        double expected = c; 
        chi2_2 += ((observed - expected) / error) * ((observed - expected) / error);
    }
    f = chi2_1 + chi2_2;
}

void fitHistogramsSimultaneously(TH1F* hist1, TH1F* hist2) {
    gHist1 = hist1;
    gHist2 = hist2;
    TMinuit minimizer(4);
    minimizer.SetFCN(combinedChi2);
    minimizer.DefineParameter(0, "c", 1.0, 0.1, 0.0, 10.0);
    minimizer.DefineParameter(1, "A", 1.0, 0.1, 0.0, 100.0);
    minimizer.DefineParameter(2, "mu", 550, 1.0, 500, 600);
    minimizer.DefineParameter(3, "sigma", 10, 0.1, 1.0, 50.0);
    minimizer.Migrad();
    
    double c, A, mu, sigma, c_err, A_err, mu_err, sigma_err;
    minimizer.GetParameter(0, c, c_err);
    minimizer.GetParameter(1, A, A_err);
    minimizer.GetParameter(2, mu, mu_err);
    minimizer.GetParameter(3, sigma, sigma_err);
    std::cout << "Fit results:" << std::endl;
    std::cout << "c = " << c << " ± " << c_err << std::endl;
    std::cout << "A = " << A << " ± " << A_err << std::endl;
    std::cout << "mu = " << mu << " ± " << mu_err << std::endl;
    std::cout << "sigma = " << sigma << " ± " << sigma_err << std::endl;

    double Nsignal = A * sigma * sqrt(2 * M_PI);
    std::cout << "Nsignal (number of events under Gaussian) = " << Nsignal << std::endl;

    TF1* func1 = new TF1("func1", "[0] + [1] * exp(-0.5 * pow((x - [2]) / [3], 2.0))", 500, 600);
    func1->SetParameters(c, A, mu, sigma);
    hist1->Fit(func1, "R");

    TF1* func2 = new TF1("func2", "[0]", 500, 600);
    func2->SetParameter(0, c);
    hist2->Fit(func2, "R");
}


void analyzeData() {
    std::vector<double> data1, data2;
    loadData("data_1.dat", data1);
    loadData("data_2.dat", data2);
    TH1F* hist1 = createHistogram(data1, "hist1");
    TH1F* hist2 = createHistogram(data2, "hist2");
    fitHistogramsSimultaneously(hist1, hist2);
    TCanvas* canvas1 = new TCanvas("canvas1", "Histogram 1", 800, 600);
    hist1->GetXaxis()->SetTitle("Units");
    hist1->GetYaxis()->SetTitle("Counts");
    hist1->Draw();
    canvas1->SaveAs("10_hist1_fit.pdf");
    TCanvas* canvas2 = new TCanvas("canvas2", "Histogram 2", 800, 600);
    hist2->GetXaxis()->SetTitle("Units");
    hist2->GetYaxis()->SetTitle("Counts");
    hist2->Draw();
    canvas2->SaveAs("10_hist2_fit.pdf");
    delete canvas1;
    delete canvas2;
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    analyzeData();
    app.Run();
    return 0;
}
/*Fit results:
c = 8.55389 ± 0.243445
A = 41.4977 ± 1.82611
mu = 548.987 ± 0.430794
sigma = 9.95201 ± 0.358879
Nsignal (number of events under Gaussian) = 1035.2
*/
