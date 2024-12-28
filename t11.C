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

double funct1(double* x, double* par) {
    double gaussian = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2.0));
    double background = par[3];
    return gaussian + background;
}

double funct2(double* x, double* par) {
    double background = par[3];
    return background;
}

void combinedLogLikelihood(int& npar, double* grad, double& f, double* par, int flag) {
    double chisq = 0.0;

    for (int i = 1; i <= gHist1->GetNbinsX(); ++i) {
        double x = gHist1->GetBinCenter(i);
        double observed = gHist1->GetBinContent(i);
        double model = funct1(&x, par);
        if (model > 0 && observed >= 0) {
            double delta = TMath::Poisson(observed, model);
            chisq += -2.0 * log(delta);
        }
    }

    for (int i = 1; i <= gHist2->GetNbinsX(); ++i) {
        double x = gHist2->GetBinCenter(i);
        double observed = gHist2->GetBinContent(i);
        double model = funct2(&x, par);
        if (model > 0 && observed >= 0) {
            double delta = TMath::Poisson(observed, model);
            chisq += -2.0 * log(delta);
        }
    }

    f = chisq;
}

void fitHistogramsSimultaneously(TH1F* hist1, TH1F* hist2) {
    gHist1 = hist1;
    gHist2 = hist2;
    TMinuit minimizer(4);
    minimizer.SetFCN(combinedLogLikelihood);
    minimizer.DefineParameter(0, "Amplitude", 1.0, 0.1, 0.0, 100.0);
    minimizer.DefineParameter(1, "Mean", 550.0, 1.0, 500.0, 600.0);
    minimizer.DefineParameter(2, "Sigma", 10.0, 0.1, 1.0, 50.0);
    minimizer.DefineParameter(3, "Background", 1.0, 0.1, 0.0, 10.0);
    minimizer.Migrad();

    double amplitude, mean, sigma, background;
    double amplitude_err, mean_err, sigma_err, background_err;
    minimizer.GetParameter(0, amplitude, amplitude_err);
    minimizer.GetParameter(1, mean, mean_err);
    minimizer.GetParameter(2, sigma, sigma_err);
    minimizer.GetParameter(3, background, background_err);

    std::cout << "Fit results:" << std::endl;
    std::cout << "Amplitude = " << amplitude << " ± " << amplitude_err << std::endl;
    std::cout << "Mean = " << mean << " ± " << mean_err << std::endl;
    std::cout << "Sigma = " << sigma << " ± " << sigma_err << std::endl;
    std::cout << "Background = " << background << " ± " << background_err << std::endl;

    double Nsignal = amplitude * sigma * sqrt(2 * M_PI);
    std::cout << "Nsignal (number of events under Gaussian) = " << Nsignal << std::endl;

    TF1* func1 = new TF1("func1", "[0] * exp(-0.5 * pow((x - [1]) / [2], 2.0)) + [3]", 500, 600);
    func1->SetParameters(amplitude, mean, sigma, background);
    hist1->Fit(func1, "R");

    TF1* func2 = new TF1("func2", "[0]", 500, 600); 
    func2->SetParameter(0, background); 
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
    canvas1->SaveAs("11_hist1_fit.pdf");

    TCanvas* canvas2 = new TCanvas("canvas2", "Histogram 2", 800, 600);
    hist2->GetXaxis()->SetTitle("Units");
    hist2->GetYaxis()->SetTitle("Counts");
    hist2->Draw();
    canvas2->SaveAs("11_hist2_fit.pdf");

    delete canvas1;
    delete canvas2;
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    analyzeData();
    app.Run();
    return 0;
}
