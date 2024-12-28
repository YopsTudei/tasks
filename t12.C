#include <iostream>
#include <fstream>
#include <vector>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>

std::vector<double> loadData(const char* filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        return {};
    }
    std::vector<double> data;
    double value;
    while (file >> value) {
        data.push_back(value);
    }
    std::cout << "Loaded " << data.size() << " values from " << filename << std::endl;
    return data;
}

TH1D* createHistogram(const std::vector<double>& data, const char* name, int bins, double min, double max) {
    TH1D* hist = new TH1D(name, "Spectrum", bins, min, max);
    for (double value : data) {
        if (value >= min && value <= max) {
            hist->Fill(value);
        } else {
            std::cerr << "Warning: Value " << value << " is outside the histogram range [" << min << ", " << max << "]" << std::endl;
        }
    }
    return hist;
}

double normalizedBreitWigner(double* x, double* par) {
    // par[1] = M (масса резонанса)
    // par[2] = G (ширина резонанса)
    double M = par[1];
    double G = par[2];
    double y = sqrt(M * M * (M * M + 1)); 
    double K = (2 * sqrt(2) * M * G * y) / (TMath::Pi() * sqrt(M * M + y)); 
    double denominator = (x[0] * x[0] - M * M) * (x[0] * x[0] - M * M) + M * M * G * G;
    return K / denominator;
}

/*
TF1* createFitFunction() {
    // Модель Брайта-Вигнера + Пуассон
    // f(x) = A / ((x - E)^2 + (Gamma^2 / 4)) + B * exp(-x / Gamma)
    TF1* func = new TF1("fitFunc", "[0] / ((x - [1])^2 + ([2]^2) / 4) + (1000-[0])/[3] * exp(-x / [3])", 0, 10);
    return func; 
}
*/

TF1* createFitFunction() {
    TF1* func = new TF1("fitFunc", [](double* x, double* par) {
        double bw = par[0] * normalizedBreitWigner(x, par);

        double expPart = ((1000 - par[0]) / par[3]) * exp(-x[0] / par[3]);

        return bw + expPart;
    }, 0, 10, 4);
    return func;
}

void fitHistogram(TH1D* hist, TF1* func) {
    func->SetParameters(25, 5, 3, 0.005); 
    hist->Fit(func, "R");
}

void plotResults(TH1D* hist) {
    TCanvas* canvas = new TCanvas("canvas", "Fit Results", 800, 600);
    hist->GetXaxis()->SetTitle("X");
    hist->GetYaxis()->SetTitle("Counts");
    hist->Draw();
    canvas->SaveAs("12_fit_results.pdf");
}
void t12() {
    std::vector<double> data = loadData("task10Nov.dat");
    if (data.empty()) {
        std::cerr << "Error: No data loaded." << std::endl;
        return;
    }
    TH1D* hist = createHistogram(data, "h", 100, 0, 10);
    TF1* fitFunc = createFitFunction();
    fitHistogram(hist, fitFunc);
    plotResults(hist);
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    std::vector<double> data = loadData("task10Nov.dat");
    TH1D* hist = createHistogram(data, "h", 100, 0, 10);
    TF1* fitFunc = createFitFunction();
    fitHistogram(hist, fitFunc);
    app.Run();
    return 0;
}
