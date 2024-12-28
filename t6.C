#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TApplication.h>

const double mKS = 497.611; 
const double mK_L = 497.614; 
const double mPi = 139.570; 

const double E_CM = 1020.0; 

double generateThetaK() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    
    double cosThetaK = dis(gen);
    return acos(cosThetaK); 
}

double generatePhiK() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 2 * M_PI);
    
    return dis(gen); 
}

double generateThetaPi() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    
    double cosThetaPi = dis(gen);
    return acos(cosThetaPi); 
}

double generatePhiPi() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 2 * M_PI);
    
    return dis(gen); 
}

void generateEvent(std::vector<double>& cosThetaPi, std::vector<double>& phiPi,
                   std::vector<double>& cosThetaK, std::vector<double>& phiK,
                   std::vector<double>& length) {

    
    double thetaK = generateThetaK();
    double phiK_val = generatePhiK();
    double thetaPi = generateThetaPi();
    double phiPi_val = generatePhiPi();
    
    double pPi = E_CM / 2; 

    if (pPi < 40.0) return; 

    cosThetaPi.push_back(cos(thetaPi));
    phiPi.push_back(phiPi_val);
    cosThetaK.push_back(cos(thetaK));
    phiK.push_back(phiK_val);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(10.0, 50.0);
    length.push_back(dis(gen)); 
}

void plotDistributions(const std::vector<double>& cosThetaPi, const std::vector<double>& phiPi,
                       const std::vector<double>& cosThetaK, const std::vector<double>& phiK,
                       const std::vector<double>& length) {

    
    TH1F* h_cosThetaPi = new TH1F("h_cosThetaPi", "Distribution of cos(Theta_pi) in KS", 100, -1, 1);
    TH1F* h_phiPi = new TH1F("h_phiPi", "Distribution of Phi_pi in KS", 100, 0, 2 * M_PI);
    TH1F* h_cosThetaK = new TH1F("h_cosThetaK", "Distribution of cos(Theta_K) in lab frame", 100, -1, 1);
    TH1F* h_phiK = new TH1F("h_phiK", "Distribution of Phi_K in lab frame", 100, 0, 2 * M_PI);
    TH1F* h_length = new TH1F("h_length", "Distribution of length of flight", 100, 10, 50);

    for (size_t i = 0; i < cosThetaPi.size(); i++) {
        h_cosThetaPi->Fill(cosThetaPi[i]);
        h_phiPi->Fill(phiPi[i]);
        h_cosThetaK->Fill(cosThetaK[i]);
        h_phiK->Fill(phiK[i]);
        h_length->Fill(length[i]);
    }
    
    TCanvas* canvas = new TCanvas("canvas", "Distributions", 800, 600);
    canvas->Divide(2, 2);
    canvas->cd(1);
    h_cosThetaPi->Draw();
    canvas->cd(2);
    h_phiPi->Draw();
    canvas->cd(3);
    h_cosThetaK->Draw();
    canvas->cd(4);
    h_length->Draw(); 
    canvas->SaveAs("distributions.pdf");
}

int t6() {   
    std::vector<double> cosThetaPi, phiPi, cosThetaK, phiK, length;

    for (int i = 0; i < 10000; i++) {
        generateEvent(cosThetaPi, phiPi, cosThetaK, phiK, length);
    }

    plotDistributions(cosThetaPi, phiPi, cosThetaK, phiK, length);
    return 0;
}
