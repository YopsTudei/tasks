#include <TLorentzRotation.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <iostream>
#include <cmath>

constexpr Double_t mKS = 497.611; // MeV
constexpr Double_t mPi = 139.57039; // MeV
constexpr Double_t E_CM = 1020.0; // MeV
constexpr Double_t MIN_ENERGY_FOR_REGISTER = 40.0;
constexpr Double_t C_TAU = 2.6844; // cm
constexpr Double_t DETECTOR_RADIUS = 30.0; // cm
constexpr Double_t DETECTOR_LENGTH = 50.0; // cm

void task6() {
    TRandom3 rnd(0);

    TH1D* hKsTheta = new TH1D("hKsTheta", "Ks Theta Distribution;#theta_{Ks} [rad];Events", 100, 0, M_PI);
    TH1D* hKsPhi = new TH1D("hKsPhi", "Ks Phi Distribution;#phi_{Ks} [rad];Events", 100, 0, 2 * M_PI);
    TH1D* hLength = new TH1D("hLength", "Ks Flight Length;Length [cm];Events", 100, 0, 50);

    for (int i = 0; i < 10000; ++i) {
        Double_t KsTheta = rnd.Uniform(0.0, M_PI);
        Double_t KsPhi = rnd.Uniform(0.0, 2 * M_PI);
        Double_t KsEnergy = (E_CM - mKS) * rnd.Rndm() + mKS;
        Double_t KsMomentum = sqrt(KsEnergy * KsEnergy - mKS * mKS);

        TLorentzVector Ks(0, 0, KsMomentum, KsEnergy);

        Double_t length = -log(rnd.Rndm()) * C_TAU * KsMomentum / KsEnergy;
        if (length > DETECTOR_LENGTH || length < 0) continue;

        hKsTheta->Fill(KsTheta);
        hKsPhi->Fill(KsPhi);
        hLength->Fill(length);
    }

    TCanvas* c = new TCanvas("c", "Distributions", 1200, 800);
    c->Divide(2, 2);

    c->cd(1); hKsTheta->Draw();
    c->cd(2); hKsPhi->Draw();
    c->cd(3); hLength->Draw();
}

