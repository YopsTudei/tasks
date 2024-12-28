#include <ROOT/TSeq.hxx>
#include <TF1.h>
#include <TLorentzRotation.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>

namespace T6 {
    constexpr Double_t FULL_ENERGY = 1020;
    constexpr Double_t MASS_OF_PI = 139.57039;
    constexpr Double_t MASS_OF_K = 497.611;

    constexpr Double_t RADIUS_OF_DETECTOR = 30;
    constexpr Double_t LEHGTH_OF_DETECTOR = 50;
    const Double_t MAX_LENGTH = TMath::Sqrt(RADIUS_OF_DETECTOR * RADIUS_OF_DETECTOR + LEHGTH_OF_DETECTOR * LEHGTH_OF_DETECTOR / 4);

    constexpr Double_t MIN_ENERGY_FOR_REGISTER_OF_PI = 40;
    constexpr Double_t TWOPI = TMath::Pi();
    constexpr Double_t AVERAGE_MILEAGE = 6;
    constexpr Double_t C_TAU = 2.6844;
}

void task6(UInt_t N = 100000) {
    using namespace TMath;

    TRandom3 rnd(0);

    TH1D *histKsLength = new TH1D("histKsLength", "Ks Flight Length", 100, 0, T6::MAX_LENGTH);
    TH1D *histKsTheta = new TH1D("histKsTheta", "Ks Theta Distribution", 100, 0, Pi() / 2);
    TH1D *histKsPhi = new TH1D("histKsPhi", "Ks Phi Distribution", 100, 0, Pi());
    TH1D *histBTW = new TH1D("histBTW", "Angle Between Vectors", 100, -Pi(), Pi());

    TF1 *ThetaDistribution = new TF1("ThetaDistribution", "sin(x)*sin(x)", 0, TMath::Pi());

    for (auto i : ROOT::TSeqI(N)) {
        Double_t KsEnergy = (T6::FULL_ENERGY - T6::MASS_OF_K) * rnd.Rndm() + T6::MASS_OF_K;
        TLorentzVector Ks(0., 0., 0., T6::MASS_OF_K);

        TLorentzRotation T;

        Double_t KsMomentum = Sqrt(KsEnergy * KsEnergy - T6::MASS_OF_K * T6::MASS_OF_K);
        Double_t Phi = 2 * Pi() * rnd.Rndm();
        TF1 *thetaDist = new TF1("thetaDist", "1 - x^2", -1, 1);
        Double_t CosTheta = thetaDist->GetRandom(); 
        Double_t SinTheta = Sqrt(1 - CosTheta * CosTheta);

        Double_t Length = -Log(rnd.Rndm()) * T6::AVERAGE_MILEAGE;
        histKsLength->Fill(Length);

        Double_t bx = KsMomentum * SinTheta * Cos(Phi) / KsEnergy;
        Double_t by = KsMomentum * SinTheta * Sin(Phi) / KsEnergy;
        Double_t bz = KsMomentum * CosTheta / KsEnergy;

        T.Boost(bx, by, bz);

        Double_t PiEnergy = (T6::MASS_OF_K * T6::MASS_OF_K + T6::MASS_OF_PI * T6::MASS_OF_PI) / (2 * T6::MASS_OF_K);
        Double_t PiMomentum = Sqrt(PiEnergy * PiEnergy - T6::MASS_OF_PI * T6::MASS_OF_PI);

        Double_t PiPlusPhi = 2 * Pi() * rnd.Rndm();
        Double_t PiPlusCosTheta = 2 * rnd.Rndm() - 1;
        Double_t PiPlusSinTheta = Sqrt(1 - PiPlusCosTheta * PiPlusCosTheta);

        TLorentzVector PiPlus(
            PiMomentum * PiPlusSinTheta * Cos(PiPlusPhi),
            PiMomentum * PiPlusSinTheta * Sin(PiPlusPhi),
            PiMomentum * PiPlusCosTheta,
            PiEnergy
        );

        TLorentzVector PiMinus(
            -PiPlus.Px(),
            -PiPlus.Py(),
            -PiPlus.Pz(),
            PiEnergy
        );

        PiPlus = T * PiPlus;
        PiMinus = T * PiMinus;

        histKsPhi->Fill(Phi);
        histKsTheta->Fill(ACos(CosTheta));
        if (PiPlus.Vect().Mag() > 0 && PiMinus.Vect().Mag() > 0 && Ks.Vect().Mag() > 0) {
            double anglePlus = Ks.Vect().Angle(PiPlus.Vect());
            double angleMinus = Ks.Vect().Angle(PiMinus.Vect());

            histBTW->Fill(anglePlus);
            histBTW->Fill(angleMinus);
        } else {
            std::cerr << "Warning: One of the vectors is zero. Skipping angle calculation.\n";
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Histograms", 1500, 800);
    c1->Divide(2, 3);

    c1->cd(1);
    histKsLength->Draw();

    c1->cd(2);
    histKsTheta->Draw();

    c1->cd(3);
    histKsPhi->Draw();

    c1->cd(4);
    histBTW->Draw();
}

int main() {
    task6();
    return 0;
}
